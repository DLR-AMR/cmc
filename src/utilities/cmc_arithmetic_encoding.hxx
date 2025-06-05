#ifndef CMC_ARITHMETIC_ENCODING_HXX
#define CMC_ARITHMETIC_ENCODING_HXX

#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_arithmetic_encoding_frequency_model.hxx"
#include "utilities/cmc_iface_entropy_alphabet.hxx"
#include "utilities/cmc_entropy_coder.hxx"

#include "mpi/cmc_mpi.hxx"

#include <cstddef>
#include <memory>
#include <vector>

namespace cmc::entropy_coding::arithmetic_coding
{


/* The implementation of the arithmetic encoder follows the guide provided by
    @techreport{TRArithmetic,
    	Author = {Eric Bodden and Malte Clasen and Joachim Kneis},
    	Institution = {Sable Research Group, McGill University},
    	Month = {May},
    	Number = {2007-5},
    	Title = {Arithmetic Coding revealed - A guided tour from theory to praxis},
    	Url = {https://www.bodden.de/pubs/sable-tr-2007-5.pdf},
    	Year = {2007},
    	Bdsk-Url-1 = {https://www.bodden.de/pubs/sable-tr-2007-5.pdf}
    }
*/

constexpr uint32_t kMSBBit = 0x80000000;

constexpr uint32_t kSymbolJumpToNextByte = 0xFFFFFFFF;
constexpr uint32_t kSymbolPlaceholder = 0xEFFFFFFF;

constexpr int kAdditionalAlphabetSymbols = 4;
constexpr int kSymbolFrequencyMPITag = 0x01;

/* It works on the least 31 significant bits of an uint32_t */
const size_t kNumWorkingBits = sizeof(Uint32_t) * bit_map::kCharBit - 1;

using Uint32_t = uint32_t; 
using Uint8_t = uint8_t;

const Uint32_t kFirstQ =  0x20000000;
const Uint32_t kMid =  0x40000000;
const Uint32_t kThirdQ =  0x60000000;

class ArithmeticEncoderAlphabet : public IEntropyAlphabet
{
public:
    void InitializeSymbols(const size_t type_size_in_bytes) override
    {
        alphabet_.reserve(type_size_in_bytes * bit_map::kCharBit + 1);
    };
    void UpdateSymbol(const uint32_t symbol) override
    {
        auto alphaiter = alphabet_.find(symbol);
        if (alphaiter == alphabet_.end())
        {
            alphabet_.insert({symbol, 1});
        } else
        {
            ++(alphabet_[symbol]);
        }
    };
    Alphabet GetSymbolFrequencies()
    {
        return alphabet_;
    }

    size_t size() const {return alphabet_.size();};
    auto begin() const  {return alphabet_.begin();};
    auto end() const  {return alphabet_.end();};

private:
    Alphabet alphabet_;
};

template <typename T>
class Encoder : public IEntropyCoder
{
public:
    const int N = sizeof(T);

    Encoder() = default;

    void InitializeAlphabet(const size_t type_size_in_bytes = 2) override {
        alphabet_ = std::make_unique<ArithmeticEncoderAlphabet>();
        alphabet_->InitializeSymbols(type_size_in_bytes);
    };

    void UpdateSymbolFrequency(const uint32_t symbol) override {
        alphabet_->UpdateSymbol(symbol);
    };

    void SetupEncoding(const MPI_Comm comm) override
    {
        /* Communicate the alphabet */
        std::vector<Letter> vectorized_alphabet = this->CommunicateSymbolFrequencies(comm);
        
        /* Create the frequency model */
        frequency_model_ = std::make_unique<StaticFrequencyModel>(vectorized_alphabet);
    }

    void EncodeSymbol(const uint32_t symbol) override;
    void FinishEncoding();

    bit_vector::BitVector EncodeAlphabet() const override {return frequency_model_->EncodeAlphabet();};

    bit_map::BitMap GetEncodedBitStream() const override {return encoded_stream_;}

    void ClearEncodedBitStream() override
    {
        lower_ = 0;
        higher_ = 0x7FFFFFFF;
        step_ = 0;
        scale_ = 0;
        encoded_stream_ = bit_map::BitMap();
    };

    void Reset() override
    {
        this->ClearEncodedBitStream();
        alphabet_.reset(nullptr);
        frequency_model_.reset(nullptr);
    };

private:
    void Encode(const Uint32_t low, const Uint32_t high, const Uint32_t total);
    void ProcessCommunicatedSymbolFrequencies(Alphabet& update_alphabet, const std::vector<Letter>& received_alphabet);
    void RemovePlaceholderSymbolsFromAlphabet(std::vector<Letter>& broadcasted_alphabet);
    std::vector<Letter> CommunicateSymbolFrequencies(const MPI_Comm comm);

    /* Initial maximum range spans 31 bits -> [0x00000000, 0x7FFFFFFF) */
    Uint32_t lower_{0};
    Uint32_t higher_{0x7FFFFFFF};
    Uint32_t step_{0};
    Uint32_t scale_{0};

    bit_map::BitMap encoded_stream_;
};

template <typename T>
inline void
Encoder<T>::EncodeSymbol(const uint32_t symbol)
{
    /* Get the cumulative and total frequency counts */
    const uint32_t low_count = frequency_model_->GetCumulativeCountOfAllSymbolFrequenciesLowerThan(symbol);
    const uint32_t high_count = frequency_model_->GetCumulativeCountOfAllSymbolFrequenciesIncluding(symbol);
    const uint32_t total_count = frequency_model_->GetTotalSymbolFrequencyCount();

    /* Encode the symbol based on its interval */
    Encode(low_count, high_count, total_count);
}

template <typename T>
inline void
Encoder<T>::Encode(const Uint32_t low, const Uint32_t high, const Uint32_t total)
{
    cmc_assert(total < (1 << 29));

    /* Calculate a uniform step size for the range. This potentially prevents overflows in a multiplication of the range and \a high (\see https://www.sable.mcgill.ca/~ebodde/pubs/sable-tr-2007-5.pdf) */
    step_ = (higher_ - lower_ + 1) / total;

    /* Adjust the open upper bound of the interval */
    higher_ = lower_ + step_ * high - 1;

    /* Adjust the lower bound of the interval */
    lower_ = lower_ + step_ * low;

    /* Apply the E1/E2 scaling. In case the current interval lies completely wihtin one half (either smaller or greater than "kMid",
     * we will not leave this half anymore which is why we can already output that information and remove it from further consideration) */
    while (higher_ < kMid || lower_ >= kMid)
    {
        if (higher_ < kMid)
        {
            /* Apply E1 scaling */
            encoded_stream_.AppendUnsetBit();
            /* Adjust lower and upper bound by shifting the bounds */
            lower_ <<= 1;
            higher_ = (higher_ << 1) + 1;

            /* Potentially, perform E3 scaling, when a carry needs to be added */
            for (; scale_ > 0; --scale_)
            {
                encoded_stream_.AppendSetBit();
            }
        } else if (lower_ >= kMid)
        {
            /* Apply E2 scaling */
            encoded_stream_.AppendSetBit();

            /* Adjust lower and upper bound by shifting the reduced bounds */
            lower_ = (lower_ - kMid) << 1;
            higher_ = ((higher_ - kMid) << 1) + 1;

            /* Potentially, perform E3 scaling, when a carry needs to be added */
            for (; scale_ > 0; --scale_)
            {
                encoded_stream_.AppendUnsetBit();
            }
        }
    }

    /* Check for E3 scaling which applies when the bounds converges around the center of the interval,
     * but it is not yet determinable in which the half the range will fall. The correct bits are added
     * as soon as the next E1 or E2 scaling is applied. */
    while (kFirstQ <= lower_ && higher_ < kThirdQ)
    {
        /* Count the necessary E3 scalings */
        ++scale_;

        /* Adjust the bounds to move out of the E3 case */
        lower_ = (lower_ - kFirstQ) << 1;
        higher_ = ((higher_ - kFirstQ) << 1) + 1; 
    }
}

template <typename T>
inline void
Encoder<T>::FinishEncoding()
{
    if (lower_ < kFirstQ)
    {
        /* Indicate the case that lower inerval bound is below the first quarter */
        encoded_stream_.AppendUnsetBit();
        /* Perform the pending E3 scaling */
        for (Uint32_t i = 0; i < scale_ + 1; ++i)
        {
            encoded_stream_.AppendSetBit();
        }
    } else
    {
        encoded_stream_.AppendSetBit();
    }

    /* Empace some buffer bits */
    encoded_stream_.FillCurrentByte();
    encoded_stream_.AppendUnsetBit();

    /* Enforce a minimum encoded stream length */
    if (encoded_stream_.size() <= kNumWorkingBits)
    {
        for (size_t i = encoded_stream_.size(); i <= kNumWorkingBits; ++i)
        {
            encoded_stream_.AppendUnsetBit();
        }
    }

    cmc_debug_msg("The arithmetic encoder encoded the message in ", encoded_stream_.size(), " bits (=> ", encoded_stream_.size_bytes(), " bytes).");
}

template <typename T>
inline void
Encoder<T>::ProcessCommunicatedSymbolFrequencies(Alphabet& update_alphabet, const std::vector<Letter>& received_alphabet)
{
    for (auto al_iter = received_alphabet.begin(); al_iter != received_alphabet.end(); ++al_iter)
    {
        /* Check if the symbol is already in the alphabet */
        auto sym_iter = update_alphabet.find(al_iter->symbol);

        /* Add the received frequency for the symbol */
        if (sym_iter != update_alphabet.end())
        {
            /* Add the frequency for the symbol */
            update_alphabet[al_iter->symbol] += al_iter->frequency;
        } else
        {
            /* Add the new symbol to the alphabet */
            update_alphabet[al_iter->symbol] = al_iter->frequency;
        }
    }
}

template <typename T>
inline void
Encoder<T>::RemovePlaceholderSymbolsFromAlphabet(std::vector<Letter>& broadcasted_alphabet)
{
    /* The placeholder elements are emplace at the end of the vectorized alphabet.
     * Therefore, we search for the first occurence of the placeholder and remove the remainder */
    size_t num_non_placeholder_symbols{0};

    for (size_t letter_idx = 0; letter_idx < broadcasted_alphabet.size(); ++letter_idx)
    {
        if (broadcasted_alphabet[letter_idx].symbol == kSymbolPlaceholder)
        {
            /* In case we encounter a placeholder value, we break the loop and will reduce the vector */
            break;
        } else
        {
            ++num_non_placeholder_symbols;
        }
    }

    /* Resize the vector to the non-placeholder values */
    broadcasted_alphabet.resize(num_non_placeholder_symbols);
}

/**
 * @brief The collected symbol frequencies are exchanged, such that each process holds the same frequencies.
 * This is important, because the decompression can be performed with a different amount of processes/distributions.
 * 
 * @param comm The communicator to use in order to exhcange the symbol frequencies 
 */
template <typename T>
inline std::vector<Letter>
Encoder<T>::CommunicateSymbolFrequencies(const MPI_Comm comm)
{
#if CMC_ENABLE_MPI

    cmc_debug_msg("The symbol frequencies of the entropy coding alphabet will be exchanged.");

    /* Define the maximum permitted alphabet size */
    const int kMaxAlphabetSize = N * bit_map::kCharBit + kAdditionalAlphabetSymbols;

    /* Get the size of the communicator */
    int comm_size{1};
    const int ret_val_comm_size = MPI_Comm_size(comm, &comm_size);
    MPICheckError(ret_val_comm_size);

    /* Define the properties of the custom 'Letter' data type */
    const int num_fields = 2;
    int array_of_blocklengths[] = {1,1};
    MPI_Aint array_of_displacements[num_fields];
    array_of_displacements[0] = offsetof(Letter, symbol);
    array_of_displacements[1] = offsetof(Letter, frequency);
    MPI_Datatype array_of_types[] = {MPI_INT, MPI_INT};

    /* Create the actual MPI data type */
    MPI_Datatype MPI_Letter;
    const int ret_val_letter_struct = MPI_Type_create_struct(num_fields, array_of_blocklengths, array_of_displacements,
                                                             array_of_types, &MPI_Letter); 
    MPICheckError(ret_val_letter_struct);
    const int ret_val_type_commit = MPI_Type_commit(&MPI_Letter);
    MPICheckError(ret_val_type_commit); 

    /* Get the rank of the process within the communicator */
    int rank{0};
    const int ret_val_rank = MPI_Comm_rank(comm, &rank);
    MPICheckError(ret_val_rank);

    /* Define the process to collect the symbol frequencies */
    const int root_rank = 0;

    /* Cast the base pointer to the actual implementation to construct the needed vectorized alphabet */
    const ArithmeticEncoderAlphabet* const  alphabet = dynamic_cast<ArithmeticEncoderAlphabet*>(alphabet_.get());
    
    /* We create a vector of Letters which will be filled and sent */
    std::vector<Letter> vectorized_alphabet;

    /* And we create a new alphabt which will hold be updated frequencies */
    Alphabet updated_alphabet = alphabet->GetSymbolFrequencies();

    /* Declare an MPI_Request object for the send */
    MPI_Request send_req;

    /* All processes send the symbol frequencies to the root process */
    if (rank != root_rank)
    {
        if (alphabet->size() > 0)
        {
            /* We only need to send data, when there is data */
            vectorized_alphabet.reserve(alphabet->size());
            
            /* Convert the symbols and frequencies */
            for (auto al_iter = alphabet->begin(); al_iter != alphabet->end(); ++al_iter)
            {
                vectorized_alphabet.emplace_back(Letter{al_iter->first, al_iter->second});
            }

            /* Send the frequencies to the root rank */
            const int ret_val_send = MPI_Isend(vectorized_alphabet.data(), vectorized_alphabet.size(), MPI_Letter,
                                              root_rank, kSymbolFrequencyMPITag, comm, &send_req);
            MPICheckError(ret_val_send);
        }
    } else
    {
        /* The root rank adds the additional symbol indicating process boundaries */
        updated_alphabet[kSymbolJumpToNextByte] = comm_size;
    }

    /* Define a barrier after all send operations in order to be sure, whether all frequencies have been sent */
    MPI_Request barrier_req;
    const int ret_val_barrier = MPI_Ibarrier(comm, &barrier_req);
    MPICheckError(ret_val_barrier);

    /* Allocate a vector for the new alphabet which will be filled by root rank and broadcasted afterwards */
    std::vector<Letter> broadcasted_alphabet(kMaxAlphabetSize, Letter{kSymbolPlaceholder, 0});

    /* On the root rank, we start to receive and process the messages */
    if (rank == root_rank)
    {
        /* We allocate a vector with the upper bound of symbols for 64-bit types */
        std::vector<Letter> received_frequencies(kMaxAlphabetSize);

        bool continue_receiving = true;

        while(continue_receiving)
        {
            int is_message_waiting{0};
            MPI_Status status;

            /* Check if there is a message to receive */
            const int ret_val = MPI_Iprobe(MPI_ANY_SOURCE, kSymbolFrequencyMPITag, comm, &is_message_waiting, &status);
            MPICheckError(ret_val);

            /* If there is a message, we receive it */
            if (is_message_waiting)
            {
                /* Get the count of symbol-frequencies */
                int count{0};
                const int ret_val_count = MPI_Get_count(&status, MPI_Letter, &count);
                MPICheckError(ret_val_count);

                /* Resize the receiption vector */
                received_frequencies.resize(count);

                /* Receive the actual message */
                const int ret_val_recv = MPI_Recv(received_frequencies.data(), count, MPI_Letter, status.MPI_SOURCE, status.MPI_TAG, comm, MPI_STATUS_IGNORE);
                MPICheckError(ret_val_recv);

                /* Process the symbol frequencies */
                ProcessCommunicatedSymbolFrequencies(updated_alphabet, received_frequencies);
            }

            /* Check if the barrier has been reached*/
            int is_barrier_complete;
            const int ret_val_barrier_check = MPI_Test(&barrier_req, &is_barrier_complete, MPI_STATUS_IGNORE);
            MPICheckError(ret_val_barrier_check);

            /* Update the loop flag, we continue until no more messages are queued */
            continue_receiving = is_message_waiting || !(is_barrier_complete);
        }

        /* Check whether the maximum permitted alphabet size is not violated */
        if (updated_alphabet.size() > static_cast<size_t>(kMaxAlphabetSize))
        {
            /* The broadcast operation is limited to the given size */
            cmc_err_msg("The alphabet size (= the amount of symbols) is larger than the maximum permitted alphabet size. Therefore, the alphabet cannot be communicated with the utilized MPI-Broadcast operation.");
        }

        /* After all messages have been received, setup the vectorized alphabet which will be broadcasted */
        int idx = 0;
        for (auto al_iter = updated_alphabet.begin(); al_iter != updated_alphabet.end(); ++al_iter, ++idx)
        {
            /* Convert the symbols and frequencies */
            broadcasted_alphabet[idx] = Letter{al_iter->first, al_iter->second};
        }
    }

    /* After all messages have been transferred and the root rank has updated the symbol frequencies,
     * it broadcasts the new "dictionary", such that each process uses the sam statistical model */
    const int ret_val_bcast = MPI_Bcast(broadcasted_alphabet.data(), kMaxAlphabetSize, MPI_Letter, root_rank, comm);
    MPICheckError(ret_val_bcast);

    /* We remove the palacehodler values from the alphabet whihc were needed by the broadcast */
    RemovePlaceholderSymbolsFromAlphabet(broadcasted_alphabet);

    cmc_debug_msg("The symbol frequencies of the entropy coding alphabet have been successfully exchanged.");

    return broadcasted_alphabet;
#else
    cmc_warn_msg("MPI-Communication of the entropy alphabet is called although cmc is not linked against MPI. However, the serial function call is perfomed.");

    /* Cast the base pointer to the actual implementation to construct the needed vectorized alphabet */
    const ArithmeticEncoderAlphabet* const  alphabet = dynamic_cast<ArithmeticEncoderAlphabet*>(alphabet_.get());
    
    /* The static frequency model takes a vector of Letters as input */
    std::vector<Letter> vectorized_alphabet;
    vectorized_alphabet.reserve(alphabet->size());
    
    /* Convert the symbols and frequencies */
    for (auto al_iter = alphabet->begin(); al_iter != alphabet->end(); ++al_iter)
    {
        vectorized_alphabet.emplace_back(Letter{al_iter->first, al_iter->second});
    }

    /* Additionally add the end of process symbol for completeness */
    vectorized_alphabet.emplace_back(Letter{kSymbolJumpToNextByte, 1});

    return vectorized_alphabet;
#endif
}


class Decoder
{
public:

    Decoder() = delete;
    Decoder(std::unique_ptr<IACModel> frequency_model)
    : frequency_model_{std::move(frequency_model)}{};
    Decoder(std::unique_ptr<IACModel> frequency_model, bit_map::BitMapView encoded_stream)
    : frequency_model_{std::move(frequency_model)}, encoded_stream_{std::move(encoded_stream)}{
        /* Fill the buffer initially */
        for (size_t i = kNumWorkingBits; i > 0; --i)
        {

            buffer_ |= (GetNextBit() << (i - 1));
        }
    };

    Uint32_t DecodeNextSymbol();

    void Reset(bit_map::BitMapView new_encoded_stream)
    {
        lower_ = 0;
        higher_ = 0x7FFFFFFF;
        step_ = 0;
        scale_ = 0;
        buffer_ = 0;
        encoded_stream_ = new_encoded_stream;

        /* Fill the buffer initially with the new stream */
        for (size_t i = kNumWorkingBits; i > 0; --i)
        {
            buffer_ |= (GetNextBit() << (i - 1));
        }
    }

    void ResetAfterProcessBoundary()
    {
        lower_ = 0;
        higher_ = 0x7FFFFFFF;
        step_ = 0;
        scale_ = 0;
        buffer_ = 0;
        encoded_stream_.MoveToNextByte();
        
        /* Fill the buffer initially with the new stream */
        for (size_t i = kNumWorkingBits; i > 0; --i)
        {
            buffer_ |= (GetNextBit() << (i - 1));
        }
    }
private:
    inline Uint32_t GetNextBit() {return (encoded_stream_.GetNextBit() == true ? Uint32_t{1} : Uint32_t{0});};
    void AdvanceDecoder(const Uint32_t low, const Uint32_t high, const Uint32_t total);

    /* Initial maximum range spans 31 bits -> [0x00000000, 0x7FFFFFFF) */
    Uint32_t lower_{0};
    Uint32_t higher_{0x7FFFFFFF};
    Uint32_t step_{0};
    Uint32_t scale_{0};

    Uint32_t buffer_{0};
    std::unique_ptr<IACModel> frequency_model_;
    bit_map::BitMapView encoded_stream_;
};

inline Uint32_t
Decoder::DecodeNextSymbol()
{
    /* Get the total frequency count of all letters */
    const uint32_t total_count = frequency_model_->GetTotalSymbolFrequencyCount();

    /* Compute the uniform split of the number range */
    step_ = (higher_ - lower_ + 1) / total_count;

    /* Get the symbol from the buffer */
    const Uint32_t decoded_symbol_value = (buffer_ - lower_) / step_;

    /* Get the actual symbol from the alphabet of the model */
    const uint32_t symbol = frequency_model_->GetSymbolFromCumulativeFrequency(decoded_symbol_value);

    /* Get the cumulative and total frequency counts */
    const uint32_t low_count = frequency_model_->GetCumulativeCountOfAllSymbolFrequenciesLowerThan(symbol);
    const uint32_t high_count = frequency_model_->GetCumulativeCountOfAllSymbolFrequenciesIncluding(symbol);

    /* Advance the decoder */
    AdvanceDecoder(low_count, high_count, total_count);

    return symbol;
}

inline void
Decoder::AdvanceDecoder(const Uint32_t low, const Uint32_t high, const Uint32_t total)
{
    /* Update the open upper bound of the interval */
    higher_ = lower_ + step_ * high - 1;

    /* Update the lower bound of the interval */
    lower_ = lower_ + step_ * low;

    /* Reverse the E1/E2 scalings */
    while (higher_ < kMid || lower_ >= kMid)
    {
        if (higher_ < kMid)
        {
            lower_ <<= 1;
            higher_ = (higher_ << 1) + 1;
            buffer_ = (buffer_ << 1) | GetNextBit();
        } else if (lower_ >= kMid)
        {
            lower_ = (lower_ - kMid) << 1;
            higher_ = ((higher_ - kMid) << 1) + 1;
            buffer_ = ((buffer_ - kMid) << 1) | GetNextBit();
        }

        /* Reset the E3 scaling counter */
        scale_ = 0;
    }

    /* Account for the E3 scalings */
    while (lower_ >= kFirstQ && higher_ < kThirdQ)
    {
        ++scale_;
        lower_ = (lower_ - kFirstQ) << 1;
        higher_ = ((higher_ - kFirstQ) << 1) + 1;
        buffer_ = ((buffer_ - kFirstQ) << 1) | GetNextBit();
    }
}

}

#endif /* !CMC_ARITHMETIC_ENCODING_HXX */
