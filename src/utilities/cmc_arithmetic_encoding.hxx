#ifndef CMC_ARITHMETIC_ENCODING_HXX
#define CMC_ARITHMETIC_ENCODING_HXX

#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_arithmetic_encoding_frequency_model.hxx"
#include "utilities/cmc_entropy_alphabet.hxx"
#include "utilities/cmc_entropy_coder.hxx"

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
    Alphabet GetSymbolFrequencies() override
    {
        return alphabet_;
    }

    size_t size() const {return alphabet_.size();};
    auto begin() const  {return alphabet_.begin();};
    auto end() const  {return alphabet_.end();};

private:
    Alphabet alphabet_;
};

class Encoder : public IEntropyCoder
{
public:
    Encoder(){
        alphabet_ = std::make_unique<ArithmeticEncoderAlphabet>();
    };

    void InitializeAlphabet(const size_t type_size_in_bytes = 2) override {
        alphabet_->InitializeSymbols(type_size_in_bytes);
    };

    void UpdateSymbolFrequency(const uint32_t symbol) override {
        alphabet_->UpdateSymbol(symbol);
    };

    void SetupEncoding() override
    {
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
        /* Create the frequency model */
        frequency_model_ = std::make_unique<StaticFrequencyModel>(vectorized_alphabet);
    }

    void EncodeSymbol(const uint32_t symbol) override;
    void FinishEncoding();

    bit_vector::BitVector EncodeAlphabet() const override {return frequency_model_->EncodeAlphabet();};

    bit_map::BitMap GetEncodedBitStream() const override {return encoded_stream_;}

    void ClearEncodedBitStream() override {
        lower_ = 0;
        higher_ = 0x7FFFFFFF;
        step_ = 0;
        scale_ = 0;
        encoded_stream_ = bit_map::BitMap();
        };

private:
    void Encode(const Uint32_t low, const Uint32_t high, const Uint32_t total);
    /* Initial maximum range spans 31 bits -> [0x00000000, 0x7FFFFFFF) */
    Uint32_t lower_{0};
    Uint32_t higher_{0x7FFFFFFF};
    Uint32_t step_{0};
    Uint32_t scale_{0};

    bit_map::BitMap encoded_stream_;
};


inline void
Encoder::EncodeSymbol(const uint32_t symbol)
{
    /* Get the cumulative and total frequency counts */
    const uint32_t low_count = frequency_model_->GetCumulativeCountOfAllSymbolFrequenciesLowerThan(symbol);
    const uint32_t high_count = frequency_model_->GetCumulativeCountOfAllSymbolFrequenciesIncluding(symbol);
    const uint32_t total_count = frequency_model_->GetTotalSymbolFrequencyCount();

    /* Encode the symbol based on its interval */
    Encode(low_count, high_count, total_count);
}

inline void
Encoder::Encode(const Uint32_t low, const Uint32_t high, const Uint32_t total)
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

inline void
Encoder::FinishEncoding()
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

    /* Enforce a minimum encoded stream length */
    if (encoded_stream_.size() < kNumWorkingBits)
    {
        for (size_t i = encoded_stream_.size(); i < kNumWorkingBits; ++i)
        {
            encoded_stream_.AppendUnsetBit();
        }
    }

    cmc_debug_msg("The arithmetic encoder encoded the message in ", encoded_stream_.size(), " bits (=> ", encoded_stream_.size_bytes(), " bytes).");
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

inline uint32_t
GetLeadingZeroCount(const uint32_t encoded_lzc)
{
    if(kMSBBit > encoded_lzc)
    {
        return encoded_lzc;
    } else
    {
        return encoded_lzc - kMSBBit;
    }
}

enum ResidualOperation {IntegerAddition, IntegerSubtraction};

inline
std::pair<ResidualOperation, uint32_t>
DecodeLZC(const uint32_t encoded_lzc)
{
    if(kMSBBit > encoded_lzc)
    {
        return std::make_pair(ResidualOperation::IntegerAddition, encoded_lzc);
    } else
    {
        return std::make_pair(ResidualOperation::IntegerSubtraction, encoded_lzc - kMSBBit);
    }
}

}

#endif /* !CMC_ARITHMETIC_ENCODING_HXX */
