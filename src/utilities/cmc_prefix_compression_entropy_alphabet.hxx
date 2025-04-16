#ifndef CMC_PREFIX_COMPRESSION_ENTROPY_ALPHABET_HXX
#define CMC_PREFIX_COMPRESSION_ENTROPY_ALPHABET_HXX

#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_iface_entropy_alphabet.hxx"

#include <vector>

namespace cmc::entropy_coding::arithmetic_coding
{

template <typename T>
constexpr static size_t GetPrefixAlphabetCount()
{
    return (sizeof(T) * bit_map::kCharBit + 2);
}

/* This alphebet is in particular tailored to the byte compression techniques and their symbols (in this case: prefix lengths)
 * The encoding of other symbols won't work. Therefore, the default ArithmeticEncoderAlphabet is needed */
template <typename T>
class PrefixCompressionAlphabet : public IByteCompressionEntropyAlphabet
{
public:
    size_t GetAlphabetSize() const override {return GetPrefixAlphabetCount<T>();};

    void InitializeSymbols([[maybe_unused]] const size_t type_size_in_bytes) override
    {
        cmc_assert(type_size_in_bytes == sizeof(T));
        alphabet_ = std::vector<uint32_t>(GetAlphabetSize(), 0);

        /* We fill each frequency count with a zero except the symbol indicating process boundaries */
        UpdateSymbol(kByteCompressionSymbolJumpToNextByte);
    };

    void UpdateSymbol(const uint32_t symbol) override
    {
        ++(alphabet_[TransformSymbolToArrayIndex(symbol)]);
    };
    
    const std::vector<uint32_t>&  GetSymbolFrequencies() const override
    {
        return alphabet_;
    }

    std::vector<uint32_t> CommunicateSymbolFrequencies(const MPI_Comm comm) override;

    bool DoesSymbolExistInAlphabet(const uint32_t symbol) override;
    uint32_t TransformSymbolToArrayIndex(const uint32_t symbol) override;
    uint32_t RevertArrayIndexToSymbol(const uint32_t index) override;

    size_t size() const {return alphabet_.size();};
    auto begin() const  {return alphabet_.begin();};
    auto end() const  {return alphabet_.end();};

private:
    std::vector<uint32_t> alphabet_;
};


template <typename T>
bool
PrefixCompressionAlphabet<T>::DoesSymbolExistInAlphabet(const uint32_t symbol)
{
    if (symbol >= GetAlphabetSize() && symbol != kByteCompressionSymbolJumpToNextByte)
    {
        return false;
    }
    
    return true;
}

template <typename T>
uint32_t
PrefixCompressionAlphabet<T>::TransformSymbolToArrayIndex(const uint32_t symbol)
{
    cmc_assert(DoesSymbolExistInAlphabet(symbol));
    
    if (not CheckIfJumpToNextByteIndicatorBitIsSet(symbol))
    {
        /* The symbol corresponds to the array index */
        return symbol;
    } else
    {
        /* Return the last accessible entry of the array */
        return GetAlphabetSize() - 1;
    }
}

template <typename T>
uint32_t
PrefixCompressionAlphabet<T>::RevertArrayIndexToSymbol(const uint32_t index)
{
    cmc_assert(index < static_cast<uint32_t>(GetAlphabetSize()));
    
    if (index < GetAlphabetSize() - 1)
    {
        /* If the index directly corresponds to the symbol */
        return index;
    } else
    {
        /* If it is the last index, we return the kByteCompressionSymbolJumpToNextByte */
        return kByteCompressionSymbolJumpToNextByte;
    }
}

/**
 * @brief The collected symbol frequencies are exchanged, such that each process holds the same frequencies.
 * This is important, because the decompression can be performed with a different amount of processes/distributions.
 * 
 * @param comm The communicator to use in order to exhcange the symbol frequencies 
 */
template <typename T>
std::vector<uint32_t>
PrefixCompressionAlphabet<T>::CommunicateSymbolFrequencies(const MPI_Comm comm)
{
#if CMC_ENABLE_MPI

    cmc_debug_msg("The symbol frequencies of the entropy coding alphabet will be exchanged.");
    
    /* Get the local symbol frequencies */
    const std::vector<uint32_t>& local_symbol_frequencies = alphabet_;

    /* Allocate the global frequencies */
    std::vector<uint32_t> global_symbol_frequencies(GetAlphabetSize());

    /* Exchange the symbol frequencies */
    const int ret_val = MPI_Allreduce(local_symbol_frequencies.data(), global_symbol_frequencies.data(), static_cast<int>(GetAlphabetSize()), MPI_UINT32_T, MPI_SUM, comm);
    MPICheckError(ret_val);

    cmc_debug_msg("The symbol frequencies of the entropy coding alphabet have been successfully exchanged.");

    return global_symbol_frequencies;
#else
    cmc_warn_msg("MPI-Communication of the entropy alphabet is called although cmc is not linked against MPI. However, the serial function call is perfomed.");

    return alphabet_;
#endif
}

}


#endif /* !CMC_PREFIX_COMPRESSION_ENTROPY_ALPHABET_HXX */
