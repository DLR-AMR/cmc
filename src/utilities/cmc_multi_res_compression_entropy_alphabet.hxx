#ifndef CMC_MULTI_RES_COMPRESSION_ENTROPY_ALPHABET_HXX
#define CMC_MULTI_RES_COMPRESSION_ENTROPY_ALPHABET_HXX

#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_iface_entropy_alphabet.hxx"

#include <vector>

namespace cmc::entropy_coding::arithmetic_coding
{

constexpr uint32_t kByteCompressionSignumBit = 0x40000000;

constexpr inline bool
CheckIfCompressionSignumBitIsSet(const uint32_t symbol)
{
    return ((symbol >> (sizeof(uint32_t) * bit_map::kCharBit - 2)) & uint32_t{1});
}

template <typename T>
constexpr static size_t GetMultiResAlphabetCount()
{
    return 2 * (sizeof(T) * bit_map::kCharBit + 1) + 1;
}

/* This alphebet is in particular tailored to the byte compression techniques and their symbols (LZC, prefix lengths)
 * The encoding of other symbols won't work. Therefore, the default ArithmeticEncoderAlphabet is needed */
template <typename T>
class MultiResCompressionAlphabet : public IByteCompressionEntropyAlphabet
{
public:
    size_t GetAlphabetSize() const override {return GetMultiResAlphabetCount<T>();};

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
    
    const std::vector<uint32_t>& GetSymbolFrequencies() const override
    {
        return alphabet_;
    }

#ifdef CMC_ENABLE_MPI
    std::vector<uint32_t> CommunicateSymbolFrequencies(const MPI_Comm comm) override;
#endif

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
MultiResCompressionAlphabet<T>::DoesSymbolExistInAlphabet(const uint32_t symbol)
{
    if (symbol < kByteCompressionSignumBit)
    {
        if (symbol <= sizeof(T) * bit_map::kCharBit)
        {
            return true;
        } else 
        {
            return false;
        }
    } else
    {
        if (symbol == kByteCompressionSymbolJumpToNextByte)
        {
            return true;
        }
        if (symbol >= kByteCompressionSignumBit && symbol <= kByteCompressionSignumBit + sizeof(T) * bit_map::kCharBit)
        {
            return true;
        } else
        {
            return false;
        }
    }
}

template <typename T>
uint32_t
MultiResCompressionAlphabet<T>::TransformSymbolToArrayIndex(const uint32_t symbol)
{
    /* Define the offset if the signum is set */
    constexpr uint32_t offset = sizeof(T) * bit_map::kCharBit + 1;
    
    if (CheckIfCompressionSignumBitIsSet(symbol))
    {
        /* Return the index for the symbol with a set signum bit */
        return (symbol - kByteCompressionSignumBit) + offset;
    } else if (CheckIfJumpToNextByteIndicatorBitIsSet(symbol))
    {
        /* Return the last accessible entry of the array*/
        return (GetAlphabetSize() - 1);
    } else
    {
        return symbol;
    }
}

template <typename T>
uint32_t
MultiResCompressionAlphabet<T>::RevertArrayIndexToSymbol(const uint32_t index)
{
    cmc_assert(index < static_cast<uint32_t>(GetAlphabetSize()));
    
    /* Define the offset if the signum is set */
    constexpr uint32_t offset = sizeof(T) * bit_map::kCharBit + 1;

    if (index < offset)
    {
        /* If the index is smaller than the offset, the index equals the symbol */
        return index;
    } else
    {
        if (index == GetAlphabetSize() - 1)
        {
            /* If it is the last index, we return the kByteCompressionSymbolJumpToNextByte */
            return kByteCompressionSymbolJumpToNextByte;
        } else
        {
            /* In case the index is greater/equal to the offset, we need to apply de-offset the value and apply te signum */
            return kByteCompressionSignumBit + (index - offset);
        }
    }
}

#ifdef CMC_ENABLE_MPI
/**
 * @brief The collected symbol frequencies are exchanged, such that each process holds the same frequencies.
 * This is important, because the decompression can be performed with a different amount of processes/distributions.
 * 
 * @param comm The communicator to use in order to exhcange the symbol frequencies 
 */
template <typename T>
std::vector<uint32_t>
MultiResCompressionAlphabet<T>::CommunicateSymbolFrequencies(const MPI_Comm comm)
{
#ifdef CMC_ENABLE_MPI

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

#endif

}


#endif /* !CMC_MULTI_RES_COMPRESSION_ENTROPY_ALPHABET_HXX */
