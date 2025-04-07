#ifndef CMC_BYTE_COMPRESSION_ENTROPY_ALPHABET_HXX
#define CMC_BYTE_COMPRESSION_ENTROPY_ALPHABET_HXX

#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_iface_entropy_alphabet.hxx"

#include <array>

namespace cmc::entropy_coding::arithmetic_coding
{

constexpr uint32_t kByteCompressionSignumBit = 0x40000000;
constexpr uint32_t kByteCompressionSymbolJumpToNextByte = 0x80000000;


constexpr inline bool
CheckIfCompressionSignumBitIsSet(const uint32_t symbol)
{
    return ((symbol >> (sizeof(uint32_t) * bit_map::kCharBit - 2)) & uint32_t{1});
}

constexpr inline bool
CheckIfJumpToNextByteIndicatorBitIsSet(const uint32_t symbol)
{
    return ((symbol >> (sizeof(uint32_t) * bit_map::kCharBit - 1)) & uint32_t{1});
}

template <typename T>
constexpr bool
DoesSymbolExistsInAlphabet(const uint32_t symbol)
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

/* This alphebet is in particular tailored to the byte compression techniques and their symbols (LZC, prefix lengths)
 * The encoding of other symbols won't work. Therefore, the default ArithmeticEncoderAlphabet is needed */
template <typename T>
class ByteCompressionAlphabet : public IEntropyAlphabet
{
public:
    constexpr static size_t N = 2 * (sizeof(T) * bit_map::kCharBit + 1) + 1;

    void InitializeSymbols([[maybe_unused]] const size_t type_size_in_bytes) override
    {
        cmc_assert(type_size_in_bytes == sizeof(T));
        alphabet_.fill(uint32_t{0});
        /* We fill each frequency count with a zero except the symbol indicating process boundaries */
        UpdateSymbol(kByteCompressionSymbolJumpToNextByte);
    };

    void UpdateSymbol(const uint32_t symbol) override
    {
        ++(alphabet_[TransformSymbolToArrayIndex<T>(symbol)]);
    };
    
    const std::array<uint32_t, N>&  GetSymbolFrequencies() const
    {
        return alphabet_;
    }

    size_t size() const {return alphabet_.size();};
    auto begin() const  {return alphabet_.begin();};
    auto end() const  {return alphabet_.end();};

private:
    std::array<uint32_t, N> alphabet_;
};


template <typename T>
constexpr uint32_t
TransformSymbolToArrayIndex(const uint32_t symbol)
{
    /* Define the offset if the signum is set */
    constexpr uint32_t offset = sizeof(T) * bit_map::kCharBit + 1;
    
    if (CheckIfCompressionSignumBitIsSet(symbol))
    {
        /* Return the index for the symbols with a set signum bit */
        return (symbol - kByteCompressionSignumBit) + offset;
    } else if (CheckIfJumpToNextByteIndicatorBitIsSet(symbol))
    {
        /* Return the last accessible entry of the array*/
        return (ByteCompressionAlphabet<T>::N - 1);
    } else
    {
        return symbol;
    }
}

template <typename T>
constexpr uint32_t
RevertArrayIndexToSymbol(const uint32_t index)
{
    cmc_assert(index < static_cast<uint32_t>(ByteCompressionAlphabet<T>::N));
    
    /* Define the offset if the signum is set */
    constexpr uint32_t offset = sizeof(T) * bit_map::kCharBit + 1;

    if (index < offset)
    {
        /* If the index is smaller than the offset, the index equals the symbol */
        return index;
    } else
    {
        if (index == ByteCompressionAlphabet<T>::N - 1)
        {
            /* If it is the last index, we return the kByteCompressionSymbolJumpToNextByte */
            return kByteCompressionSymbolJumpToNextByte;
        } else
        {
            /* In case the index is greater/eual to the offset, we need to apply de-offset the value and apply te signum */
            return kByteCompressionSignumBit + (index - offset);
        }
    }
}


}


#endif /* !CMC_BYTE_COMPRESSION_ENTROPY_ALPHABET_HXX */
