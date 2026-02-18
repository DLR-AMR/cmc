#ifndef CMC_PAR_MULTI_RES_EXTRACTION_UTIL_HXX
#define CMC_PAR_MULTI_RES_EXTRACTION_UTIL_HXX

#include "utilities/cmc_bits_vector.hxx"

#include <array>
#include <limits>

namespace cmc::par::lossless::multi_res
{

template<typename T>
concept ArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T>);

template<typename T>
concept OneByteArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T> && (sizeof(T) == 1));

template<typename T>
concept TwoByteArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T> && (sizeof(T) == 2));

template<typename T>
concept FourByteArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T> && (sizeof(T) == 4));

template<typename T>
concept EightByteArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T> && (sizeof(T) == 8));

template<typename T>
concept UnsignedIntegerType = (std::is_unsigned_v<T> && std::is_integral_v<T>);

using SymbolType = uint8_t;

inline constexpr int32_t kMaxPresentElementLevelUnknown = -1;

inline constexpr int32_t kNumMaxChildrenElements = 8;

constexpr inline uint8_t kResidualSignumIndication = 0x80;

template<ArithmeticType T>
constexpr inline uint8_t kProcessEndSymbol = kResidualSignumIndication + sizeof(T) * cmc::bits::kCharBit;

constexpr inline
template<UnsignedIntegerType T>
constexpr SymbolType
CreateEntropySymbol(const bool is_approx_greater, const T residual)
{
    /* Check if the residual is zero, in this case we only store a single entropy code for +/- 0 */
    if (cmc::bits::GetLZC(residual) == sizeof(T) * cmc::bits::kCharBit) [[unlikely]]
    {
        return static_cast<SymbolType>(cmc::bits::GetLZC(residual));
    }

    return ((SymbolType{is_approx_greater} << 7) | static_cast<SymbolType>(cmc::bits::GetLZC(residual)));
}

template<ArithmeticType T>
constexpr inline int
MapEntropySymbolToArrayIndex(const uint8_t entropy_symbol)
{
    return (entropy_symbol & SymbolType{0x7F}) + (entropy_symbol >> 7) * (sizeof(T) * cmc::bits::kCharBit + 1);
}

template<ArithmeticType T>
constexpr inline SymbolType
MapArrayIndexToEntropySymbol(int array_idx)
{
    return static_cast<SymbolType>(array_idx + (array_idx > sizeof(T) * cmc::bits::kCharBit ? kResidualSignumIndication - sizeof(T) * cmc::bits::kCharBit - 1 : 0));
}

template<ArithmeticType T>
constexpr inline int
GetNumEntropySymbols()
{
    return 2 * (sizeof(T) * cmc::bits::kCharBit + 1);
}

constexpr inline SymbolType
GetLZCFromEntropySymbol(const SymbolType entropy_symbol)
{
    return (entropy_symbol & SymbolType{0x7F});
}

template<ArithmeticType T>
constexpr inline void
AddProcessEndSymbol(std::array<uint64_t>& entropy_symbols_frequency, const uint64_t num_local_proc_end_symbols)
{
    /* We store the process-end-symbol in the last array entry */
    entropy_symbols_frequency[sizeof(T) * cmc::bits::kCharBit + kResidualSignumIndication] = num_local_proc_end_symbols;
}

template<OneByteType T>
struct LevelEncodingData
{
    uint8_t num_elements;
    std::array<SymbolType, kNumMaxChildrenElements> entropy_symbols;
    std::array<uint8_t, kNumMaxChildrenElements> residuals;
};

template<TwoByteType T>
struct LevelEncodingData
{
    uint16_t num_elements;
    std::array<SymbolType, kNumMaxChildrenElements> entropy_symbols;
    std::array<uint16_t, kNumMaxChildrenElements> residuals;
};

template<FourByteType T>
struct LevelEncodingData
{
    uint32_t num_elements;
    std::array<SymbolType, kNumMaxChildrenElements> entropy_symbols;
    std::array<uint32_t, kNumMaxChildrenElements> residuals;
};

template<EightByteType T>
struct LevelEncodingData
{
    uint64_t num_elements;
    std::array<SymbolType, kNumMaxChildrenElements> entropy_symbols;
    std::array<uint64_t, kNumMaxChildrenElements> residuals;
};

template<ArithmeticType T>
constexpr T
ComputeArithmeticMean(const T* values, const int num_values)
{
    T sum = 0;
    for (int idx{0}; idx < num_values; ++idx)
    {
        sum += values[idx];
    }
    return sum / values.size();
}

template<ArithmeticType T>
constexpr T
ComputeMidRange(const T* values, const int num_values)
{
    cmc_assert(values.size() >= 2);

    T min = std::numeric_limits<T>::max();
    /* Find the minimum */
    for (int idx{0}; idx < num_values; ++idx)
    {
        if (min > values[idx])
        {
            min = values[idx];
        }
    }

    T max = std::numeric_limits<T>::lowest();
    /* Find the maximum */
    for (int idx{0}; idx < num_values; ++idx)
    {
        if (max < values[idx])
        {
            max = values[idx];
        }
    }

    /* Compute the mid-range*/
    return ((max / 2) + (min / 2));
}

}

#endif /* !CMC_PAR_MULTI_RES_EXTRACTION_UTIL_HXX */
