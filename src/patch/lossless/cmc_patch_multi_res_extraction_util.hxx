#ifndef CMC_PATCH_LOSSLESS_CMC_PATCH_MULTI_RES_EXTRACTION_UTIL_HXX
#define CMC_PATCH_LOSSLESS_CMC_PATCH_MULTI_RES_EXTRACTION_UTIL_HXX

#include "cmc.hxx"
#include "utilities/cmc_bits.hxx"
#include "utilities/cmc_bits_vector.hxx"
#include "utilities/cmc_huffman_coder.hxx"
#include "utilities/cmc_bits_stream_decoder.hxx"

#include <cstdint>
#include <array>
#include <limits>
#include <execution>

namespace cmc::serial::patch::lossless::multi_res
{

template<int32_t DIM>
concept Dimension = (DIM >= 1 && DIM <= 4);

constexpr int kDimReductionFactor = 2;

using SizeType = uint64_t;
constexpr SizeType kNumCharsVariableName = 256;

template<int32_t DIM>
constexpr int kNumMaxChildrenElements;
template<>
constexpr int kNumMaxChildrenElements<1> = kDimReductionFactor;
template<>
constexpr int kNumMaxChildrenElements<2> = kDimReductionFactor * kDimReductionFactor;
template<>
constexpr int kNumMaxChildrenElements<3> = kDimReductionFactor * kDimReductionFactor * kDimReductionFactor;
template<>
constexpr int kNumMaxChildrenElements<4> = kDimReductionFactor * kDimReductionFactor * kDimReductionFactor * kDimReductionFactor;

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

using OneByteResidualType = uint8_t;

using TwoByteResidualType = uint16_t;

using FourByteResidualType = uint32_t;

using EightByteResidualType = uint64_t;

using SymbolType = uint8_t;

template<typename T, int32_t DIM>
struct PatchEncodingData;

template<typename T, int32_t DIM>
requires Dimension<DIM> && OneByteArithmeticType<T>
struct PatchEncodingData<T, DIM>
{
    T coarse_value{};
    std::array<OneByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
    std::array<SymbolType, kNumMaxChildrenElements<DIM>> entropy_symbols{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && TwoByteArithmeticType<T>
struct PatchEncodingData<T, DIM>
{
    T coarse_value{};
    std::array<TwoByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
    std::array<SymbolType, kNumMaxChildrenElements<DIM>> entropy_symbols{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && FourByteArithmeticType<T>
struct PatchEncodingData<T, DIM>
{
    T coarse_value{};
    std::array<FourByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
    std::array<SymbolType, kNumMaxChildrenElements<DIM>> entropy_symbols{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && EightByteArithmeticType<T>
struct PatchEncodingData<T, DIM>
{
    T coarse_value{};
    std::array<EightByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
    std::array<SymbolType, kNumMaxChildrenElements<DIM>> entropy_symbols{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM>
struct Residuals;

template<typename T, int32_t DIM>
requires Dimension<DIM> && OneByteArithmeticType<T>
struct Residuals<T, DIM>
{
    std::array<OneByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && TwoByteArithmeticType<T>
struct Residuals<T, DIM>
{
    std::array<TwoByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && FourByteArithmeticType<T>
struct Residuals<T, DIM>
{
    std::array<FourByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && EightByteArithmeticType<T>
struct Residuals<T, DIM>
{
    std::array<EightByteResidualType, kNumMaxChildrenElements<DIM>> residuals{};
};

template<OneByteArithmeticType T>
constexpr inline
uint8_t
TransformToUInteger(const T value)
{
    return std::bit_cast<uint8_t>(value);
}

template<TwoByteArithmeticType T>
constexpr inline
uint16_t
TransformToUInteger(const T value)
{
    return std::bit_cast<uint16_t>(value);
}

template<FourByteArithmeticType T>
constexpr inline
uint32_t
TransformToUInteger(const T value)
{
    return std::bit_cast<uint32_t>(value);
}

template<EightByteArithmeticType T>
constexpr inline
uint64_t
TransformToUInteger(const T value)
{
    return std::bit_cast<uint64_t>(value);
}

using SymbolType = uint8_t;

constexpr inline SymbolType kResidualSignumIndication = 0x80;

template<UnsignedIntegerType T>
constexpr inline SymbolType
CreateEntropySymbol(const bool is_approx_greater, const T residual)
{
    /* Check if the residual is zero, in this case we only store a single entropy code for +/- 0 */
    if (cmc::bits::GetLZC(residual) == sizeof(T) * cmc::bits::kCharBit) [[unlikely]]
    {
        return static_cast<SymbolType>(cmc::bits::GetLZC(residual));
    }

    return ((SymbolType{is_approx_greater} << 7) | static_cast<SymbolType>(cmc::bits::GetLZC(residual)));
}

constexpr inline bool
IsApproximationGreater(const SymbolType entropy_symbol)
{
    return (entropy_symbol >> 7);
}

template<ArithmeticType T>
constexpr inline int
MapEntropySymbolToArrayIndex(const SymbolType entropy_symbol)
{
    return (entropy_symbol & SymbolType{0x7F}) + (entropy_symbol >> 7) * (sizeof(T) * cmc::bits::kCharBit + 1);
}

template<ArithmeticType T>
constexpr inline SymbolType
MapArrayIndexToEntropySymbol(int array_idx)
{
    return static_cast<SymbolType>(array_idx + (array_idx > static_cast<int>(sizeof(T) * cmc::bits::kCharBit) ? kResidualSignumIndication - sizeof(T) * cmc::bits::kCharBit - 1 : 0));
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


template<ArithmeticType T, int32_t N>
constexpr inline T
ComputeArithmeticMean(const std::array<T, N>& values)
{
    T sum = static_cast<T>(0);
    for (int idx{0}; idx < N; ++idx)
    {
        sum += values[idx];
    }
    return sum / static_cast<T>(N);
}

template<ArithmeticType T, int32_t N>
constexpr inline T
ComputeMidRange(const std::array<T, N>& values)
{
    static_assert(N >= 2);

    T min = std::numeric_limits<T>::max();
    /* Find the minimum */
    for (int idx{0}; idx < N; ++idx)
    {
        if (min > values[idx])
        {
            min = values[idx];
        }
    }

    T max = std::numeric_limits<T>::lowest();
    /* Find the maximum */
    for (int idx{0}; idx < N; ++idx)
    {
        if (max < values[idx])
        {
            max = values[idx];
        }
    }

    /* Compute the mid-range */
    return ((max / 2) + (min / 2));
}

template<ArithmeticType T>
inline T
ComputeArithmeticMean(const std::span<T> values)
{
    cmc_assert(values.size() >= 2);

    T sum = static_cast<T>(0);
    for (size_t idx{0}; idx < values.size(); ++idx)
    {
        sum += values[idx];
    }
    return sum / static_cast<T>(values.size());
}

template<ArithmeticType T>
inline T
ComputeMidRange(const std::span<T> values)
{
    cmc_assert(values.size() >= 2);

    T min{std::numeric_limits<T>::max()};
    /* Find the minimum */
    for (size_t idx{0}; idx < values.size(); ++idx)
    {
        if (min > values[idx])
        {
            min = values[idx];
        }
    }

    T max{std::numeric_limits<T>::lowest()};
    /* Find the maximum */
    for (size_t idx{0}; idx < values.size(); ++idx)
    {
        if (max < values[idx])
        {
            max = values[idx];
        }
    }

    /* Compute the mid-range */
    return ((max / 2) + (min / 2));
}

template <int DIM>
requires Dimension<DIM>
inline int
ComputeNumCompressionLevels(const int max_dimension_length)
{
    cmc_assert(max_dimension_length > 0);
    if (max_dimension_length == 1)
    {
        return 0;
    }

    int exp = 1;
    for (int i{1}; i < max_dimension_length; ++i)
    {
        exp *= kDimReductionFactor;
        if (exp >= max_dimension_length)
        {
            return i;
        }
    }

    return 0;
}



}

#endif /* !CMC_PATCH_LOSSLESS_CMC_PATCH_MULTI_RES_EXTRACTION_UTIL_HXX */
