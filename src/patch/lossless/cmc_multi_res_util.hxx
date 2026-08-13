#ifndef CMC_PATCH_LOSSLESS_MULTI_RES_UTIL_HXX
#define CMC_PATCH_LOSSLESS_MULTI_RES_UTIL_HXX

#include "utilities/cmc_bits.hxx"
#include "utilities/cmc_bits_vector.hxx"
#include "utilities/cmc_huffman_coder.hxx"
#include "utilities/cmc_bits_stream_decoder.hxx"

#include <array>
#include <limits>
#include <execution>
#include <cmath>
#include <bitset>

namespace cmc::patch::lossless::multi_res
{

/* Switch to minimize the compression error for the fast multi res compression */
constexpr bool kTryToCorrectMeanFastCompression = true;

template<int32_t DIM>
concept Dimension = (DIM >= 1 && DIM <= 4);

/* The data reduction per dimension */    
constexpr int32_t kDIMReductionFactor = 2;

/* Set the number of adjacent data points that will be coarsened given a certain dimensionality */
template<int32_t DIM>
constexpr int kPackSize;
template<>
constexpr int kPackSize<1> = kDIMReductionFactor;
template<>
constexpr int kPackSize<2> = kDIMReductionFactor * kDIMReductionFactor;
template<>
constexpr int kPackSize<3> = kDIMReductionFactor * kDIMReductionFactor * kDIMReductionFactor;
template<>
constexpr int kPackSize<4> = kDIMReductionFactor * kDIMReductionFactor * kDIMReductionFactor * kDIMReductionFactor;

constexpr int k1DNumControlValues = 3;
constexpr int k1DNumPredicitionValues = 2;
constexpr int k2DNumControlValues = 5;
constexpr int k2DNumPredicitionValues = 4;
constexpr int k3DNumControlValues = 7;
constexpr int k3DNumPredicitionValues = 8;
constexpr int k4DNumControlValues = 9;
constexpr int k4DNumPredicitionValues = 16;

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

template<typename T>
concept ResidualType = (std::is_unsigned_v<T> && std::is_integral_v<T>);

using OneByteResidualType = uint8_t;

using TwoByteResidualType = uint16_t;

using FourByteResidualType = uint32_t;

using EightByteResidualType = uint64_t;

template<typename T>
concept FloatType = (std::is_floating_point_v<T>);

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

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline int32_t
GetMaximumDimensionLength(const std::array<int32_t, DIM>& dimension_lengths)
{
    /* Get the largest dimension length */
    const auto max_dim_iter = std::max_element(dimension_lengths.begin(), dimension_lengths.end());

    if (max_dim_iter == dimension_lengths.end())
    {
        cmc_err_msg("A maximum dimension length could not be retrieved from the dimension-length vector!");
    }

    return *max_dim_iter;
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
constexpr int32_t 
GetNumCompressionIterations(const int32_t max_dim_length)
{
    return static_cast<int32_t>(std::ceil(std::log2(max_dim_length)));
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
constexpr int32_t 
GetNumCompressionIterations(const std::array<int32_t, DIM>& dimension_lengths)
{
    const int32_t max_dim_length = GetMaximumDimensionLength<T, DIM>(dimension_lengths);
    return GetNumCompressionIterations<T, DIM>(max_dim_length);
}

/** Entropy Symbol Definitions **/

using SymbolType = uint8_t;

inline constexpr int32_t kMaxPresentElementLevelUnknown = -1;

constexpr inline SymbolType kResidualSignumIndication = 0x80;

template<ArithmeticType T>
constexpr inline SymbolType kProcessEndSymbol = kResidualSignumIndication + sizeof(T) * cmc::bits::kCharBit;

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

template<ArithmeticType T>
constexpr inline void
AddProcessEndSymbol(std::array<uint64_t, GetNumEntropySymbols<T>()>& entropy_symbols_frequency, const uint64_t num_local_proc_end_symbols)
{
    /* We store the process-end-symbol in the last array entry */
    entropy_symbols_frequency[GetNumEntropySymbols<T>() - 1] = num_local_proc_end_symbols;
}

/** END of Entropy Symbol Definitions **/

template<typename T, int32_t DIM>
struct PatchEncoding;

template<typename T, int32_t DIM>
requires Dimension<DIM> && OneByteArithmeticType<T>
struct PatchEncoding<T, DIM>
{
    uint8_t num_elements;
    std::array<SymbolType, kPackSize<DIM>> entropy_symbols{};
    std::array<OneByteResidualType, kPackSize<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && TwoByteArithmeticType<T>
struct PatchEncoding<T, DIM>
{
    uint16_t num_elements;
    std::array<SymbolType, kPackSize<DIM>> entropy_symbols{};
    std::array<TwoByteResidualType, kPackSize<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && FourByteArithmeticType<T>
struct PatchEncoding<T, DIM>
{
    uint32_t num_elements;
    std::array<SymbolType, kPackSize<DIM>> entropy_symbols{};
    std::array<FourByteResidualType, kPackSize<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && EightByteArithmeticType<T>
struct PatchEncoding<T, DIM>
{
    uint64_t num_elements;
    std::array<SymbolType, kPackSize<DIM>> entropy_symbols{};
    std::array<EightByteResidualType, kPackSize<DIM>> residuals{};
};

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
ComputeArithmeticMean(const std::span<const T> values)
{
    cmc_assert(values.size() >= 2);

    T sum = static_cast<T>(0);
    for (unsigned idx{0}; idx < values.size(); ++idx)
    {
        sum += values[idx];
    }
    return sum / static_cast<T>(values.size());
}

template<ArithmeticType T>
inline T
ComputeMidRange(const std::span<const T> values)
{
    cmc_assert(values.size() >= 2);

    T min{std::numeric_limits<T>::max()};
    /* Find the minimum */
    for (unsigned idx{0}; idx < values.size(); ++idx)
    {
        if (min > values[idx])
        {
            min = values[idx];
        }
    }

    T max{std::numeric_limits<T>::lowest()};
    /* Find the maximum */
    for (unsigned idx{0}; idx < values.size(); ++idx)
    {
        if (max < values[idx])
        {
            max = values[idx];
        }
    }

    /* Compute the mid-range */
    return ((max / 2) + (min / 2));
}

template<ArithmeticType T, int32_t N>
inline std::array<T, N+2>
CreatePredictors(const std::array<T, N>& init_data, const int num_elems)
{
    /* Utilize each individual value as a predictor */
    std::array<T, N + 2> predictors;
    std::copy_n(init_data.begin(), num_elems, predictors.begin());

    /* Create a view on the data */
    const std::span<const T> data_span(init_data.data(), num_elems);

    /* Append the arithmetic mean as a predictor */
    predictors[num_elems] = ComputeArithmeticMean<T>(data_span);

    /* Append the mid-range as a predictor */
    predictors[num_elems + 1] = ComputeMidRange<T>(data_span);

    return predictors;
}

inline std::pair<long double, long double>
Fast2Sum(const long double x, const long double y)
{
    //cmc_assert(std::abs(x) > std::abs(y));
    const long double sum = x + y;
    const long double ya = sum - x;
    const long double corr = y - ya;
    return std::make_pair(sum, corr);
}

inline std::pair<long double, long double>
TwoSum(const long double x, const long double y)
{
    const long double sum = x + y;
    const long double xa = sum - y;
    const long double ya = sum - xa;
    const long double corrx = x - xa;
    const long double corry = y - ya;
    const long double corr = corrx + corry;
    return std::make_pair(sum, corr);
}

template<ArithmeticType T, int32_t N>
inline long double
KahanTwoSum(const std::array<T, N>& values, const int num_vals)
{
    static_assert(N >= 1);
    cmc_assert(num_vals >= 1);
    long double sum{0.0};
    long double correction{0.0};
    
    for (int val_idx{0}; val_idx < num_vals; ++val_idx)
    {
        long double next_val_corrected = static_cast<long double>(values[val_idx]) + correction;
        auto [new_sum, new_corr] = TwoSum(sum, next_val_corrected);
        sum = new_sum;
        correction = new_corr;
    }

    return static_cast<T>(sum);
}

template<FloatType T, int32_t N>
inline T
NeumaierKahanBabuskaSum(const std::array<T, N>& values, const int num_vals)
{
    static_assert(N >= 1);
    cmc_assert(num_vals >= 1);
    long double sum{0.0};
    long double correction{0.0};
    
    for (int val_idx{0}; val_idx < num_vals; ++val_idx)
    {
        volatile long double sn = sum + static_cast<long double>(values[val_idx]);
        
        if (std::abs(sum) >= std::abs(values[val_idx]))
        {
            volatile long double corr = (sum - sn) + static_cast<long double>(values[val_idx]);
            correction += corr;
        } else
        {
            volatile long double corr = (static_cast<long double>(values[val_idx]) - sn) + sum;
            correction += corr;
        }

        sum = sn;
    }

    const long double sum_corrected = sum + correction;

    return static_cast<T>(sum_corrected);
}


template<FloatType T>
inline T
NeumaierKahanBabuskaSum(const std::vector<T>& values)
{
    cmc_assert(values.size() >= 1);
    long double sum{0.0};
    long double correction{0.0};
    
    for (const T& val : values)
    {
        volatile long double sn = sum + static_cast<long double>(val);
        
        if (std::abs(sum) >= std::abs(val))
        {
            volatile long double corr = (sum - sn) + static_cast<long double>(val);
            correction += corr;
        } else
        {
            volatile long double corr = (static_cast<long double>(val) - sn) + sum;
            correction += corr;
        }

        sum = sn;
    }

    const long double sum_corrected = sum + correction;

    return static_cast<T>(sum_corrected);
}

template<ArithmeticType T, int32_t N>
inline T
CreateArithmeticMeanPredictor(const std::array<T, N>& init_data, const int num_elems)
{
    const T sum = NeumaierKahanBabuskaSum<T, N>(init_data, num_elems);
    const T mean = sum / static_cast<T>(num_elems);
    return mean;
}

constexpr uint32_t kCheckFloatSign = 0x80000000U;
constexpr uint32_t kAllFloatBitsSet = 0xFFFFFFFFU;

constexpr uint64_t kCheckDoubleSign = 0x8000000000000000ULL;
constexpr uint64_t kAllDoubleBitsSet = 0xFFFFFFFFFFFFFFFFULL;

/* Assign a monotonically ordering to all floats */
inline uint32_t
OrderFloat(const float value)
{
    const uint32_t u_value = std::bit_cast<uint32_t>(value);

    if (u_value & kCheckFloatSign)
    {
        return u_value ^ kAllFloatBitsSet;
    } else
    {
        return u_value ^ kCheckFloatSign;
    }
}

/* Reverse the ordering to the actual float */
inline float
ReFloatOrder(const uint32_t order)
{
    if (order & kCheckFloatSign)
    {
        return std::bit_cast<float>(order ^ kCheckFloatSign);
    } else
    {
        return std::bit_cast<float>(order ^ kAllFloatBitsSet);
    }
}

inline uint64_t
OrderFloat(const double value)
{
    const uint64_t u_value = std::bit_cast<uint64_t>(value);

    if (u_value & kCheckDoubleSign)
    {
        return u_value ^ kAllDoubleBitsSet;
    } else
    {
        return u_value ^ kCheckDoubleSign;
    }
}

/* Reverse the ordering to the actual float */
inline double
ReFloatOrder(const uint64_t order)
{
    if (order & kCheckDoubleSign)
    {
        return std::bit_cast<double>(order ^ kCheckDoubleSign);
    } else
    {
        return std::bit_cast<double>(order ^ kAllDoubleBitsSet);
    }
}

template <FloatType T>
inline T
ComputeFMAForImplicitValue(const T N, const T mean, const T partial_sum)
{
    return std::fma(N, mean, -partial_sum);
}

inline float
FindMatchingMeanValue(const float N, const float partial_sum, const float value_to_match)
{
    /* Perform a binary search such that the values is matched */
    uint32_t lbound = OrderFloat(std::numeric_limits<float>::lowest());
    uint32_t ubound = OrderFloat(std::numeric_limits<float>::max());

    /* Define the match we are looking for */
    const uint32_t match = OrderFloat(value_to_match);

    while (lbound < ubound)
    {
        /* Compute the mid value from the interval */
        const uint32_t mid = lbound + static_cast<uint64_t>((ubound - lbound)) / 2;

        /* Get the order from the mid value */
        const float fmid = ReFloatOrder(mid);

        /* Compute the corresponding target value */
        const float fmid_target = ComputeFMAForImplicitValue<float>(N, fmid, partial_sum);

        /* Get the order from the fmid value */
        const uint32_t order_fmid_target = OrderFloat(fmid_target);

        /* Check in which interval the search continues */
        if (order_fmid_target < match)
        {
            lbound = mid + 1;
        } else
        {
            ubound = mid;
        }
    }

    /* Get the output value fromt he binary search */
    const float eval_mean = ReFloatOrder(lbound);

    /* Compute the matching value */
    #ifdef CMC_ENABLE_DEBUG
    [[maybe_unused]] const float eval_match = ComputeFMAForImplicitValue<float>(N, eval_mean, partial_sum);
    cmc_debug_msg("Evaluated match: ", eval_match, ", ", std::bitset<32>(std::bit_cast<uint32_t>(eval_match)));
    cmc_debug_msg(" Value to match: ", value_to_match, ", ", std::bitset<32>(std::bit_cast<uint32_t>(value_to_match)));
    #endif

    return eval_mean;
}

inline double
FindMatchingMeanValue(const double N, const double partial_sum, const double value_to_match)
{
    /* Perform a binary search such that the values is matched */
    uint64_t lbound = OrderFloat(std::numeric_limits<double>::lowest());
    uint64_t ubound = OrderFloat(std::numeric_limits<double>::max());

    /* Define the match we are looking for */
    const uint64_t match = OrderFloat(static_cast<double>(value_to_match));

    while (lbound < ubound)
    {
        /* Compute the mid value from the interval */
        const uint64_t mid = lbound + static_cast<uint64_t>((ubound - lbound)) / 2;

        /* Get the order from the mid value */
        const double fmid = ReFloatOrder(mid);

        /* Compute the corresponding target value */
        const double fmid_target = ComputeFMAForImplicitValue<double>(N, fmid, partial_sum);

        /* Get the order from the fmid value */
        const uint64_t order_fmid_target = OrderFloat(fmid_target);

        /* Check in which interval the search continues */
        if (order_fmid_target < match)
        {
            lbound = mid + 1;
        } else
        {
            ubound = mid;
        }
    }

    /* Get the output value fromt he binary search */
    const double eval_mean = ReFloatOrder(lbound);

    /* Compute the matching value */
    #ifdef CMC_ENABLE_DEBUG
    [[maybe_unused]] const double eval_match = ComputeFMAForImplicitValue<double>(N, eval_mean, partial_sum);
    cmc_debug_msg("Evaluated match: ", eval_match, ", ", std::bitset<64>(std::bit_cast<uint64_t>(eval_match)));
    cmc_debug_msg(" Value to match: ", value_to_match, ", ", std::bitset<64>(std::bit_cast<uint64_t>(value_to_match)));
    #endif

    return eval_mean;
}

template<FloatType T, int32_t N>
inline T
CreateMatchingMeanPredictor(const std::array<T, N>& init_data, const int num_elems)
{
    cmc_assert(num_elems > 1);
    const T num_vals = static_cast<T>(num_elems);
    const auto d_partial_sum = NeumaierKahanBabuskaSum<T, N>(init_data, num_elems - 1);
    const T partial_sum = static_cast<T>(d_partial_sum);
    const T value_to_match = init_data[num_elems - 1];
    /* Create the mean predictor */
    const T matching_mean = FindMatchingMeanValue(num_vals, partial_sum, value_to_match);
    return matching_mean;
}

template<FloatType T>
inline T
ComputeImplicitValue(const T mean_value, const std::vector<T>& values)
{
    cmc_assert(values.size() > 1ULL);
    const T num_vals = static_cast<T>(values.size() + 1);
    const auto d_partial_sum = NeumaierKahanBabuskaSum<T>(values);
    const T partial_sum = static_cast<T>(d_partial_sum);
    const T implicit_val = ComputeFMAForImplicitValue<T>(num_vals, mean_value, partial_sum);
    return implicit_val;
}

template<typename T, int32_t DIM>
requires Dimension<DIM>
struct Residuals;

template<typename T, int32_t DIM>
requires Dimension<DIM> && OneByteArithmeticType<T>
struct Residuals<T, DIM>
{
    std::array<OneByteResidualType, kPackSize<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && TwoByteArithmeticType<T>
struct Residuals<T, DIM>
{
    std::array<TwoByteResidualType, kPackSize<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && FourByteArithmeticType<T>
struct Residuals<T, DIM>
{
    std::array<FourByteResidualType, kPackSize<DIM>> residuals{};
};

template<typename T, int32_t DIM>
requires Dimension<DIM> && EightByteArithmeticType<T>
struct Residuals<T, DIM>
{
    std::array<EightByteResidualType, kPackSize<DIM>> residuals{};
};

template <typename T>
inline T
GetValue(const std::vector<T>& data, const int time, const int lev, const int lat, const int lon, const int kLonLength, const int kLatLength, const int kLevLength, [[maybe_unused]] const int kTimeLength)
{
    cmc_assert(time * (kLevLength * kLatLength * kLonLength) + lev * (kLatLength * kLonLength) + lat * kLonLength + lon < data.size());
    return data[time * (kLevLength * kLatLength * kLonLength) + lev * (kLatLength * kLonLength) + lat * kLonLength + lon];
}

template <typename T>
inline void
SetValue(std::vector<T>& data, const T& value, const int time, const int lev, const int lat, const int lon, const int kLonLength, const int kLatLength,  const int kLevLength, [[maybe_unused]] const int kTimeLength)
{
    cmc_assert(time * (kLevLength * kLatLength * kLonLength) + lev * (kLatLength * kLonLength) + lat * kLonLength + lon < data.size());
    data[time * (kLevLength * kLatLength * kLonLength) + lev * (kLatLength * kLonLength) + lat * kLonLength + lon] = value;
}

template <typename T>
inline T
GetValue(const std::vector<T>& data, const int lev, const int lat, const int lon, const int kLonLength, const int kLatLength, [[maybe_unused]] const int kLevLength)
{
    cmc_assert(lev * (kLatLength * kLonLength) + lat * kLonLength + lon < data.size());
    return data[lev * (kLatLength * kLonLength) + lat * kLonLength + lon];
}

template <typename T>
inline void
SetValue(std::vector<T>& data, const T& value, const int lev, const int lat, const int lon, const int kLonLength, const int kLatLength, [[maybe_unused]] const int kLevLength)
{
    cmc_assert(lev * (kLatLength * kLonLength) + lat * kLonLength + lon < data.size());
    data[lev * (kLatLength * kLonLength) + lat * kLonLength + lon] = value;
}

template <typename T>
inline T
GetValue(const std::vector<T>& data, const int lat, const int lon, const int kLonLength, [[maybe_unused]] const int kLatLength)
{
    cmc_assert(lat * kLonLength + lon < data.size());
    return data[lat * kLonLength + lon];
}

template <typename T>
inline void
SetValue(std::vector<T>& data, const T& value, const int lat, const int lon, const int kLonLength, [[maybe_unused]] const int kLatLength)
{
    cmc_assert(lat * kLonLength + lon < data.size());
    data[lat * kLonLength + lon] = value;
}

template <typename T>
inline T
GetValue(const std::vector<T>& data, const int lon, [[maybe_unused]] const int kLonLength)
{
    cmc_assert(lon < kLonLength);
    cmc_assert(lon < data.size());
    return data[lon];
}

template <typename T>
inline void
SetValue(std::vector<T>& data, const T& value, const int lon, [[maybe_unused]] const int kLonLength)
{
    cmc_assert(lon < kLonLength);
    cmc_assert(lon < data.size());
    data[lon] = value;
}


}

#endif /* !CMC_PATCH_LOSSLESS_MULTI_RES_UTIL_HXX */
