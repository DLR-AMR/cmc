#ifndef CMC_PATCH_LOSSY_MULTI_RES_UTIL_HXX
#define CMC_PATCH_LOSSY_MULTI_RES_UTIL_HXX

#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "amr/lossy/cmc_par_multi_res_interpolation_util.hxx"

#include <cmath>
#include <string>
#include <span>
#include <filesystem>
#include <array>
#include <vector>
#include <algorithm>
#include <execution>

namespace cmc::patch::lossy::multi_res
{

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

using OneByteResidualType = uint8_t;

using TwoByteResidualType = uint16_t;

using FourByteResidualType = uint32_t;

using EightByteResidualType = uint64_t;

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
using SymbolType = uint16_t;

constexpr inline SymbolType kPredictionWithinBound = 0;

constexpr inline int kResidualMaxDeviationFactor = 255;
constexpr inline int kCheckResidualMaxDeviationFactor = kResidualMaxDeviationFactor;
static_assert(kResidualMaxDeviationFactor % 2 != 0 && kResidualMaxDeviationFactor > 0 && kResidualMaxDeviationFactor < std::numeric_limits<SymbolType>::max() - 2);

constexpr inline int kEntropySymbolBinSignSwitch = (kResidualMaxDeviationFactor - 1) / 2;

constexpr inline SymbolType kFlagUnpredictable = kResidualMaxDeviationFactor;

constexpr inline SymbolType kProcessEndSymbol = kResidualMaxDeviationFactor + 1;

inline SymbolType
CreateEntropySymbolFromQuantizationBin(const SymbolType bin, const bool is_prediction_greater)
{
    cmc_assert(bin > static_cast<SymbolType>(0) && bin <= static_cast<SymbolType>(kEntropySymbolBinSignSwitch));
    return (is_prediction_greater ? kEntropySymbolBinSignSwitch + bin : bin);
}

inline std::pair<bool, int>
GetQuantizationBinFromEntropySymbol(const SymbolType symbol)
{
    cmc_assert(symbol != kProcessEndSymbol && symbol != kFlagUnpredictable);
    const bool is_prediction_greater = (symbol > kEntropySymbolBinSignSwitch);
    const int bin = (is_prediction_greater ? symbol - kEntropySymbolBinSignSwitch : symbol);
    return std::make_pair(is_prediction_greater, bin);
}

template<ArithmeticType T>
constexpr inline int
MapEntropySymbolToArrayIndex(const SymbolType entropy_symbol)
{
    return static_cast<int>(entropy_symbol);
}

template<ArithmeticType T>
constexpr inline SymbolType
MapArrayIndexToEntropySymbol(int array_idx)
{
    return static_cast<SymbolType>(array_idx);
}

constexpr inline int
GetNumEntropySymbols()
{
    return kResidualMaxDeviationFactor + 2;
}

constexpr inline void
AddProcessEndSymbol(std::array<uint64_t, GetNumEntropySymbols()>& entropy_symbols_frequency, const uint64_t num_local_proc_end_symbols)
{
    /* We store the process-end-symbol in the last array entry */
    entropy_symbols_frequency[GetNumEntropySymbols() - 1] = num_local_proc_end_symbols;
}
/** END OF Entropy Symbol Definitions **/

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
struct PatchEncoding
{
    std::array<SymbolType, kPackSize<DIM>> quantization_bins{};
    std::array<T, kPackSize<DIM>> unpredictable_values{};
    int32_t num_elements;
};

template<ArithmeticType T>
inline T
GetAbsResidual(const T& value1, const T& value2)
{
    if constexpr (std::is_signed_v<T>)
    {
        return std::fabs(value2 - value1);
    } else
    {
        return (value1 > value2 ? value1 - value2 : value2 - value1);
    }
}

inline SymbolType
ComputeBin(const float abs_residual, const float abs_permitted_error)
{
    return static_cast<SymbolType>(std::floor(abs_residual / (2.0 * abs_permitted_error) + 0.5));
}

template<ArithmeticType T>
inline T
DequantizeValue(const T prediction, const float permitted_abs_error, const bool is_prediction_greater, const SymbolType bin)
{
    return static_cast<T>(prediction + (is_prediction_greater ? -2.0 : +2.0) * permitted_abs_error * bin);
}

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

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline std::array<T, k4DNumControlValues>
GetFaceControlValues(const std::vector<T>& data,
                     const int32_t time_idx, const int32_t lev_idx, const int32_t lat_idx, const int32_t lon_idx, const std::array<int32_t, DIM>& dim_lengths)
{
    static_assert(DIM == 4);
    static_assert(k4DNumControlValues == 9);

    constexpr int32_t kTimeID = 0;
    constexpr int32_t kLevID = 1;
    constexpr int32_t kLatID = 2;
    constexpr int32_t kLonID = 3;

    const T patch_control_value = GetValue<T>(data, time_idx, lev_idx, lat_idx, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID], dim_lengths[kTimeID]);

    std::array<T, k4DNumControlValues> control_values;
    control_values.fill(patch_control_value);

    if (lon_idx > 0) [[likely]]
    {
        control_values[1] = GetValue<T>(data, time_idx, lev_idx, lat_idx, lon_idx - 1, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID], dim_lengths[kTimeID]);
    }
    if (lon_idx < dim_lengths[kLonID] - 1) [[likely]]
    {
        control_values[2] = GetValue<T>(data, time_idx, lev_idx, lat_idx, lon_idx + 1, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID], dim_lengths[kTimeID]);
    }

    if (lat_idx > 0) [[likely]]
    {
        control_values[3] = GetValue<T>(data, time_idx, lev_idx, lat_idx - 1, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID], dim_lengths[kTimeID]);
    }
    if (lat_idx < dim_lengths[kLatID] - 1) [[likely]]
    {
        control_values[4] = GetValue<T>(data, time_idx, lev_idx, lat_idx + 1, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID], dim_lengths[kTimeID]);
    }

    if (lev_idx > 0) [[likely]]
    {
        control_values[5] = GetValue<T>(data, time_idx, lev_idx - 1, lat_idx, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID], dim_lengths[kTimeID]);
    }
    if (lev_idx < dim_lengths[kLevID] - 1) [[likely]]
    {
        control_values[6] = GetValue<T>(data, time_idx, lev_idx + 1, lat_idx, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID], dim_lengths[kTimeID]);
    }

    if (time_idx > 0) [[likely]]
    {
        control_values[5] = GetValue<T>(data, time_idx - 1, lev_idx, lat_idx, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID], dim_lengths[kTimeID]);
    }
    if (time_idx < dim_lengths[kTimeID] - 1) [[likely]]
    {
        control_values[6] = GetValue<T>(data, time_idx + 1, lev_idx, lat_idx, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID], dim_lengths[kTimeID]);
    }

    return control_values;
}


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline std::array<T, k3DNumControlValues>
GetFaceControlValues(const std::vector<T>& data,
                     const int32_t lev_idx, const int32_t lat_idx, const int32_t lon_idx, const std::array<int32_t, DIM>& dim_lengths)
{
    static_assert(DIM == 3);
    static_assert(k3DNumControlValues == 7);

    constexpr int32_t kLevID = 0;
    constexpr int32_t kLatID = 1;
    constexpr int32_t kLonID = 2;

    const T patch_control_value = GetValue<T>(data, lev_idx, lat_idx, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID]);
    
    std::array<T, k3DNumControlValues> control_values;
    control_values.fill(patch_control_value);

    if (lon_idx > 0) [[likely]]
    {
        control_values[1] = GetValue<T>(data, lev_idx, lat_idx, lon_idx - 1, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID]);
    }
    if (lon_idx < dim_lengths[kLonID] - 1) [[likely]]
    {
        control_values[2] = GetValue<T>(data, lev_idx, lat_idx, lon_idx + 1, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID]);
    }

    if (lat_idx > 0) [[likely]]
    {
        control_values[3] = GetValue<T>(data, lev_idx, lat_idx - 1, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID]);
    }
    if (lat_idx < dim_lengths[kLatID] - 1) [[likely]]
    {
        control_values[4] = GetValue<T>(data, lev_idx, lat_idx + 1, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID]);
    }

    if (lev_idx > 0) [[likely]]
    {
        control_values[5] = GetValue<T>(data, lev_idx - 1, lat_idx, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID]);
    }
    if (lev_idx < dim_lengths[kLevID] - 1) [[likely]]
    {
        control_values[6] = GetValue<T>(data, lev_idx + 1, lat_idx, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID]);
    }

    return control_values;
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline std::array<T, k2DNumControlValues>
GetFaceControlValues(const std::vector<T>& data,
                     const int32_t lat_idx, const int32_t lon_idx, const std::array<int32_t, DIM>& dim_lengths)
{
    static_assert(DIM == 2);
    static_assert(k2DNumControlValues == 5);

    constexpr int32_t kLatID = 0;
    constexpr int32_t kLonID = 1;

    const T patch_control_value = GetValue<T>(data, lat_idx, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID]);
    
    std::array<T, k2DNumControlValues> control_values;
    control_values.fill(patch_control_value);

    if (lon_idx > 0) [[likely]]
    {
        control_values[1] = GetValue<T>(data, lat_idx, lon_idx - 1, dim_lengths[kLonID], dim_lengths[kLatID]);
    }
    if (lon_idx < dim_lengths[kLonID] - 1) [[likely]]
    {
        control_values[2] = GetValue<T>(data, lat_idx, lon_idx + 1, dim_lengths[kLonID], dim_lengths[kLatID]);
    }

    if (lat_idx > 0) [[likely]]
    {
        control_values[3] = GetValue<T>(data, lat_idx - 1, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID]);
    }
    if (lat_idx < dim_lengths[kLatID] - 1) [[likely]]
    {
        control_values[4] = GetValue<T>(data, lat_idx + 1, lon_idx, dim_lengths[kLonID], dim_lengths[kLatID]);
    }

    return control_values;
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline std::array<T, k1DNumControlValues>
GetFaceControlValues(const std::vector<T>& data,
                      const int32_t lon_idx, const std::array<int32_t, DIM>& dim_lengths)
{
    static_assert(DIM == 1);
    static_assert(k1DNumControlValues == 3);

    constexpr int32_t kLonID = 0;

    const T patch_control_value = GetValue<T>(data, lon_idx, dim_lengths[kLonID]);
    
    std::array<T, k1DNumControlValues> control_values;
    control_values.fill(patch_control_value);

    if (lon_idx > 0) [[likely]]
    {
        control_values[1] = GetValue<T>(data, lon_idx - 1, dim_lengths[kLonID]);
    }
    if (lon_idx < dim_lengths[kLonID] - 1) [[likely]]
    {
        control_values[2] = GetValue<T>(data, lon_idx + 1, dim_lengths[kLonID]);
    }

    return control_values;
}

namespace idw
{
/* 2D-IDW interpolation */
namespace two_dimensional
{
constexpr int32_t kExpDist = 8;
constexpr float ExpDist(const double dist)
{
    double exp_dist = dist;
    for (int32_t i{1}; i < kExpDist; ++i)
    {
        exp_dist *= dist;
    }
    return static_cast<float>(exp_dist);
}
//TODO: The distances account for the neighboring element's midpoint and not the actual first child's midpoint
constexpr float pred1_d0 = 1.0f / ExpDist(std::sqrt(2.0));
constexpr float pred1_d1 = 1.0f / ExpDist(std::sqrt(26.0));
constexpr float pred1_d2 = 1.0f / ExpDist(std::sqrt(10.0));
constexpr float pred1_d3 = 1.0f / ExpDist(std::sqrt(10.0));
constexpr float pred1_d4 = 1.0f / ExpDist(std::sqrt(26.0));
constexpr float sum_pred1 = pred1_d0 + pred1_d1 + pred1_d2 + pred1_d3 + pred1_d4;

constexpr float pred2_d0 = 1.0f / ExpDist(std::sqrt(2.0));
constexpr float pred2_d1 = 1.0f / ExpDist(std::sqrt(10.0));
constexpr float pred2_d2 = 1.0f / ExpDist(std::sqrt(26.0));
constexpr float pred2_d3 = 1.0f / ExpDist(std::sqrt(26.0));
constexpr float pred2_d4 = 1.0f / ExpDist(std::sqrt(10.0));
constexpr float sum_pred2 = pred2_d0 + pred2_d1 + pred2_d2 + pred2_d3 + pred2_d4;

constexpr float pred3_d0 = 1.0f / ExpDist(std::sqrt(8.0));
constexpr float pred3_d1 = 1.0f / ExpDist(std::sqrt(26.0));
constexpr float pred3_d2 = 1.0f / ExpDist(std::sqrt(10.0));
constexpr float pred3_d3 = 1.0f / ExpDist(std::sqrt(26.0));
constexpr float pred3_d4 = 1.0f / ExpDist(std::sqrt(10.0));
constexpr float sum_pred3 = pred3_d0 + pred3_d1 + pred3_d2 + pred3_d3 + pred3_d4;
}
/* Define the prediction matrix */
constexpr std::array<std::array<float, k2DNumControlValues>, k2DNumPredicitionValues> PredictionMatrix2DIDW{{
    {1.0,0.0,0.0,0.0,0.0},
    {two_dimensional::pred1_d0 / two_dimensional::sum_pred1, two_dimensional::pred1_d1 / two_dimensional::sum_pred1, two_dimensional::pred1_d2 / two_dimensional::sum_pred1, two_dimensional::pred1_d3 / two_dimensional::sum_pred1, two_dimensional::pred1_d4 / two_dimensional::sum_pred1},
    {two_dimensional::pred2_d0 / two_dimensional::sum_pred2, two_dimensional::pred2_d1 / two_dimensional::sum_pred2, two_dimensional::pred2_d2 / two_dimensional::sum_pred2, two_dimensional::pred2_d3 / two_dimensional::sum_pred2, two_dimensional::pred2_d4 / two_dimensional::sum_pred2},
    {two_dimensional::pred3_d0 / two_dimensional::sum_pred3, two_dimensional::pred3_d1 / two_dimensional::sum_pred3, two_dimensional::pred3_d2 / two_dimensional::sum_pred3, two_dimensional::pred3_d3 / two_dimensional::sum_pred3, two_dimensional::pred3_d4 / two_dimensional::sum_pred3},
    }};

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
constexpr inline 
std::array<T, k2DNumPredicitionValues>
Perform2DIDWPrediction(const std::array<T, k2DNumControlValues>& control_values)
{
    static_assert(DIM == 2);

    const std::array<T, k2DNumPredicitionValues> prediction{
        control_values[0],
        static_cast<T>(control_values[0] * PredictionMatrix2DIDW[1][0] + control_values[1] * PredictionMatrix2DIDW[1][1] + control_values[2] * PredictionMatrix2DIDW[1][2] + control_values[3] * PredictionMatrix2DIDW[1][3] + control_values[4] * PredictionMatrix2DIDW[1][4]),
        static_cast<T>(control_values[0] * PredictionMatrix2DIDW[2][0] + control_values[1] * PredictionMatrix2DIDW[2][1] + control_values[2] * PredictionMatrix2DIDW[2][2] + control_values[3] * PredictionMatrix2DIDW[2][3] + control_values[4] * PredictionMatrix2DIDW[2][4]),
        static_cast<T>(control_values[0] * PredictionMatrix2DIDW[3][0] + control_values[1] * PredictionMatrix2DIDW[3][1] + control_values[2] * PredictionMatrix2DIDW[3][2] + control_values[3] * PredictionMatrix2DIDW[3][3] + control_values[4] * PredictionMatrix2DIDW[3][4]),
    };

    return prediction;
}



/* 3D-IDW interpolation */
namespace three_dimensional
{
constexpr int32_t kExpDist = 8;
constexpr float ExpDist(const double dist)
{
    double exp_dist = dist;
    for (int32_t i{1}; i < kExpDist; ++i)
    {
        exp_dist *= dist;
    }
    return static_cast<float>(exp_dist);
}
constexpr float pred1_d0 = 1.0f / ExpDist(2.0);
constexpr float pred1_d1 = 1.0f / ExpDist(6.0);
constexpr float pred1_d2 = 1.0f / ExpDist(2.0);
constexpr float pred1_d3 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float pred1_d4 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float pred1_d5 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float pred1_d6 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float sum_pred1 = pred1_d0 + pred1_d1 + pred1_d2 + pred1_d3 + pred1_d4 + pred1_d5 + pred1_d6;

constexpr float pred2_d0 = 1.0f / ExpDist(2.0);
constexpr float pred2_d1 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float pred2_d2 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float pred2_d3 = 1.0f / ExpDist(6.0);
constexpr float pred2_d4 = 1.0f / ExpDist(2.0);
constexpr float pred2_d5 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float pred2_d6 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float sum_pred2 = pred2_d0 + pred2_d1 + pred2_d2 + pred2_d3 + pred2_d4 + pred2_d5 + pred2_d6;

constexpr float pred3_d0 = 1.0f / ExpDist(std::sqrt(8.0));
constexpr float pred3_d1 = 1.0f / ExpDist(std::sqrt(26.0));
constexpr float pred3_d2 = 1.0f / ExpDist(std::sqrt(8.0));
constexpr float pred3_d3 = 1.0f / ExpDist(std::sqrt(26.0));
constexpr float pred3_d4 = 1.0f / ExpDist(std::sqrt(8.0));
constexpr float pred3_d5 = 1.0f / ExpDist(std::sqrt(24.0));
constexpr float pred3_d6 = 1.0f / ExpDist(std::sqrt(24.0));
constexpr float sum_pred3 = pred3_d0 + pred3_d1 + pred3_d2 + pred3_d3 + pred3_d4 + pred3_d5 + pred3_d6;

constexpr float pred4_d0 = 1.0f / ExpDist(2.0);
constexpr float pred4_d1 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float pred4_d2 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float pred4_d3 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float pred4_d4 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float pred4_d5 = 1.0f / ExpDist(6.0);
constexpr float pred4_d6 = 1.0f / ExpDist(2.0);
constexpr float sum_pred4 = pred4_d0 + pred4_d1 + pred4_d2 + pred4_d3 + pred4_d4 + pred4_d5 + pred4_d6;

constexpr float pred5_d0 = 1.0f / ExpDist(std::sqrt(8.0));
constexpr float pred5_d1 = 1.0f / ExpDist(std::sqrt(40.0));
constexpr float pred5_d2 = 1.0f / ExpDist(std::sqrt(8.0));
constexpr float pred5_d3 = 1.0f / ExpDist(std::sqrt(24.0));
constexpr float pred5_d4 = 1.0f / ExpDist(std::sqrt(24.0));
constexpr float pred5_d5 = 1.0f / ExpDist(std::sqrt(40.0));
constexpr float pred5_d6 = 1.0f / ExpDist(std::sqrt(8.0));
constexpr float sum_pred5 = pred5_d0 + pred5_d1 + pred5_d2 + pred5_d3 + pred5_d4 + pred5_d5 + pred5_d6;

constexpr float pred6_d0 = 1.0f / ExpDist(std::sqrt(8.0));
constexpr float pred6_d1 = 1.0f / ExpDist(std::sqrt(24.0));
constexpr float pred6_d2 = 1.0f / ExpDist(std::sqrt(24.0));
constexpr float pred6_d3 = 1.0f / ExpDist(std::sqrt(40.0));
constexpr float pred6_d4 = 1.0f / ExpDist(std::sqrt(8.0));
constexpr float pred6_d5 = 1.0f / ExpDist(std::sqrt(20.0));
constexpr float pred6_d6 = 1.0f / ExpDist(std::sqrt(8.0));
constexpr float sum_pred6 = pred6_d0 + pred6_d1 + pred6_d2 + pred6_d3 + pred6_d4 + pred6_d5 + pred6_d6;

constexpr float pred7_d0 = 1.0f / ExpDist(std::sqrt(12.0));
constexpr float pred7_d1 = 1.0f / ExpDist(std::sqrt(44.0));
constexpr float pred7_d2 = 1.0f / ExpDist(std::sqrt(12.0));
constexpr float pred7_d3 = 1.0f / ExpDist(std::sqrt(44.0));
constexpr float pred7_d4 = 1.0f / ExpDist(std::sqrt(12.0));
constexpr float pred7_d5 = 1.0f / ExpDist(std::sqrt(24.0));
constexpr float pred7_d6 = 1.0f / ExpDist(std::sqrt(12.0));
constexpr float sum_pred7 = pred7_d0 + pred7_d1 + pred7_d2 + pred7_d3 + pred7_d4 + pred7_d5 + pred7_d6;
}

/* Define the prediction matrix */
constexpr std::array<std::array<float, k3DNumControlValues>, k3DNumPredicitionValues> PredictionMatrix3DIDW{{
    {1.0,0.0,0.0,0.0,0.0,0.0,0.0},
    {three_dimensional::pred1_d0 / three_dimensional::sum_pred1, three_dimensional::pred1_d1 / three_dimensional::sum_pred1, three_dimensional::pred1_d2 / three_dimensional::sum_pred1, three_dimensional::pred1_d3 / three_dimensional::sum_pred1, three_dimensional::pred1_d4 / three_dimensional::sum_pred1, three_dimensional::pred1_d5 / three_dimensional::sum_pred1, three_dimensional::pred1_d6 / three_dimensional::sum_pred1},
    {three_dimensional::pred2_d0 / three_dimensional::sum_pred2, three_dimensional::pred2_d1 / three_dimensional::sum_pred2, three_dimensional::pred2_d2 / three_dimensional::sum_pred2, three_dimensional::pred2_d3 / three_dimensional::sum_pred2, three_dimensional::pred2_d4 / three_dimensional::sum_pred2, three_dimensional::pred2_d5 / three_dimensional::sum_pred2, three_dimensional::pred2_d6 / three_dimensional::sum_pred2},
    {three_dimensional::pred3_d0 / three_dimensional::sum_pred3, three_dimensional::pred3_d1 / three_dimensional::sum_pred3, three_dimensional::pred3_d2 / three_dimensional::sum_pred3, three_dimensional::pred3_d3 / three_dimensional::sum_pred3, three_dimensional::pred3_d4 / three_dimensional::sum_pred3, three_dimensional::pred3_d5 / three_dimensional::sum_pred3, three_dimensional::pred3_d6 / three_dimensional::sum_pred3},
    {three_dimensional::pred4_d0 / three_dimensional::sum_pred4, three_dimensional::pred4_d1 / three_dimensional::sum_pred4, three_dimensional::pred4_d2 / three_dimensional::sum_pred4, three_dimensional::pred4_d3 / three_dimensional::sum_pred4, three_dimensional::pred4_d4 / three_dimensional::sum_pred4, three_dimensional::pred4_d5 / three_dimensional::sum_pred4, three_dimensional::pred4_d6 / three_dimensional::sum_pred4},
    {three_dimensional::pred5_d0 / three_dimensional::sum_pred5, three_dimensional::pred5_d1 / three_dimensional::sum_pred5, three_dimensional::pred5_d2 / three_dimensional::sum_pred5, three_dimensional::pred5_d3 / three_dimensional::sum_pred5, three_dimensional::pred5_d4 / three_dimensional::sum_pred5, three_dimensional::pred5_d5 / three_dimensional::sum_pred5, three_dimensional::pred5_d6 / three_dimensional::sum_pred5},
    {three_dimensional::pred6_d0 / three_dimensional::sum_pred6, three_dimensional::pred6_d1 / three_dimensional::sum_pred6, three_dimensional::pred6_d2 / three_dimensional::sum_pred6, three_dimensional::pred6_d3 / three_dimensional::sum_pred6, three_dimensional::pred6_d4 / three_dimensional::sum_pred6, three_dimensional::pred6_d5 / three_dimensional::sum_pred6, three_dimensional::pred6_d6 / three_dimensional::sum_pred6},
    {three_dimensional::pred7_d0 / three_dimensional::sum_pred7, three_dimensional::pred7_d1 / three_dimensional::sum_pred7, three_dimensional::pred7_d2 / three_dimensional::sum_pred7, three_dimensional::pred7_d3 / three_dimensional::sum_pred7, three_dimensional::pred7_d4 / three_dimensional::sum_pred7, three_dimensional::pred7_d5 / three_dimensional::sum_pred7, three_dimensional::pred7_d6 / three_dimensional::sum_pred7},
    }};

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
constexpr inline 
std::array<T, k3DNumPredicitionValues>
Perform3DIDWPrediction(const std::array<T, k3DNumControlValues>& control_values)
{
    static_assert(DIM == 3);

    const std::array<T, k3DNumPredicitionValues> prediction{
        control_values[0],
        static_cast<T>(control_values[0] * PredictionMatrix3DIDW[1][0] + control_values[1] * PredictionMatrix3DIDW[1][1] + control_values[2] * PredictionMatrix3DIDW[1][2] + control_values[3] * PredictionMatrix3DIDW[1][3] + control_values[4] * PredictionMatrix3DIDW[1][4] + control_values[5] * PredictionMatrix3DIDW[1][5] + control_values[6] * PredictionMatrix3DIDW[1][6]),
        static_cast<T>(control_values[0] * PredictionMatrix3DIDW[2][0] + control_values[1] * PredictionMatrix3DIDW[2][1] + control_values[2] * PredictionMatrix3DIDW[2][2] + control_values[3] * PredictionMatrix3DIDW[2][3] + control_values[4] * PredictionMatrix3DIDW[2][4] + control_values[5] * PredictionMatrix3DIDW[2][5] + control_values[6] * PredictionMatrix3DIDW[2][6]),
        static_cast<T>(control_values[0] * PredictionMatrix3DIDW[3][0] + control_values[1] * PredictionMatrix3DIDW[3][1] + control_values[2] * PredictionMatrix3DIDW[3][2] + control_values[3] * PredictionMatrix3DIDW[3][3] + control_values[4] * PredictionMatrix3DIDW[3][4] + control_values[5] * PredictionMatrix3DIDW[3][5] + control_values[6] * PredictionMatrix3DIDW[3][6]),
        static_cast<T>(control_values[0] * PredictionMatrix3DIDW[4][0] + control_values[1] * PredictionMatrix3DIDW[4][1] + control_values[2] * PredictionMatrix3DIDW[4][2] + control_values[3] * PredictionMatrix3DIDW[4][3] + control_values[4] * PredictionMatrix3DIDW[4][4] + control_values[5] * PredictionMatrix3DIDW[4][5] + control_values[6] * PredictionMatrix3DIDW[4][6]),
        static_cast<T>(control_values[0] * PredictionMatrix3DIDW[5][0] + control_values[1] * PredictionMatrix3DIDW[5][1] + control_values[2] * PredictionMatrix3DIDW[5][2] + control_values[3] * PredictionMatrix3DIDW[5][3] + control_values[4] * PredictionMatrix3DIDW[5][4] + control_values[5] * PredictionMatrix3DIDW[5][5] + control_values[6] * PredictionMatrix3DIDW[5][6]),
        static_cast<T>(control_values[0] * PredictionMatrix3DIDW[6][0] + control_values[1] * PredictionMatrix3DIDW[6][1] + control_values[2] * PredictionMatrix3DIDW[6][2] + control_values[3] * PredictionMatrix3DIDW[6][3] + control_values[4] * PredictionMatrix3DIDW[6][4] + control_values[5] * PredictionMatrix3DIDW[6][5] + control_values[6] * PredictionMatrix3DIDW[6][6]),
        static_cast<T>(control_values[0] * PredictionMatrix3DIDW[7][0] + control_values[1] * PredictionMatrix3DIDW[7][1] + control_values[2] * PredictionMatrix3DIDW[7][2] + control_values[3] * PredictionMatrix3DIDW[7][3] + control_values[4] * PredictionMatrix3DIDW[7][4] + control_values[5] * PredictionMatrix3DIDW[7][5] + control_values[6] * PredictionMatrix3DIDW[7][6]),
    };

    return prediction;
}

/* 4D-IDW interpolation */
template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
constexpr inline 
std::array<T, k4DNumPredicitionValues>
Perform4DIDWPrediction(const std::array<T, k4DNumControlValues>& control_values)
{
    static_assert(DIM == 4);
    static_assert(false, "4D IDW Prediction not yet implemented!");
    return std::array<T, k4DNumPredicitionValues>{};
}

/* 1D-IDW interpolation */
template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
constexpr inline 
std::array<T, k1DNumPredicitionValues>
Perform1DIDWPrediction(const std::array<T, k1DNumControlValues>& control_values)
{
    static_assert(DIM == 1);
    static_assert(false, "1D IDW Prediction not yet implemented!");
    return std::array<T, k1DNumPredicitionValues>{};
}

}

namespace rbf
{

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
constexpr 
inline std::array<T, k2DNumPredicitionValues>
Perform2DRBFPrediction(const std::array<T, k2DNumControlValues>& control_values)
{
    static_assert(std::is_floating_point_v<T>, "The RBF Prediction is currently only applicable to floating point data!");
    
    /* Get the inverse sytem matrix */
    constexpr std::array<std::array<T, cmc::par::lossy::rbf::util::kNumQuadControlPoints>, cmc::par::lossy::rbf::util::kNumQuadControlPoints> M_inv = cmc::par::lossy::rbf::util::GetQuadMatrixInverse<T>();

    /* Compute the weights for the quad prediction */
    const std::array<T, cmc::par::lossy::rbf::util::kNumQuadControlPoints> weights{
        M_inv[0][0] * control_values[0] + M_inv[0][1] * control_values[1] + M_inv[0][2] * control_values[2] + M_inv[0][3] * control_values[3] + M_inv[0][4] * control_values[4],
        M_inv[1][0] * control_values[0] + M_inv[1][1] * control_values[1] + M_inv[1][2] * control_values[2] + M_inv[1][3] * control_values[3] + M_inv[1][4] * control_values[4],
        M_inv[2][0] * control_values[0] + M_inv[2][1] * control_values[1] + M_inv[2][2] * control_values[2] + M_inv[2][3] * control_values[3] + M_inv[2][4] * control_values[4],
        M_inv[3][0] * control_values[0] + M_inv[3][1] * control_values[1] + M_inv[3][2] * control_values[2] + M_inv[3][3] * control_values[3] + M_inv[3][4] * control_values[4],
        M_inv[4][0] * control_values[0] + M_inv[4][1] * control_values[1] + M_inv[4][2] * control_values[2] + M_inv[4][3] * control_values[3] + M_inv[4][4] * control_values[4]
    };

    /* Get the evaluation coordinates */
    constexpr std::array<std::array<T, cmc::par::lossy::rbf::util::kNumQuadControlPoints>, cmc::par::lossy::rbf::util::kNumQuadPredictionPoints> eval_coords = cmc::par::lossy::rbf::util::GetQuadPredictionCoordsEvaluation<T>();

    /* Perform the prediction */
    const std::array<T, k2DNumPredicitionValues> prediction{
        control_values[0],
        weights[0] * eval_coords[0][0] + weights[1] * eval_coords[0][1] + weights[2] * eval_coords[0][2] + weights[3] * eval_coords[0][3] + weights[4] * eval_coords[0][4],
        weights[0] * eval_coords[1][0] + weights[1] * eval_coords[1][1] + weights[2] * eval_coords[1][2] + weights[3] * eval_coords[1][3] + weights[4] * eval_coords[1][4],
        weights[0] * eval_coords[2][0] + weights[1] * eval_coords[2][1] + weights[2] * eval_coords[2][2] + weights[3] * eval_coords[2][3] + weights[4] * eval_coords[2][4]
    };

    return prediction;
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
constexpr 
inline std::array<T, k3DNumPredicitionValues>
Perform3DRBFPrediction(const std::array<T, k3DNumControlValues>& control_values)
{
    static_assert(std::is_floating_point_v<T>, "The RBF Prediction is currently only applicable to floating point data!");

    /* Get the inverse sytem matrix */
    constexpr std::array<std::array<T, cmc::par::lossy::rbf::util::kNumHexControlPoints>, cmc::par::lossy::rbf::util::kNumHexControlPoints> M_inv = cmc::par::lossy::rbf::util::GetHexMatrixInverse<T>();

    /* Compute the weights for the quad prediction */
    const std::array<T, cmc::par::lossy::rbf::util::kNumHexControlPoints> weights{
        M_inv[0][0] * control_values[0] + M_inv[0][1] * control_values[1] + M_inv[0][2] * control_values[2] + M_inv[0][3] * control_values[3] + M_inv[0][4] * control_values[4] + M_inv[0][5] * control_values[5] + M_inv[0][6] * control_values[6],
        M_inv[1][0] * control_values[0] + M_inv[1][1] * control_values[1] + M_inv[1][2] * control_values[2] + M_inv[1][3] * control_values[3] + M_inv[1][4] * control_values[4] + M_inv[1][5] * control_values[5] + M_inv[1][6] * control_values[6],
        M_inv[2][0] * control_values[0] + M_inv[2][1] * control_values[1] + M_inv[2][2] * control_values[2] + M_inv[2][3] * control_values[3] + M_inv[2][4] * control_values[4] + M_inv[2][5] * control_values[5] + M_inv[2][6] * control_values[6],
        M_inv[3][0] * control_values[0] + M_inv[3][1] * control_values[1] + M_inv[3][2] * control_values[2] + M_inv[3][3] * control_values[3] + M_inv[3][4] * control_values[4] + M_inv[3][5] * control_values[5] + M_inv[3][6] * control_values[6],
        M_inv[4][0] * control_values[0] + M_inv[4][1] * control_values[1] + M_inv[4][2] * control_values[2] + M_inv[4][3] * control_values[3] + M_inv[4][4] * control_values[4] + M_inv[4][5] * control_values[5] + M_inv[4][6] * control_values[6],
        M_inv[5][0] * control_values[0] + M_inv[5][1] * control_values[1] + M_inv[5][2] * control_values[2] + M_inv[5][3] * control_values[3] + M_inv[5][4] * control_values[4] + M_inv[5][5] * control_values[5] + M_inv[5][6] * control_values[6],
        M_inv[6][0] * control_values[0] + M_inv[6][1] * control_values[1] + M_inv[6][2] * control_values[2] + M_inv[6][3] * control_values[3] + M_inv[6][4] * control_values[4] + M_inv[6][5] * control_values[5] + M_inv[6][6] * control_values[6]
    };

    /* Get the evaluation coordinates */
    constexpr std::array<std::array<T, cmc::par::lossy::rbf::util::kNumHexControlPoints>, cmc::par::lossy::rbf::util::kNumHexPredictionPoints> eval_coords = cmc::par::lossy::rbf::util::GetHexPredictionCoordsEvaluation<T>();

    /* Perform the prediction */
    const std::array<T, k3DNumPredicitionValues> prediction{
        control_values[0],
        weights[0] * eval_coords[0][0] + weights[1] * eval_coords[0][1] + weights[2] * eval_coords[0][2] + weights[3] * eval_coords[0][3] + weights[4] * eval_coords[0][4] + weights[5] * eval_coords[0][5] + weights[6] * eval_coords[0][6],
        weights[0] * eval_coords[1][0] + weights[1] * eval_coords[1][1] + weights[2] * eval_coords[1][2] + weights[3] * eval_coords[1][3] + weights[4] * eval_coords[1][4] + weights[5] * eval_coords[1][5] + weights[6] * eval_coords[1][6],
        weights[0] * eval_coords[2][0] + weights[1] * eval_coords[2][1] + weights[2] * eval_coords[2][2] + weights[3] * eval_coords[2][3] + weights[4] * eval_coords[2][4] + weights[5] * eval_coords[2][5] + weights[6] * eval_coords[2][6],
        weights[0] * eval_coords[3][0] + weights[1] * eval_coords[3][1] + weights[2] * eval_coords[3][2] + weights[3] * eval_coords[3][3] + weights[4] * eval_coords[3][4] + weights[5] * eval_coords[3][5] + weights[6] * eval_coords[3][6],
        weights[0] * eval_coords[4][0] + weights[1] * eval_coords[4][1] + weights[2] * eval_coords[4][2] + weights[3] * eval_coords[4][3] + weights[4] * eval_coords[4][4] + weights[5] * eval_coords[4][5] + weights[6] * eval_coords[4][6],
        weights[0] * eval_coords[5][0] + weights[1] * eval_coords[5][1] + weights[2] * eval_coords[5][2] + weights[3] * eval_coords[5][3] + weights[4] * eval_coords[5][4] + weights[5] * eval_coords[5][5] + weights[6] * eval_coords[5][6],
        weights[0] * eval_coords[6][0] + weights[1] * eval_coords[6][1] + weights[2] * eval_coords[6][2] + weights[3] * eval_coords[6][3] + weights[4] * eval_coords[6][4] + weights[5] * eval_coords[6][5] + weights[6] * eval_coords[6][6],
    };

    return prediction;
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
constexpr 
inline std::array<T, k1DNumPredicitionValues>
Perform1DRBFPrediction(const std::array<T, k1DNumControlValues>& control_values)
{
    static_assert(false, "The RBF Prediction for 1D is currently not implemented!");
    return std::array<T, k1DNumPredicitionValues>();
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
constexpr 
inline std::array<T, k4DNumPredicitionValues>
Perform4DRBFPrediction(const std::array<T, k4DNumControlValues>& control_values)
{
    static_assert(false, "The RBF Prediction for 4D is currently not implemented!");
    return std::array<T, k4DNumPredicitionValues>();
}

}

}

#endif /* !CMC_PATCH_LOSSY_MULTI_RES_UTIL_HXX */
