#ifndef CMC_T8_ADAPT_TRACK_INACCURACY_FORWARD_HXX
#define CMC_T8_ADAPT_TRACK_INACCURACY_FORWARD_HXX
/**
 * @file cmc_t8_adapt_track_inaccuracy_forward.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_vector_view.hxx"

namespace cmc
{

enum TrackingOption {TrackFullInaccuracy = 0, TrackMinimalWorkingInaccuracy};

constexpr size_t kInvalidSizeHintForInaccuracyContainer = 0;

template<typename T>
using ComputeInaccuracy = std::vector<double>(*)(const VectorView<T>& values, const T& interpolated_value, const std::vector<double>& previous_absolute_max_deviation, const T& missing_value);

template<typename T>
using ComputeSingleInaccuracy = double(*)(const T value, const T interpolated_value, const T& missing_value);

template<class T>
class InaccuracyComputer;

template<typename T>
auto
ComputeRelativeDeviation(const VectorView<T>& values, const T& nominal_value, const double previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, std::vector<double>>;


template<typename T>
auto
ComputeRelativeDeviation(const VectorView<T>& values, const T& nominal_value, const double previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, std::vector<double>>;


template<typename T>
auto
ComputeAbsoluteDeviation(const VectorView<T>& values, const T& nominal_value, const double previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, std::vector<double>>;

template<typename T>
auto
ComputeAbsoluteDeviation(const VectorView<T>& values, const T& nominal_value, const double previous_absolute_max_deviation, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, std::vector<double>>;


template<typename T>
auto
ComputeSingleRelativeDeviation(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>;

template<typename T>
auto
ComputeSingleRelativeDeviation(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>;


template<typename T>
auto
ComputeSingleAbsoluteDeviation(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>;

template<typename T>
auto
ComputeSingleAbsoluteDeviation(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>;

class InaccuracyContainer;

/* Track the inaccurcy of all elements throughout the whole adaptation.
 * This allows (for example) for additional lossy comrpession steps after the AMR lossy compression */
class FullInaccuracyTracker;

/* Track the inaccurcay of the last adaptation step only. This saves storage, but after the compression,
 * we are not able to check the deviation of each element. 
 */
class MinimalInaccuracyTracker;

}

#endif
