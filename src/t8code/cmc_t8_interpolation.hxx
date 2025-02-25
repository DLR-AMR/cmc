#ifndef CMC_T8_INTERPOLATION_HXX
#define CMC_T8_INTERPOLATION_HXX
/**
 * @file cmc_t8_interpolation.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_vector_view.hxx"
#include "utilities/cmc_log_functions.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest.h>
#endif

#include <type_traits>

namespace cmc
{

constexpr bool kProhibitCoarseningOfMissingValues = false;

template <typename T>
using InterpolateSkipMissingValues = T(*)(const VectorView<T>& values, const t8_forest_t old_forest,const int ltree_id, const t8_locidx_t lelement_id, const int num_elements, const T missing_value);

template <typename T>
using Interpolate = T(*)(const VectorView<T>& values, const t8_forest_t old_forest, const int ltree_id, const t8_locidx_t lelement_id, const int num_elements);

template<class T>
class InterpolationFunctionalSkipMissingValues
{
public:
    InterpolationFunctionalSkipMissingValues() = delete;
    InterpolationFunctionalSkipMissingValues(const InterpolateSkipMissingValues<T> interpolation_function)
    : interpolate_{interpolation_function}{};

    ~InterpolationFunctionalSkipMissingValues(){};
    
    InterpolationFunctionalSkipMissingValues(const InterpolationFunctionalSkipMissingValues& other) = default;
    InterpolationFunctionalSkipMissingValues& operator=(const InterpolationFunctionalSkipMissingValues& other) = default;
    InterpolationFunctionalSkipMissingValues(InterpolationFunctionalSkipMissingValues&& other) = default;
    InterpolationFunctionalSkipMissingValues& operator=(InterpolationFunctionalSkipMissingValues&& other) = default;

    T operator()(const VectorView<T>& values, const t8_forest_t old_forest, const int ltree_id, const t8_locidx_t lelement_id, const int num_elements, const T missing_value) const
    {
        return interpolate_(values, old_forest, ltree_id, lelement_id, num_elements, missing_value);
    };

    void SetInterpolation(const InterpolateSkipMissingValues<T> interpolation_function) {interpolate_ = interpolation_function;};
private:
    InterpolateSkipMissingValues<T> interpolate_;
};


template<class T>
class InterpolationFunctional
{
public:
    InterpolationFunctional() = delete;
    InterpolationFunctional(const Interpolate<T> interpolation_function)
    : interpolate_{interpolation_function}{};

    ~InterpolationFunctional(){};
    
    InterpolationFunctional(const InterpolationFunctional& other) = default;
    InterpolationFunctional& operator=(const InterpolationFunctional& other) = default;
    InterpolationFunctional(InterpolationFunctional&& other) = default;
    InterpolationFunctional& operator=(InterpolationFunctional&& other) = default;

    T operator()(const VectorView<T>& values, const t8_forest_t old_forest, const int ltree_id, const t8_locidx_t lelement_id, const int num_elements, const T missing_value) const
    {
        return interpolate_(values, old_forest, ltree_id, lelement_id, num_elements, missing_value);
    };

    void SetInterpolation(const Interpolate<T> interpolation_function) {interpolate_ = interpolation_function;};
private:
    Interpolate<T> interpolate_;
};

template <typename T>
T InterpolateToArithmeticMeanSkipMissingValues(const VectorView<T>& values, [[maybe_unused]] const t8_forest_t old_forest, [[maybe_unused]] const int ltree_id, [[maybe_unused]] const t8_locidx_t lelement_id, [[maybe_unused]] const int num_elements, T missing_value)
{
    cmc_assert(std::is_arithmetic_v<T>);
    cmc_assert(!values.empty());

    double number_non_missing_values = 0;
    double sum = 0;

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        if (!ApproxCompare(*iter, missing_value))
        {
            sum += static_cast<double>(*iter);
            ++number_non_missing_values;
        }
    }

    const bool perform_interpolation = (kProhibitCoarseningOfMissingValues ? number_non_missing_values == num_elements : number_non_missing_values > 0);

    if (perform_interpolation)
    {
        return static_cast<T>(sum / number_non_missing_values);
    } else
    {
        return missing_value;
    }
}

template <typename T>
T InterpolateToMidRangeSkipMissingValues(const VectorView<T>& values, [[maybe_unused]] const t8_forest_t old_forest, [[maybe_unused]] const int ltree_id, [[maybe_unused]] const t8_locidx_t lelement_id, [[maybe_unused]] const int num_elements, T missing_value)
{
    cmc_assert(std::is_arithmetic_v<T>);
    cmc_assert(!values.empty());

    T min = std::numeric_limits<T>::max();
    T max = std::numeric_limits<T>::min();

    bool min_found = false;
    bool max_found = false;

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        if (!ApproxCompare(*iter, missing_value))
        {
            if (min > *iter)
            {
                min = *iter;
                min_found = true;
            }
            if (max < *iter)
            {
                max = *iter;
                max_found = true;
            }
        }
    }

    if (min_found && max_found)
    {
        return ((max + min) / 2);
    } else if (min_found && !max_found)
    {
        return min;
    } else if (!min_found && max_found)
    {
        return max;
    } else
    {
        return missing_value;
    }
}


template <typename T>
T InterpolateToArithmeticMean(const VectorView<T>& values, [[maybe_unused]] const t8_forest_t old_forest, [[maybe_unused]] const int ltree_id, [[maybe_unused]] const t8_locidx_t lelement_id, [[maybe_unused]] const int num_elements)
{
    cmc_assert(std::is_arithmetic_v<T>);
    cmc_assert(!values.empty());

    double sum = 0;

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        sum += static_cast<double>(*iter);
    }

    return static_cast<T>(sum / static_cast<double>(values.size()));
}

template <typename T>
T InterpolateToMidRange(const VectorView<T>& values, [[maybe_unused]] const t8_forest_t old_forest, [[maybe_unused]] const int ltree_id, [[maybe_unused]] const t8_locidx_t lelement_id, [[maybe_unused]] const int num_elements)
{
    cmc_assert(std::is_arithmetic_v<T>);
    cmc_assert(!values.empty());
    cmc_assert(values.size() >= 2);

    T min = std::numeric_limits<T>::max();
    T max = std::numeric_limits<T>::min();

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        if (min > *iter)
        {
            min = *iter;
        }
        if (max < *iter)
        {
            max = *iter;
        }
    }

    return ((max / 2) + (min / 2));
}

}

#endif /* !CMC_T8_INTERPOLATION_HXX */
