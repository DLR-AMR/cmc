#ifndef CMC_T8_INTERPOLATION_HXX
#define CMC_T8_INTERPOLATION_HXX
/**
 * @file cmc_t8_interpolation.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_span.hxx"
#include "utilities/cmc_log_functions.h"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include "t8_forest.h"
#endif

#include <type_traits>

namespace cmc
{

constexpr bool kProhibitCoarseningOfMissingValues = false;

template <typename T>
using Interpolate2 = T(*)(const VectorView<T>& values, const t8_forest_t old_forest, const t8_locidx_t lelement_id, const int num_elements, const T missing_value);

template<class T>
class InterpolationFunctional2
{
public:
    InterpolationFunctional2() = delete;
    InterpolationFunctional2(const Interpolate2<T> interpolation_function)
    : interpolate_{interpolation_function}{};

    ~InterpolationFunctional2(){};
    
    InterpolationFunctional2(const InterpolationFunctional2& other) = default;
    InterpolationFunctional2& operator=(const InterpolationFunctional2& other) = default;
    InterpolationFunctional2(InterpolationFunctional2&& other) = default;
    InterpolationFunctional2& operator=(InterpolationFunctional2&& other) = default;

    T operator()(const VectorView<T>& values, const t8_forest_t old_forest, const t8_locidx_t lelement_id, const int num_elements, const T missing_value) const
    {
        return interpolate_(values, old_forest, lelement_id, num_elements, missing_value);
    };

    void SetInterpolation(const Interpolate2<T> interpolation_function) {interpolate_ = interpolation_function;};
private:
    Interpolate2<T> interpolate_;
};

template <typename T>
T InterpoalteToArithmeticMean(const VectorView<T>& values, [[maybe_unused]] const t8_forest_t old_forest, [[maybe_unused]] const t8_locidx_t lelement_id, const int num_elements, T missing_value)
{
    cmc_assert(std::is_arithmetic_v<T>);
    cmc_assert(!values.empty());

    T number_non_missing_values = 0;
    T sum = 0;

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        if (!ApproxCompare(*iter, missing_value))
        {
            sum += *iter;
            ++number_non_missing_values;
        }
    }

    const bool perform_interpolation = (kProhibitCoarseningOfMissingValues ? number_non_missing_values == num_elements : number_non_missing_values > 0);

    if (perform_interpolation)
    {
        return sum / number_non_missing_values;
    } else
    {
        return missing_value;
    }
}

}

#endif /* !CMC_T8_INTERPOLATION_HXX */
