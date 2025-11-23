#ifndef CMC_INTERPOLATION_FN_HXX
#define CMC_INTERPOLATION_FN_HXX

#include "utilities/cmc_vector_view.hxx"

#include <numeric>
#include <type_traits>
#include <vector>

namespace cmc
{

template <typename T>
T InterpolateToArithmeticMean(const VectorView<T>& values)
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
T InterpolateToMidRange(const VectorView<T>& values)
{
    cmc_assert(std::is_arithmetic_v<T>);
    cmc_assert(!values.empty());
    if (values.size() == 1)
    {
        return values.front();
    }

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

template <typename T>
T InterpolateToArithmeticMean(const std::vector<T>& values)
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
T InterpolateToMidRange(const std::vector<T>& values)
{
    cmc_assert(std::is_arithmetic_v<T>);
    cmc_assert(!values.empty());
    if (values.size() == 1)
    {
        return values.front();
    }
    
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

#endif /* !CMC_INTERPOLATION_FN_HXX */
