#ifndef LOSSY_CMC_MULTI_RES_EXTRACTION_RESIDUAL_COMPUTATION_HXX
#define LOSSY_CMC_MULTI_RES_EXTRACTION_RESIDUAL_COMPUTATION_HXX

#include "utilities/cmc_vector_view.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "t8code/cmc_t8_adapt_track_inaccuracy.hxx"

#include <type_traits>
#include <cstring>
#include <cstddef>
#include <utility>

namespace cmc::lossy::multi_res
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template <typename T>
auto
GetResidualsLZC(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint32_t), int>
{
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint32_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint32_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    /* Compute the unsigned residual and determine the LZC */
    const uint32_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);
    return diff_compression_val.GetNumberLeadingZeros();
}

template <typename T>
auto
GetResidualsLZC(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint16_t), int>
{
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint16_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint16_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    /* Compute the unsigned residual and determine the LZC */
    const uint16_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);
    return diff_compression_val.GetNumberLeadingZeros();
}

template <typename T>
auto
GetResidualsLZC(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint8_t), int>
{
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint8_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint8_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    /* Compute the unsigned residual and determine the LZC */
    const uint8_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);
    return diff_compression_val.GetNumberLeadingZeros();
}

template <typename T>
auto
GetResidualsLZC(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint64_t), int>
{
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint64_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint64_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    /* Compute the unsigned residual and determine the LZC */
    const uint64_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);
    return diff_compression_val.GetNumberLeadingZeros();
}

template <typename T>
auto
GetCumulativeResidualsLZC(const T& approximation, const VectorView<CompressionValue<T>>& real_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint32_t), int>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Get an unsigned representation of the approximation value */
    uint32_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = real_values.begin(); val_iter != real_values.end(); ++val_iter)
    {
        /* Get an unsigned representation of the real value */
        uint32_t unsigned_real_value;
        std::memcpy(&unsigned_real_value, val_iter->GetMemoryForReading().data(), N);

        /* Compute the unsigned residual and determine the LZC */
        const uint32_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();
    }

    return cumulative_lzc;
}

template <typename T>
auto
GetCumulativeResidualsLZC(const T& approximation, const VectorView<CompressionValue<T>>& real_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint16_t), int>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Get an unsigned representation of the approximation value */
    uint16_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = real_values.begin(); val_iter != real_values.end(); ++val_iter)
    {
        /* Get an unsigned representation of the real value */
        uint16_t unsigned_real_value;
        std::memcpy(&unsigned_real_value, val_iter->GetMemoryForReading().data(), N);

        /* Compute the unsigned residual and determine the LZC */
        const uint16_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();
    }

    return cumulative_lzc;
}

template <typename T>
auto
GetCumulativeResidualsLZC(const T& approximation, const VectorView<CompressionValue<T>>& real_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint8_t), int>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Get an unsigned representation of the approximation value */
    uint8_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = real_values.begin(); val_iter != real_values.end(); ++val_iter)
    {
        /* Get an unsigned representation of the real value */
        uint8_t unsigned_real_value;
        std::memcpy(&unsigned_real_value, val_iter->GetMemoryForReading().data(), N);

        /* Compute the unsigned residual and determine the LZC */
        const uint8_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();
    }

    return cumulative_lzc;
}

template <typename T>
auto
GetCumulativeResidualsLZC(const T& approximation, const VectorView<CompressionValue<T>>& real_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint64_t), int>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Get an unsigned representation of the approximation value */
    uint64_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = real_values.begin(); val_iter != real_values.end(); ++val_iter)
    {
        /* Get an unsigned representation of the real value */
        uint64_t unsigned_real_value;
        std::memcpy(&unsigned_real_value, val_iter->GetMemoryForReading().data(), N);

        /* Compute the unsigned residual and determine the LZC */
        const uint64_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();
    }

    return cumulative_lzc;
}


template <typename T>
auto
ComputeResidual(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint32_t), std::pair<bool, CompressionValue<T>>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint32_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint32_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    const bool is_approximation_greater = (unsigned_approximation >= unsigned_real_value);

    /* Compute the residual */
    const uint32_t diff = (unsigned_approximation >= unsigned_real_value ? unsigned_approximation - unsigned_real_value : unsigned_real_value - unsigned_approximation);

    /* Convert the difference to a CompressionValue */
    const CompressionValue<T> diff_compression_val(diff);

    /* Return the flag and the residual */
    return std::make_pair(is_approximation_greater, diff_compression_val);
}


template <typename T>
auto
ComputeResidual(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint8_t), std::pair<bool, CompressionValue<T>>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint8_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint8_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    const bool is_approximation_greater = (unsigned_approximation >= unsigned_real_value);

    /* Compute the residual */
    const uint8_t diff = (unsigned_approximation >= unsigned_real_value ? unsigned_approximation - unsigned_real_value : unsigned_real_value - unsigned_approximation);

    /* Convert the difference to a CompressionValue */
    const CompressionValue<T> diff_compression_val(diff);

    /* Return the flag and the residual */
    return std::make_pair(is_approximation_greater, diff_compression_val);
}

template <typename T>
auto
ComputeResidual(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint16_t), std::pair<bool, CompressionValue<T>>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint16_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint16_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    const bool is_approximation_greater = (unsigned_approximation >= unsigned_real_value);

    /* Compute the residual */
    const uint16_t diff = (unsigned_approximation >= unsigned_real_value ? unsigned_approximation - unsigned_real_value : unsigned_real_value - unsigned_approximation);

    /* Convert the difference to a CompressionValue */
    const CompressionValue<T> diff_compression_val(diff);

    /* Return the flag and the residual */
    return std::make_pair(is_approximation_greater, diff_compression_val);
}

template <typename T>
auto
ComputeResidual(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint64_t), std::pair<bool, CompressionValue<T>>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint64_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint64_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    const bool is_approximation_greater = (unsigned_approximation >= unsigned_real_value);

    /* Compute the residual */
    const uint64_t diff = (unsigned_approximation >= unsigned_real_value ? unsigned_approximation - unsigned_real_value : unsigned_real_value - unsigned_approximation);

    /* Convert the difference to a CompressionValue */
    const CompressionValue<T> diff_compression_val(diff);

    /* Return the flag and the residual */
    return std::make_pair(is_approximation_greater, diff_compression_val);
}


/////////////////////////////////////////////////

struct ErrorCompliance
{
    ErrorCompliance() = delete;
    ErrorCompliance(const bool is_error_threshold_fulfilled, const double max_error)
    : is_error_threshold_satisfied{is_error_threshold_fulfilled}, max_introduced_error{max_error}{};

    const bool is_error_threshold_satisfied;
    const double max_introduced_error;
};

template<class T>
ErrorCompliance
IsValueErrorCompliant(const std::vector<PermittedError>& permitted_errors, const T initial_value, const T nominal_value)
{
    cmc_assert(permitted_errors.size() > 0 && permitted_errors.size() <= 2);
    /* Get the absolute inaccuracy for all values (the absolute error is needed in any case) */
    const double abs_inaccuracy = ComputeSingleAbsoluteDeviation(initial_value, nominal_value);

    for (auto pe_iter = permitted_errors.begin(); pe_iter != permitted_errors.end(); ++pe_iter)
    {
        switch (pe_iter->criterion)
        {
            case CompressionCriterion::AbsoluteErrorThreshold:
            {
                /* Check if it is compliant with the permitted error */
                if (abs_inaccuracy > pe_iter->error)
                {
                    return ErrorCompliance(false, 0.0);
                }
            }
            break;
            case CompressionCriterion::RelativeErrorThreshold:
            {
                /* Get the relative inaccuracy for all values */
                const double rel_inaccuracy = ComputeSingleRelativeDeviation(initial_value, nominal_value);
                /* Check if it is compliant with the permitted error */
                if (rel_inaccuracy > pe_iter->error)
                {
                    return ErrorCompliance(false, 0.0);
                }
            }
            break;
                default:
                cmc_err_msg("The error specifications hold an unrecognized criterion.");
                return ErrorCompliance(false, 0.0);
        }
    }

    /* If the funciton reaches this point, the interpolated value complies with the permitted errors */
    return ErrorCompliance(true, abs_inaccuracy);
}


template<class T>
ErrorCompliance
IsValueErrorCompliantRegardingPreviousDeviations(const std::vector<PermittedError>& permitted_errors, const T initial_value, const T nominal_value, const double previous_abs_deviation) 
{
    cmc_assert(permitted_errors.size() > 0 && permitted_errors.size() <= 2);

    /* Get the absolute inaccuracy for all values (the absolute error is needed in any case) */
    const double abs_inaccuracy = ComputeAbsoluteDeviation(initial_value, nominal_value, previous_abs_deviation);

    for (auto pe_iter = permitted_errors.begin(); pe_iter != permitted_errors.end(); ++pe_iter)
    {
        switch (pe_iter->criterion)
        {
            case CompressionCriterion::AbsoluteErrorThreshold:
            {
                /* Check if it is compliant with the permitted error */
                if (abs_inaccuracy > pe_iter->error)
                {
                    return ErrorCompliance(false, 0.0);
                }
            }
            break;
            case CompressionCriterion::RelativeErrorThreshold:
            {
                /* Get the relative inaccuracy for all values */
                const double rel_inaccuracy = ComputeRelativeDeviation(initial_value, nominal_value, previous_abs_deviation);
                /* Check if it is compliant with the permitted error */
                if (rel_inaccuracy > pe_iter->error)
                {
                    return ErrorCompliance(false, 0.0);
                }
            }
            break;
                default:
                cmc_err_msg("The error specifications hold an unrecognized criterion.");
                return ErrorCompliance(false, 0.0);
        }
    }

    /* If the function reaches this point, the interpolated value complies with the permitted errors */
    return ErrorCompliance(true, previous_abs_deviation + abs_inaccuracy);
}


template<typename T>
std::pair<CompressionValue<T>, double>
GetMaximumTailClearedResidual(const std::vector<PermittedError>& permitted_errors, const double previous_abs_deviation, const T initial_value, const T approximation, const CompressionValue<T>& initial_residual, const bool is_approximation_greater)
{
    bool is_clearing_progressing = true;

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT - static_cast<int>(initial_residual.GetFrontBit());

    const CompressionValue<T> approximation_value(approximation);
    CompressionValue<T> residual = initial_residual;
    double inaccuracy{previous_abs_deviation};

    while (is_clearing_progressing && iteration_count < max_iteration_count)
    {
        /* Save the previous residual */
        const CompressionValue<T> save_previous_residual = residual;

        /* Clear the next set bit from the tail */
        residual.ClearNextSetBitFromTail();

        CompressionValue<T> residual_applied_value = approximation_value;

        if (is_approximation_greater)
        {
            /* In case the approximation is greater than the initial value, we need to subtract the residual */
            residual_applied_value.PerformIntegerSubtraction(residual);
        } else
        {
            /* In case the approximation is smaller than the initial value, we need to add the residual */
            residual_applied_value.PerformIntegerAddition(residual);
        }

        /* Reinterpret the approximation with the trimmed residual */
        const T trimmed_approximation_value = residual_applied_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = IsValueErrorCompliantRegardingPreviousDeviations<T>(permitted_errors, initial_value, trimmed_approximation_value, previous_abs_deviation);
    
        if (!error_evaluation.is_error_threshold_satisfied)
        {
            /* Revert the changes to the value */
            residual = save_previous_residual;
            is_clearing_progressing = false;
        } else
        {
            /* Store the compliant introduced absolute error */
            inaccuracy = error_evaluation.max_introduced_error;
        }

        ++iteration_count;
    }

    return std::make_pair(residual, inaccuracy);
}


template<typename T>
std::pair<CompressionValue<T>, double>
GetMaximumTailToggledResidual(const std::vector<PermittedError>& permitted_errors, const double previous_abs_deviation, const T initial_value, const T approximation, const CompressionValue<T>& initial_residual, const bool is_approximation_greater)
{
    bool is_toogling_progressing = true;

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT - static_cast<int>(initial_residual.GetFrontBit());

    const CompressionValue<T> approximation_value(approximation);
    CompressionValue<T> residual = initial_residual;
    double inaccuracy{previous_abs_deviation};

    while (is_toogling_progressing && iteration_count < max_iteration_count)
    {
        /* Save the previous residual */
        const CompressionValue<T> save_previous_residual = residual;

        /* Clear the next set bit from the tail */
        residual.ToggleTailUntilNextUnsetBit();

        CompressionValue<T> residual_applied_value = approximation_value;

        if (is_approximation_greater)
        {
            /* In case the approximation is greater than the initial value, we need to subtract the residual */
            residual_applied_value.PerformIntegerSubtraction(residual);
        } else
        {
            /* In case the approximation is smaller than the initial value, we need to add the residual */
            residual_applied_value.PerformIntegerAddition(residual);
        }

        /* Reinterpret the approximation with the trimmed residual */
        const T trimmed_approximation_value = residual_applied_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = IsValueErrorCompliantRegardingPreviousDeviations<T>(permitted_errors, initial_value, trimmed_approximation_value, previous_abs_deviation);
    
        if (!error_evaluation.is_error_threshold_satisfied)
        {
            /* Revert the changes to the value */
            residual = save_previous_residual;
            is_toogling_progressing = false;
        } else
        {
            /* Store the compliant introduced absolute error */
            inaccuracy = error_evaluation.max_introduced_error;
        }

        ++iteration_count;
    }

    return std::make_pair(residual, inaccuracy);
}


template<typename T>
std::pair<CompressionValue<T>, double>
GetMaximumTailClearedResidual(const std::vector<PermittedError>& permitted_errors, const T initial_value, const T approximation, const CompressionValue<T>& initial_residual, const bool is_approximation_greater)
{
    bool is_clearing_progressing = true;

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT - static_cast<int>(initial_residual.GetFrontBit());

    const CompressionValue<T> approximation_value(approximation);
    CompressionValue<T> residual = initial_residual;
    double inaccuracy{0.0};

    while (is_clearing_progressing && iteration_count < max_iteration_count)
    {
        /* Save the previous residual */
        const CompressionValue<T> save_previous_residual = residual;

        /* Clear the next set bit from the tail */
        residual.ClearNextSetBitFromTail();

        CompressionValue<T> residual_applied_value = approximation_value;

        if (is_approximation_greater)
        {
            /* In case the approximation is greater than the initial value, we need to subtract the residual */
            residual_applied_value.PerformIntegerSubtraction(residual);
        } else
        {
            /* In case the approximation is smaller than the initial value, we need to add the residual */
            residual_applied_value.PerformIntegerAddition(residual);
        }

        /* Reinterpret the approximation with the trimmed residual */
        const T trimmed_approximation_value = residual_applied_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = IsValueErrorCompliant<T>(permitted_errors, initial_value, trimmed_approximation_value);
    
        if (!error_evaluation.is_error_threshold_satisfied)
        {
            /* Revert the changes to the value */
            residual = save_previous_residual;
            is_clearing_progressing = false;
        } else
        {
            /* Store the compliant introduced absolute error */
            inaccuracy = error_evaluation.max_introduced_error;
        }

        ++iteration_count;
    }

    return std::make_pair(residual, inaccuracy);
}


template<typename T>
std::pair<CompressionValue<T>, double>
GetMaximumTailToggledResidual(const std::vector<PermittedError>& permitted_errors, const T initial_value, const T approximation, const CompressionValue<T>& initial_residual, const bool is_approximation_greater)
{
    bool is_toogling_progressing = true;

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT - static_cast<int>(initial_residual.GetFrontBit());

    const CompressionValue<T> approximation_value(approximation);
    CompressionValue<T> residual = initial_residual;
    double inaccuracy{0.0};

    while (is_toogling_progressing && iteration_count < max_iteration_count && residual.GetTailBit() < sizeof(T) * CHAR_BIT)
    {
        /* Save the previous residual */
        const CompressionValue<T> save_previous_residual = residual;

        /* Clear the next set bit from the tail */
        residual.ToggleTailUntilNextUnsetBit();

        CompressionValue<T> residual_applied_value = approximation_value;

        if (is_approximation_greater)
        {
            /* In case the approximation is greater than the initial value, we need to subtract the residual */
            residual_applied_value.PerformIntegerSubtraction(residual);
        } else
        {
            /* In case the approximation is smaller than the initial value, we need to add the residual */
            residual_applied_value.PerformIntegerAddition(residual);
        }

        /* Reinterpret the approximation with the trimmed residual */
        const T trimmed_approximation_value = residual_applied_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = IsValueErrorCompliant<T>(permitted_errors, initial_value, trimmed_approximation_value);
    
        if (!error_evaluation.is_error_threshold_satisfied)
        {
            /* Revert the changes to the value */
            residual = save_previous_residual;
            is_toogling_progressing = false;
        } else
        {
            /* Store the compliant introduced absolute error */
            inaccuracy = error_evaluation.max_introduced_error;
        }

        ++iteration_count;
    }

    return std::make_pair(residual, inaccuracy);
}

static int count_zero_res = 0;

template <typename T>
std::tuple<bool, CompressionValue<T>, double>
ComputeMinimalResidual(const T& approximation, const CompressionValue<T>& real_value, const std::vector<PermittedError>& permitted_errors, const double previous_absolute_error)
{
    /* Get the initial value within the correct type */
    const T initial_value = real_value.template ReinterpretDataAs<T>();

    /* Check whether the residual is smaller than the permitted error */
    const ErrorCompliance check_zero_residual = IsValueErrorCompliantRegardingPreviousDeviations<T>(permitted_errors, initial_value, approximation, previous_absolute_error);

    /* In case the residual is smaller than the permitted error */
    if (check_zero_residual.is_error_threshold_satisfied)
    {
        ++count_zero_res;
        /* In this case the residual is zero */
        return std::make_tuple(false, CompressionValue<T>(), check_zero_residual.max_introduced_error);
    }

    /* In case it the permitted error is smaller than the residual */

    /* Compute the initial residual */
    auto [is_approximation_greater, residual] = ComputeResidual<T>(approximation, real_value);

    /* We trim the residual by the means of clearing and toggeling the least significant bits from the residual */
    auto [tail_cleared_residual, tail_cleared_inaccurcy] = GetMaximumTailClearedResidual(permitted_errors, previous_absolute_error, initial_value, approximation, residual, is_approximation_greater);
    auto [tail_toggled_residual, tail_toggled_inaccuracy] = GetMaximumTailToggledResidual(permitted_errors, previous_absolute_error, initial_value, approximation, residual, is_approximation_greater);

    /* Check which approach leads to less significant bits */
    if (tail_cleared_residual.GetTailBit() >= tail_toggled_residual.GetTailBit())
    {
        /* If the tail cleared residual has less significant bits */
        return std::make_tuple(is_approximation_greater, tail_cleared_residual, tail_cleared_inaccurcy);
    } else
    {
        /* If the tail toggled residual has less significant bits */
        return std::make_tuple(is_approximation_greater, tail_toggled_residual, tail_toggled_inaccuracy);
    }
}


template <typename T>
std::tuple<bool, CompressionValue<T>, double>
ComputeMinimalResidual(const T& approximation, const CompressionValue<T>& real_value, const std::vector<PermittedError>& permitted_errors)
{
    /* Get the initial value within the correct type */
    const T initial_value = real_value.template ReinterpretDataAs<T>();

    /* Check whether the residual is smaller than the permitted error */
    const ErrorCompliance check_zero_residual = IsValueErrorCompliant<T>(permitted_errors, initial_value, approximation);

    /* In case the residual is smaller than the permitted error */
    if (check_zero_residual.is_error_threshold_satisfied)
    {
        /* In this case the residual is zero */
        return std::make_tuple(false, CompressionValue<T>(), check_zero_residual.max_introduced_error);
    }

    /* In case it the permitted error is smaller than the residual */

    /* Compute the initial residual */
    auto [is_approximation_greater, residual] = ComputeResidual<T>(approximation, real_value);

    /* We trim the residual by the means of clearing and toggeling the least significant bits from the residual */
    auto [tail_cleared_residual, tail_cleared_inaccurcy] = GetMaximumTailClearedResidual(permitted_errors, initial_value, approximation, residual, is_approximation_greater);
    auto [tail_toggled_residual, tail_toggled_inaccuracy] = GetMaximumTailToggledResidual(permitted_errors, initial_value, approximation, residual, is_approximation_greater);

    /* Check which approach lead to less significant bits */
    if (tail_cleared_residual.GetTailBit() >= tail_toggled_residual.GetTailBit())
    {
        /* If the tail cleared residual has less significant bits */
        return std::make_tuple(is_approximation_greater, tail_cleared_residual, tail_cleared_inaccurcy);
    } else
    {
        /* If the tail toggled residual has less significant bits */
        return std::make_tuple(is_approximation_greater, tail_toggled_residual, tail_toggled_inaccuracy);
    }
}

}

#endif /* !LOSSY_CMC_MULTI_RES_EXTRACTION_RESIDUAL_COMPUTATION_HXX */
