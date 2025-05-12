#ifndef CMC_VARIABLE_UTILITIES_HXX
#define CMC_VARIABLE_UTILITIES_HXX

#include "t8code/cmc_t8_data_variables_forward.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "t8code/cmc_t8_interpolation.hxx"
#include "t8code/cmc_t8_adapt_track_inaccuracy.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "utilities/cmc_error_domain.hxx"

#include <vector>

namespace cmc
{

struct ErrorCompliance
{
    ErrorCompliance() = delete;
    ErrorCompliance(const bool is_error_threshold_fulfilled, const double max_error)
    : is_error_threshold_satisfied{is_error_threshold_fulfilled}, max_introduced_error{max_error}{};

    const bool is_error_threshold_satisfied;
    const double max_introduced_error;
};

template<class T>
class VariableUtilities
{
public:    
    VariableUtilities() = default;
    ~VariableUtilities() = default;

    VariableUtilities(const VariableUtilities& other)
     : interpolate_{other.interpolate_},
       is_inaccuracy_storage_set_{other.is_inaccuracy_storage_set_},
       tracking_option_{other.tracking_option_} {
        	if (other.inaccuracy_storage_ != nullptr)
            {
                inaccuracy_storage_ = std::unique_ptr<InaccuracyContainer>(other.inaccuracy_storage_->clone());
            }
       };
    VariableUtilities(VariableUtilities&& other)
     : interpolate_{std::move(other.interpolate_)},
       is_inaccuracy_storage_set_{std::move(other.is_inaccuracy_storage_set_)},
       tracking_option_{std::move(other.tracking_option_)} {
            if (other.inaccuracy_storage_ != nullptr)
            {
                inaccuracy_storage_ = std::move(other.inaccuracy_storage_);
            }
       };

    VariableUtilities& operator=(const VariableUtilities& other);
    VariableUtilities& operator=(VariableUtilities&& other);

    void SetInterpolation(Interpolate<T> interpolate_function) {interpolate_ = interpolate_function;};
    InterpolationFunctional<T>& GetInterpolation() {return interpolate_;};

    void SetInterpolation(InterpolateSkipMissingValues<T> interpolate_function) {interpolate_skip_missing_values_ = interpolate_function;};
    InterpolateSkipMissingValues<T>& GetInterpolationSkipMissingValues() {return interpolate_skip_missing_values_;};
    
    T Interpolation(const VectorView<T>& values, const t8_forest_t previous_mesh, const int tree_id, const t8_locidx_t start_index, const int num_element) const;
    T Interpolation(const VectorView<T>& values, const t8_forest_t previous_mesh, const int tree_id, const t8_locidx_t start_index, const int num_elements, const T missing_value) const;
    

    ErrorCompliance IsCoarseningErrorCompliantSkipMissingValues(const std::vector<PermittedError>& permitted_errors, const VectorView<T>& previous_values, const std::vector<double>& previous_deviations, const T interpolated_value, const T missing_value) const;
    ErrorCompliance IsCoarseningErrorCompliant(const std::vector<PermittedError>& permitted_errors, const VectorView<T>& previous_values, const std::vector<double>& previous_deviations, const T interpolated_value) const;

    //std::vector<ErrorCompliance> AreAlternativeValuesErrorCompliant(const std::vector<T>& alternative_values, const std::vector<PermittedError>& permitted_errors, const VectorView<T>& previous_values, const std::vector<double>& previous_deviations, const T missing_value) const;
    //ErrorCompliance IsValueErrorCompliant(const std::vector<PermittedError>& permitted_errors, const T initial_value, const T nominal_value, const T& missing_value) const;
    //ErrorCompliance IsValueErrorCompliantRegardingPreviousDeviations(const std::vector<PermittedError>& permitted_errors, const T initial_value, const T nominal_value, const double previous_abs_deviation, const T& missing_value) const;
    //ErrorCompliance IsCoarseningErrorCompliantRegardingInitialData(const std::vector<PermittedError>& permitted_errors, const std::vector<T>& initial_data, const std::vector<DomainIndex>& initial_elem_indices, const T interpolated_value, const T missing_value) const;
    
    void SetInaccuracyStorage(const TrackingOption& tracking_option) {tracking_option_ = tracking_option;};
    void SetUpInaccuracyStorage(const size_t size_hint = kInvalidSizeHintForInaccuracyContainer);
    bool IsInaccuracyStorageAllocated() const {return is_inaccuracy_storage_set_;};

    void StoreInaccuracy(const int index, const double error) {inaccuracy_storage_->StoreInaccuracy(index, error);};

    void TransferPreviousDeviations(const int start_index_previous_values, const int num_previous_values) {inaccuracy_storage_->TransferPreviousDeviations(start_index_previous_values, num_previous_values);};
    void TransferPreviousDeviation(const int start_index_previous_values) {inaccuracy_storage_->TransferPreviousDeviation(start_index_previous_values);};

    void AllocateDeviationStorage(const t8_locidx_t num_elements) {inaccuracy_storage_->AllocateDeviationStorage(num_elements);};
    std::vector<double> GetPreviousDeviations(const int start_index, const int num_elements) const;
    double GetPreviousDeviation(const int index) const;
    void SwitchDeviations() {inaccuracy_storage_->SwitchDeviations();};

    void RepartitionInaccuracyData(t8_forest_t initial_forest, t8_forest_t partitioned_forest);

    //friend class TransformerInputToCompressionVariable;
    //friend class TransformerCompressionToOutputVariable;

private:
    const InaccuracyComputerSkipMissingValues<T>& CompressionCriterionToInaccuracyComputerSkipMissingValues(const CompressionCriterion criterion) const;

    Interpolate<T> interpolate_{InterpolateToArithmeticMean};
    InterpolateSkipMissingValues<T> interpolate_skip_missing_values_{InterpolateToArithmeticMeanSkipMissingValues};

    InaccuracyComputerSkipMissingValues<T> compute_absolute_inaccuracy_skip_missing_values_{ComputeAbsoluteDeviationSkipMissingValues, ComputeSingleAbsoluteDeviationSkipMissingValues};
    InaccuracyComputerSkipMissingValues<T> compute_relative_inaccuracy_skip_missing_values_{ComputeRelativeDeviationSkipMissingValues, ComputeSingleRelativeDeviationSkipMissingValues};

    InaccuracyComputer<T> compute_absolute_inaccuracy_{ComputeAbsoluteDeviation, ComputeSingleAbsoluteDeviation};
    InaccuracyComputer<T> compute_relative_inaccuracy_{ComputeRelativeDeviation, ComputeSingleRelativeDeviation};

    bool is_inaccuracy_storage_set_{false};
    TrackingOption tracking_option_{TrackingOption::TrackFullInaccuracy};
    std::unique_ptr<InaccuracyContainer> inaccuracy_storage_{nullptr};
};


template<class T>
VariableUtilities<T>&
VariableUtilities<T>::operator=(const VariableUtilities& other)
{
    interpolate_ = other.interpolate_;
    is_inaccuracy_storage_set_ = other.is_inaccuracy_storage_set_;
    tracking_option_ = other.tracking_option_;
    inaccuracy_storage_.reset(other.inaccuracy_storage_->clone());
    return *this;
}

template<class T>
VariableUtilities<T>& 
VariableUtilities<T>::operator=(VariableUtilities&& other)
{
    interpolate_ = std::move(other.interpolate_);
    is_inaccuracy_storage_set_ = std::move(other.is_inaccuracy_storage_set_);
    tracking_option_ = std::move(other.tracking_option_);
    inaccuracy_storage_ = std::move(other.inaccuracy_storage_);
    return *this;
}

template<class T>
void
VariableUtilities<T>::SetUpInaccuracyStorage(const size_t size_hint)
{
    if (is_inaccuracy_storage_set_)
    {
        cmc_warn_msg("The inaccuracy container has to be set once. It cannot change throughout the compression!");
    } else
    {
        switch (tracking_option_)
        {
            case TrackingOption::TrackFullInaccuracy:
                inaccuracy_storage_ = std::make_unique<FullInaccuracyTracker>(size_hint);
            break;
            case TrackingOption::TrackMinimalWorkingInaccuracy:
                inaccuracy_storage_ = std::make_unique<MinimalInaccuracyTracker>(size_hint);
            break;
            default:
                cmc_err_msg("An invalid tracking option has been supplied for the InaccuracyContainer.");
        }
        is_inaccuracy_storage_set_ = true;
    }
}


template<class T>
T
VariableUtilities<T>::Interpolation(const VectorView<T>& values, const t8_forest_t previous_mesh, const int tree_id, const t8_locidx_t start_index, const int num_elements, const T missing_value) const
{
    return interpolate_skip_missing_values_(values, previous_mesh, tree_id, start_index, num_elements, missing_value);
}

template<class T>
T
VariableUtilities<T>::Interpolation(const VectorView<T>& values, const t8_forest_t previous_mesh, const int tree_id, const t8_locidx_t start_index, const int num_elements) const
{
    return interpolate_(values, previous_mesh, tree_id, start_index, num_elements);
}

template<class T>
ErrorCompliance 
VariableUtilities<T>::IsCoarseningErrorCompliantSkipMissingValues(const std::vector<PermittedError>& permitted_errors, const VectorView<T>& previous_values, const std::vector<double>& previous_deviations, const T interpolated_value, const T missing_value) const
{
    double max_introduced_error = 0.0;

    /* Get the absolute inaccuracy for all values (the absolute error is needed in any case) */
    const std::vector<double> abs_inaccuracy = compute_absolute_inaccuracy_skip_missing_values_(previous_values, interpolated_value, previous_deviations, missing_value);

    for (auto pe_iter = permitted_errors.begin(); pe_iter != permitted_errors.end(); ++pe_iter)
    {
        switch (pe_iter->criterion)
        {
            case CompressionCriterion::AbsoluteErrorThreshold:
            {
                /* Check if it is compliant with the permitted error */
                for (auto iacc_iter = abs_inaccuracy.begin(); iacc_iter != abs_inaccuracy.end(); ++iacc_iter)
                {
                    if (*iacc_iter > pe_iter->error)
                    {
                        return ErrorCompliance(false, 0.0);
                    } else if (*iacc_iter > max_introduced_error)
                    {
                        max_introduced_error = *iacc_iter;
                    }
                }
            }
            break;
            case CompressionCriterion::RelativeErrorThreshold:
            {
                /* Get the relative inaccuracy for all values */
                const std::vector<double> rel_inaccuracy = compute_relative_inaccuracy_skip_missing_values_(previous_values, interpolated_value, previous_deviations, missing_value);
                /* Check if it is compliant with the permitted error */
                int index = 0;
                for (auto iacc_iter = rel_inaccuracy.begin(); iacc_iter != rel_inaccuracy.end(); ++iacc_iter, ++index)
                {
                    if (*iacc_iter > pe_iter->error)
                    {
                        return ErrorCompliance(false, 0.0);
                    } else if (abs_inaccuracy[index] > max_introduced_error)
                    {
                        max_introduced_error = abs_inaccuracy[index];
                    }
                }
            }
            break;
                default:
                cmc_err_msg("The error specifications hold an unrecognized criterion.");
                return ErrorCompliance(false, 0.0);
        }
    }

    /* If the funciton reaches this point, the interpolated value complies with the permitted errors */
    return ErrorCompliance(true, max_introduced_error);
}


template<class T>
ErrorCompliance 
VariableUtilities<T>::IsCoarseningErrorCompliant(const std::vector<PermittedError>& permitted_errors, const VectorView<T>& previous_values, const std::vector<double>& previous_deviations, const T interpolated_value) const
{
    double max_introduced_error = 0.0;

    /* Get the absolute inaccuracy for all values (the absolute error is needed in any case) */
    const std::vector<double> abs_inaccuracy = compute_absolute_inaccuracy_(previous_values, interpolated_value, previous_deviations);

    for (auto pe_iter = permitted_errors.begin(); pe_iter != permitted_errors.end(); ++pe_iter)
    {
        switch (pe_iter->criterion)
        {
            case CompressionCriterion::AbsoluteErrorThreshold:
            {
                /* Check if it is compliant with the permitted error */
                for (auto iacc_iter = abs_inaccuracy.begin(); iacc_iter != abs_inaccuracy.end(); ++iacc_iter)
                {
                    if (*iacc_iter > pe_iter->error)
                    {
                        return ErrorCompliance(false, 0.0);
                    } else if (*iacc_iter > max_introduced_error)
                    {
                        max_introduced_error = *iacc_iter;
                    }
                }
            }
            break;
            case CompressionCriterion::RelativeErrorThreshold:
            {
                /* Get the relative inaccuracy for all values */
                const std::vector<double> rel_inaccuracy = compute_relative_inaccuracy_(previous_values, interpolated_value, previous_deviations);
                /* Check if it is compliant with the permitted error */
                int index = 0;
                for (auto iacc_iter = rel_inaccuracy.begin(); iacc_iter != rel_inaccuracy.end(); ++iacc_iter, ++index)
                {
                    if (*iacc_iter > pe_iter->error)
                    {
                        return ErrorCompliance(false, 0.0);
                    } else if (abs_inaccuracy[index] > max_introduced_error)
                    {
                        max_introduced_error = abs_inaccuracy[index];
                    }
                }
            }
            break;
                default:
                cmc_err_msg("The error specifications hold an unrecognized criterion.");
                return ErrorCompliance(false, 0.0);
        }
    }

    /* If the funciton reaches this point, the interpolated value complies with the permitted errors */
    return ErrorCompliance(true, max_introduced_error);
}

template<typename T>
std::vector<double> 
VariableUtilities<T>::GetPreviousDeviations(const int start_index, const int num_elements) const
{
    cmc_assert(inaccuracy_storage_ != nullptr);
    return inaccuracy_storage_->GetInaccuracyForRange(start_index, num_elements);
}

template<typename T>
double
VariableUtilities<T>::GetPreviousDeviation(const int index) const
{
    return inaccuracy_storage_->GetInaccuracy(index);
}

template<typename T>
const InaccuracyComputerSkipMissingValues<T>&
VariableUtilities<T>::CompressionCriterionToInaccuracyComputerSkipMissingValues(const CompressionCriterion criterion) const
{
    switch (criterion)
    {
        case CompressionCriterion::RelativeErrorThreshold:
            return compute_relative_inaccuracy_skip_missing_values_;
        break;
        case CompressionCriterion::AbsoluteErrorThreshold:
            return compute_absolute_inaccuracy_skip_missing_values_;
        break;
        default:
            cmc_err_msg("The supplied compression criterion does not correspond to an inaccuracy computation funciton.");
            return compute_absolute_inaccuracy_skip_missing_values_;
    }
}

template<typename T>
void
VariableUtilities<T>::RepartitionInaccuracyData(t8_forest_t initial_forest, t8_forest_t partitioned_forest)
{
    cmc_assert(is_inaccuracy_storage_set_);
    cmc_assert(inaccuracy_storage_ != nullptr);

    /* Repartition the tracked deviations compliant to the new partitioned forest */
    inaccuracy_storage_->RepartitionDeviations(initial_forest, partitioned_forest);
}


#if 0



template<class T>
ErrorCompliance 
VariableUtilities<T>::IsCoarseningErrorCompliantRegardingInitialData(const std::vector<PermittedError>& permitted_errors, const std::vector<T>& initial_data, const std::vector<DomainIndex>& initial_elem_indices, const T interpolated_value, const T missing_value) const 
{
    double max_introduced_error = 0.0;

    for (auto pe_iter = permitted_errors.begin(); pe_iter != permitted_errors.end(); ++pe_iter)
    {
        switch (pe_iter->criterion)
        {
            case CompressionCriterion::AbsoluteErrorThreshold:
            {
                for (size_t idx = 0; idx < initial_elem_indices.size(); ++idx)
                {
                    /* Compute the absolute deviation from the interpolated value to the initial value */
                    const double abs_error = ComputeSingleAbsoluteDeviation(initial_data[initial_elem_indices[idx]], interpolated_value, missing_value);
                    if (abs_error > pe_iter->error)
                    {
                        return ErrorCompliance(false, 0.0);
                    } else if (abs_error > max_introduced_error)
                    {
                        max_introduced_error = abs_error;
                    }
                }
            }
            break;
            case CompressionCriterion::RelativeErrorThreshold:
            {
                for (size_t idx = 0; idx < initial_elem_indices.size(); ++idx)
                {
                    /* Compute the absolute deviation from the interpolated value to the initial value */
                    const double rel_error = ComputeSingleRelativeDeviation(initial_data[initial_elem_indices[idx]], interpolated_value, missing_value);
                    if (rel_error > pe_iter->error)
                    {
                        return ErrorCompliance(false, 0.0);
                    } else if (rel_error > max_introduced_error)
                    {
                        max_introduced_error = rel_error;
                    }
                }
            }
            break;
                default:
                cmc_err_msg("The error specifications hold an unrecognized criterion.");
                return ErrorCompliance(false, 0.0);
        }
    }

    /* If the funciton reaches this point, the interpolated value complies with the permitted errors */
    return ErrorCompliance(true, max_introduced_error);
}

template<class T>
ErrorCompliance
VariableUtilities<T>::IsValueErrorCompliant(const std::vector<PermittedError>& permitted_errors, const T initial_value, const T nominal_value, const T& missing_value) const
{
    /* Get the absolute inaccuracy for all values (the absolute error is needed in any case) */
    const double abs_inaccuracy = compute_absolute_inaccuracy_(initial_value, nominal_value, missing_value);

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
                const double rel_inaccuracy = compute_relative_inaccuracy_(initial_value, nominal_value, missing_value);
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
VariableUtilities<T>::IsValueErrorCompliantRegardingPreviousDeviations(const std::vector<PermittedError>& permitted_errors, const T initial_value, const T nominal_value, const double previous_abs_deviation, const T& missing_value) const
{
    /* Get the absolute inaccuracy for all values (the absolute error is needed in any case) */
    const double abs_inaccuracy = ComputeAbsoluteDeviation(initial_value, nominal_value, previous_abs_deviation, missing_value);

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
                const double rel_inaccuracy = ComputeRelativeDeviation(initial_value, nominal_value, previous_abs_deviation, missing_value);
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

template<class T>
std::vector<ErrorCompliance> 
VariableUtilities<T>::AreAlternativeValuesErrorCompliant(const std::vector<T>& alternative_values, const std::vector<PermittedError>& permitted_errors, const VectorView<T>& previous_values, const std::vector<double>& previous_deviations, const T missing_value) const
{
    std::vector<ErrorCompliance> error_compliances;
    error_compliances.reserve(alternative_values.size());

    for (auto val_iter = alternative_values.begin(); val_iter != alternative_values.end(); ++val_iter)
    {
        error_compliances.push_back(IsCoarseningErrorCompliant(permitted_errors, previous_values, previous_deviations, *val_iter, missing_value));
    }

    return error_compliances;
}



#endif

}


#endif /* !CMC_VARIABLE_UTILITIES_HXX */
