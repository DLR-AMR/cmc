#ifndef CMC_T8_VARIABLES_HXX
#define CMC_T8_VARIABLES_HXX

#include "t8code/cmc_t8_data_variables_forward.hxx"
#include "utilities/cmc_variable_transformer_forward.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "t8code/cmc_t8_interpolation.hxx"
#include "t8code/cmc_t8_adapt_track_inaccuracy.hxx"
#include "utilities/cmc_compression_settings.hxx"

#ifdef CMC_WITH_NETCDF
#include "netcdf/cmc_netcdf.hxx"
#endif

#include <cstdlib>
#include <memory>
#include <type_traits>
#include <vector>

namespace cmc
{


constexpr int kGlobalContextInformationNotGiven = std::numeric_limits<int>::lowest();

struct PermittedError
{
    PermittedError() = delete;
    PermittedError(const CompressionCriterion etype, const double permitted_error)
    : criterion{etype}, error{permitted_error}{};

    const CompressionCriterion criterion{CompressionCriterion::CriterionUndefined};
    const double error{0.0};
};

struct ErrorCompliance
{
    ErrorCompliance() = delete;
    ErrorCompliance(const bool is_error_threshold_fulfilled, const double max_error)
    : is_error_threshold_satisfied{is_error_threshold_fulfilled}, max_introduced_error{max_error}{};

    const bool is_error_threshold_satisfied;
    const double max_introduced_error;
};

class VariableErrorDomains
{
public:
    VariableErrorDomains() = default;
    VariableErrorDomains(const double domain_specific_error, const GeoDomain& domain, const CompressionCriterion criterion)
    : specific_error_{domain_specific_error}, domain_{domain}, criterion_{criterion}{};
    ~VariableErrorDomains() = default;

    VariableErrorDomains(const VariableErrorDomains& other) = default;
    VariableErrorDomains& operator=(const VariableErrorDomains& other) = default;
    VariableErrorDomains(VariableErrorDomains&& other) = default;
    VariableErrorDomains& operator=(VariableErrorDomains&& other) = default;

    double GetError() const { return specific_error_;};
    const GeoDomain& GetDomain() const {return domain_;};
    CompressionCriterion GetCriterion() const {return criterion_;};
private:
    double specific_error_;
    GeoDomain domain_;
    CompressionCriterion criterion_;
};

template<class T>
class VariableAttributes
{
public:
    VariableAttributes() = default;
    VariableAttributes(const T missing_value, const DataLayout initial_layout, const DataLayout pre_compression_layout,
                       const int global_context_information)
    : missing_value_{missing_value}, initial_data_layout_{initial_layout}, pre_compression_layout_{pre_compression_layout},
      global_context_information_{global_context_information}{};
    ~VariableAttributes(){};

    VariableAttributes(const VariableAttributes& other) = default;
    VariableAttributes& operator=(const VariableAttributes& other) = default;
    VariableAttributes(VariableAttributes&& other) = default;
    VariableAttributes& operator=(VariableAttributes&& other) = default;

    void SetMissingValue(const T missing_value) {missing_value_ = missing_value;};
    T GetMissingValue() const {return missing_value_;};

    DataLayout GetInitialDataLayout() const {return initial_data_layout_;};

    DataLayout GetPreCompressionLayout() const {return pre_compression_layout_;};

    Dimension HasSplitDimension() const {return has_split_dimension_;};

    int GetGlobalContextInformation() const {return global_context_information_;};

    void SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const;

    friend class TransformerInputToCompressionVariable;
    friend class TransformerCompressionToOutputVariable;

private:
    T missing_value_{std::numeric_limits<T>::lowest()};
    CmcUniversalType add_offset_{static_cast<double>(0.0)};
    CmcUniversalType scale_factor_{static_cast<double>(1.0)};
    bool is_scaling_and_offset_applied_{true};
    DataLayout initial_data_layout_{DataLayout::LayoutUndefined};
    DataLayout pre_compression_layout_{DataLayout::LayoutUndefined};
    Dimension has_split_dimension_{Dimension::DimensionUndefined};
    int global_context_information_{kGlobalContextInformationNotGiven}; //!< A variable for describing an optional additional relation (the meaning of this varibale may dependent on the context)
};

template<class T>
class VariableUtilities
{
public:
    using ed_iterator = std::vector<VariableErrorDomains>::iterator;
    using ed_const_iterator = std::vector<VariableErrorDomains>::const_iterator;
    
    VariableUtilities() = default;
    ~VariableUtilities() = default;

    VariableUtilities(const VariableUtilities& other)
     : interpolate_{other.interpolate_},
       error_domains_(other.error_domains_),
       is_inaccuracy_storage_set_{other.is_inaccuracy_storage_set_},
       tracking_option_{other.tracking_option_} {
        	if (other.inaccuracy_storage_ != nullptr)
            {
                inaccuracy_storage_ = std::unique_ptr<InaccuracyContainer>(other.inaccuracy_storage_->clone());
            }
       };
    VariableUtilities(VariableUtilities&& other)
     : interpolate_{std::move(other.interpolate_)},
       error_domains_(std::move(other.error_domains_)),
       is_inaccuracy_storage_set_{std::move(other.is_inaccuracy_storage_set_)},
       tracking_option_{std::move(other.tracking_option_)} {
            if (other.inaccuracy_storage_ != nullptr)
            {
                inaccuracy_storage_ = std::move(other.inaccuracy_storage_);
            }
       };

    VariableUtilities& operator=(const VariableUtilities& other);
    VariableUtilities& operator=(VariableUtilities&& other);

    void SetInterpolation(Interpolate2<T> interpolate_function) {interpolate_ = interpolate_function;};
    InterpolationFunctional2<T>& GetInterpolation() {return interpolate_;};

    void SetUpVariableErrorDomains(const CompressionSpecifications& variable_specifications, const GeoDomain& variable_global_domain);

    T Interpolate(const VectorView<T>& values, const t8_forest_t previous_mesh, const t8_locidx_t start_index, const int num_elements, const T missing_value) const;
    ErrorCompliance IsCoarseningErrorCompliant(const std::vector<PermittedError>& permitted_errors, const VectorView<T>& previous_values, const std::vector<double>& previous_deviations, const T interpolated_value, const T missing_value) const;
    std::vector<ErrorCompliance> AreAlternativeValuesErrorCompliant(const std::vector<T>& alternative_values, const std::vector<PermittedError>& permitted_errors, const VectorView<T>& previous_values, const std::vector<double>& previous_deviations, const T missing_value) const;
    ErrorCompliance IsValueErrorCompliant(const std::vector<PermittedError>& permitted_errors, const T initial_value, const T nominal_value, const T& missing_value) const;

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

    ed_iterator GetErrorDomainsBegin() { return error_domains_.begin(); };
    ed_iterator GetErrorDomainsEnd() { return error_domains_.end(); };

    ed_const_iterator GetErrorDomainsBegin() const { return error_domains_.begin(); };
    ed_const_iterator GetErrorDomainsEnd() const { return error_domains_.end(); };

    ed_const_iterator GetErrorDomainsCBegin() const { return error_domains_.cbegin(); };
    ed_const_iterator GetErrorDomainsCEnd() const { return error_domains_.cend(); };

    friend class TransformerInputToCompressionVariable;
    friend class TransformerCompressionToOutputVariable;

private:
    const InaccuracyComputer<T>& CompressionCriterionToInaccuracyComputer(const CompressionCriterion criterion) const;

    InterpolationFunctional2<T> interpolate_{InterpoalteToArithmeticMean};

    InaccuracyComputer<T> compute_absolute_inaccuracy_{ComputeAbsoluteDeviation, ComputeSingleAbsoluteDeviation};
    InaccuracyComputer<T> compute_relative_inaccuracy_{ComputeRelativeDeviation, ComputeSingleRelativeDeviation};
    
    std::vector<VariableErrorDomains> error_domains_;

    bool is_inaccuracy_storage_set_{false};
    TrackingOption tracking_option_{TrackingOption::TrackFullInaccuracy};
    std::unique_ptr<InaccuracyContainer> inaccuracy_storage_{nullptr};
};



/** VARIABLE_ATTRIBUTES<T> MEMBER FUNCITONS **/
template<class T>
void
VariableAttributes<T>::SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const
{
    #ifdef CMC_WITH_NETCDF
    int err = nc_put_att(ncid, var_id, "mv", nc_type, 1, static_cast<const void*>(&missing_value_));
    NcCheckError(err);
    #endif
}


/** VARIABLE_UTILITIES<T> MEMBER FUNCITONS **/
template<class T>
VariableUtilities<T>&
VariableUtilities<T>::operator=(const VariableUtilities& other)
{
    interpolate_ = other.interpolate_;
    error_domains_ = other.error_domains_;
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
    error_domains_ = std::move(other.error_domains_);
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

/* Create the variable error domains for the variable. The general criterion considered on the whole domain is always put first within the 
vector of error domains */
template<class T>
void
VariableUtilities<T>::SetUpVariableErrorDomains(const CompressionSpecifications& variable_specifications, const GeoDomain& variable_global_domain)
{
    /* Create the general error domain */
    error_domains_.emplace_back(variable_specifications.general_max_error, variable_global_domain, variable_specifications.general_compression_criterion);
    
    /* Create all specific error domains */
    for (auto ed_iter = variable_specifications.specific_error_domains.begin(); ed_iter != variable_specifications.specific_error_domains.end(); ++ed_iter)
    {
        error_domains_.emplace_back(ed_iter->allowed_error, ed_iter->domain, ed_iter->compression_criterion);
    }
};

template<class T>
T
VariableUtilities<T>::Interpolate(const VectorView<T>& values, const t8_forest_t previous_mesh, const t8_locidx_t start_index, const int num_elements, const T missing_value) const
{
    return interpolate_(values, previous_mesh, start_index, num_elements, missing_value);
};

template<class T>
ErrorCompliance 
VariableUtilities<T>::IsCoarseningErrorCompliant(const std::vector<PermittedError>& permitted_errors, const VectorView<T>& previous_values, const std::vector<double>& previous_deviations, const T interpolated_value, const T missing_value) const
{
    double max_introduced_error = 0.0;

    /* Get the absolute inaccuracy for all values (the absolute error is needed in any case) */
    const std::vector<double> abs_inaccuracy = compute_absolute_inaccuracy_(previous_values, interpolated_value, previous_deviations, missing_value);

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
                const std::vector<double> rel_inaccuracy = compute_relative_inaccuracy_(previous_values, interpolated_value, previous_deviations, missing_value);
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
const InaccuracyComputer<T>&
VariableUtilities<T>::CompressionCriterionToInaccuracyComputer(const CompressionCriterion criterion) const
{
    switch (criterion)
    {
        case CompressionCriterion::RelativeErrorThreshold:
            return compute_relative_inaccuracy_;
        break;
        case CompressionCriterion::AbsoluteErrorThreshold:
            return compute_absolute_inaccuracy_;
        break;
        default:
            cmc_err_msg("The supplied compression criterion does not correspond to an inaccuracy computation funciton.");
            return compute_absolute_inaccuracy_;
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

}

#endif
