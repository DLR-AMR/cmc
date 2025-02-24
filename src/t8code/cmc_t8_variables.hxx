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

/** VARIABLE_ATTRIBUTES<T> MEMBER FUNCITONS **/
template<class T>
void
VariableAttributes<T>::SetMissingValueInNCFile(const int ncid, const int var_id, const int nc_type) const
{
    #ifdef CMC_WITH_NETCDF
    int err = nc_put_att(ncid, var_id, "mv", nc_type, 1, static_cast<const void*>(&missing_value_));
    nc::CheckError(err);
    #endif
}

}

#endif
