#ifndef CMC_EMBEDDED_VARIABLE_ATTRIBUTES_HXX
#define CMC_EMBEDDED_VARIABLE_ATTRIBUTES_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"

#include <limits>

namespace cmc
{

constexpr int kGlobalContextInformationNotGiven = std::numeric_limits<int>::lowest();

template<class T>
class VariableAttributes
{
public:
    VariableAttributes() = default;
    VariableAttributes(const GeoDomain& global_domain, const T missing_value, const DataLayout initial_layout, const DataLayout pre_compression_layout,
                       const int global_context_information)
    :global_domain_{global_domain}, missing_value_{missing_value}, initial_data_layout_{initial_layout}, pre_compression_layout_{pre_compression_layout},
      global_context_information_{global_context_information}{};
    ~VariableAttributes(){};

    VariableAttributes(const VariableAttributes& other) = default;
    VariableAttributes& operator=(const VariableAttributes& other) = default;
    VariableAttributes(VariableAttributes&& other) = default;
    VariableAttributes& operator=(VariableAttributes&& other) = default;

    const GeoDomain& GetGlobalDomain() const {return global_domain_;};
    void SetGlobalDomain(const GeoDomain& global_domain) {global_domain_ = global_domain;};
    void SetGlobalDomain(GeoDomain&& global_domain) {global_domain_ = std::move(global_domain);};

    void SetMissingValue(const T missing_value) {missing_value_ = missing_value;};
    T GetMissingValue() const {return missing_value_;};

    DataLayout GetInitialDataLayout() const {return initial_data_layout_;};

    DataLayout GetPreCompressionLayout() const {return pre_compression_layout_;};

    Dimension HasSplitDimension() const {return has_split_dimension_;};

    int GetGlobalContextInformation() const {return global_context_information_;};

private:
    GeoDomain global_domain_;
    T missing_value_{std::numeric_limits<T>::lowest()};
    CmcUniversalType add_offset_{static_cast<double>(0.0)};
    CmcUniversalType scale_factor_{static_cast<double>(1.0)};
    bool is_scaling_and_offset_applied_{true};
    DataLayout initial_data_layout_{DataLayout::LayoutUndefined};
    DataLayout pre_compression_layout_{DataLayout::LayoutUndefined};
    Dimension has_split_dimension_{Dimension::DimensionUndefined};
    int global_context_information_{kGlobalContextInformationNotGiven}; //!< A variable for describing an optional additional relation (the meaning of this varibale may dependent on the context)
};



}

#endif /* !CMC_EMBEDDED_VARIABLE_ATTRIBUTES_HXX */