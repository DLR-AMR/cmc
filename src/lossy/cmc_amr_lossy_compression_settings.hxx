#ifndef CMC_AMR_LOSSY_COMPRESSION_SETTINGS_HXX
#define CMC_AMR_LOSSY_COMPRESSION_SETTINGS_HXX

/** @file cmc_amr_lossy_compressor.hxx
 */
#include "utilities/cmc_utilities.hxx"

#include <array>
#include <vector>

namespace cmc {

enum CompressionCriterion {CriterionUndefined, RelativeErrorThreshold, AbsoluteErrorThreshold, DomainExclusion, _CombinedCriteria};

enum CompressionMode {CompressionModeUndefined, OneForAll, OneForOne};

constexpr int kErrorCriterionHoldsForAllVariables = INT_MIN;

struct CompressionSpecifications
{
    CompressionCriterion compression_criterion;
    int variable_id{kErrorCriterionHoldsForAllVariables};
    double max_error{0.0};
};

struct DimensionInterval
{
    int start_index{0};
    int end_index{0};
    Dimension dim{Dimension::DimensionUndefined};
};

struct CertainErrorDomain
{
    std::vector<DimensionInterval> domain;
    double allowed_error{0.0};
};

struct SplitVariable
{
    int variable_id;
    DataLayout preferred_data_layout{DataLayout::LayoutUndefined};
};

class CompressionSettings
{
public:
    CompressionSettings(){};
    ~CompressionSettings(){};
    
    CompressionSettings(const CompressionSettings& other) = default;
    CompressionSettings& operator=(const CompressionSettings& other) = default;
    CompressionSettings(CompressionSettings&& other) = default;
    CompressionSettings& operator=(CompressionSettings&& other) = default;

    void SetRelativeErrorCriterion(const double max_error, const int variable_id = kErrorCriterionHoldsForAllVariables);
    void SetAbsoluteErrorCriterion(const double max_error, const int variable_id = kErrorCriterionHoldsForAllVariables);
    void SetDomainExclusionCriterion(std::vector<DimensionInterval> domain_to_exclude_from_compression);

    void SplitVariableByDimension(SplitVariable variable_to_split);

private:
    std::vector<CompressionSpecifications> compression_criteria_per_variable;
    std::vector<CertainErrorDomain> error_domains;
    std::vector<SplitVariable> variables_to_split;
};

}

#endif /* !CMC_AMR_LOSSY_COMPRESSION_SETTINGS_HXX */
