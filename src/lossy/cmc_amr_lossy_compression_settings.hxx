#ifndef CMC_AMR_LOSSY_COMPRESSION_SETTINGS_HXX
#define CMC_AMR_LOSSY_COMPRESSION_SETTINGS_HXX

/** @file cmc_amr_lossy_compression_settings.hxx
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
    constexpr CompressionSpecifications(const CompressionCriterion criterion, int var_id, double allowed_max_error)
    : compression_criterion{criterion}, variable_id{var_id}, max_error{allowed_max_error}{};
    CompressionCriterion compression_criterion;
    int variable_id{kErrorCriterionHoldsForAllVariables};
    double max_error{0.0};
};

struct CertainErrorDomain
{
    constexpr CertainErrorDomain(double error, const GeoDomain& domain_specs)
    : allowed_error{error}, domain{domain_specs}{};
    CertainErrorDomain(double error, GeoDomain&& domain_specs)
    : allowed_error{error}, domain{std::move(domain_specs)}{};

    double allowed_error{0.0};
    GeoDomain domain;
};

struct SplitVariable
{
    constexpr SplitVariable(int var_id, DataLayout preferred_layout)
    : variable_id{var_id}, preferred_data_layout{preferred_layout}{};
    int variable_id{-1};
    DataLayout preferred_data_layout{DataLayout::LayoutUndefined};
};

class CompressionSettings
{
public:
    CompressionSettings() = default;
    
    CompressionSettings(const CompressionSettings& other) = default;
    CompressionSettings& operator=(const CompressionSettings& other) = default;
    CompressionSettings(CompressionSettings&& other) = default;
    CompressionSettings& operator=(CompressionSettings&& other) = default;

    ~CompressionSettings() = default;

    /* Set an error criterion (either a relative or absolute criterion) for a single or all variables */
    void SetRelativeErrorCriterion(const double max_error, const int variable_id = kErrorCriterionHoldsForAllVariables);
    void SetAbsoluteErrorCriterion(const double max_error, const int variable_id = kErrorCriterionHoldsForAllVariables);
    
    /* Set an exclusive error for a specified domain, an error of 0.0 excludes the domain from the compression */
    void SetCertainErrorForDomain(const double error, const GeoDomain& domain);
    void SetCertainErrorForDomain(const double error, GeoDomain&& domain);
    void SetCertainErrorForDomain(const CertainErrorDomain& error_domain);
    void SetCertainErrorForDomain(CertainErrorDomain&& error_domain);

    /* Split a variable at a specific dimension (e.g. 3D variable "Lat x Lon x Height" may be split into #Height 2D variables defined on the domain "Lat x Lon") */
    void SplitVariableByDimension(const SplitVariable& variable_to_split);
    void SplitVariableByDimension(SplitVariable&& variable_to_split);

private:
    std::vector<CompressionSpecifications> criterion_per_variable_;
    std::vector<CertainErrorDomain> error_domains_;
    std::vector<SplitVariable> variables_to_split_;
};

}

#endif /* !CMC_AMR_LOSSY_COMPRESSION_SETTINGS_HXX */
