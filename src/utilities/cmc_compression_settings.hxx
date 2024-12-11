#ifndef CMC_COMPRESSION_SETTINGS_HXX
#define CMC_COMPRESSION_SETTINGS_HXX

/** @file cmc_compression_settings.hxx
 */
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"

namespace cmc {

enum CompressionCriterion {CriterionUndefined, RelativeErrorThreshold, AbsoluteErrorThreshold, DomainExclusion, _CombinedCriteria};

enum CompressionMode {CompressionModeUndefined, OneForAll, OneForOne};

constexpr int kErrorCriterionHoldsForAllVariables = std::numeric_limits<int>::lowest();
constexpr int kSplitAllVariables = std::numeric_limits<int>::lowest();

struct CertainErrorDomain
{
    constexpr CertainErrorDomain(const CompressionCriterion criterion, const double error, const GeoDomain& domain_specs)
    : compression_criterion{criterion},allowed_error{error}, domain{domain_specs}{};
    constexpr CertainErrorDomain(const CompressionCriterion criterion, const double error, GeoDomain&& domain_specs)
    : compression_criterion{criterion}, allowed_error{error}, domain{std::move(domain_specs)}{};

    CompressionCriterion compression_criterion{CompressionCriterion::CriterionUndefined};
    double allowed_error{0.0};
    GeoDomain domain;
};

struct CompressionSpecifications
{
    CompressionSpecifications(const int var_id)
    : variable_id{var_id}{};
    CompressionSpecifications(const int var_id, const CompressionCriterion criterion, const double allowed_max_error)
    : variable_id{var_id}, general_compression_criterion{criterion}, general_max_error{allowed_max_error}{};

    int GetVariableID() const {return variable_id;};
    
    int variable_id{kErrorCriterionHoldsForAllVariables};
    CompressionCriterion general_compression_criterion;
    double general_max_error{0.0};
    std::vector<CertainErrorDomain> specific_error_domains;
};

struct SplitVariable
{
    constexpr SplitVariable(int var_id, const Dimension dimension_to_split)
    : variable_id{var_id}, split_dimension{dimension_to_split}{};

    Dimension GetSplitDimension() const {return split_dimension;};
    int GetVariableID() const {return variable_id;};

    int variable_id{-1};
    Dimension split_dimension{Dimension::DimensionUndefined};
};

class CompressionSettings
{
public:
    using sv_iterator = std::vector<SplitVariable>::iterator;
    using sv_const_iterator = std::vector<SplitVariable>::const_iterator;

    using cs_iterator = std::vector<CompressionSpecifications>::iterator;
    using cs_const_iterator = std::vector<CompressionSpecifications>::const_iterator;

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
    void SetCertainErrorForDomain(const CompressionCriterion criterion, double error, const GeoDomain& domain, const int variable_id = kErrorCriterionHoldsForAllVariables);
    void SetCertainErrorForDomain(const CompressionCriterion criterion, double error, GeoDomain&& domain, const int variable_id = kErrorCriterionHoldsForAllVariables);
    void SetCertainErrorForDomain(const CertainErrorDomain& error_domain, const int variable_id = kErrorCriterionHoldsForAllVariables);
    void SetCertainErrorForDomain(CertainErrorDomain&& error_domain, const int variable_id = kErrorCriterionHoldsForAllVariables);

    /* Split a variable at a specific dimension (e.g. 3D variable "Lat x Lon x Height" may be split into #Height 2D variables defined on the domain "Lat x Lon") */
    void SplitVariableByDimension(const SplitVariable& variable_to_split);
    void SplitVariableByDimension(SplitVariable&& variable_to_split);
    void SetSplitVariables(const std::vector<SplitVariable>& split_variables);
    void SetSplitVariables(std::vector<SplitVariable>&& split_variables);
    sv_iterator GetSplitVariablesBegin() { return variables_to_split_.begin(); };
    sv_iterator GetSplitVariablesEnd() { return variables_to_split_.end(); };
    sv_const_iterator GetSplitVariablesBegin() const { return variables_to_split_.begin(); };
    sv_const_iterator GetSplitVariablesEnd() const { return variables_to_split_.end(); };
    sv_const_iterator GetSplitVariablesCBegin() const { return variables_to_split_.cbegin(); };
    sv_const_iterator GetSplitVariablesCEnd() const { return variables_to_split_.cend(); };

    cs_iterator GetSpecificationsBegin() { return specifications_.begin(); };
    cs_iterator GetSpecificationsEnd() { return specifications_.end(); };
    cs_const_iterator GetSpecificationsBegin() const { return specifications_.begin(); };
    cs_const_iterator GetSpecificationsEnd() const { return specifications_.end(); };
    cs_const_iterator GetSpecificationsCBegin() const { return specifications_.cbegin(); };
    cs_const_iterator GetSpecificationsCEnd() const { return specifications_.cend(); };

    bool AreThereVariablesToSplit() const;
    bool AreTheSettingsValid() const;

private:
    void SetGeneralErrorCriterion(const CompressionCriterion, const double, const int);

    std::vector<CompressionSpecifications> specifications_;
    std::vector<SplitVariable> variables_to_split_;
};

}

#endif /* !CMC_COMPRESSION_SETTINGS_HXX */
