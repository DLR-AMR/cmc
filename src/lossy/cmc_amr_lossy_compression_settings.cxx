/**
 * @file cmc_amr_lossy_compression_settings.cxx
 */

#include "lossy/cmc_amr_lossy_compression_settings.hxx"

inline void
CompressionSettings::SetRelativeErrorCriterion(const double max_error, const int variable_id)
{
    criterion_per_variable_.emplace_back(CompressionSpecifications(
        CompressionCriterion::RelativeErrorThreshold,
        variable_id,
        max_error
    ));
}

inline void
CompressionSettings::SetAbsoluteErrorCriterion(const double max_error, const int variable_id)
{
    criterion_per_variable_.emplace_back(CompressionSpecifications(
        CompressionCriterion::AbsoluteErrorThreshold,
        variable_id,
        max_error
    ));
}

inline void
CompressionSettings::SetCertainErrorForDomain(const double error, const GeoDomain& domain)
{
    error_domains_.emplace_back(CertainErrorDomain(error, domain));
}

inline void
CompressionSettings::SetCertainErrorForDomain(const double error, GeoDomain&& domain)
{
    error_domains_.emplace_back(CertainErrorDomain(error, std::move(domain)));
}

inline void
CompressionSettings::SetCertainErrorForDomain(const CertainErrorDomain& error_domain)
{
    error_domains_.push_back(error_domain);
}

inline void
CompressionSettings::SetCertainErrorForDomain(CertainErrorDomain&& error_domain)
{
    error_domains_.emplace_back(std::move(error_domain));
}

inline void
CompressionSettings::SplitVariableByDimension(const SplitVariable& variable_to_split)
{
    variables_to_split_.push_back(variable_to_split);
}

inline void
CompressionSettings::SplitVariableByDimension(SplitVariable&& variable_to_split)
{
    variables_to_split_.emplace_back(std::move(variable_to_split));
}

