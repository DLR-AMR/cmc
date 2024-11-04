/**
 * @file cmc_compression_settings.cxx
 */

#include "utilities/cmc_compression_settings.hxx"

namespace cmc
{

static
std::vector<CompressionSpecifications>::iterator
FindVariableInSpecifications(const int var_id, std::vector<CompressionSpecifications>& specifications)
{
    for (auto iter = specifications.begin(); iter != specifications.end(); ++iter)
    {
        if (iter->variable_id == var_id)
        {
            return iter;
        }
    }

    return specifications.end();
}

[[maybe_unused]] static
std::vector<SplitVariable>::iterator
FindVariableInSplitVariables(const int var_id, std::vector<SplitVariable>& split_vars)
{
    for (auto iter = split_vars.begin(); iter != split_vars.end(); ++iter)
    {
        if (iter->variable_id == var_id)
        {
            return iter;
        }
    }

    return split_vars.end();
}

void
CompressionSettings::SetGeneralErrorCriterion(const CompressionCriterion criterion, const double max_error, const int variable_id)
{
    auto var_iter = FindVariableInSpecifications(variable_id, specifications_);

    if (var_iter != specifications_.end())
    {
        /* There are already some specifications for this variable */
        var_iter->general_compression_criterion = criterion;
        var_iter->general_max_error = max_error;
    } else
    {
        /* Up until now, no specifications have been set for the variable */
        specifications_.emplace_back(variable_id, criterion, max_error);
    }
}

void
CompressionSettings::SetRelativeErrorCriterion(const double max_error, const int variable_id)
{
    SetGeneralErrorCriterion(CompressionCriterion::RelativeErrorThreshold, max_error, variable_id);
}

void
CompressionSettings::SetAbsoluteErrorCriterion(const double max_error, const int variable_id)
{
    SetGeneralErrorCriterion(CompressionCriterion::AbsoluteErrorThreshold, max_error, variable_id);
}

void
CompressionSettings::SetCertainErrorForDomain(const CompressionCriterion criterion, const double error, const GeoDomain& domain, const int variable_id)
{
    auto var_iter = FindVariableInSpecifications(variable_id, specifications_);

    if (var_iter != specifications_.end())
    {
        var_iter->specific_error_domains.emplace_back(criterion, error, domain);
    } else
    {
        specifications_.emplace_back(variable_id);
        specifications_.back().specific_error_domains.emplace_back(criterion, error, domain);
    }
}

void
CompressionSettings::SetCertainErrorForDomain(const CompressionCriterion criterion, const double error, GeoDomain&& domain, const int variable_id)
{
    auto var_iter = FindVariableInSpecifications(variable_id, specifications_);

    if (var_iter != specifications_.end())
    {
        var_iter->specific_error_domains.emplace_back(criterion, error, std::move(domain));
    } else
    {
        specifications_.emplace_back(variable_id);
        specifications_.back().specific_error_domains.emplace_back(criterion, error, std::move(domain));
    }
}

void
CompressionSettings::SetCertainErrorForDomain(const CertainErrorDomain& error_domain, const int variable_id)
{
   auto var_iter = FindVariableInSpecifications(variable_id, specifications_);

    if (var_iter != specifications_.end())
    {
        var_iter->specific_error_domains.push_back(error_domain);
    } else
    {
        specifications_.emplace_back(variable_id);
        specifications_.back().specific_error_domains.push_back(error_domain);
    }
}

void
CompressionSettings::SetCertainErrorForDomain(CertainErrorDomain&& error_domain, const int variable_id)
{
    auto var_iter = FindVariableInSpecifications(variable_id, specifications_);

    if (var_iter != specifications_.end())
    {
        var_iter->specific_error_domains.push_back(std::move(error_domain));
    } else
    {
        specifications_.emplace_back(variable_id);
        specifications_.back().specific_error_domains.push_back(std::move(error_domain));
    }
}

void
CompressionSettings::SplitVariableByDimension(const SplitVariable& variable_to_split)
{
    cmc_assert(FindVariableInSplitVariables(variable_to_split.variable_id, variables_to_split_) == variables_to_split_.end());
    variables_to_split_.push_back(variable_to_split);
}

void
CompressionSettings::SplitVariableByDimension(SplitVariable&& variable_to_split)
{
    cmc_assert(FindVariableInSplitVariables(variable_to_split.variable_id, variables_to_split_) == variables_to_split_.end());
    variables_to_split_.emplace_back(std::move(variable_to_split));
}

bool 
CompressionSettings::AreThereVariablesToSplit() const
{
    return (variables_to_split_.size() > 0 ? true : false);
}

bool
CompressionSettings::AreTheSettingsValid() const
{
    bool is_criterion_for_all_vars_present = false;

    for (auto iter = specifications_.begin(); iter != specifications_.end(); ++iter)
    {
        if (is_criterion_for_all_vars_present && iter->variable_id == kErrorCriterionHoldsForAllVariables)
        {
            /* Only one criterion is possible for all variables */
            return false;
        } else if (iter->variable_id == kErrorCriterionHoldsForAllVariables)
        {
            is_criterion_for_all_vars_present = true;
        }
    }

    for (auto iter = specifications_.begin(); iter != specifications_.end(); ++iter)
    {
        for (auto ed_iter = iter->specific_error_domains.begin(); ed_iter != iter->specific_error_domains.end(); ++ed_iter)
        {
            if (ed_iter->compression_criterion == CompressionCriterion::CriterionUndefined)
            {
                return false;
            }
            if (!(ed_iter->domain.IsValid()))
            {
                return false;
            }
        }
    }

    for (auto iter = variables_to_split_.begin(); iter != variables_to_split_.end(); ++iter)
    {
        if (iter->split_dimension == Dimension::DimensionUndefined)
        {
            return false;
        }
    }

    return true;
}

}
