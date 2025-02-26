#ifndef CMC_COMPRESSION_SETTINGS_HXX
#define CMC_COMPRESSION_SETTINGS_HXX

/** @file cmc_compression_settings.hxx
 */
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_error_domain.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#include <vector>

namespace cmc {

class CompressionSettings
{
public:
    void SetGeneralErrorCriterion(const CompressionCriterion criterion, const double permitted_error)
    {
        error_domains_.emplace_back(PermittedError{criterion, permitted_error}, error_domain_fn::GeneralErrorCriterion);
    }
    void SetErrorDomain(const ErrorDomain& error_domain){error_domains_.push_back(error_domain);};
    void SetErrorDomain(ErrorDomain&& error_domain){error_domains_.push_back(std::move(error_domain));};

    std::vector<PermittedError>
    FindRestrictingErrors(t8_forest_t forest, int tree_id, int lelement_id,
                          t8_eclass_scheme_c* ts, const int num_elements, const t8_element_t* elements[]) const;

private:
    std::vector<ErrorDomain> error_domains_;
};


std::vector<PermittedError>
CompressionSettings::FindRestrictingErrors(t8_forest_t forest, int tree_id, int lelement_id,
                                           t8_eclass_scheme_c* ts, const int num_elements, const t8_element_t* elements[]) const
{
    bool is_absolute_error_found = false;
    double min_absolute_error = std::numeric_limits<double>::max();
    bool is_relative_error_found = false;
    double min_relative_error = std::numeric_limits<double>::max();

    /* Iterate over all domains and find the minimum applicable error criterion */
    for (auto error_domain_iter = error_domains_.begin(); error_domain_iter != error_domains_.end(); ++error_domain_iter)
    {
        /* Check if any of the considered elements is within the domain of this error criterion */
        const bool is_any_element_in_domain = error_domain_iter->IsAnyElementWithinDomain(forest, tree_id,lelement_id, ts, num_elements, elements);

        if (is_any_element_in_domain == true)
        {
            const CompressionCriterion type = error_domain_iter->GetCompressionCriterion();
            const double error = error_domain_iter->GetErrorThreshold();

            /* Check whether an absolute or relative criterion is defined for this area */
           if (type == CompressionCriterion::AbsoluteErrorThreshold)
           {
                is_absolute_error_found = true;
                if (min_absolute_error > error)
                {
                    min_absolute_error = error;
                }
           } else if (type == CompressionCriterion::RelativeErrorThreshold)
           {
                is_relative_error_found = true;
                if (min_relative_error > error)
                {
                    min_relative_error = error;
                }
           }
        }
    }

    /* Accumualate all restraining compression criteria */
    std::vector<PermittedError> permitted_errors;
    permitted_errors.reserve(2);

    if (is_absolute_error_found)
    {
        permitted_errors.emplace_back(CompressionCriterion::AbsoluteErrorThreshold, min_absolute_error);
    }

    if (is_relative_error_found)
    {
        permitted_errors.emplace_back(CompressionCriterion::RelativeErrorThreshold, min_relative_error);
    }

    return permitted_errors;
}

}

#endif /* !CMC_COMPRESSION_SETTINGS_HXX */
