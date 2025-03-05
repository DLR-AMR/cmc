#ifndef CMC_ERROR_DOMAIN_HXX
#define CMC_ERROR_DOMAIN_HXX

#include "t8code/cmc_t8_mesh.hxx"
#include "utilities/cmc_error_domain_functions.hxx"

#include <vector>
#include <utility>

namespace cmc
{

enum CompressionCriterion {CriterionUndefined, RelativeErrorThreshold, AbsoluteErrorThreshold};

struct PermittedError
{
    PermittedError() = delete;
    PermittedError(const CompressionCriterion etype, const double permitted_error)
    : criterion{etype}, error{permitted_error}{};

    const CompressionCriterion criterion{CompressionCriterion::CriterionUndefined};
    const double error{0.0};
};

class ErrorDomain
{
public:
    ErrorDomain() = delete;

    ErrorDomain(const PermittedError& error_criterion, const IsElementInDomainFn* is_inside_domain_check)
    : error_criterion_{error_criterion}, is_inside_domain_check_{is_inside_domain_check}{};
    
    inline bool
    IsElementWithinDomain(t8_forest_t forest, int tree_id, int lelement_id,
                          t8_eclass_scheme_c* ts, const t8_element_t* element) const
    {
        return is_inside_domain_check_(forest, tree_id, lelement_id, ts, element);
    }
    
    inline bool
    IsAnyElementWithinDomain(t8_forest_t forest, int tree_id, int lelement_id,
                             t8_eclass_scheme_c* ts, const int num_elements, const t8_element_t* elements[]) const
    {
        for (int idx = 0; idx < num_elements; ++idx)
        {
            const bool is_inside = is_inside_domain_check_(forest, tree_id, lelement_id + idx, ts, elements[idx]);

            if (is_inside == true)
            {
                return true;
            }
        }
        return false;
    }

    inline bool
    AreAllElementsWithinDomain(t8_forest_t forest, int tree_id, int lelement_id,
                              t8_eclass_scheme_c* ts, const int num_elements, const t8_element_t* elements[]) const
    {
        for (int idx = 0; idx < num_elements; ++idx)
        {
            const bool is_inside = is_inside_domain_check_(forest, tree_id, lelement_id + idx, ts, elements[idx]);

            if (is_inside == false)
            {
                return false;
            }
        }
        return true;
    }

    inline PermittedError GetPermittedError() const {return error_criterion_;};
    
    inline CompressionCriterion GetCompressionCriterion() const {return error_criterion_.criterion;};

    inline double GetErrorThreshold() const {return error_criterion_.error;};

private:
    const PermittedError error_criterion_;
    const IsElementInDomainFn* is_inside_domain_check_;
};

}

#endif /* !CMC_ERROR_DOMAIN_HXX */
