#ifndef CMC_ERROR_DOMAIN_FUNCTIONS_HXX
#define CMC_ERROR_DOMAIN_FUNCTIONS_HXX

#include "t8code/cmc_t8_mesh.hxx"

namespace cmc
{

/* Typedef for a general error domain function determining
 * whether a specific element is considered as "inside" the domain */
typedef bool IsElementInDomainFn(t8_forest_t forest, int tree_id, int lelement_id,
                                 const t8_scheme_c* ts, const t8_element_t * element);

namespace error_domain_fn
{

/* An error domain function resembling a global error criterion that holds for all elements */
inline bool
GeneralErrorCriterion([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] int tree_id, [[maybe_unused]] int first_lelement_id,
                      [[maybe_unused]] const t8_scheme_c* ts, [[maybe_unused]] const t8_element_t * element)
{
    return true;
}

}

}

#endif /* !CMC_ERROR_DOMAIN_FUNCTIONS_HXX */
