#ifndef CMC_ERROR_DOMAIN_FUNCTIONS_HXX
#define CMC_ERROR_DOMAIN_FUNCTIONS_HXX

#include "t8code/cmc_t8_mesh.hxx"

namespace cmc
{

/* Typedef for a general error domain function determining
 * whether a specific element is considered as "inside" the domain */
typedef bool IsElementInDomainFn(t8_forest_t forest, const t8_locidx_t which_tree, const t8_eclass_t tree_class, const t8_locidx_t lelement_id,
                                 const t8_scheme_c* ts, const t8_element_t* element);
namespace error_domain_fn
{

/* An error domain function resembling a global error criterion that holds for all elements */
inline bool
GeneralErrorCriterion([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] const t8_locidx_t which_tree, [[maybe_unused]] const t8_eclass_t tree_class,
                      [[maybe_unused]] const t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme_c* ts, [[maybe_unused]] const t8_element_t* element)
{
    return true;
}

}

}

#endif /* !CMC_ERROR_DOMAIN_FUNCTIONS_HXX */
