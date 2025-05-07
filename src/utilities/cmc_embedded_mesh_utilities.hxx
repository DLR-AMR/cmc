#ifndef CMC_EMBEDDED_MESH_UTILITIES_HXX
#define CMC_EMBEDDED_MESH_UTILITIES_HXX

#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_log_functions.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx> 
#include <t8_forest/t8_forest_iterate.h> 
#include <t8_forest/t8_forest_partition.h>
#include <t8_schemes/t8_default/t8_default.hxx>

#endif

#include <array>
#include <tuple>

namespace cmc
{

struct AdaptDataInitialEmbeddedMesh
{
public:
    AdaptDataInitialEmbeddedMesh() = delete;
    AdaptDataInitialEmbeddedMesh(const GeoDomain& domain, const int initial_refinement_lvl, const DataLayout layout)
    : global_domain{domain}, initial_refinement_level{initial_refinement_lvl}, initial_layout{layout}{};

    const GeoDomain& global_domain;
    const int initial_refinement_level;
    const DataLayout initial_layout;
};

inline t8_locidx_t
RefineToInitialEmbeddedMesh (t8_forest_t forest,
                             [[maybe_unused]] t8_forest_t forest_from,
                             t8_locidx_t which_tree,
                             const t8_eclass_t tree_class,
                             t8_locidx_t lelement_id,
                             const t8_scheme_c * ts,
                             [[maybe_unused]] const int is_family,
                             const int num_elements,
                             t8_element_t * elements[])
{
    AdaptDataInitialEmbeddedMesh* adapt_data = static_cast<AdaptDataInitialEmbeddedMesh*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Check if the element is already on the initial refinement level (or if it is still coarser) */
    if (ts->element_get_level(tree_class, elements[0]) >= adapt_data->initial_refinement_level)
    {
        /* If the element's level is already on the initial refinement level the refinement process stops */
        return 0;
    }

    /* If the element is inside the global domain, it will be refined until the intial refinement level is reached */
    if (IsMeshElementWithinGlobalDomain(tree_class, elements[0], ts, adapt_data->global_domain, adapt_data->initial_refinement_level, adapt_data->initial_layout))
    {
        return 1;
    } else
    {
        return 0;
    }
}

inline t8_eclass_t
EmbeddedMeshDimensionToElementClass(const int dimensionality)
{
    cmc_assert(dimensionality == 2 || dimensionality == 3);
    return (dimensionality == 3 ? T8_ECLASS_HEX : T8_ECLASS_QUAD);
}

inline int
CalculateInitialEmbeddedRefinementLevel(const GeoDomain& global_domain)
{
    const size_t max_elem_per_direction = global_domain.GetLargestDimensionLength();
    /* Calculate the induced initial refinement level needed in order to build an enclosing mesh */
    return static_cast<int>(std::ceil(std::log2(max_elem_per_direction) + std::numeric_limits<double>::epsilon()));
}

inline 
void ValidateInitialRefinementLevelForDimensionality([[maybe_unused]] const int initial_refinement_level, [[maybe_unused]] const int dimensionality)
{
    #ifdef CMC_ENABLE_DEBUG
    cmc_assert(dimensionality == 2 || dimensionality == 3);
    
    const int maximum_possible_refinement_level = (dimensionality == 2 ? P4EST_MAXLEVEL : P8EST_MAXLEVEL);
    
    if (initial_refinement_level < 0 ||
        initial_refinement_level > maximum_possible_refinement_level)
    {
        cmc_err_msg("The corresponding refinement level is not within the range of an computationally posiible refinement level.");
    }

    #endif
}

inline std::tuple<t8_forest_t, int, int>
BuildInitialEmbeddedMesh(const GeoDomain& domain, const DataLayout initial_layout, MPI_Comm comm)
{
    cmc_debug_msg("The intial embedded mesh will be constructed");
    /* Get the dimensionality of the domain on which the variable is defined */
    const int dimensionality = domain.GetDimensionality();
    cmc_debug_msg("Comm is: ", comm);
    /* Create the cmesh; either a quad or hex tree based on the dimension */
    t8_cmesh_t cmesh = t8_cmesh_new_hypercube(EmbeddedMeshDimensionToElementClass(dimensionality), comm, 0, 0, 1);

    /* Determine the initial refinement level for embedding the data */
    const int initial_refinement_level = CalculateInitialEmbeddedRefinementLevel(domain);
    ValidateInitialRefinementLevelForDimensionality(initial_refinement_level, dimensionality);

    /* Construct a forest from the cmesh */
    t8_forest_t initial_forest;
    t8_forest_init(&initial_forest);
    t8_forest_set_cmesh(initial_forest, cmesh, comm);
    t8_forest_set_scheme(initial_forest, t8_scheme_new_default());
    t8_forest_set_level(initial_forest, 0);
    t8_forest_commit(initial_forest);

    AdaptDataInitialEmbeddedMesh adapt_data(domain, initial_refinement_level, initial_layout);

    t8_forest_t adapted_forest;
    for (int adaptation_steps = 0; adaptation_steps <= initial_refinement_level; ++adaptation_steps)
    {
        t8_forest_init(&adapted_forest);
        t8_forest_set_adapt(adapted_forest, initial_forest, RefineToInitialEmbeddedMesh, 0);
        const int set_partition_for_coarsening = 0; //TODO change to one later
        t8_forest_set_partition(adapted_forest, NULL, set_partition_for_coarsening);
        t8_forest_set_user_data(adapted_forest, static_cast<void*>(&adapt_data));
        t8_forest_commit(adapted_forest);

        initial_forest = adapted_forest;
    }

    cmc_debug_msg("The intial embedded mesh (of dimensionality ", dimensionality, " and with a maximum initial refinement level of ", initial_refinement_level, ") has been constructed");

    return std::make_tuple(initial_forest, initial_refinement_level, dimensionality); 
}



}




#endif /* !CMC_EMBEDDED_MESH_UTILITIES_HXX */
