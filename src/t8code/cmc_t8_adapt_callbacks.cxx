#include "t8code/cmc_t8_adapt_callbacks.h"
#include "utilities/cmc_utilities.hxx"
#include "t8code/cmc_t8_adapt.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_data_variables.hxx"

#ifdef CMC_WITH_T8CODE
namespace cmc
{

constexpr bool kCoarsenOutsideOfDomainBoundary = true;

static bool
IsCoarsenedFromTheBeginningOnwards(t8_eclass_scheme_c * ts, t8_element_t * element, const int initial_refinement_level, const int count_adaptation_step)
{
    /* Check if the current element level resembles a coarsening process that started at the beginning of the compression */
    if (t8_element_level(ts, element) != initial_refinement_level - count_adaptation_step)
    {
        /* The element has not been coarsened before */
        return false;
    } else
    {
        /* The element has been coarsened from the start onwards */
        return true;
    }
}


/**
 * @brief Function determining whether the element (respectively the family of elements) complies with the supplied relative error threshold
 * 
 * @param adapt_data The adapt_data from the forest which will be adapted
 * @param num_elements The amount of elements passed to the current call of the adapt fuction
 * @param elements Pointer to the corresponding element(s of the family) 
 * @return true If the elements fulfill the relative error threshold and could be coarsened
 * @return false If the elements does not fulfill the error threshold and therefore cannot be coarsened 
 */
t8_locidx_t
PerformAdaptiveCoarseningOneForOne (t8_forest_t forest,
                                    [[maybe_unused]] t8_forest_t forest_from,
                                    [[maybe_unused]] int which_tree,
                                    int lelement_id_,
                                    t8_eclass_scheme_c * ts,
                                    const int is_family,
                                    const int num_elements,
                                    t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    const int lelement_id = t8_forest_get_tree_element_offset(forest_from, which_tree) + lelement_id_;
    /* Get the adapt data from the forest */
    AdaptData* adapt_data = static_cast<AdaptData*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Get the current comrpession variable */
    Var& compression_variable = adapt_data->GetCurrentCompressionVariable();

    if (is_family == 0 ||
        !IsCoarsenedFromTheBeginningOnwards(ts, elements[0], adapt_data->GetInitialRefinementLevelOfMesh(), adapt_data->GetAdaptationStepCount()))
    {
        compression_variable.LeaveElementUnchanged(lelement_id);
        return kLeaveElementUnchanged;
    }

    if constexpr (!kCoarsenOutsideOfDomainBoundary)
    {
        for (int j = 0; j < num_elements; ++j)
        {
            if (!IsMeshElementWithinGeoDomain(elements[j], ts, compression_variable.GetGlobalDomain(), 12, compression_variable.GetInitialDataLayout()))
            {
                compression_variable.LeaveElementUnchanged(lelement_id);
                return kLeaveElementUnchanged;
            }
        }
    }

    const std::vector<PermittedError> permitted_errors = compression_variable.GetPermittedError(num_elements, elements, ts);
    
    const ErrorCompliance evaluation = compression_variable.EvaluateCoarsening(permitted_errors, forest, lelement_id, num_elements);

    if (evaluation.is_error_threshold_satisfied)
    {
        compression_variable.ApplyInterpolation(lelement_id, evaluation);
        return kCoarsenElements;
    } else
    {
        compression_variable.PopInterpolation();
        compression_variable.LeaveElementUnchanged(lelement_id);
        return kLeaveElementUnchanged;
    }

    #else
    return kLeaveElementUnchanged;
    #endif

}


/**
 * @brief Function determining whether the element (respectively the family of elements) complies with the supplied relative error threshold
 * 
 * @param adapt_data The adapt_data from the forest which will be adapted
 * @param num_elements The amount of elements passed to the current call of the adapt fuction
 * @param elements Pointer to the corresponding element(s of the family) 
 * @return true If the elements fulfill the relative error threshold and could be coarsened
 * @return false If the elements does not fulfill the error threshold and therefore cannot be coarsened 
 */
//TODO: implement
t8_locidx_t
PerformAdaptiveCoarseningOneForAll ([[maybe_unused]] t8_forest_t forest,
                                    [[maybe_unused]] t8_forest_t forest_from,
                                    [[maybe_unused]] int which_tree,
                                    [[maybe_unused]] int lelement_id,
                                    [[maybe_unused]] t8_eclass_scheme_c * ts,
                                    [[maybe_unused]] const int is_family,
                                    [[maybe_unused]] const int num_elements,
                                    [[maybe_unused]] t8_element_t * elements[])
{
    //TODO: implement
    return 0;
}



t8_locidx_t
RefineToInitialMesh (t8_forest_t forest,
                     [[maybe_unused]] t8_forest_t forest_from,
                     [[maybe_unused]] t8_locidx_t which_tree,
                     [[maybe_unused]] t8_locidx_t lelement_id,
                     t8_eclass_scheme_c * ts,
                     [[maybe_unused]] const int is_family,
                     [[maybe_unused]] const int num_elements,
                     t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE

    AdaptDataInitialMesh* adapt_data = static_cast<AdaptDataInitialMesh*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);
    
    /* Check if the element is already on the initial refinement level (or if it is still coarser) */
    if(t8_element_level(ts, elements[0]) >= adapt_data->initial_refinement_level)
    {
        /* If the element's level is already on the initial refinement level the refinement process stops */
        return 0;
    }

    /* If the element is inside the global domain, it will be refined until the intial refinement level is reached */
    if (IsMeshElementWithinGeoDomain(elements[0], ts, adapt_data->global_domain, adapt_data->initial_refinement_level, adapt_data->initial_layout))
    {
        return 1;
    } else
    {
        return 0;
    }

    #else
    return CMC_ERR;
    #endif
}



static void
PushRefinementBit(RefinementBits* rd, const bool refined)
{
    if (rd->current_bit_position == CHAR_BIT)
    {
        rd->current_bit_position = 0;
        rd->refinement_indicator.push_back(uint8_t{0});
    }
    if (refined)
    {
        /* Set a one for this element indicating that this element was refined further */
        rd->refinement_indicator.back() |= (uint8_t{1} << rd->current_bit_position);
    } else
    {
        /* Set a zero for this element indicating that ot wont be refined any further */
        rd->refinement_indicator.back() &= ~(uint8_t{1} << rd->current_bit_position);
    }
    ++(rd->current_bit_position);
}


t8_locidx_t
FindRefinementBits (t8_forest_t forest,
                    [[maybe_unused]] t8_forest_t forest_from,
                    [[maybe_unused]] int which_tree,
                    [[maybe_unused]] int lelement_id,
                    [[maybe_unused]] t8_eclass_scheme_c * ts,
                    const int is_family,
                    [[maybe_unused]] const int num_elements,
                    [[maybe_unused]] t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    /* Get the adapt data from the forest */
    RefinementBits* adapt_data = static_cast<RefinementBits*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /** Only a family of elements can be coarsened.
     * If only a single element is passed to the callback, the callback function can be returned immediately.
    */
    if (is_family == 0)
    {
        PushRefinementBit(adapt_data, false);
        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return kLeaveElementUnchanged;
    } else
    {
        PushRefinementBit(adapt_data, true);
        return kCoarsenElements;
    }

    #else
    return CMC_ERR;
    #endif
}


t8_locidx_t
ApplyRefinementBits (t8_forest_t forest,
                     [[maybe_unused]] t8_forest_t forest_from,
                     [[maybe_unused]] int which_tree,
                     [[maybe_unused]] int lelement_id,
                     [[maybe_unused]] t8_eclass_scheme_c * ts,
                     [[maybe_unused]] const int is_family,
                     [[maybe_unused]] const int num_elements,
                     [[maybe_unused]] t8_element_t * elements[])
{
    /* Get the adapt data from the forest */
    DecompressionRefinementBits* adapt_data = static_cast<DecompressionRefinementBits*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    if (adapt_data->bit_position == 8)
    {
        ++(adapt_data->byte_position);
        adapt_data->bit_position = 0;
    }

    const uint8_t current_byte = adapt_data->encoded_refinements[adapt_data->byte_position];

    /* Check if the currently considered bit is set */
    if ((current_byte >> adapt_data->bit_position) & (uint8_t) 1)
    {
        /* If it is set, the element will be refined */
        ++(adapt_data->bit_position);
        return kRefineElement;
    } else
    {
        /* If it is not set, the element remains unchanged */
        ++(adapt_data->bit_position);
        return kLeaveElementUnchanged;
    }
}

t8_locidx_t
FindPrefixBitsEGU (t8_forest_t forest,
                    [[maybe_unused]] t8_forest_t forest_from,
                    [[maybe_unused]] int which_tree,
                    int lelement_id,
                    [[maybe_unused]] t8_eclass_scheme_c * ts,
                    const int is_family,
                    const int num_elements,
                    [[maybe_unused]] t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    /* Get the adapt data from the forest */
    PrefixAdaptDataEGU* adapt_data = static_cast<PrefixAdaptDataEGU*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /** Only a family of elements can be coarsened.
     * If only a single element is passed to the callback, the callback function can be returned immediately.
    */
    if (is_family == 0)
    {
        adapt_data->LeaveCoarsePrefixUnchangedEGU(lelement_id);

        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return kLeaveElementUnchanged;
    } else
    {
        adapt_data->EvaluateCommonPrefixEGU(lelement_id, num_elements);
        /* If there is a common prefix */
        return kCoarsenElements;
    }

    #else
    return CMC_ERR;
    #endif
}


}
#endif
