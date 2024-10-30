#include "t8code/cmc_t8_adapt_callbacks.h"
#include "utilities/cmc_utilities.hxx"
#include "t8code/cmc_t8_adapt.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_data_variables.hxx"

#include "utilities/cmc_log_functions.hxx"

#ifdef CMC_WITH_T8CODE
namespace cmc
{

constexpr bool kCoarsenOutsideOfDomainBoundary = true;

static bool
IsCoarsenedFromTheBeginningOnwards(t8_eclass_scheme_c * ts, const t8_element_t * element, const int initial_refinement_level, const int count_adaptation_step)
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
            if (!IsMeshElementWithinGlobalDomain(elements[j], ts, compression_variable.GetGlobalDomain(), adapt_data->GetInitialRefinementLevelOfMesh(), compression_variable.GetInitialDataLayout()))
            {
                compression_variable.LeaveElementUnchanged(lelement_id);
                return kLeaveElementUnchanged;
            }
        }
    }
    std::vector<const t8_element_t*> const_elements;
    const_elements.reserve(num_elements);
    std::copy_n(elements, num_elements, std::back_inserter(const_elements));

    const std::vector<PermittedError> permitted_errors = compression_variable.GetPermittedError(num_elements, const_elements.data(), ts);
    
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
    if (IsMeshElementWithinGlobalDomain(elements[0], ts, adapt_data->global_domain, adapt_data->initial_refinement_level, adapt_data->initial_layout))
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
ExtractCommonPrefixes (t8_forest_t forest,
                       [[maybe_unused]] t8_forest_t forest_from,
                       [[maybe_unused]] int which_tree,
                       int lelement_id,
                       t8_eclass_scheme_c * ts,
                       const int is_family,
                       const int num_elements,
                       t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    /* Get the adapt data from the forest */
    PrefixAdaptData* adapt_data = static_cast<PrefixAdaptData*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* A prefix can only be extarcted from a family of elements. If a single element is passed to the adaptation funciton,
     * it is stored as prefix and has to wait until it's siblings are present within the mesh in order to then perform
     * a common prefix evaluation and extraction. */
    if (is_family == 0)
    {
        /* Since, only a single element has been passed, we are dragging it's value along and try to coarsen it's corresponding family later */
        adapt_data->LeaveCoarsePrefixUnchanged(lelement_id);

        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return kLeaveElementUnchanged;
    } else
    {
        //std::vector<int> non_missig_val_element_ids;
        //non_missig_val_element_ids.reserve(num_elements);
        //
        //const int end_elem_id = lelement_id + num_elements;
        //
        //for (int elem_id = lelement_id; elem_id < end_elem_id; ++elem_id)
        //{
        //    //cmc_debug_msg("Elem Id: ", elem_id, ", ref lvl: ", adapt_data->GetInitialRefinementLevelOfMesh(), ", layout: ", adapt_data->GetInitialDataLayout());
        //    if (IsMeshElementWithinGlobalDomain(elements[elem_id - lelement_id], ts, adapt_data->GetGlobalDomain(), adapt_data->GetInitialRefinementLevelOfMesh(), adapt_data->GetInitialDataLayout()))
        //    {
        //        non_missig_val_element_ids.push_back(elem_id);
        //    }
        //}
        //adapt_data->ExtractCommonPrefix(non_missig_val_element_ids);
        
        /* Extract a common prefix of the family */
        adapt_data->ExtractCommonPrefix(lelement_id, num_elements);

        /* Therefore, we need to coarsen the mesh at this point */
        return kCoarsenElements;
    }

    #else
    return CMC_ERR;
    #endif
}

t8_locidx_t
DecompressPrefixEncoding (t8_forest_t forest,
                          t8_forest_t forest_from,
                          int which_tree,
                          int lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements,
                          t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    /* Get the adapt data from the forest */
    DecompressPrefixAdaptData* adapt_data = static_cast<DecompressPrefixAdaptData*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    bool refine_element{false};

    /* Check if a prefix is given for this element */
    if (adapt_data->IsPrefixGiven())
    {
        /* Get the next prefix and it's length */
        auto [prefix, prefix_length] = adapt_data->GetNextPrefix();

        //cmc_debug_msg("elem id: ", lelement_id, " mit prefix laenge: ", prefix_length);

        /* Check if a refinement is given as well and append the prefix accordingly */
        if (adapt_data->IsRefinementGiven())
        {
            adapt_data->ApplyPrefixAndRefine(lelement_id, prefix, prefix_length);
            refine_element = true;
        } else
        {
            adapt_data->ApplyPrefixAndLeaveElementUnchanged(lelement_id, prefix, prefix_length);
        }
    } else
    {
        /* If no prefix is given, we check whether the element will be refined */
        if (adapt_data->IsRefinementGiven())
        {
            adapt_data->Refine(lelement_id);
            refine_element = true;
        } else
        {
            adapt_data->LeaveElementUnchanged(lelement_id);
        }
    }


    return refine_element;

    #else
    return CMC_ERR;
    #endif
}

t8_locidx_t
DecompressSuffixEncoding (t8_forest_t forest,
                          t8_forest_t forest_from,
                          int which_tree,
                          int lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements,
                          t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    /* Get the adapt data from the forest */
    DecompressPrefixAdaptData* adapt_data = static_cast<DecompressPrefixAdaptData*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Check if a prefix is given for this element */
    if (adapt_data->IsPrefixGiven())
    {
        /* Get the next prefix and it's length */
        auto [prefix, prefix_length] = adapt_data->GetNextPrefix();

        //cmc_debug_msg("suff: elem id: ", lelement_id, " mit prefix laenge: ", prefix_length); 
        adapt_data->ApplyPrefixAndLeaveElementUnchanged(lelement_id, prefix, prefix_length);
        
    } else
    {
        adapt_data->LeaveElementUnchanged(lelement_id);
    }


    return 0;

    #else
    return CMC_ERR;
    #endif
}

}
#endif
