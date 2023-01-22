#include "t8code/cmc_t8_adapt_callbacks.h"
#include "t8code/cmc_t8code_data.hxx"
#include "cmc_t8code_geo_mesh.h"

#ifdef CMC_WITH_T8CODE

t8_locidx_t
cmc_t8_adapt_coarsen_geo_mesh_callback (t8_forest_t forest,
                                        t8_forest_t forest_from,
                                        t8_locidx_t which_tree,
                                        t8_locidx_t lelement_id,
                                        t8_eclass_scheme_c * ts,
                                        const int is_family,
                                        const int num_elements,
                                        t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    cmc_t8_adapt_data_t adapt_data = static_cast<cmc_t8_adapt_data_t>(t8_forest_get_user_data(forest));
    
    cmc_assert(adapt_data != NULL);
    cmc_assert(ts != NULL);
    
    /* Check if a family of elements is considered */
    if (is_family)
    {
        /* Check the location of the element */
        for(int elem_id{0}; elem_id < num_elements; ++elem_id)
        {
            if (cmc_t8_elem_inside_geo_mesh(elements[elem_id], ts, *(adapt_data->t8_data), adapt_data->current_var_id) != 0)
            {
                /* The family of this element cannot be coarsened since this element lies within the "lat x lon x lev"-mesh */
                return 0;
            }
        }
        /* If none of the elements lies within the "lat x lon x lev"-mesh, this family of elements can be coarsened */
        return -1;
    }
    /* Otherwise if no family of elements is considerd, but just a single element, it will remain unchanged within the mesh */
    return 0;

    #else
    return CMC_ERR;
    #endif
}

t8_locidx_t
cmc_t8_adapt_callback_refine_to_initial_lvl (t8_forest_t forest,
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
    cmc_t8_adapt_data_t adapt_data = static_cast<cmc_t8_adapt_data_t>(t8_forest_get_user_data(forest));
    
    /* Check which compression mode is used */
    if (adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
        adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
    {
        /* If a 'One for All' compression mode is chosen */
        /* If the element is inside the domain of the geo mesh and its current refinement level is smaller than the initial refinement level, it will be refined */
        if (cmc_t8_elem_inside_geo_mesh(elements[0], ts, *(adapt_data->t8_data), 0) != 0 &&
            t8_element_level(ts, elements[0]) < adapt_data->t8_data->assets->initial_refinement_lvl)
        {
            /* Refine the element */
            return 1;
        } else
        {
            /* Otherwise the element stays the same */
            return 0;
        }
    }
    else
    {
        /* If a 'One for One' compression mode is chosen */
        /* If the element is inside the domain of the geo mesh and its current refinement level is smaller than the initial refinement level, it will be refined */
        if (cmc_t8_elem_inside_geo_mesh(elements[0], ts, *(adapt_data->t8_data), adapt_data->current_var_id) != 0 &&
            t8_element_level(ts, elements[0]) < adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl)
        {
            /* Refine the element */
            return 1;
        } else
        {
            /* Otherwise the element stays the same */
            return 0;
        }
    }
    #endif
}

t8_locidx_t
cmc_t8_adapt_callback_coarsen_exclude_area(t8_forest_t forest,
                                           t8_forest_t forest_from,
                                           int which_tree,
                                           int lelement_id,
                                           t8_eclass_scheme_c * ts,
                                           const int is_family,
                                           const int num_elements,
                                           t8_element_t * elements[]) 
{
    #ifdef CMC_WITH_T8CODE
    /* Only a family of elements can be coarsened.
     * If only a single element is passed to the callback, the callback function can be returned immediately */
    if (is_family == 0)
    {
        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    } else
    {
        /* Get the adapt data from the forest */
        cmc_t8_adapt_data_t adapt_data = static_cast<cmc_t8_adapt_data_t>(t8_forest_get_user_data(forest));

        cmc_assert(adapt_data != NULL);

        /* Only perform two coarsening steps */
        if (adapt_data->adapt_step >= 1)
        {
            return 0;
        }

        /* Check if at least one element of the family lies within the defined geo-spatial domain */
        for (int elem_id{0}; elem_id < num_elements; ++elem_id)
        {
            /* Example check for the domain of South-America */
            if (cmc_t8_elem_inside_geo_domain(elements[elem_id], ts, *(adapt_data->t8_data), adapt_data->current_var_id, 13, 45, 112, 132, 0, 0) == 1)
            {
                /* Do not coarsen elements which lay inside the given domain */
                return 0;
            }
        }

        /* If all checks have passed, coarsen the elements' family */
        return -1;
    }
    #else
    return CMC_ERR;
    #endif
}

t8_locidx_t
cmc_t8_adapt_callback_coarsen_error_threshold (t8_forest_t forest,
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
    cmc_t8_adapt_data_t adapt_data = static_cast<cmc_t8_adapt_data_t>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Only a family of elements can be coarsened.
     * If only a single element is passed to the callback, the callback function can be returned immediately */
    if (is_family == 0)
    {
        ++(adapt_data->_counter);
        ++(adapt_data->_counter_nxt_lvl);
        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    }
    /* Only try to coarsen elements which has been coarsened in the step before (Other is not possible and this check prevents misleading calculations that it is possible nonetheless */
    const int initial_ref_lvl{(adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
                               adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D) ?
                               adapt_data->t8_data->assets->initial_refinement_lvl :
                               adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl};
    if (t8_element_level(ts, elements[0]) != initial_ref_lvl - adapt_data->adapt_step)
    {
        ++(adapt_data->_counter);
        ++(adapt_data->_counter_nxt_lvl);
        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    }
    /* Do not coarsen overlapping boundary elements / element-families) */
    for (int elem_id{0}; elem_id < num_elements; ++elem_id)
    {
        /* Check if the element is inside the (lat x lon x lev) mesh */
        /* Since all vars are defined on the same domain, just pass the var_id of the first variable to the function */
        if (cmc_t8_elem_inside_geo_mesh(elements[elem_id], ts, *(adapt_data->t8_data), 0) == 0)
        {
            ++(adapt_data->_counter);
            ++(adapt_data->_counter_nxt_lvl);
            /* The element stays the same, it will not be refined and it cannot be coarsened */
            return 0;
        }
    }
    /* Now it is strictly sepereated between the compression modes */
    /* This is due to the fact, that in a 'One for All' compression mode, the threshold needs to be fullfilled for every varibale
     * whereas in a 'One for One' compression mode the threshold only needs to hold for a single variable itself */
    /* Check which mode is used */
    if (adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
        adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
    {
        /* In case a 'One for All' compression mode is used */
        /* Check the element families mean value, and calculate the deviation introduced by this coarsening */
        /* If it complies with the error boundary, coarsening will be performed */
        const size_t num_elems{static_cast<size_t>(pow(t8_element_num_children(ts, elements[0]), (adapt_data->t8_data->assets->initial_refinement_lvl + 1 - t8_element_level(ts, elements[0]))))};
        const size_t starting_elem_id{static_cast<size_t>((t8_element_level(ts, elements[0]) == adapt_data->t8_data->assets->initial_refinement_lvl) ?
                                    lelement_id : (adapt_data->initial_ref_lvl_ids[0])[adapt_data->_counter])};

        /* Save the mean value of the variable */
        std::vector<cmc_universal_type_t> mean;
        mean.reserve(adapt_data->t8_data->vars.size());

        /* Iterte over variables and check if the mean fullfills the deviation defined by the max error */
        for (size_t var_id{0}; var_id < adapt_data->t8_data->vars.size(); ++var_id)
        {
            /* Calculate the mean value */
            mean.push_back((adapt_data->t8_data->initial_data[var_id]).mean_value_same_type_wo_missing_vals(starting_elem_id, starting_elem_id + num_elems -1, adapt_data->t8_data->vars[var_id]->var->missing_value));
            /* Check if the given error threshold would be fullfilled by the coarsening */
            if (!(adapt_data->t8_data->initial_data[var_id].check_range_fullfills_deviation_threshold_from_value(starting_elem_id, starting_elem_id + num_elems -1, adapt_data->t8_data->settings.max_err, mean[var_id], adapt_data->t8_data->vars[var_id]->var->missing_value)))
            {
                /* If the threshold is not fullfilled */
                ++(adapt_data->_counter);
                ++(adapt_data->_counter_nxt_lvl);

                /* Do not coarsen the element's family */
                return 0;
            }
        }

        /* If all variables fullfill the threshold, we are able to coarsen the family */
        /* Save the first element id of the range/domain the coarse element covers (in dependecy of the initial data/mesh) */
        /* Due to the SFC all children lay next to each other in memory */ 
        if (t8_element_level(ts, elements[0]) == adapt_data->t8_data->geo_data->initial_refinement_lvl)
        {
            (adapt_data->initial_ref_lvl_ids[0])[adapt_data->_counter_nxt_lvl] = lelement_id;
        } else
        {
            (adapt_data->initial_ref_lvl_ids[0])[adapt_data->_counter_nxt_lvl] = (adapt_data->initial_ref_lvl_ids[0])[adapt_data->_counter];
        }

        /* Save the calculatd mean value for the interpolation afterwards */
        for (size_t var_id{0}; var_id < adapt_data->t8_data->vars.size(); ++var_id)
        {
            (*(adapt_data->adapted_data))[var_id].assign_value(adapt_data->coarsening_counter, mean[var_id]);
        }
        /* Increment the pointers */
        ++(adapt_data->coarsening_counter);
        adapt_data->_counter += num_elements;
        ++(adapt_data->_counter_nxt_lvl);

        /* Coarsen the element's family */
        return -1;
    }
    else {
        /* In case a 'One for One' compression mode is used */
        /* Check the element families mean value, and calculate the deviation introduced by this coarsening */
        /* If it complies with the error boundary, coarsening will be performed */
        const size_t num_elems{static_cast<size_t>(pow(t8_element_num_children(ts, elements[0]), (adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl + 1 - t8_element_level(ts, elements[0]))))};
        const size_t starting_elem_id{static_cast<size_t>((t8_element_level(ts, elements[0]) == adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl) ?
                                    lelement_id : (adapt_data->initial_ref_lvl_ids[adapt_data->current_var_id])[adapt_data->_counter])};

        /* Calculate the mean value of all initial data points (in the coarse element's area) except missing_values */
        //TODO: Is it possible to have missing values here? Or is this case excluded by the check above if all elements are inside the geo mesh (guess depends on the variable - from netCDf or from MESSy for example)
        cmc_universal_type_t mean = adapt_data->t8_data->initial_data[adapt_data->current_var_id].mean_value_same_type_wo_missing_vals(starting_elem_id, starting_elem_id + num_elems -1, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);

        /* Check if the given error threshold would be fullfilled by the coarsening */
        if (!(adapt_data->t8_data->initial_data[adapt_data->current_var_id].check_range_fullfills_deviation_threshold_from_value(starting_elem_id, starting_elem_id + num_elems -1, adapt_data->t8_data->settings.max_err, mean, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value)))
        {
            /* If the threshold is not fullfilled */
            ++(adapt_data->_counter);
            ++(adapt_data->_counter_nxt_lvl);
            /* Do not coarsen the element's family */
            return 0;
        }
        
        /* Save the id of the first element of the initial mesh element */
        /* Due to the SFC all children lay next to each other in memory */ 
        if (t8_element_level(ts, elements[0]) == adapt_data->t8_data->geo_data->initial_refinement_lvl)
        {
            (adapt_data->initial_ref_lvl_ids[adapt_data->current_var_id])[adapt_data->_counter_nxt_lvl] = lelement_id;
        } else
        {
            (adapt_data->initial_ref_lvl_ids[adapt_data->current_var_id])[adapt_data->_counter_nxt_lvl] = (adapt_data->initial_ref_lvl_ids[adapt_data->current_var_id])[adapt_data->_counter];
        }

        /* Save the calculatd mean value for the interpolation afterwards */
        (*(adapt_data->adapted_data))[0].assign_value(adapt_data->coarsening_counter, mean);

        /* Increment the pointers */
        ++(adapt_data->coarsening_counter);
        adapt_data->_counter += num_elements;
        ++(adapt_data->_counter_nxt_lvl);
        /* Coarsen the element's family */
        return -1;
    }
    
    #else
    return CMC_ERR;
    #endif
}

#endif
