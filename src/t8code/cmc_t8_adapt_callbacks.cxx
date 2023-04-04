#include "t8code/cmc_t8_adapt_callbacks.h"
#include "t8code/cmc_t8code_data.hxx"
#include "cmc_t8code_geo_mesh.h"

#ifdef CMC_WITH_T8CODE

/** Static function blocks which may be used for the creation of adaptation criteria **/

/**
 * @brief Function determining whether the element (respectively the family of elements) lies completly inside the geo mesh
 * 
 * @param adapt_data The adapt_data from the forest which will be adapted
 * @param num_elements The amount of elements passed to the current call of the adapt fuction
 * @param elements Pointer to the corresponding element(s of the family) 
 * @return true If all elements are inside the geo mesh
 * @return false If at least one element does not reside inside the geo mesh
 */
static inline
bool cmc_t8_elements_all_inside_geo_mesh(cmc_t8_adapt_data_t adapt_data, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    /* Set the current_var_id based on the compression mode (in a 'One for All'-compression, current var_id is smaller thn zero, otherwise in a 'One for One'-compression it resembles the current intern id of the variable) */
    const int var_id = (adapt_data->current_var_id < 0 ? 0 : adapt_data->current_var_id);

    /* Check the location of the element */
    for(int elem_id{0}; elem_id < num_elements; ++elem_id)
    {
        if (cmc_t8_elem_inside_geo_mesh(elements[elem_id], ts, *(adapt_data->t8_data), var_id) == 0)
        {
            /* If at least one element does not lay within the geo mesh */
            return false;
        }
    }
    /* If all of the elements lie within the "lat x lon x lev"-mesh */
    return true;
    #endif
}

/**
 * @brief Function determining whether the element (respectively the family of elements) lies completly outside of the geo mesh
 * 
 * @param adapt_data The adapt_data from the forest which will be adapted
 * @param num_elements The amount of elements passed to the current call of the adapt fuction
 * @param elements Pointer to the corresponding element(s of the family) 
 * @return true If all elements are outside of the geo mesh
 * @return false If at least one element is inside of the geo mesh
 */
static inline
bool cmc_t8_elements_all_outside_geo_mesh(cmc_t8_adapt_data_t adapt_data, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    /* Set the current_var_id based on the compression mode (in a 'One for All'-compression, current var_id is smaller thn zero, otherwise in a 'One for One'-compression it resembles the current intern id of the variable) */
    const int var_id = (adapt_data->current_var_id < 0 ? 0 : adapt_data->current_var_id);

    /* Check the location of the element */
    for(int elem_id{0}; elem_id < num_elements; ++elem_id)
    {
        if (cmc_t8_elem_inside_geo_mesh(elements[elem_id], ts, *(adapt_data->t8_data), var_id) != 0)
        {
            /* If at least one element of the family lies within the "lat x lon x lev"-mesh */
            return false;
        }
    }
    /* If all of the elements lay outside of the "lat x lon x lev"-mesh */
    return true;
    #endif
}

/**
 * @brief Function determining whether the element (respectively the family of elements) lies/lay completly inside the geo domain specified by the 'exclude_area'-compression-criterion
 * 
 * @param adapt_data The adapt_data from the forest which will be adapted
 * @param num_elements The amount of elements passed to the current call of the adapt fuction
 * @param elements Pointer to the corresponding element(s of the family) 
 * @return true If all elements are inside the specified domain 
 * @return false If at leat one of the elements is outside of the specidifed domain
 */
static inline
bool cmc_t8_elements_all_inside_geo_domain(cmc_t8_adapt_data_t adapt_data, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    /* Set the current_var_id based on the compression mode (in a 'One for All'-compression, current var_id is smaller thn zero, otherwise in a 'One for One'-compression it resembles the current intern id of the variable) */
    const int var_id = (adapt_data->current_var_id < 0 ? 0 : adapt_data->current_var_id);

    /* Check if at least one element of the family lies within the excluded geo-spatial domain */
    for (int elem_id{0}; elem_id < num_elements; ++elem_id)
    {
        /* Check if the element is inside the given domain (which will be excluded by the compression) */
        if (cmc_t8_elem_inside_geo_domain(elements[elem_id], ts, *(adapt_data->t8_data), var_id, adapt_data->t8_data->settings.exclude_area_start_indices[CMC_COORD_IDS::CMC_LAT], adapt_data->t8_data->settings.exclude_area_end_indices[CMC_COORD_IDS::CMC_LAT],
                                          adapt_data->t8_data->settings.exclude_area_start_indices[CMC_COORD_IDS::CMC_LON], adapt_data->t8_data->settings.exclude_area_end_indices[CMC_COORD_IDS::CMC_LON],
                                          adapt_data->t8_data->settings.exclude_area_start_indices[CMC_COORD_IDS::CMC_LEV], adapt_data->t8_data->settings.exclude_area_end_indices[CMC_COORD_IDS::CMC_LEV]) == 0)
        {
            /* If at least one element lies outside the given domain */
            return false;
        }
    }

    /* If all elements lay inside of the domain */
    return true;
    #endif
}

/**
 * @brief Function determining whether the element (respectively the family of elements) lies/lay completly outside the geo domain specified by the 'exclude_area'-compression-criterion
 * 
 * @param adapt_data The adapt_data from the forest which will be adapted
 * @param num_elements The amount of elements passed to the current call of the adapt fuction
 * @param elements Pointer to the corresponding element(s of the family) 
 * @return true If all elements are outside the specified domain 
 * @return false If at leat one of the elements is inside of the specidifed domain
 */
static inline
bool cmc_t8_elements_all_outside_geo_domain(cmc_t8_adapt_data_t adapt_data, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE

    /* Set the current_var_id based on the compression mode (in a 'One for All'-compression, current var_id is smaller thn zero, otherwise in a 'One for One'-compression it resembles the current intern id of the variable) */
    const int var_id = (adapt_data->current_var_id < 0 ? 0 : adapt_data->current_var_id);

    /* Check if at least one element of the family lies within the excluded geo-spatial domain */
    for (int elem_id{0}; elem_id < num_elements; ++elem_id)
    {
        /* Check if the element is inside the given domain (which will be excluded by the compression) */
        if (cmc_t8_elem_inside_geo_domain(elements[elem_id], ts, *(adapt_data->t8_data), var_id, adapt_data->t8_data->settings.exclude_area_start_indices[CMC_COORD_IDS::CMC_LAT], adapt_data->t8_data->settings.exclude_area_end_indices[CMC_COORD_IDS::CMC_LAT],
                                          adapt_data->t8_data->settings.exclude_area_start_indices[CMC_COORD_IDS::CMC_LON], adapt_data->t8_data->settings.exclude_area_end_indices[CMC_COORD_IDS::CMC_LON],
                                          adapt_data->t8_data->settings.exclude_area_start_indices[CMC_COORD_IDS::CMC_LEV], adapt_data->t8_data->settings.exclude_area_end_indices[CMC_COORD_IDS::CMC_LEV]) == 1)
        {
            /* If at least one element lies inside the given domain */
            return false;
        }
    }

    /* If all elements lay outside of the domain */
    return true;
    #endif
}

static inline
bool cmc_t8_element_level_coarser_than_initial_refinement_level(cmc_t8_adapt_data_t adapt_data, t8_eclass_scheme_c * ts, t8_element_t* element)
{
    #ifdef CMC_WITH_T8CODE
    const int initial_ref_lvl{(adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
                               adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D) ?
                               adapt_data->t8_data->assets->initial_refinement_lvl :
                               adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl};
    
    /* Check if the eleemnt's current level is smaller than the initial refinement level */
    if (t8_element_level(ts, element) < initial_ref_lvl - adapt_data->adapt_step)
    {
        /* The element's level is smaller/coarser */
        return true;
    } else
    {
        /* The element's level is the same or finer than the initial refinement level */
        return false;
    } 
    #endif
} 

static inline
bool cmc_t8_elements_coarsened_from_the_beginning_onwards(cmc_t8_adapt_data_t adapt_data, t8_eclass_scheme_c * ts, t8_element_t* element)
{
    #ifdef CMC_WITH_T8CODE

    const int initial_ref_lvl{(adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
                               adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D) ?
                               adapt_data->t8_data->assets->initial_refinement_lvl :
                               adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl};
    if (t8_element_level(ts, element) != initial_ref_lvl - adapt_data->adapt_step)
    {
        /* The element has not been coarsened before */
        return false;
    } else
    {
        /* The element has been coarsened from the start onwards (especially during the last adaptation) */
        return true;
    }
    #endif
}


static inline
void cmc_t8_update_counters_for_rel_error_threshold(cmc_t8_adapt_data_t adapt_data)
{
    #ifdef CMC_WITH_T8CODE
    ++(adapt_data->_counter);
    ++(adapt_data->_counter_nxt_lvl);
    #endif
}

/**
 * @brief Function determining whether the element (respectively the family of elements) lies completly inside the geo domain specified by the 'exclude_area'-compression-criterion
 * 
 * @param adapt_data The adapt_data from the forest which will be adapted
 * @param num_elements The amount of elements passed to the current call of the adapt fuction
 * @param elements Pointer to the corresponding element(s of the family) 
 * @return true If the elements fulfill the relative error threshold and could be coarsened
 * @return false If the elements does not fulfill the error threshold and therefore cannot be coarsened 
 */
static inline
bool cmc_t8_elements_fulfill_rel_error_threshold(cmc_t8_adapt_data_t adapt_data, t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    
    /* It is strictly sepereated between the compression modes */
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

                /* The elements does not fulfill the threshold */
                return false;
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
        return true;
    }
    else {
        /* In case of a 'One for One' compression mode is used */
        /* Check the element families mean value, and calculate the deviation introduced by this coarsening */
        /* If it complies with the error boundary, coarsening will be performed */
        const size_t num_elems{static_cast<size_t>(pow(t8_element_num_children(ts, elements[0]), (adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl + 1 - t8_element_level(ts, elements[0]))))};
        const size_t starting_elem_id{static_cast<size_t>((t8_element_level(ts, elements[0]) == adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl) ?
                                    lelement_id : (adapt_data->initial_ref_lvl_ids[adapt_data->current_var_id])[adapt_data->_counter])};

        /* Calculate the mean value of all initial data points (in the coarse element's area) except missing_values */
        //TODO: Is it possible to have missing values here? Or is this case excluded by the checks beforehand if all elements are inside the geo mesh (guess depends on the variable - from netCDf or from MESSy for example)
        cmc_universal_type_t mean = adapt_data->t8_data->initial_data[adapt_data->current_var_id].mean_value_same_type_wo_missing_vals(starting_elem_id, starting_elem_id + num_elems -1, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);

        /* Check if the given error threshold would be fulfilled by the coarsening */
        if (!(adapt_data->t8_data->initial_data[adapt_data->current_var_id].check_range_fullfills_deviation_threshold_from_value(starting_elem_id, starting_elem_id + num_elems -1, adapt_data->t8_data->settings.max_err, mean, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value)))
        {
            /* If the threshold is not fullfilled */
            ++(adapt_data->_counter);
            ++(adapt_data->_counter_nxt_lvl);
            /* The elements do not fulfill the relative error threshold */
            return false;
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
        return true;
    }
    #endif
}


/** End of staticfunction building blocks **/

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

    /* Only a family of elements can be coarsened.
     * If only a single element is passed to the callback, the callback function can be returned immediately */
    if (is_family == 0)
    {
        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    }

    cmc_t8_adapt_data_t adapt_data = static_cast<cmc_t8_adapt_data_t>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != NULL);
    cmc_assert(ts != NULL);
    
    /* Check if all elements are outisde of the geo mesh */
    if (cmc_t8_elements_all_outside_geo_mesh(adapt_data, ts, num_elements, elements))
    {
        /* If all elements are outside of the geo mesh, we are able to coarse them in order to obtain our initial mesh for the compression */
        return -1;
    } else
    {
        /* If not all elements are outside the geo mesh, the elements stay the same */
        return 0;
    }

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
    
    if (!cmc_t8_elements_all_inside_geo_mesh(adapt_data, ts, 1, elements))
    {
        /* The element is not inside the geo mesh, therefore we cannot refine it */
        return 0;
    } 

    /* Check if the element is already on the initial refinement level (or if it is still coarser) */
    if(cmc_t8_element_level_coarser_than_initial_refinement_level(adapt_data, ts, elements[0]))
    {
        /* If the element's level is still coarser than the initial refinement level, we can refine the element */
        return 1;
    } else
    {
        /* If we have already reached the initial refinement level, the element stays the same */
        return 0;
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
    }

    /* Get the adapt data from the forest */
    cmc_t8_adapt_data_t adapt_data = static_cast<cmc_t8_adapt_data_t>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != NULL);
    
    /* Check if the family of elements is inside the geo mesh */
    if (!cmc_t8_elements_all_inside_geo_mesh(adapt_data, ts, num_elements, elements))
    {
        /* If not all elements of the family are inside the geo mesh, the elements stay the same and we can return immediately */
        return 0;
    }

    /* Check if the family of elements lay outside the specified domain */
    if (cmc_t8_elements_all_outside_geo_domain(adapt_data, ts, num_elements, elements))
    {
        /* If all elements are outside of the domain, we can coarsen them */
        return -1;
    } else
    {
        /* If some elements of the family lay within the domain, we cannot coarsen them (the family of elements stays the same) */
        return 0;
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

    /** Only a family of elements can be coarsened.
     * If only a single element is passed to the callback, the callback function can be returned immediately. 
    */
    if (is_family == 0)
    {
        /* Due to the nature of the relative error criterion, we need to keep track which elements are coarsened.
         * Therefore, when this criterion is used, we need to update the counters in the 'adapt_data' every time we prematurely return without coarsening
         */
        cmc_t8_update_counters_for_rel_error_threshold(adapt_data);

        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    }

    /** Furthermore, if the elements are not inside the geo mesh, we cannot coarsen them. 
     * And elements (resp. family of elements) that have not been coarsened before (resp. from the start of coarsening onwards) cannot be coarsened
     * since this would discard artefacts/high variability of the data in a small scale (which were initially present). So once we have decided to no coarsen the elements, the decision is absolute
     * (This keeps 'details' in the data)
     */
    if (!cmc_t8_elements_coarsened_from_the_beginning_onwards(adapt_data, ts, elements[0]) ||
        !cmc_t8_elements_all_inside_geo_mesh(adapt_data, ts, num_elements, elements))
    {
        /* Update the counters within the adapt data */
        cmc_t8_update_counters_for_rel_error_threshold(adapt_data);

        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    }

    /* Families of elements that have passed the previous checks are now conducted if they fulfill the given relative error threshold */
    if (cmc_t8_elements_fulfill_rel_error_threshold(adapt_data, lelement_id, ts, num_elements, elements))
    {
        /* If the elements family fulfills the relative error threshold, we are able to coarsen the elements */
        return -1;
    } else
    {
        /* If the family of elements do not fulfill the relative error threshold, the element stays the same */
        return 0;
    }

    #else
    return CMC_ERR;
    #endif
}


t8_locidx_t
cmc_t8_adapt_callback_coarsen_combined_criteria (t8_forest_t forest,
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

    /** Only a family of elements can be coarsened.
     * If only a single element is passed to the callback, the callback function can be returned immediately. 
    */
    if (is_family == 0)
    {
        /* Due to the nature of the relative error criterion, we need to keep track which elements are coarsened.
         * Therefore, when this criterion is used, we need to update the counters in the 'adapt_data' every time we prematurely return without coarsening
         */
        cmc_t8_update_counters_for_rel_error_threshold(adapt_data);
        
        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    }

    /** Furthermore, if the elements are not inside the geo mesh, we cannot coarsen them. 
     * And elements (resp. family of elements) that have not been coarsened before (resp. from the start of coarsening onwards) cannot be coarsened
     * since this would discard artefacts/high variability of the data in a small scale (which were initially present). So once we have decided to no coarsen the elements, the decision is absolute
     * (This keeps 'details' in the data).
     * Moreover the exclude-area criterion is checked.
     */
    if (!cmc_t8_elements_coarsened_from_the_beginning_onwards(adapt_data, ts, elements[0]) ||
        !cmc_t8_elements_all_inside_geo_mesh(adapt_data, ts, num_elements, elements) ||
        !cmc_t8_elements_all_outside_geo_domain(adapt_data, ts, num_elements, elements))
    {
        /* Update the counters within the adapt data */
        cmc_t8_update_counters_for_rel_error_threshold(adapt_data);

        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    }

    /* Families of elements that have passed the previous checks are now conducted if they fulfill the given relative error threshold */
    if (cmc_t8_elements_fulfill_rel_error_threshold(adapt_data, lelement_id, ts, num_elements, elements))
    {
        /* If the elements family fulfills the relative error threshold, we are able to coarsen the elements */
        return -1;
    } else
    {
        /* If the family of elements do not fulfill the relative error threshold, the element stays the same */
        return 0;
    }

    #else
    return CMC_ERR;
    #endif
}

#endif
