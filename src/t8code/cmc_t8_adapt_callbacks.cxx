#include "t8code/cmc_t8_adapt_callbacks.h"
#include "t8code/cmc_t8code_data.hxx"
#include "cmc_t8code_geo_mesh.h"

#ifdef CMC_WITH_T8CODE
namespace cmc
{




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
cmc_compress_by_adaptive_coarsening (t8_forest_t forest,
                                     t8_forest_t forest_from,
                                     int which_tree,
                                     int lelement_id,
                                     t8_eclass_scheme_c * ts,
                                     const int is_family,
                                     const int num_elements,
                                     t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE

    /** Only a family of elements can be coarsened.
     * If only a single element is passed to the callback, the callback function can be returned immediately. 
    */
    if (is_family == 0)
    {
        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    }

    /* Get the adapt data from the forest */
    cmc_t8_adapt_data_t adapt_data = static_cast<cmc_t8_adapt_data_t>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    //Check geo-location

    //
    /* Check whether the error tolerance is still satisfied by the potential coarsening */
    const bool perform_coarsening = adapt_data->deviations_.CheckInaccuracy();

    if (perform_coarsening)
    {
        /* Coarsen the family of elements */
        return -1;
    } else
    {
        /* The element('s family) stays the same - no adaptation is performed */
        return 0;
    }
    #endif

}


















/* Static function building blocks */
/**
 * @brief This function calculates (an approximation) of the deviation up to the initial data and decides based on the supplied error threshold whether the element's family shall be coarsened or not
 * 
 * @param adapt_data A pointer to the current @struct cmc_t8_adapt_data holding the status and data of the adaptation
 * @param current_forest The current forest which is built from the adaptation
 * @param lelement_id the first local element id of it's family members
 * @param num_elements The number of elements within this family
 * @return true The element's family complies with the relative error threshold and shall be coarsened 
 * @return false The element's family does not comply with the relative error threshold and therefore is not allowed to be coarsened
 *
 * @note In order to keep track of the deviations it is necessary that if this function returns true, the family strictly has to be coarsened.
 *       In case this functions returns false, the element's family is not allowed to be coarsened. Otherwise this criterion will break!
 */
static bool
cmc_t8_elements_comply_to_relative_error_threshold(AdaptData& adapt_data, t8_forest_t current_forest, const int which_tree, const int lelement_id, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE

    /* Iterate over all variables which need to be considered */
    for (auto var)
    {
        /* Calculate the inerpolation result if this coarsening would have happen */
        cmc_universal_type_t interpolation_result = cmc_t8_adapt_evaluate_interpolation_function(adapt_data, current_forest, which_tree, lelement_id, ts, num_elements, elements);

        /* Since this will be (eventually) the first coarsening step, we can just calculate the normal relative deviation */
        const std::vector<double> single_deviations = adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->calculate_relative_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);

        /* Get the maximum deviation */
        auto max_elem = std::max_element(single_deviations.begin(), single_deviations.end());
        const double max_dev = (max_elem != single_deviations.end() ? *max_elem : 0.0);

        /* Check if the deviation complies with the error threshold */
        if (max_dev <= adapt_data->t8_data->settings.max_err)
        {
            /* Store the global element id of the first element */
            adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

            /* Calculate the maximal difference range for the values */
            const std::vector<double> abs_deviation = adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);
            auto max_elem_abs_dev = std::max_element(abs_deviation.begin(), abs_deviation.end());
            double max_dev_abs = (max_elem_abs_dev != single_deviations.end() ? *max_elem_abs_dev : 0.0);

            adapt_data->associated_max_deviations_new.front().push_back(max_dev_abs);

            /* Increment the coarsening counter */
            ++(adapt_data->coarsening_counter);

            /* We are going to save the already interpolated value when the coarsening is indeed applied */
            adapt_data->interpolation_data->interpolated_data.front().push_back(std::make_pair(lelement_id, interpolation_result));

            /* If so, we are able to coarsen the element's family */
            return true;
        } else
        {
            /* Since the threshold is not fullfilled, we cannot coarsen the family */
            return false;
        }

        /////////////////////////////////

        /* Declare a vector storing the adjusted deviations */
            std::vector<double> adjusted_deviations;
            adjusted_deviations.reserve(num_elements);
            std::vector<double> previous_max_range;
            previous_max_range.reserve(num_elements);

            /* Check if the error threshold is fulfilled for each element of this family */
            for (int i = 0; i < num_elements; ++i)
            {
                /* Get the previous maximum deviation */
                auto iter_last_max_dev_elem_id = std::lower_bound(adapt_data->associated_deviations_gelement_id.begin(), adapt_data->associated_deviations_gelement_id.end(), static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id + i));
                int pos = std::distance(adapt_data->associated_deviations_gelement_id.begin(), iter_last_max_dev_elem_id);

                /* If the element was not found (it was maybe a missing_value dummy element) and therefore, we do not coarsen at the moment */
                if (iter_last_max_dev_elem_id == adapt_data->associated_deviations_gelement_id.end())
                {
                    return false;
                }

                /* Save the previous maximum range */
                previous_max_range.push_back(adapt_data->associated_max_deviations.front()[pos]);

                /* Calculate the maximum deviation taking the the previous deviation into account as well */
                adjusted_deviations.push_back(calculate_two_step_relative_max_deviation(adapt_data->associated_max_deviations.front()[pos], adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->operator[](static_cast<size_t>(lelement_id + i)), interpolation_result));

                /* Check whether the error threshold is violated */
                if (adjusted_deviations.back() > adapt_data->t8_data->settings.max_err)
                {
                    /* In this case, we can stop and cannot coarsen the element's family */
                    return false;
                }
            }

            /* If we reach this line, the error threshold has not been violated and we can indeed coarsen the element's family */
            /* Save the global element id of the first element */
            adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

            /* Calculate the maximal difference range for the values */
            std::vector<double> abs_deviation = adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);
            /* Add the maximum range from previous coarsening steps */
            for (int i = 0; i < num_elements; ++i)
            {
                abs_deviation[i] += previous_max_range[i];
            }

            /* Calculate the maximal value range within this element's family central to the interpolated value */
            auto max_elem_abs_dev = std::max_element(abs_deviation.begin(), abs_deviation.end());
            double max_dev_abs = (max_elem_abs_dev != abs_deviation.end() ? *max_elem_abs_dev : 0.0);

            adapt_data->associated_max_deviations_new.front().push_back(max_dev_abs);

            /* Increment the coarsening counter */
            ++(adapt_data->coarsening_counter);

            /* We are going to save the already interpolated value when the coarsening is indeed applied */
            adapt_data->interpolation_data->interpolated_data.front().push_back(std::make_pair(lelement_id, interpolation_result));

            return true;



    }







    /* We distinguish between the compression mode, because dependent on it, we check has to be perfomed for all variables or just for a single variable */
    if (adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
    {
        /* On the finest level we can just calculate the relative deviation, but on coarser levels we need to take previous coarsening steps into account */
        if (ts->t8_element_level(elements[0]) != adapt_data->t8_data->assets->initial_refinement_lvl)
        {
            /* Iterate over all variables */
            for (size_t var_iter = 0; var_iter < adapt_data->t8_data->vars.size(); ++var_iter)
            {
                /* We need to store the current variable id for the calculation of the potential interpolation */
                adapt_data->interpolation_data->current_var_id = var_iter;
                /* Calculate the value resulting from the coarsening (if we would apply the coarsening) */
                cmc_universal_type_t interpolation_result = cmc_t8_adapt_evaluate_interpolation_function(adapt_data, current_forest, which_tree, lelement_id, ts, num_elements, elements);
               
                /* In this case, we are already up higher in the hierachy and need to take previous deviations into account */

                /* Declare a vector storing the adjusted deviations */
                std::vector<double> adjusted_deviations;
                adjusted_deviations.reserve(num_elements);
                std::vector<double> previous_max_range;
                previous_max_range.reserve(num_elements);

                /* Check if the error threshold is fulfilled for each element of this family */
                for (int i = 0; i < num_elements; ++i)
                {
                    /* Get the previous maximum deviation */
                    auto iter_last_max_dev_elem_id = std::lower_bound(adapt_data->associated_deviations_gelement_id.begin(), adapt_data->associated_deviations_gelement_id.end(), static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id + i));
                    int pos = std::distance(adapt_data->associated_deviations_gelement_id.begin(), iter_last_max_dev_elem_id);

                    /* If the element was not found (it was maybe a missing_value dummy element) and therefore, we do not coarsen at the moment */
                    if (iter_last_max_dev_elem_id == adapt_data->associated_deviations_gelement_id.end())
                    {
                        /* If we have to stop early, we need to pop back the last deviations and interpolated values */
                        for (size_t rm_var_iter = 0; rm_var_iter < var_iter; ++rm_var_iter)
                        {
                            adapt_data->associated_max_deviations_new[rm_var_iter].pop_back();
                            adapt_data->interpolation_data->interpolated_data[rm_var_iter].pop_back();
                        }
                        /* Do not coarsen the family */
                        return false;
                    }

                    /* Save the previous maximum range */
                    previous_max_range.push_back(adapt_data->associated_max_deviations[var_iter][pos]);

                    /* Calculate the maximum deviation taking the previous deviations into account as well */
                    adjusted_deviations.push_back(calculate_two_step_relative_max_deviation(adapt_data->associated_max_deviations[var_iter][pos], adapt_data->t8_data->vars[var_iter]->var->data->operator[](static_cast<size_t>(lelement_id + i)), interpolation_result));

                    /* Check whether the error threshold is violated */
                    if (adjusted_deviations.back() > adapt_data->t8_data->settings.max_err)
                    {
                        /* If we have to stop early, we need to pop back the last deviations and interpolated values */
                        for (size_t rm_var_iter = 0; rm_var_iter < var_iter; ++rm_var_iter)
                        {
                            adapt_data->associated_max_deviations_new[rm_var_iter].pop_back();
                            adapt_data->interpolation_data->interpolated_data[rm_var_iter].pop_back();
                        }

                        /* In this case, we can stop and cannot coarsen the element's family */
                        return false;
                    }
                }

                /* Calculate the maximal difference range for the values */
                std::vector<double> abs_deviation = adapt_data->t8_data->vars[var_iter]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[var_iter]->var->missing_value);
                /* Add the maximum range from previous coarsening steps */
                for (int i = 0; i < num_elements; ++i)
                {
                    abs_deviation[i] += previous_max_range[i];
                }

                /* Calculate the maximal value range within this element's family central to the interpolated value */
                auto max_elem_abs_dev = std::max_element(abs_deviation.begin(), abs_deviation.end());
                double max_dev_abs = (max_elem_abs_dev != abs_deviation.end() ? *max_elem_abs_dev : 0.0);

                adapt_data->associated_max_deviations_new[var_iter].push_back(max_dev_abs);

                /* We are going to save the already interpolated value when the coarsening is likely to be applied */
                adapt_data->interpolation_data->interpolated_data[var_iter].push_back(std::make_pair(lelement_id, interpolation_result));
            }

            /* Save the first global element id for this coarsening process */
            adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

            /* Increment the coarsening counter */
            ++(adapt_data->coarsening_counter);

            /* Coarsen the element's family */
            return true;
        } else
        {
            /* Iterate over all variables */
            for (size_t var_iter = 0; var_iter < adapt_data->t8_data->vars.size(); ++var_iter)
            {
                /* We need to store the current variable id for the calculation of the potential interpolation */
                adapt_data->interpolation_data->current_var_id = var_iter;
                /* Calculate the value resulting from the coarsening (if we would apply the coarsening) */
                cmc_universal_type_t interpolation_result = cmc_t8_adapt_evaluate_interpolation_function(adapt_data, current_forest, which_tree, lelement_id, ts, num_elements, elements);
                /* Since this will be (eventually) the first coarsening step, we can just calculate the normal relative deviation */
                const std::vector<double> single_deviations = adapt_data->t8_data->vars[var_iter]->var->data->calculate_relative_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[var_iter]->var->missing_value);
                
                /* Get the maximum deviation */
                auto max_elem = std::max_element(single_deviations.begin(), single_deviations.end());
                const double max_dev = (max_elem != single_deviations.end() ? *max_elem : 0.0);

                /* Check if the deviation complies with the error threshold */
                if (max_dev <= adapt_data->t8_data->settings.max_err)
                {
                    /* Calculate the maximal difference range for the values */
                    const std::vector<double> abs_deviation = adapt_data->t8_data->vars[var_iter]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[var_iter]->var->missing_value);
                    auto max_elem_abs_dev = std::max_element(abs_deviation.begin(), abs_deviation.end());
                    double max_dev_abs = (max_elem_abs_dev != single_deviations.end() ? *max_elem_abs_dev : 0.0);

                    adapt_data->associated_max_deviations_new[var_iter].push_back(max_dev_abs);
                    
                    /* We are going to save the already interpolated value when the coarsening is likely to be applied */
                    adapt_data->interpolation_data->interpolated_data[var_iter].push_back(std::make_pair(lelement_id, interpolation_result));
                } else
                {
                    /* If we have to stop early, we need to pop back the last deviations */
                    for (size_t rm_var_iter = 0; rm_var_iter < var_iter; ++rm_var_iter)
                    {
                        adapt_data->associated_max_deviations_new[rm_var_iter].pop_back();
                        adapt_data->interpolation_data->interpolated_data[rm_var_iter].pop_back();
                    }

                    /* Since the threshold is not fullfilled, we cannot coarsen the family */
                    return false;
                }
            }

            /* Save the first global element id for this coarsening process */
            adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

            /* When we have reached this point, the threshold is fullfilled for each varibale */
            /* Increment the coarsening counter */
            ++(adapt_data->coarsening_counter);

            /* If so, we are able to coarsen the element's family */
            return true;
        }
    } else {
        /* In case a 'One for One' compression mode is used */

        /* Calculate the inerpolation result if this coarsening would have happen */
        cmc_universal_type_t interpolation_result = cmc_t8_adapt_evaluate_interpolation_function(adapt_data, current_forest, which_tree, lelement_id, ts, num_elements, elements);

        /* On the finest level we can just calculate the relative deviation, but on coarser levels we need to take previous coarsening steps into account */
        if (ts->t8_element_level(elements[0]) != adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl)
        {
            /* In this case, we are already up higher in the hierachy and need to take previous deviations into account as well */

            /* Declare a vector storing the adjusted deviations */
            std::vector<double> adjusted_deviations;
            adjusted_deviations.reserve(num_elements);
            std::vector<double> previous_max_range;
            previous_max_range.reserve(num_elements);

            /* Check if the error threshold is fulfilled for each element of this family */
            for (int i = 0; i < num_elements; ++i)
            {
                /* Get the previous maximum deviation */
                auto iter_last_max_dev_elem_id = std::lower_bound(adapt_data->associated_deviations_gelement_id.begin(), adapt_data->associated_deviations_gelement_id.end(), static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id + i));
                int pos = std::distance(adapt_data->associated_deviations_gelement_id.begin(), iter_last_max_dev_elem_id);

                /* If the element was not found (it was maybe a missing_value dummy element) and therefore, we do not coarsen at the moment */
                if (iter_last_max_dev_elem_id == adapt_data->associated_deviations_gelement_id.end())
                {
                    return false;
                }

                /* Save the previous maximum range */
                previous_max_range.push_back(adapt_data->associated_max_deviations.front()[pos]);

                /* Calculate the maximum deviation taking the the previous deviation into account as well */
                adjusted_deviations.push_back(calculate_two_step_relative_max_deviation(adapt_data->associated_max_deviations.front()[pos], adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->operator[](static_cast<size_t>(lelement_id + i)), interpolation_result));

                /* Check whether the error threshold is violated */
                if (adjusted_deviations.back() > adapt_data->t8_data->settings.max_err)
                {
                    /* In this case, we can stop and cannot coarsen the element's family */
                    return false;
                }
            }

            /* If we reach this line, the error threshold has not been violated and we can indeed coarsen the element's family */
            /* Save the global element id of the first element */
            adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

            /* Calculate the maximal difference range for the values */
            std::vector<double> abs_deviation = adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);
            /* Add the maximum range from previous coarsening steps */
            for (int i = 0; i < num_elements; ++i)
            {
                abs_deviation[i] += previous_max_range[i];
            }

            /* Calculate the maximal value range within this element's family central to the interpolated value */
            auto max_elem_abs_dev = std::max_element(abs_deviation.begin(), abs_deviation.end());
            double max_dev_abs = (max_elem_abs_dev != abs_deviation.end() ? *max_elem_abs_dev : 0.0);

            adapt_data->associated_max_deviations_new.front().push_back(max_dev_abs);

            /* Increment the coarsening counter */
            ++(adapt_data->coarsening_counter);

            /* We are going to save the already interpolated value when the coarsening is indeed applied */
            adapt_data->interpolation_data->interpolated_data.front().push_back(std::make_pair(lelement_id, interpolation_result));

            return true;
        } else
        {
            /* Since this will be (eventually) the first coarsening step, we can just calculate the normal relative deviation */
            const std::vector<double> single_deviations = adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->calculate_relative_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);

            /* Get the maximum deviation */
            auto max_elem = std::max_element(single_deviations.begin(), single_deviations.end());
            const double max_dev = (max_elem != single_deviations.end() ? *max_elem : 0.0);

            /* Check if the deviation complies with the error threshold */
            if (max_dev <= adapt_data->t8_data->settings.max_err)
            {
                /* Store the global element id of the first element */
                adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

                /* Calculate the maximal difference range for the values */
                const std::vector<double> abs_deviation = adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);
                auto max_elem_abs_dev = std::max_element(abs_deviation.begin(), abs_deviation.end());
                double max_dev_abs = (max_elem_abs_dev != single_deviations.end() ? *max_elem_abs_dev : 0.0);

                adapt_data->associated_max_deviations_new.front().push_back(max_dev_abs);

                /* Increment the coarsening counter */
                ++(adapt_data->coarsening_counter);

                /* We are going to save the already interpolated value when the coarsening is indeed applied */
                adapt_data->interpolation_data->interpolated_data.front().push_back(std::make_pair(lelement_id, interpolation_result));

                /* If so, we are able to coarsen the element's family */
                return true;
            } else
            {
                /* Since the threshold is not fullfilled, we cannot coarsen the family */
                return false;
            }
        }
    }
    #endif
}





/* Coarsening Criteria */
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
cmc_t8_adapt_callback_coarsen_based_on_relative_error_threshold (t8_forest_t forest,
                                                                 t8_forest_t forest_from,
                                                                 int which_tree,
                                                                 int lelement_id,
                                                                 t8_eclass_scheme_c * ts,
                                                                 const int is_family,
                                                                 const int num_elements,
                                                                 t8_element_t * elements[]) 
{
    /* Get the adapt data from the forest */
    AdaptData* adapt_data = static_cast<AdaptData*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /** Only a family of elements can be coarsened.
     * If only a single element is passed to the callback, the callback function can be returned immediately. 
    */
    if (is_family == 0)
    {
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
        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    }

}




}
#endif










////////////////////////
////////////////////////
////// Old functions
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
    /* Get the initial refinement level of the data based on the compression criterion */
    const int initial_ref_lvl{(adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL) ?
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
    /* Get the initial refinement level of the data based on the compression criterion */
    const int initial_ref_lvl{(adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL) ?
                               adapt_data->t8_data->assets->initial_refinement_lvl :
                               adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl};
    /* Check if the current element level resemles a coarsening process of those elements from the beginning of the compression onwards */
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

/**
 * @brief This functions preliminary evalutes the interpolation function (e.g. in order to check if a coarsening step would comply with the error threshold)
 * 
 * @param adapt_data The adapt_data holding all information of the current adaptation as well as a pointer to the interpolation struct respectively the interpolation function
 * @param current_forest The forest which is currently constructed during the adaptation 
 * @param which_tree The current tree in whihc the elements are evaluated
 * @param lelement_id The first element id of the element which are about to be evaluated
 * @param ts The element scheme
 * @param num_elements The size of the family of elements
 * @param elements The array of pointeers of elements of the considered family
 * @return cmc_universal_type_t The interpolated value with whom the element' values will be (potentially) replaced
 *
 * @note Since this function preliminary evaluates the interpolation, the @var current_forest is not allowed to be used for the computation duirng the interpolation
 */
static inline
cmc_universal_type_t
cmc_t8_adapt_evaluate_interpolation_function(cmc_t8_adapt_data_t adapt_data, t8_forest_t current_forest, const int which_tree, const int lelement_id, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    /* We need to set the interpolation data struct to the old_forest in order to access the struct during the interpolation call */
    t8_forest_set_user_data(adapt_data->forest_begin, adapt_data->interpolation_data);

    /* We calculate the potential value for a coarsening step */
    return adapt_data->interpolation_data->interpolate(adapt_data->forest_begin, current_forest, which_tree, ts, -1, num_elements, lelement_id, 1, lelement_id - (adapt_data->coarsening_counter) * num_elements);
    #endif
}

/**
 * @brief This function's only application is within the serial relative error threshold criterion
 * 
 * @param adapt_data The current adapt_data resembling the status of the adaptation
 */
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
bool cmc_t8_elements_fulfill_rel_error_threshold_serial_with_initial_data_present(cmc_t8_adapt_data_t adapt_data, t8_locidx_t lelement_id, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE
    
    /* It is strictly sepereated between the compression modes */
    /* This is due to the fact, that in a 'One for All' compression mode, the threshold needs to be fullfilled for every varibale
     * whereas in a 'One for One' compression mode the threshold only needs to hold for a single variable itself */
    /* Check which mode is used */
    if (adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
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

/**
 * @brief This function calculates (an approximation) of the deviation up to the initial data and decides based on the supplied error threshold whether the element's family shall be coarsened or not
 * 
 * @param adapt_data A pointer to the current @struct cmc_t8_adapt_data holding the status and data of the adaptation
 * @param current_forest The current forest which is built from the adaptation
 * @param lelement_id the first local element id of it's family members
 * @param num_elements The number of elements within this family
 * @return true The element's family complies with the relative error threshold and shall be coarsened 
 * @return false The element's family does not comply with the relative error threshold and therefore is not allowed to be coarsened
 *
 * @note In order to keep track of the deviations it is necessary that if this function returns true, the family strictly has to be coarsened.
 *       In case this functions returns false, the element's family is not allowed to be coarsened. Otherwise this criterion will break!
 */
static inline
bool
cmc_t8_elements_comply_to_relative_error_threshold(cmc_t8_adapt_data_t adapt_data, t8_forest_t current_forest, const int which_tree, const int lelement_id, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE

    /* We distinguish between the compression mode, because dependent on it, we check has to be perfomed for all variables or just for a single variable */
    if (adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
    {
        /* On the finest level we can just calculate the relative deviation, but on coarser levels we need to take previous coarsening steps into account */
        if (ts->t8_element_level(elements[0]) != adapt_data->t8_data->assets->initial_refinement_lvl)
        {
            /* Iterate over all variables */
            for (size_t var_iter = 0; var_iter < adapt_data->t8_data->vars.size(); ++var_iter)
            {
                /* We need to store the current variable id for the calculation of the potential interpolation */
                adapt_data->interpolation_data->current_var_id = var_iter;
                /* Calculate the value resulting from the coarsening (if we would apply the coarsening) */
                cmc_universal_type_t interpolation_result = cmc_t8_adapt_evaluate_interpolation_function(adapt_data, current_forest, which_tree, lelement_id, ts, num_elements, elements);
               
                /* In this case, we are already up higher in the hierachy and need to take previous deviations into account */

                /* Declare a vector storing the adjusted deviations */
                std::vector<double> adjusted_deviations;
                adjusted_deviations.reserve(num_elements);
                std::vector<double> previous_max_range;
                previous_max_range.reserve(num_elements);

                /* Check if the error threshold is fulfilled for each element of this family */
                for (int i = 0; i < num_elements; ++i)
                {
                    /* Get the previous maximum deviation */
                    auto iter_last_max_dev_elem_id = std::lower_bound(adapt_data->associated_deviations_gelement_id.begin(), adapt_data->associated_deviations_gelement_id.end(), static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id + i));
                    int pos = std::distance(adapt_data->associated_deviations_gelement_id.begin(), iter_last_max_dev_elem_id);

                    /* If the element was not found (it was maybe a missing_value dummy element) and therefore, we do not coarsen at the moment */
                    if (iter_last_max_dev_elem_id == adapt_data->associated_deviations_gelement_id.end())
                    {
                        /* If we have to stop early, we need to pop back the last deviations and interpolated values */
                        for (size_t rm_var_iter = 0; rm_var_iter < var_iter; ++rm_var_iter)
                        {
                            adapt_data->associated_max_deviations_new[rm_var_iter].pop_back();
                            adapt_data->interpolation_data->interpolated_data[rm_var_iter].pop_back();
                        }
                        /* Do not coarsen the family */
                        return false;
                    }

                    /* Save the previous maximum range */
                    previous_max_range.push_back(adapt_data->associated_max_deviations[var_iter][pos]);

                    /* Calculate the maximum deviation taking the previous deviations into account as well */
                    adjusted_deviations.push_back(calculate_two_step_relative_max_deviation(adapt_data->associated_max_deviations[var_iter][pos], adapt_data->t8_data->vars[var_iter]->var->data->operator[](static_cast<size_t>(lelement_id + i)), interpolation_result));

                    /* Check whether the error threshold is violated */
                    if (adjusted_deviations.back() > adapt_data->t8_data->settings.max_err)
                    {
                        /* If we have to stop early, we need to pop back the last deviations and interpolated values */
                        for (size_t rm_var_iter = 0; rm_var_iter < var_iter; ++rm_var_iter)
                        {
                            adapt_data->associated_max_deviations_new[rm_var_iter].pop_back();
                            adapt_data->interpolation_data->interpolated_data[rm_var_iter].pop_back();
                        }

                        /* In this case, we can stop and cannot coarsen the element's family */
                        return false;
                    }
                }

                /* Calculate the maximal difference range for the values */
                std::vector<double> abs_deviation = adapt_data->t8_data->vars[var_iter]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[var_iter]->var->missing_value);
                /* Add the maximum range from previous coarsening steps */
                for (int i = 0; i < num_elements; ++i)
                {
                    abs_deviation[i] += previous_max_range[i];
                }

                /* Calculate the maximal value range within this element's family central to the interpolated value */
                auto max_elem_abs_dev = std::max_element(abs_deviation.begin(), abs_deviation.end());
                double max_dev_abs = (max_elem_abs_dev != abs_deviation.end() ? *max_elem_abs_dev : 0.0);

                adapt_data->associated_max_deviations_new[var_iter].push_back(max_dev_abs);

                /* We are going to save the already interpolated value when the coarsening is likely to be applied */
                adapt_data->interpolation_data->interpolated_data[var_iter].push_back(std::make_pair(lelement_id, interpolation_result));
            }

            /* Save the first global element id for this coarsening process */
            adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

            /* Increment the coarsening counter */
            ++(adapt_data->coarsening_counter);

            /* Coarsen the element's family */
            return true;
        } else
        {
            /* Iterate over all variables */
            for (size_t var_iter = 0; var_iter < adapt_data->t8_data->vars.size(); ++var_iter)
            {
                /* We need to store the current variable id for the calculation of the potential interpolation */
                adapt_data->interpolation_data->current_var_id = var_iter;
                /* Calculate the value resulting from the coarsening (if we would apply the coarsening) */
                cmc_universal_type_t interpolation_result = cmc_t8_adapt_evaluate_interpolation_function(adapt_data, current_forest, which_tree, lelement_id, ts, num_elements, elements);
                /* Since this will be (eventually) the first coarsening step, we can just calculate the normal relative deviation */
                const std::vector<double> single_deviations = adapt_data->t8_data->vars[var_iter]->var->data->calculate_relative_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[var_iter]->var->missing_value);
                
                /* Get the maximum deviation */
                auto max_elem = std::max_element(single_deviations.begin(), single_deviations.end());
                const double max_dev = (max_elem != single_deviations.end() ? *max_elem : 0.0);

                /* Check if the deviation complies with the error threshold */
                if (max_dev <= adapt_data->t8_data->settings.max_err)
                {
                    /* Calculate the maximal difference range for the values */
                    const std::vector<double> abs_deviation = adapt_data->t8_data->vars[var_iter]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[var_iter]->var->missing_value);
                    auto max_elem_abs_dev = std::max_element(abs_deviation.begin(), abs_deviation.end());
                    double max_dev_abs = (max_elem_abs_dev != single_deviations.end() ? *max_elem_abs_dev : 0.0);

                    adapt_data->associated_max_deviations_new[var_iter].push_back(max_dev_abs);
                    
                    /* We are going to save the already interpolated value when the coarsening is likely to be applied */
                    adapt_data->interpolation_data->interpolated_data[var_iter].push_back(std::make_pair(lelement_id, interpolation_result));
                } else
                {
                    /* If we have to stop early, we need to pop back the last deviations */
                    for (size_t rm_var_iter = 0; rm_var_iter < var_iter; ++rm_var_iter)
                    {
                        adapt_data->associated_max_deviations_new[rm_var_iter].pop_back();
                        adapt_data->interpolation_data->interpolated_data[rm_var_iter].pop_back();
                    }

                    /* Since the threshold is not fullfilled, we cannot coarsen the family */
                    return false;
                }
            }

            /* Save the first global element id for this coarsening process */
            adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

            /* When we have reached this point, the threshold is fullfilled for each varibale */
            /* Increment the coarsening counter */
            ++(adapt_data->coarsening_counter);

            /* If so, we are able to coarsen the element's family */
            return true;
        }
    } else {
        /* In case a 'One for One' compression mode is used */

        /* Calculate the inerpolation result if this coarsening would have happen */
        cmc_universal_type_t interpolation_result = cmc_t8_adapt_evaluate_interpolation_function(adapt_data, current_forest, which_tree, lelement_id, ts, num_elements, elements);

        /* On the finest level we can just calculate the relative deviation, but on coarser levels we need to take previous coarsening steps into account */
        if (ts->t8_element_level(elements[0]) != adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl)
        {
            /* In this case, we are already up higher in the hierachy and need to take previous deviations into account as well */

            /* Declare a vector storing the adjusted deviations */
            std::vector<double> adjusted_deviations;
            adjusted_deviations.reserve(num_elements);
            std::vector<double> previous_max_range;
            previous_max_range.reserve(num_elements);

            /* Check if the error threshold is fulfilled for each element of this family */
            for (int i = 0; i < num_elements; ++i)
            {
                /* Get the previous maximum deviation */
                auto iter_last_max_dev_elem_id = std::lower_bound(adapt_data->associated_deviations_gelement_id.begin(), adapt_data->associated_deviations_gelement_id.end(), static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id + i));
                int pos = std::distance(adapt_data->associated_deviations_gelement_id.begin(), iter_last_max_dev_elem_id);

                /* If the element was not found (it was maybe a missing_value dummy element) and therefore, we do not coarsen at the moment */
                if (iter_last_max_dev_elem_id == adapt_data->associated_deviations_gelement_id.end())
                {
                    return false;
                }

                /* Save the previous maximum range */
                previous_max_range.push_back(adapt_data->associated_max_deviations.front()[pos]);

                /* Calculate the maximum deviation taking the the previous deviation into account as well */
                adjusted_deviations.push_back(calculate_two_step_relative_max_deviation(adapt_data->associated_max_deviations.front()[pos], adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->operator[](static_cast<size_t>(lelement_id + i)), interpolation_result));

                /* Check whether the error threshold is violated */
                if (adjusted_deviations.back() > adapt_data->t8_data->settings.max_err)
                {
                    /* In this case, we can stop and cannot coarsen the element's family */
                    return false;
                }
            }

            /* If we reach this line, the error threshold has not been violated and we can indeed coarsen the element's family */
            /* Save the global element id of the first element */
            adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

            /* Calculate the maximal difference range for the values */
            std::vector<double> abs_deviation = adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);
            /* Add the maximum range from previous coarsening steps */
            for (int i = 0; i < num_elements; ++i)
            {
                abs_deviation[i] += previous_max_range[i];
            }

            /* Calculate the maximal value range within this element's family central to the interpolated value */
            auto max_elem_abs_dev = std::max_element(abs_deviation.begin(), abs_deviation.end());
            double max_dev_abs = (max_elem_abs_dev != abs_deviation.end() ? *max_elem_abs_dev : 0.0);

            adapt_data->associated_max_deviations_new.front().push_back(max_dev_abs);

            /* Increment the coarsening counter */
            ++(adapt_data->coarsening_counter);

            /* We are going to save the already interpolated value when the coarsening is indeed applied */
            adapt_data->interpolation_data->interpolated_data.front().push_back(std::make_pair(lelement_id, interpolation_result));

            return true;
        } else
        {
            /* Since this will be (eventually) the first coarsening step, we can just calculate the normal relative deviation */
            const std::vector<double> single_deviations = adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->calculate_relative_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);

            /* Get the maximum deviation */
            auto max_elem = std::max_element(single_deviations.begin(), single_deviations.end());
            const double max_dev = (max_elem != single_deviations.end() ? *max_elem : 0.0);

            /* Check if the deviation complies with the error threshold */
            if (max_dev <= adapt_data->t8_data->settings.max_err)
            {
                /* Store the global element id of the first element */
                adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

                /* Calculate the maximal difference range for the values */
                const std::vector<double> abs_deviation = adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);
                auto max_elem_abs_dev = std::max_element(abs_deviation.begin(), abs_deviation.end());
                double max_dev_abs = (max_elem_abs_dev != single_deviations.end() ? *max_elem_abs_dev : 0.0);

                adapt_data->associated_max_deviations_new.front().push_back(max_dev_abs);

                /* Increment the coarsening counter */
                ++(adapt_data->coarsening_counter);

                /* We are going to save the already interpolated value when the coarsening is indeed applied */
                adapt_data->interpolation_data->interpolated_data.front().push_back(std::make_pair(lelement_id, interpolation_result));

                /* If so, we are able to coarsen the element's family */
                return true;
            } else
            {
                /* Since the threshold is not fullfilled, we cannot coarsen the family */
                return false;
            }
        }
    }
    #endif
}


/**
 * @brief This function calculates (an approximation) of the deviation up to the initial data and decides based on the supplied error threshold whether the element's family shall be coarsened or not
 * 
 * @param adapt_data A pointer to the current @struct cmc_t8_adapt_data holding the status and data of the adaptation
 * @param current_forest The current forest which is built from the adaptation
 * @param lelement_id the first local element id of it's family members
 * @param num_elements The number of elements within this family
 * @return true The element's family complies with the relative error threshold and shall be coarsened 
 * @return false The element's family does not comply with the relative error threshold and therefore is not allowed to be coarsened
 *
 * @note In order to keep track of the deviations it is necessary that if this function returns true, the family strictly has to be coarsened.
 *       In case this functions returns false, the element's family is not allowed to be coarsened. Otherwise this criterion will break!
 */
static inline
bool
cmc_t8_elements_comply_to_absolute_error_threshold(cmc_t8_adapt_data_t adapt_data, t8_forest_t current_forest, const int which_tree, const int lelement_id, t8_eclass_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    #ifdef CMC_WITH_T8CODE

    /* We distinguish between the compression mode, because dependent on it, we check has to be perfomed for all variables or just for a single variable */
    if (adapt_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
    {
        /* On the finest level we can just calculate the relative deviation, but on coarser levels we need to take previous coarsening steps into account */
        if (ts->t8_element_level(elements[0]) != adapt_data->t8_data->assets->initial_refinement_lvl)
        {
            /* Iterate over all variables */
            for (size_t var_iter = 0; var_iter < adapt_data->t8_data->vars.size(); ++var_iter)
            {
                /* We need to store the current variable id for the calculation of the potential interpolation */
                adapt_data->interpolation_data->current_var_id = var_iter;
                /* Calculate the value resulting from the coarsening (if we would apply the coarsening) */
                cmc_universal_type_t interpolation_result = cmc_t8_adapt_evaluate_interpolation_function(adapt_data, current_forest, which_tree, lelement_id, ts, num_elements, elements);
               
                /* In this case, we are already up higher in the hierachy and need to take previous deviations into account */

                /* Declare a vector storing the adjusted deviations */
                std::vector<double> adjusted_deviations;
                adjusted_deviations.reserve(num_elements);

                /* Check if the error threshold is fulfilled for each element of this family */
                for (int i = 0; i < num_elements; ++i)
                {
                    /* Get the previous maximum deviation */
                    auto iter_last_max_dev_elem_id = std::lower_bound(adapt_data->associated_deviations_gelement_id.begin(), adapt_data->associated_deviations_gelement_id.end(), static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id + i));
                    int pos = std::distance(adapt_data->associated_deviations_gelement_id.begin(), iter_last_max_dev_elem_id);

                    /* If the element was not found (it was maybe a missing_value dummy element) and therefore, we do not coarsen at the moment */
                    if (iter_last_max_dev_elem_id == adapt_data->associated_deviations_gelement_id.end())
                    {
                        /* If we have to stop early, we need to pop back the last deviations and interpolated values */
                        for (size_t rm_var_iter = 0; rm_var_iter < var_iter; ++rm_var_iter)
                        {
                            adapt_data->associated_max_deviations_new[rm_var_iter].pop_back();
                            adapt_data->interpolation_data->interpolated_data[rm_var_iter].pop_back();
                        }
                        /* Do not coarsen the family */
                        return false;
                    }
                    /* Calculate the maximum deviation taking the previous deviations into account as well */
                    adjusted_deviations.push_back(calculate_two_step_absolute_max_deviation(adapt_data->associated_max_deviations[var_iter][pos], adapt_data->t8_data->vars[var_iter]->var->data->operator[](static_cast<size_t>(lelement_id + i)), interpolation_result));

                    /* Check whether the error threshold is violated */
                    if (adjusted_deviations.back() > adapt_data->t8_data->settings.max_err)
                    {
                        /* If we have to stop early, we need to pop back the last deviations and interpolated values */
                        for (size_t rm_var_iter = 0; rm_var_iter < var_iter; ++rm_var_iter)
                        {
                            adapt_data->associated_max_deviations_new[rm_var_iter].pop_back();
                            adapt_data->interpolation_data->interpolated_data[rm_var_iter].pop_back();
                        }

                        /* In this case, we can stop and cannot coarsen the element's family */
                        return false;
                    }
                }

                /* Get the maximum value of the deviations */
                auto max_elem_adjusted = std::max_element(adjusted_deviations.begin(), adjusted_deviations.end());
                double max_dev_new_adjusted = (max_elem_adjusted != adjusted_deviations.end() ? *max_elem_adjusted : 0.0);

                /* If we reach this line, the error threshold has not been violated and we can indeed coarsen the element's family */
                /* We need to store the maximum deviation for further adaptation steps */
                adapt_data->associated_max_deviations_new[var_iter].push_back(max_dev_new_adjusted);

                /* We are going to save the already interpolated value when the coarsening is likely to be applied */
                adapt_data->interpolation_data->interpolated_data[var_iter].push_back(std::make_pair(lelement_id, interpolation_result));
            }

            /* Save the first global element id for this coarsening process */
            adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

            /* Increment the coarsening counter */
            ++(adapt_data->coarsening_counter);

            /* Coarsen the element's family */
            return true;
        } else
        {
            /* Iterate over all variables */
            for (size_t var_iter = 0; var_iter < adapt_data->t8_data->vars.size(); ++var_iter)
            {
                /* We need to store the current variable id for the calculation of the potential interpolation */
                adapt_data->interpolation_data->current_var_id = var_iter;
                /* Calculate the value resulting from the coarsening (if we would apply the coarsening) */
                cmc_universal_type_t interpolation_result = cmc_t8_adapt_evaluate_interpolation_function(adapt_data, current_forest, which_tree, lelement_id, ts, num_elements, elements);
                /* Since this will be (eventually) the first coarsening step, we can just calculate the normal relative deviation */
                const std::vector<double> single_deviations = adapt_data->t8_data->vars[var_iter]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[var_iter]->var->missing_value);
                
                /* Get the maximum deviation */
                auto max_elem = std::max_element(single_deviations.begin(), single_deviations.end());
                double max_dev = (max_elem != single_deviations.end() ? *max_elem : 0.0);

                /* Check if the deviation complies with the error threshold */
                if (max_dev <= adapt_data->t8_data->settings.max_err)
                {
                    /* We need to store the maximum deviation for further adaptation steps */
                    adapt_data->associated_max_deviations_new[var_iter].push_back(max_dev);

                    /* We are going to save the already interpolated value when the coarsening is likely to be applied */
                    adapt_data->interpolation_data->interpolated_data[var_iter].push_back(std::make_pair(lelement_id, interpolation_result));
                } else
                {
                    /* If we have to stop early, we need to pop back the last deviations */
                    for (size_t rm_var_iter = 0; rm_var_iter < var_iter; ++rm_var_iter)
                    {
                        adapt_data->associated_max_deviations_new[rm_var_iter].pop_back();
                        adapt_data->interpolation_data->interpolated_data[rm_var_iter].pop_back();
                    }

                    /* Since the threshold is not fullfilled, we cannot coarsen the family */
                    return false;
                }
            }

            /* Save the first global element id for this coarsening process */
            adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

            /* When we have reached this point, the threshold is fullfilled for each varibale */
            /* Increment the coarsening counter */
            ++(adapt_data->coarsening_counter);

            /* If so, we are able to coarsen the element's family */
            return true;
        }
    } else {
        /* In case a 'One for One' compression mode is used */

        /* Calculate the inerpolation result if this coarsening would have happen */
        cmc_universal_type_t interpolation_result = cmc_t8_adapt_evaluate_interpolation_function(adapt_data, current_forest, which_tree, lelement_id, ts, num_elements, elements);

        /* On the finest level we can just calculate the relative deviation, but on coarser levels we need to take previous coarsening steps into account */
        if (ts->t8_element_level(elements[0]) != adapt_data->t8_data->vars[adapt_data->current_var_id]->assets->initial_refinement_lvl)
        {
            /* In this case, we are already up higher in the hierachy and need to take previous deviations into account as well */

            /* Declare a vector storing the adjusted deviations */
            std::vector<double> adjusted_deviations;
            adjusted_deviations.reserve(num_elements);

            /* Check if the error threshold is fulfilled for each element of this family */
            for (int i = 0; i < num_elements; ++i)
            {
                /* Get the previous maximum deviation */
                auto iter_last_max_dev_elem_id = std::lower_bound(adapt_data->associated_deviations_gelement_id.begin(), adapt_data->associated_deviations_gelement_id.end(), static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id + i));
                int pos = std::distance(adapt_data->associated_deviations_gelement_id.begin(), iter_last_max_dev_elem_id);

                /* If the element was not found (it was maybe a missing_value dummy element) and therefore, we do not coarsen at the moment */
                if (iter_last_max_dev_elem_id == adapt_data->associated_deviations_gelement_id.end())
                {
                    return false;
                }

                /* Calculate the maximum deviation taking the the previous deviation into account as well */
                adjusted_deviations.push_back(calculate_two_step_absolute_max_deviation(adapt_data->associated_max_deviations.front()[pos], adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->operator[](static_cast<size_t>(lelement_id + i)), interpolation_result));

                /* Check whether the error threshold is violated */
                if (adjusted_deviations.back() > adapt_data->t8_data->settings.max_err)
                {
                    /* In this case, we can stop and cannot coarsen the element's family */
                    return false;
                }
            }

            /* If we reach this line, the error threshold has not been violated and we can indeed coarsen the element's family */
            /* Save the global element id of the first element */
            adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));
            
            /* Get the maximum value of the deviations */
            auto max_elem_adjusted = std::max_element(adjusted_deviations.begin(), adjusted_deviations.end());
            double max_dev_new_adjusted = (max_elem_adjusted != adjusted_deviations.end() ? *max_elem_adjusted : 0.0);

            /* We need to store the maximum deviation for further adaptation steps */
            adapt_data->associated_max_deviations_new.front().push_back(max_dev_new_adjusted);

            /* Increment the coarsening counter */
            ++(adapt_data->coarsening_counter);

            /* We are going to save the already interpolated value when the coarsening is indeed applied */
            adapt_data->interpolation_data->interpolated_data.front().push_back(std::make_pair(lelement_id, interpolation_result));

            return true;
        } else
        {
            /* Since this will be (eventually) the first coarsening step, we can just calculate the normal relative deviation */
            const std::vector<double> single_deviations = adapt_data->t8_data->vars[adapt_data->current_var_id]->var->data->calculate_absolute_deviations_w_missing_values(lelement_id, lelement_id + num_elements -1, interpolation_result, adapt_data->t8_data->vars[adapt_data->current_var_id]->var->missing_value);

            /* Get the maximum deviation */
            auto max_elem = std::max_element(single_deviations.begin(), single_deviations.end());
            double max_dev = (max_elem != single_deviations.end() ? *max_elem : 0.0);

            /* Check if the deviation complies with the error threshold */
            if (max_dev <= adapt_data->t8_data->settings.max_err)
            {
                /* Store the global element id of the first element */
                adapt_data->associated_deviations_gelement_id_new.push_back(static_cast<uint64_t>(adapt_data->first_global_elem_id + lelement_id));

                /* We need to store the maximum deviation for further adaptation steps */
                adapt_data->associated_max_deviations_new.front().push_back(max_dev);

                /* Increment the coarsening counter */
                ++(adapt_data->coarsening_counter);

                /* We are going to save the already interpolated value when the coarsening is indeed applied */
                adapt_data->interpolation_data->interpolated_data.front().push_back(std::make_pair(lelement_id, interpolation_result));

                /* If so, we are able to coarsen the element's family */
                return true;
            } else
            {
                /* Since the threshold is not fullfilled, we cannot coarsen the family */
                return false;
            }
        }
    }
    #endif
}

/** End of static function building blocks **/

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

/**
 * @brief This function is a variation of the relative error criterion, but it does only work in serial and when the initial data (prior to the compression is present).
 *        The memory demand is quite high for this function because the initial data of all variables has to be kept. However, this enables the maximum possible compression for a relative error criterion.
 * @note This function is deprecated and should not be used, the parallel counterpart is to be prefered. (@see cmc_t8_adapt_callback_coarsen_based_on_relative_error_threshold(...))
 * 
 */
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
    if (cmc_t8_elements_fulfill_rel_error_threshold_serial_with_initial_data_present(adapt_data, lelement_id, ts, num_elements, elements))
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
    if (cmc_t8_elements_fulfill_rel_error_threshold_serial_with_initial_data_present(adapt_data, lelement_id, ts, num_elements, elements))
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
cmc_t8_adapt_callback_coarsen_based_on_relative_error_threshold (t8_forest_t forest,
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
        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    }
    
    /* Check if the element's family complies with the relative error threshold and therefore, will be coarsened */
    if (cmc_t8_elements_comply_to_relative_error_threshold(adapt_data, forest, which_tree, lelement_id, ts, num_elements, elements))
    {
        /* The family will be coarsened */
        return -1;
    } else
    {
        /* The element's family does not fullfill the relative error threshold and therefore, cannot be coarsened */
        return 0;
    }

    #else
    return CMC_ERR;
    #endif
}


/**
 * @brief Function determining whether the element (respectively the family of elements) complies with the supplied absolute error threshold
 * 
 * @param adapt_data The adapt_data from the forest which will be adapted
 * @param num_elements The amount of elements passed to the current call of the adapt fuction
 * @param elements Pointer to the corresponding element(s of the family) 
 * @return true If the elements fulfill the absolute error threshold and could be coarsened
 * @return false If the elements does not fulfill the error threshold and therefore cannot be coarsened 
 */
t8_locidx_t
cmc_t8_adapt_callback_coarsen_based_on_absolute_error_threshold (t8_forest_t forest,
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
        /* The element stays the same, it will not be refined and it cannot be coarsened */
        return 0;
    }
    
    /* Check if the element's family complies with the absolute error threshold and therefore, will be coarsened */
    if (cmc_t8_elements_comply_to_absolute_error_threshold(adapt_data, forest, which_tree, lelement_id, ts, num_elements, elements))
    {
        /* The family will be coarsened */
        return -1;
    } else
    {
        /* The element's family does not fullfill the absolute error threshold and therefore, cannot be coarsened */
        return 0;
    }

    #else
    return CMC_ERR;
    #endif
}

#endif
