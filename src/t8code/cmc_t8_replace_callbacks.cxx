#include "t8code/cmc_t8_replace_callbacks.h"
#include "t8code/cmc_t8code_data.hxx"
#include "utilities/cmc_container.h"
#include "utilities/cmc_log_functions.h"

#ifdef CMC_WITH_T8CODE

/**
 * @note A interpolation function which shall be used in combination with a relative error criterion is not allowed to use the @var forest_new
 * during the computation of the interpolaed value. This is due to the fact that the interpolation function is already called within the adaptation
 * in order to check for compliance with the the relative error threshold and at that stage the @var forest_new ist not yet committed/fully built.
 */


void
cmc_t8_general_interpolation_during_compression (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                                 t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                                 int num_incoming, t8_locidx_t first_incoming)
{
    #ifdef CMC_WITH_T8CODE
    /* Get the interpolation data */
    cmc_t8_interpolation_data_t interpolation_data = static_cast<cmc_t8_interpolation_data_t>(t8_forest_get_user_data(forest_old));
    
    /** During compression, it is assumed that there is no refinement **/

    /* Check which compression mode is used */
    if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
    {
        /* In case a 'One for One' compression is used */
        if (refine == -1)
        {
            /* The element's family will be coarsened */
            /* In case a 'One for All' compression is chosen */
            for (size_t var_id{0}; var_id < interpolation_data->t8_data->vars.size(); ++var_id)
            {
                /* Save the current id of the variable within the interpolation struct */
                interpolation_data->current_var_id = var_id;
                /* Call the actual replace function */
                cmc_universal_type_t coarsened_value = (*(interpolation_data->interpolate))(forest_old, forest_new, which_tree, ts, refine, num_outgoing, first_outgoing, num_incoming, first_incoming);
                /* Assign the value to the correct position in the data_new */
                interpolation_data->t8_data->vars[var_id]->var->data_new->assign_value(static_cast<size_t>(first_incoming),
                                                                                                           coarsened_value);
            }
        } else if (refine == 0)
        {
            /* The element stays the same, therefore, we can copy the value directly for each variable */
            for (size_t var_id{0}; var_id < interpolation_data->t8_data->vars.size(); ++var_id)
            {
                interpolation_data->t8_data->vars[var_id]->var->data_new->assign_value(static_cast<size_t>(first_incoming),
                                                                                       interpolation_data->t8_data->vars[var_id]->var->data->operator[](first_outgoing));
            }
        } else
        {
            cmc_err_msg("A refinement during a compression step should have not happened.");
        }
    } else
    {
        /* In case a 'One for One' compression is used */
        if (refine == -1)
        {
            /* The element's family will be coarsened */
            /* Call the actual replace function */
            cmc_universal_type_t coarsened_value = (interpolation_data->interpolate)(forest_old, forest_new, which_tree, ts, refine, num_outgoing, first_outgoing, num_incoming, first_incoming);
            /* Assign the value to the correct position in the data_new */
            interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->data_new->assign_value(static_cast<size_t>(first_incoming),
                                                                                                               coarsened_value);
        } else if (refine == 0)
        {
            /* The element stays the same, therefore, we can copy the value directly */
            interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->data_new->assign_value(static_cast<size_t>(first_incoming),
                                                                                                               interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->data->operator[](first_outgoing));
        } else
        {
            cmc_err_msg("A refinement during a compression step should have not happened.");
        }
    }
    #endif
}

/* The actual interpolation function for a standard arithmetic mean during a coarsening step */
static
cmc_universal_type_t
cmc_t8_compression_interpolation_std_mean(t8_forest_t forest_old, t8_forest_t t8_forest_new, t8_locidx_t which_tree,
                                          t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                          int num_incoming, t8_locidx_t first_incoming)
{
    #ifdef CMC_WITH_T8CODE
    /* Get the interpolation data */
    cmc_t8_interpolation_data_t interpolation_data = static_cast<cmc_t8_interpolation_data_t>(t8_forest_get_user_data(forest_old));

    if (refine == -1)
    {
        /* In case a coarsening has happened, we can calculate the standard mean */
        return interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->data->mean_value_same_type_wo_missing_vals(static_cast<size_t>(first_outgoing),
                                                                                                          static_cast<size_t>(first_outgoing + num_outgoing -1),
                                                                                                          interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->missing_value);
    } else if (refine == 0)
    {
        /* If the element stays the same, we can return the previous value */
        return interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->data->operator[](first_outgoing);
    } else
    {
        cmc_err_msg("A refinement during a compression step should have not happened.");
        return cmc_universal_type_t(static_cast<int>(CMC_ERR));
    }

    #endif
}

/* This a wrapper function returns a functor containing the actual interpolation function when called
 * This due to the fact that a interpolation function should return a cmc_universal_type_t which is a std::variant which does not work with the C-compiler
*/
cmc_t8_forest_interpolate_t
cmc_t8_compression_interpolation_standard_mean()
{
    #ifdef CMC_WITH_T8CODE
    /* Returns a functor holding the actual interpolation function */
    return new cmc_t8_forest_interpolate(cmc_t8_compression_interpolation_std_mean);
    #endif
}

/* The actual interpolation function that chooses the maximum value during a coarsening step */
static
cmc_universal_type_t
cmc_t8_compression_interpolation_maximum_function(t8_forest_t forest_old, t8_forest_t t8_forest_new, t8_locidx_t which_tree,
                                                  t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                                  int num_incoming, t8_locidx_t first_incoming)
{
    #ifdef CMC_WITH_T8CODE
    /* Get the interpolation data */
    cmc_t8_interpolation_data_t interpolation_data = static_cast<cmc_t8_interpolation_data_t>(t8_forest_get_user_data(forest_old));

    if (refine == -1)
    {
        /* In case a coarsening has happened, we can choose the maximum value of all of the elements */
        return interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->data->maximum_w_missing_vals(static_cast<size_t>(first_outgoing),
                                                                                                          static_cast<size_t>(first_outgoing + num_outgoing -1),
                                                                                                          interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->missing_value);
    } else if (refine == 0)
    {
        /* If the element stays the same, we can return the previous value */
        return interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->data->operator[](first_outgoing);
    } else
    {
        cmc_err_msg("A refinement during a compression step should have not happened.");
        return cmc_universal_type_t(static_cast<int>(CMC_ERR));
    }
    #endif
}

cmc_t8_forest_interpolate_t
cmc_t8_compression_interpolation_maximum()
{
    #ifdef CMC_WITH_T8CODE
    /* Returns a functor holding the actual interpolation function */
    return new cmc_t8_forest_interpolate(cmc_t8_compression_interpolation_maximum_function);
    #endif
}


/* The actual interpolation function that chooses the minimum value during a coarsening step */
static
cmc_universal_type_t
cmc_t8_compression_interpolation_minimum_function(t8_forest_t forest_old, t8_forest_t t8_forest_new, t8_locidx_t which_tree,
                                                  t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                                  int num_incoming, t8_locidx_t first_incoming)
{
    #ifdef CMC_WITH_T8CODE
    /* Get the interpolation data */
    cmc_t8_interpolation_data_t interpolation_data = static_cast<cmc_t8_interpolation_data_t>(t8_forest_get_user_data(forest_old));

    if (refine == -1)
    {
        /* In case a coarsening has happened, we can choose the maximum value of all of the elements */
        return interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->data->minimum_w_missing_vals(static_cast<size_t>(first_outgoing),
                                                                                                          static_cast<size_t>(first_outgoing + num_outgoing -1),
                                                                                                          interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->missing_value);
    } else if (refine == 0)
    {
        /* If the element stays the same, we can return the previous value */
        return interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->data->operator[](first_outgoing);
    } else
    {
        cmc_err_msg("A refinement during a compression step should have not happened.");
        return cmc_universal_type_t(static_cast<int>(CMC_ERR));
    }
    #endif
}

cmc_t8_forest_interpolate_t
cmc_t8_compression_interpolation_minimum()
{
    #ifdef CMC_WITH_T8CODE
    /* Returns a functor holding the actual interpolation function */
    return new cmc_t8_forest_interpolate(cmc_t8_compression_interpolation_minimum_function);
    #endif
}

void
cmc_t8_geo_data_interpolate_std_mean(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                     t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                     int num_incoming, t8_locidx_t first_incoming)
{
    /* (Currently not!) If overlapping boundaries can be coarsened, this function has to check it and set midpoints and data right */
    #ifdef CMC_WITH_T8CODE
    cmc_t8_interpolation_data_t interpolation_data = static_cast<cmc_t8_interpolation_data_t>(t8_forest_get_user_data(forest_new));
    if (refine == -1)
    {
        /* If a coarsening has been introduced */
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
        {
            /* In case a 'One for All' compression is chosen */
            for (size_t var_id{0}; var_id < interpolation_data->t8_data->vars.size(); ++var_id)
            {
                interpolation_data->t8_data->vars[var_id]->var->data_new->assign_value(static_cast<size_t>(first_incoming),
                                                                                       interpolation_data->t8_data->vars[var_id]->var->data->mean_value_same_type_wo_missing_vals(static_cast<size_t>(first_outgoing),
                                                                                                                                                                                  static_cast<size_t>(first_outgoing + num_outgoing -1),
                                                                                                                                                                                  interpolation_data->t8_data->vars[var_id]->var->missing_value));
            }
        }
        else
        {
            /* In case a 'One for One' compression is chosen */
            /* If a coarsening has been introduced */
            interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->data_new->assign_value(static_cast<size_t>(first_incoming),
                                                                                                               interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->data->mean_value_same_type_wo_missing_vals(static_cast<size_t>(first_outgoing),
                                                                                                               static_cast<size_t>(first_outgoing + num_outgoing -1),
                                                                                                               interpolation_data->t8_data->vars[interpolation_data->current_var_id]->var->missing_value));
        }
    } else if (refine == 0)
    {
        /* If the element stays the same */
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
        {
            /* In case a 'One for All' compression is chosen */
            for (size_t var_id{0}; var_id < interpolation_data->t8_data->vars.size(); ++var_id)
            {
                /* If no coarsening has happened, just copy the previous element's value */
                memcpy(static_cast<std::byte*>(interpolation_data->t8_data->vars[var_id]->get_initial_data_new_ptr()) + first_incoming * interpolation_data->t8_data->vars[var_id]->get_data_size(),
                       static_cast<std::byte*>(interpolation_data->t8_data->vars[var_id]->get_initial_data_ptr()) + first_outgoing * interpolation_data->t8_data->vars[var_id]->get_data_size(),
                       interpolation_data->t8_data->vars[var_id]->get_data_size());

            }
        }
        else
        {
            /* In case a 'One for One' compression is chosen */
            /* If no coarsening has happened, just copy the previous element's value */
            memcpy(static_cast<std::byte*>(interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_initial_data_new_ptr()) + first_incoming * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                   static_cast<std::byte*>(interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_initial_data_ptr()) + first_outgoing * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                   interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size());

        }
    } else
    {
        cmc_err_msg("A refinement should have not happend.");
    }
    #endif
}


void
cmc_t8_geo_data_interpolate_plain_copy_values(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                              t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                              int num_incoming, t8_locidx_t first_incoming)
{
    #ifdef CMC_WITH_T8CODE
    /* Get the interpolation data */
    cmc_t8_interpolation_data_t interpolation_data = static_cast<cmc_t8_interpolation_data_t>(t8_forest_get_user_data(forest_new));
    
    /* Check if elements got replaced in the new adapted forest (num_incoming > num_outgoing) */
    if (refine == 1)
    {
        /* If a refinement was introduced */
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
        {
            /* In case a 'One for All' compression is chosen */
            for (size_t var_id{0}; var_id < interpolation_data->t8_data->vars.size(); ++var_id)
            {
                /* If no coarsening has happened, just copy the previous element's value */
                for (int i{0}; i < num_incoming; ++ i)
                {
                    memcpy(static_cast<std::byte*>(interpolation_data->t8_data->vars[var_id]->get_initial_data_new_ptr()) + (first_incoming + i) * interpolation_data->t8_data->vars[var_id]->get_data_size(),
                           static_cast<std::byte*>(interpolation_data->t8_data->vars[var_id]->get_initial_data_ptr()) + first_outgoing * interpolation_data->t8_data->vars[var_id]->get_data_size(),
                           interpolation_data->t8_data->vars[var_id]->get_data_size());
                }
            }
        }
        else
        {
            /* In case a 'One for One' compression is chosen */
            /* If no coarsening has happened, just copy the previous element's value */
            for (int i{0}; i < num_incoming; ++ i)
            {
                memcpy(static_cast<std::byte*>(interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_initial_data_new_ptr()) + (first_incoming + i) * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                       static_cast<std::byte*>(interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_initial_data_ptr()) + first_outgoing * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                       interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size());
            }
        }
    } else if (refine == 0)
    {
        /* If the element stays the same */
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
        {
            /* In case a 'One for All' compression is chosen */
            for (size_t var_id{0}; var_id < interpolation_data->t8_data->vars.size(); ++var_id)
            {
                /* If no coarsening has happened, jsut copy the previous element's value */
                memcpy(static_cast<std::byte*>(interpolation_data->t8_data->vars[var_id]->get_initial_data_new_ptr()) + first_incoming * interpolation_data->t8_data->vars[var_id]->get_data_size(),
                       static_cast<std::byte*>(interpolation_data->t8_data->vars[var_id]->get_initial_data_ptr()) + first_outgoing * interpolation_data->t8_data->vars[var_id]->get_data_size(),
                       interpolation_data->t8_data->vars[var_id]->get_data_size());

            }
        }
        else
        {
            /* In case a 'One for One' compression is chosen */
            /* If no coarsening has happened, jsut copy the previous element's value */
            memcpy(static_cast<std::byte*>(interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_initial_data_new_ptr()) + first_incoming * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                   static_cast<std::byte*>(interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_initial_data_ptr()) + first_outgoing * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                   interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size());

        }
    } else
    {
        /* If coarsening was introduced */
        cmc_err_msg("Plain copying of values only is vald during refinement");
    }
    #endif
}

void
cmc_t8_geo_data_interpolate_error_threshold_adaption(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                                     t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                                     int num_incoming, t8_locidx_t first_incoming)
{
    #ifdef CMC_WITH_T8CODE
    /* Get the interpolation data */
    cmc_t8_interpolation_data_t interpolation_data = static_cast<cmc_t8_interpolation_data_t>(t8_forest_get_user_data(forest_new));
    
    if (refine == 0)
    {
        /* If the element stays the same */
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
        {
            /* In case a 'One for All' compression is chosen */
            for (size_t var_id{0}; var_id < interpolation_data->t8_data->vars.size(); ++var_id)
            {
                /* If no coarsening has happened, just copy the previous element's value */
                memcpy(static_cast<std::byte*>(interpolation_data->t8_data->vars[var_id]->get_initial_data_new_ptr()) + first_incoming * interpolation_data->t8_data->vars[var_id]->get_data_size(),
                       static_cast<std::byte*>(interpolation_data->t8_data->vars[var_id]->get_initial_data_ptr()) + first_outgoing * interpolation_data->t8_data->vars[var_id]->get_data_size(),
                       interpolation_data->t8_data->vars[var_id]->get_data_size());

            }
        }
        else
        {
            /* In case a 'One for One' compression is chosen */
            /* If no coarsening has happened, just copy the previous element's value */
            memcpy(static_cast<std::byte*>(interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_initial_data_new_ptr()) + first_incoming * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                   static_cast<std::byte*>(interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_initial_data_ptr()) + first_outgoing * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                   interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size());

        }
    } else if (refine == -1)
    {
        /* If coarsening has been introduced */
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL)
        {
            /* In case a 'One for All' compression is chosen */
            for (size_t var_id{0}; var_id < interpolation_data->t8_data->vars.size(); ++var_id)
            {
                memcpy(static_cast<std::byte*>(interpolation_data->t8_data->vars[var_id]->get_initial_data_new_ptr()) + first_incoming * interpolation_data->t8_data->vars[var_id]->get_data_size(),
                       static_cast<std::byte*>((*(interpolation_data->adapted_data))[var_id].get_initial_data_ptr()) + interpolation_data->coarsening_counter * interpolation_data->t8_data->vars[var_id]->get_data_size(),
                       interpolation_data->t8_data->vars[var_id]->get_data_size());
            }
            /* Increase the counter */
            ++(interpolation_data->coarsening_counter);
        }
        else
        {
            /* In case a 'One for One' compression is chosen */
            memcpy(static_cast<std::byte*>(interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_initial_data_new_ptr()) + first_incoming * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                   static_cast<std::byte*>((*(interpolation_data->adapted_data))[0].get_initial_data_ptr()) + interpolation_data->coarsening_counter * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                   interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size());
            /* Increase the counter */
            ++(interpolation_data->coarsening_counter);
        }
    }  else
    {
        cmc_err_msg("Refinement should have not happened.");
    }
    #endif
}

#endif
