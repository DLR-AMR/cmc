#include "t8code/cmc_t8_replace_callbacks.h"
#include "t8code/cmc_t8code_data.hxx"
#include "utilities/cmc_container.h"
#include "utilities/cmc_log_functions.h"

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
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
            interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
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
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
            interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
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
    
    /* Check if elements get replaced in the new adapted forest (num_incoming < num_outgoing) */
    if (refine == 1)
    {
        /* If a refinement was introduced */
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
            interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
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
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
            interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
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
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
            interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
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
            /* If no coarsening has happened, just copy the previous element's value */
            memcpy(static_cast<std::byte*>(interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_initial_data_new_ptr()) + first_incoming * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                   static_cast<std::byte*>(interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_initial_data_ptr()) + first_outgoing * interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size(),
                   interpolation_data->t8_data->vars[interpolation_data->current_var_id]->get_data_size());

        }
    } else if (refine == -1)
    {
        /* If coarsening has been introduced */
        if (interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_2D ||
            interpolation_data->t8_data->compression_mode == CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL_3D)
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
