#ifndef CMC_T8_REPLACE_CALLBACKS_H
#define CMC_T8_REPLACE_CALLBACKS_H

#ifdef CMC_WITH_T8CODE

#include "cmc_t8code.h"

/**
 * @brief A general interpolation function to use during the coarsening/compression step as needed by t8code.
 *        This interpolation function uses one of the below 'cmc_t8_forest_interpolation_t' functions during the interpolation.
 */
void
cmc_t8_general_interpolation_during_compression (t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                                 t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                                 int num_incoming, t8_locidx_t first_incoming);

/**
 * @brief A wrapper function for the standard arithmetic mean interpolation during the compression step. Calling this function
 *        returns a functor which provides the standard arithmetic mean function when called
 * 
 * @return cmc_t8_forest_interpolate_t Returns a functor which holds the actual arithmetic mean function
 */
cmc_t8_forest_interpolate_t
cmc_t8_compression_interpolation_standard_mean();

void
cmc_t8_geo_data_interpolate_std_mean(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                     t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                     int num_incoming, t8_locidx_t first_incoming);

void
cmc_t8_geo_data_interpolate_plain_copy_values(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                              t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                              int num_incoming, t8_locidx_t first_incoming);


void
cmc_t8_geo_data_interpolate_error_threshold_adaption(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
                                                     t8_eclass_scheme_c* ts, int refine, int num_outgoing, t8_locidx_t first_outgoing,
                                                     int num_incoming, t8_locidx_t first_incoming);

#endif                                                   
#endif /* CMC_T8_REPLACE_CALLBACKS_H */
