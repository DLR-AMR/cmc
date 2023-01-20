#ifndef CMC_T8_REPLACE_CALLBACKS_H
#define CMC_T8_REPLACE_CALLBACKS_H

#include "cmc_t8code.h"

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
                                                     
#endif /* CMC_T8_REPLACE_CALLBACKS_H */