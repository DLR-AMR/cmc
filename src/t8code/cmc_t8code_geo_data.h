#ifndef CMC_T8CODE_GEO_DATA_H
#define CMC_T8CODE_GEO_DATA_H

#include "cmc.h"
#include "utilities/cmc_util.h"
#include "cmc_t8code_data.h"

void
cmc_t8_apply_zcurve_ordering(cmc_t8_data& t8_data, const int var_id);

void
cmc_t8_apply_offset_and_scaling(cmc_t8_data_t t8_data, const int var_id);

void
cmc_geo_data_transform_3d_var_to_2d(cmc_t8_data_t t8_data, const int var_id, const DATA_LAYOUT preferred_layout);

void
cmc_t8_geo_data_set_error_criterium(cmc_t8_data_t t8_data, const double maximum_error_tolerance);

void
cmc_t8_geo_data_set_exclude_area(cmc_t8_data_t t8_data, const CMC_COORD_IDS coord_id, const cmc_universal_type_t& start_value, const cmc_universal_type_t& end_value);

void
cmc_t8_coarsen_data(cmc_t8_data_t t8_data, t8_forest_adapt_t adapt_function, t8_forest_replace_t interpolation_function);

void
cmc_t8_refine_to_initial_level(cmc_t8_data_t t8_data);

void
cmc_t8_geo_data_distribute_initially(cmc_t8_data_t t8_data);

void
cmc_t8_geo_data_distribute_and_apply_ordering(cmc_t8_data_t t8_data);

void
cmc_t8_distribute_data(cmc_t8_data_t t8_data);

#endif /* CMC_T8CODE_GEO_DATA_H */