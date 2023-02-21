#ifndef CMC_T8CODE_GEO_DATA_H
#define CMC_T8CODE_GEO_DATA_H

#include "cmc.h"
#include "cmc_t8code_data.h"

#if 0
#define CMC_APPLY_ZCURVE_TO_ALL_VARS -1
#define CMC_APPLY_OFFSET_AND_SCALING_TO_ALL_VARS -1
#define CMC_APPLY_TRANSFORMATION_3D_TO_2D_TO_ALL_VARS -1
#endif

void
cmc_t8_apply_zcurve_ordering(cmc_t8_data& t8_data, const int var_id);

void
cmc_t8_apply_offset_and_scaling(cmc_t8_data_t t8_data, const int var_id);

void
cmc_geo_data_transform_3d_var_to_2d(cmc_t8_data_t t8_data, const int var_id, const DATA_LAYOUT preferred_layout);

void
cmc_t8_geo_data_set_error_criterium(cmc_t8_data_t t8_data, const double maximim_error_tolerance);

void
cmc_t8_coarsen_data(cmc_t8_data_t t8_data, t8_forest_adapt_t adapt_function, t8_forest_replace_t interpolation_function);

void
cmc_t8_refine_to_initial_level(cmc_t8_data_t t8_data);

void
cmc_t8_distribute_data(cmc_t8_data_t t8_data);

#endif /* CMC_T8CODE_GEO_DATA_H */