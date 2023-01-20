#ifndef CMC_T8CODE_GEO_MESH_H
#define CMC_T8CODE_GEO_MESH_H

#include "cmc.h"
#include "utilities/cmc_geo_util.h"
#include "cmc_t8code_data.h"

bool
compare_geo_domain_equality_of_data_layouts(const DATA_LAYOUT first, const DATA_LAYOUT second);

void
cmc_t8_create_enclosing_geo_mesh(cmc_t8_data& t8_data);

int
cmc_t8_elem_inside_geo_mesh(const t8_element_t* element, t8_eclass_scheme_c* ts, const cmc_t8_data& t8_data, const int var_id);

int
cmc_t8_elem_inside_geo_domain(const t8_element_t* element, t8_eclass_scheme_c* ts, const cmc_t8_data& t8_data, const int var_id,
                              const int lat_min, const int lat_max, const int lon_min, const int lon_max, const int lev_min, const int lev_max);

#endif /*CMC_T8CODE_GEO_MESH_H */