#ifndef CMC_T8_ADAPT_CALLBACKS_H
#define CMC_T8_ADAPT_CALLBACKS_H

#ifdef CMC_WITH_T8CODE

#include "cmc_t8code.h"



t8_locidx_t
cmc_t8_adapt_coarsen_geo_mesh_callback (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          t8_locidx_t lelement_id,
                          t8_eclass_scheme_c * ts,
                          const int is_family,
                          const int num_elements,
                          t8_element_t * elements[]);

t8_locidx_t
cmc_t8_adapt_callback_refine_to_initial_lvl (t8_forest_t forest,
                                             t8_forest_t forest_from,
                                             int which_tree,
                                             int lelement_id,
                                             t8_eclass_scheme_c * ts,
                                             const int is_family,
                                             const int num_elements,
                                             t8_element_t * elements[]);

t8_locidx_t
cmc_t8_adapt_callback_coarsen_exclude_area(t8_forest_t forest,
                                           t8_forest_t forest_from,
                                           int which_tree,
                                           int lelement_id,
                                           t8_eclass_scheme_c * ts,
                                           const int is_family,
                                           const int num_elements,
                                           t8_element_t * elements[]);

t8_locidx_t
cmc_t8_adapt_callback_coarsen_error_threshold (t8_forest_t forest,
                                               t8_forest_t forest_from,
                                               int which_tree,
                                               int lelement_id,
                                               t8_eclass_scheme_c * ts,
                                               const int is_family,
                                               const int num_elements,
                                               t8_element_t * elements[]);

t8_locidx_t
cmc_t8_adapt_callback_coarsen_combined_criteria (t8_forest_t forest,
                                                 t8_forest_t forest_from,
                                                 int which_tree,
                                                 int lelement_id,
                                                 t8_eclass_scheme_c * ts,
                                                 const int is_family,
                                                 const int num_elements,
                                                 t8_element_t * elements[]);

t8_locidx_t
cmc_t8_adapt_callback_coarsen_error_threshold_parallel (t8_forest_t forest,
                                                        t8_forest_t forest_from,
                                                        int which_tree,
                                                        int lelement_id,
                                                        t8_eclass_scheme_c * ts,
                                                        const int is_family,
                                                        const int num_elements,
                                                        t8_element_t * elements[]);


#endif
#endif /* CMC_T8_ADAPT_CALLBACKS_H */
