#include "amr/lossy/cmc_adaptive_coarsening.hxx"

namespace cmc::lossy
{

/* Explicit template instantiations */
template class DefaultAdaptData<int8_t>;
template class DefaultAdaptData<char>;
template class DefaultAdaptData<int16_t>;
template class DefaultAdaptData<int32_t>;
template class DefaultAdaptData<float>;
template class DefaultAdaptData<double>;
template class DefaultAdaptData<uint8_t>;
template class DefaultAdaptData<uint16_t>;
template class DefaultAdaptData<uint32_t>;
template class DefaultAdaptData<int64_t>;
template class DefaultAdaptData<uint64_t>;

template class CompressionVariable<int8_t>;
template class CompressionVariable<char>;
template class CompressionVariable<int16_t>;
template class CompressionVariable<int32_t>;
template class CompressionVariable<float>;
template class CompressionVariable<double>;
template class CompressionVariable<uint8_t>;
template class CompressionVariable<uint16_t>;
template class CompressionVariable<uint32_t>;
template class CompressionVariable<int64_t>;
template class CompressionVariable<uint64_t>;

template t8_locidx_t DefaultAdaptiveCoarsening<int8_t> (t8_forest_t forest, t8_forest_t forest_from,
                                                            t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                            const t8_scheme_c * ts, const int is_family,
                                                            const int num_elements, t8_element_t * elements[]);
template t8_locidx_t DefaultAdaptiveCoarsening<char> (t8_forest_t forest, t8_forest_t forest_from,
                                                            t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                            const t8_scheme_c * ts, const int is_family,
                                                            const int num_elements, t8_element_t * elements[]);
template t8_locidx_t DefaultAdaptiveCoarsening<int16_t> (t8_forest_t forest, t8_forest_t forest_from,
                                                            t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                            const t8_scheme_c * ts, const int is_family,
                                                            const int num_elements, t8_element_t * elements[]);
template t8_locidx_t DefaultAdaptiveCoarsening<int32_t> (t8_forest_t forest, t8_forest_t forest_from,
                                                            t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                            const t8_scheme_c * ts, const int is_family,
                                                            const int num_elements, t8_element_t * elements[]);
template t8_locidx_t DefaultAdaptiveCoarsening<float> (t8_forest_t forest, t8_forest_t forest_from,
                                                            t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                            const t8_scheme_c * ts, const int is_family,
                                                            const int num_elements, t8_element_t * elements[]);
template t8_locidx_t DefaultAdaptiveCoarsening<double> (t8_forest_t forest, t8_forest_t forest_from,
                                                            t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                            const t8_scheme_c * ts, const int is_family,
                                                            const int num_elements, t8_element_t * elements[]);
template t8_locidx_t DefaultAdaptiveCoarsening<uint8_t> (t8_forest_t forest, t8_forest_t forest_from,
                                                            t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                            const t8_scheme_c * ts, const int is_family,
                                                            const int num_elements, t8_element_t * elements[]);
template t8_locidx_t DefaultAdaptiveCoarsening<uint16_t> (t8_forest_t forest, t8_forest_t forest_from,
                                                            t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                            const t8_scheme_c * ts, const int is_family,
                                                            const int num_elements, t8_element_t * elements[]);
template t8_locidx_t DefaultAdaptiveCoarsening<uint32_t> (t8_forest_t forest, t8_forest_t forest_from,
                                                            t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                            const t8_scheme_c * ts, const int is_family,
                                                            const int num_elements, t8_element_t * elements[]);
template t8_locidx_t DefaultAdaptiveCoarsening<int64_t> (t8_forest_t forest, t8_forest_t forest_from,
                                                            t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                            const t8_scheme_c * ts, const int is_family,
                                                            const int num_elements, t8_element_t * elements[]);
template t8_locidx_t DefaultAdaptiveCoarsening<uint64_t> (t8_forest_t forest, t8_forest_t forest_from,
                                                            t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                            const t8_scheme_c * ts, const int is_family,
                                                            const int num_elements, t8_element_t * elements[]);

}
