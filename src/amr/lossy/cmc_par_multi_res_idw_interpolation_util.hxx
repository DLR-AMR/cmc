#ifndef CMC_AMR_LOSSY_PAR_MULTI_RES_IDW_INTERPOLATION_UTIL_HXX
#define CMC_AMR_LOSSY_PAR_MULTI_RES_IDW_INTERPOLATION_UTIL_HXX

#include "cmc.hxx"
#include "utilities/cmc_idw_interpolation_util.hxx"
#include "amr/lossy/cmc_par_multi_res_extraction_util.hxx"

namespace cmc::par::lossy::idw::util
{

template <cmc::util::idw::ArithmeticType T, int32_t DIM>
requires cmc::util::idw::Dimension<DIM>
T
ComputePrediction(const std::vector<cmc::util::idw::PointData<T, DIM>>& points, const cmc::util::idw::Coordinate<DIM>& evaluation_coords)
{
    /* Compute the distances to the control points */
    std::vector<cmc::util::idw::FloatType> weights = cmc::util::idw::ComputeSQDistances<T, DIM>(points, evaluation_coords);

    /* Compute exponential distances */
    std::transform(std::execution::par_unseq, weights.begin(), weights.end(), weights.begin(), [](const cmc::util::idw::FloatType& weight){
        return 1.0 / cmc::util::idw::ComputeExpDistances(weight);
    });

    /* Normalize the weights */
    const T weigths_sum = std::reduce(std::execution::par_unseq, weights.begin(), weights.end());

    std::transform(std::execution::par_unseq, weights.begin(), weights.end(), weights.begin(), [&weigths_sum](const T& weight){
        return weight / weigths_sum;
    });

    /* Compute prediction */
    const T prediction = static_cast<T>(std::transform_reduce(std::execution::par_unseq, weights.begin(), weights.end(), points.cbegin(), static_cast<cmc::util::idw::FloatType>(0), std::plus<cmc::util::idw::FloatType>(),
                                               [](const cmc::util::idw::FloatType& weight, const cmc::util::idw::PointData<T, DIM>& point_data){
                                                return weight * static_cast<cmc::util::idw::FloatType>(point_data.data);
                                               }));

    return prediction;
}

template<int32_t DIM>
requires cmc::util::idw::Dimension<DIM>
std::array<cmc::util::idw::FloatType, DIM> 
FillElementChildIDWCoordinates(t8_forest_t mesh, t8_locidx_t tree_idx, const t8_element_t* elem)
{
    std::array<double, 3> coords{};
    //t8_forest_element_coordinate (mesh, tree_idx, elem, 0, coords.data());
    t8_forest_element_linear_centroid (mesh, tree_idx, elem, coords.data());
    std::array<cmc::util::idw::FloatType, DIM> coordinates;
    cmc::util::idw::for_constexpr<0, DIM, 1>([&](auto IDX)
    {
        coordinates[IDX] = static_cast<cmc::util::idw::FloatType>(coords[IDX]);
    });

    return coordinates;
}

template<int32_t DIM>
requires cmc::util::idw::Dimension<DIM>
std::array<cmc::util::idw::FloatType, DIM> 
FillElementNeighborIDWCoordinates(t8_forest_t mesh, const t8_scheme_c* scheme, const t8_eclass_t tree_class, t8_locidx_t tree_idx, const t8_element_t* elem)
{
    constexpr t8_locidx_t kConstructedChildIDX = 0;

    /* Get the first child element */
    t8_element_t* child_elem{nullptr};
    t8_element_new (scheme, tree_class, 1, &child_elem);
    t8_element_get_child (scheme, tree_class, elem, kConstructedChildIDX, child_elem);

    /* Compute the midpoint of the element */
    std::array<double, 3> coords{};
    t8_forest_element_linear_centroid (mesh, tree_idx, child_elem, coords.data());

    /* Destroy the child element */
    t8_element_destroy(scheme, tree_class, 1, &child_elem);

    /* Return an array of the correct dimensionality */
    std::array<cmc::util::idw::FloatType, DIM> coordinates;
    cmc::util::idw::for_constexpr<0, DIM, 1>([&](auto IDX)
    {
        coordinates[IDX] = static_cast<cmc::util::idw::FloatType>(coords[IDX]);
    });

    return coordinates;
}

}

#endif /* !CMC_AMR_LOSSY_PAR_MULTI_RES_EXTRACTION_UTIL_HXX */
