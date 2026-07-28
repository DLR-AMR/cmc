#ifndef CMC_AMR_LOSSY_PAR_MULTI_RES_IDW_INTERPOLATION_UTIL_HXX
#define CMC_AMR_LOSSY_PAR_MULTI_RES_IDW_INTERPOLATION_UTIL_HXX

#include "cmc.hxx"
#include "amr/lossy/cmc_par_multi_res_extraction_util.hxx"

#include <cstdint>
#include <string>
#include <type_traits>
#include <limits>
#include <utility>
#include <array>
#include <vector>
#include <climits>
#include <algorithm>
#include <numeric>
#include <concepts>

namespace cmc::par::lossy::idw::util
{

    #if 0

template<typename T>
concept FloatType = (std::is_floating_point_v<T> && std::is_arithmetic_v<T>);

constexpr int kIDWExponent = 4;

/**
 *  A constexpr for-loop construct
 */
template <auto START, auto END, auto INCREMENT, class FUNC>
constexpr void
for_constexpr(FUNC&& f)
{
    if constexpr (START < END)
    {
        f(std::integral_constant<decltype(START), START>());
        for_constexpr<START + INCREMENT, END, INCREMENT>(f);
    }
}

template <FloatType T, Dimension DIM>
struct Coordinate
{
    Coordinate(const std::array<T, DIM>& coords)
    : coordinates(coords) {}
    Coordinate(std::array<T, DIM>&& coords)
    : coordinates(std::move(coords)) {}

    std::array<T, DIM> coordinates{};
}

template <FloatType T, Dimension DIM>
struct PointData
{
    T data{};
    std::array<T, DIM> coordinates{};
}

template <FloatType T, Dimension DIM>
inline
T
ComputeSQDist(const PointData<T, DIM>& points, const Coordinate<T, DIM>& coords);

template <FloatType T>
inline
T
ComputeSQDist(const PointData<T, 1>& point, const Coordinate<T, 1>& coords)
{
    const T sq_dist = (point_coords[0] - coords.coordinates[0]) * (point_coords[0] - coords.coordinates[0]);
    return sq_dist
}

template <FloatType T>
inline
T
ComputeSQDist(const PointData<T, 2>& point, const Coordinate<T, 2>& coords)
{
    const T sq_dist = (point_coords[0] - coords.coordinates[0]) * (point_coords[0] - coords.coordinates[0]) +
                      (point_coords[1] - coords.coordinates[1]) * (point_coords[1] - coords.coordinates[1]);
    return sq_dist
}

template <FloatType T>
inline
T
ComputeSQDist(const PointData<T, 3>& point, const Coordinate<T, 3>& coords)
{
    const T sq_dist = (point_coords[0] - coords.coordinates[0]) * (point_coords[0] - coords.coordinates[0]) +
                      (point_coords[1] - coords.coordinates[1]) * (point_coords[1] - coords.coordinates[1]) +
                      (point_coords[2] - coords.coordinates[2]) * (point_coords[2] - coords.coordinates[2]);
    return sq_dist
}

template <FloatType T>
inline
T
ComputeSQDist(const PointData<T, 4>& point, const Coordinate<T, 4>& coords)
{
    const T sq_dist = (point_coords[0] - coords.coordinates[0]) * (point_coords[0] - coords.coordinates[0]) +
                      (point_coords[1] - coords.coordinates[1]) * (point_coords[1] - coords.coordinates[1]) +
                      (point_coords[2] - coords.coordinates[2]) * (point_coords[2] - coords.coordinates[2]) +
                      (point_coords[3] - coords.coordinates[3]) * (point_coords[3] - coords.coordinates[3]);
    return sq_dist
}

template <FloatType T, Dimension DIM>
inline
std::vector<T>
ComputeSQDistances(const std::vector<PointData<T, DIM>>& points, const Coordinate<T, DIM>& coords)
{
    std::vector<T> dist(points.size());

    /* Compute the distances to all points */
    std::transform(std::execution::par_unseq, points.cbegin(), points.cend(),
                   dist.begin(), [&coords](const auto& point){
                                    return ComputeSQDist<T, DIM>(point, coords);
                                  });
}

template <FloatType T>
constexpr
ComputeExpDistances(const T dist)
{
    T exp_dist = dist;
    for_constexpr<0, kIDWExponent, 1>([&exp_dist](const auto& idx){
        exp_dist *= dist;
    });
    return exp_dist;
}

template <int32_t DIM>
requires Dimension<DIM>
T
ComputePrediction(const std::vector<PointData<T, DIM>>& points, const Coordinate<T, DIM>& evaluation_coords)
{
    /* Compute the distances to the control points */
    std::vector<T> weights = ComputeSQDistances<T, DIM>(points, evaluation_coords);

    /* Compute exponential distances */
    std::transform(std::execution::par_unseq, weights.begin(), weights.end(), weights.begin(), [](T& weight){
        weight = 1.0 / ComputeExpDistances<T>(weight);
    });

    /* Normalize the weights */
    const T weigths_sum = std::reduce(std::execution::par_unseq, weights.begin(), weights.end());

    std::transform(std::execution::par_unseq, weights.begin(), weights.end(), weights.begin(), [](T& weight){
        weight = weight / weigths_sum;
    });

    /* Compute prediction */
    const T prediction = std::transform_reduce(std::execution::par_unseq, weights.cbegin(), weights.cend(), points.cbegin(), static_cast<T>(0), std::plus{},
                                               [](const T& weight, const PointData<T, DIM>& point_data){
                                                return weight * (point_data.data);
                                               });

    return prediction;
}

#else



template<typename T>
concept ArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T>);

template<int32_t DIM>
concept Dimension = (DIM >= 1 && DIM <= 4);

using FloatType = float;

constexpr int kIDWExponent = 4;

/**
 *  A constexpr for-loop construct
 */
template <auto START, auto END, auto INCREMENT, class FUNC>
constexpr void
for_constexpr(FUNC&& f)
{
    if constexpr (START < END)
    {
        f(std::integral_constant<decltype(START), START>());
        for_constexpr<START + INCREMENT, END, INCREMENT>(f);
    }
}

template <int32_t DIM>
requires Dimension<DIM>
struct Coordinate
{
    Coordinate() {coordinates.fill(0.0);}
    Coordinate(const std::array<FloatType, DIM>& coords)
    : coordinates(coords) {}
    Coordinate(std::array<FloatType, DIM>&& coords)
    : coordinates(std::move(coords)) {}
    
    std::array<FloatType, DIM> coordinates{};

    FloatType& operator[](std::size_t idx) {return coordinates[idx];}
    const FloatType& operator[](std::size_t idx) const {return coordinates[idx];}
};

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
struct PointData
{
    T data{};
    Coordinate<DIM> coordinates;
};

template <int32_t DIM>
inline
FloatType
ComputeSQDist(const Coordinate<DIM>& points, const Coordinate<DIM>& coords);

template <>
inline
FloatType
ComputeSQDist<1>(const Coordinate<1>& point_coords, const Coordinate<1>& coords)
{
    const FloatType sq_dist = (point_coords[0] - coords.coordinates[0]) * (point_coords[0] - coords.coordinates[0]);
    return sq_dist;
}

template <>
inline
FloatType
ComputeSQDist<2>(const Coordinate<2>& point_coords, const Coordinate<2>& coords)
{
    const FloatType sq_dist = (point_coords[0] - coords.coordinates[0]) * (point_coords[0] - coords.coordinates[0]) +
                      (point_coords[1] - coords.coordinates[1]) * (point_coords[1] - coords.coordinates[1]);
    return sq_dist;
}

template <>
inline
FloatType
ComputeSQDist<3>(const Coordinate<3>& point_coords, const Coordinate<3>& coords)
{
    const FloatType sq_dist = (point_coords[0] - coords.coordinates[0]) * (point_coords[0] - coords.coordinates[0]) +
                      (point_coords[1] - coords.coordinates[1]) * (point_coords[1] - coords.coordinates[1]) +
                      (point_coords[2] - coords.coordinates[2]) * (point_coords[2] - coords.coordinates[2]);
    return sq_dist;
}

template <>
inline
FloatType
ComputeSQDist<4>(const Coordinate<4>& point_coords, const Coordinate<4>& coords)
{
    const FloatType sq_dist = (point_coords[0] - coords.coordinates[0]) * (point_coords[0] - coords.coordinates[0]) +
                      (point_coords[1] - coords.coordinates[1]) * (point_coords[1] - coords.coordinates[1]) +
                      (point_coords[2] - coords.coordinates[2]) * (point_coords[2] - coords.coordinates[2]) +
                      (point_coords[3] - coords.coordinates[3]) * (point_coords[3] - coords.coordinates[3]);
    return sq_dist;
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline
std::vector<FloatType>
ComputeSQDistances(const std::vector<PointData<T, DIM>>& points, const Coordinate<DIM>& coords)
{
    std::vector<FloatType> dist(points.size());

    /* Compute the distances to all points */
    std::transform(std::execution::par_unseq, points.cbegin(), points.cend(),
                   dist.begin(), [&coords](const auto& point){
                                    return ComputeSQDist<DIM>(point.coordinates, coords);
                                  });
    return dist;
}

constexpr
inline FloatType
ComputeExpDistances(const FloatType dist)
{
    FloatType exp_dist = dist;
    for_constexpr<0, kIDWExponent, 1>([&]([[maybe_unused]] const auto& idx){
        exp_dist *= dist;
    });
    return exp_dist;
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
T
ComputePrediction(const std::vector<PointData<T, DIM>>& points, const Coordinate<DIM>& evaluation_coords)
{
    /* Compute the distances to the control points */
    std::vector<FloatType> weights = ComputeSQDistances<T, DIM>(points, evaluation_coords);

    /* Compute exponential distances */
    std::transform(std::execution::par_unseq, weights.begin(), weights.end(), weights.begin(), [](const FloatType& weight){
        return 1.0 / ComputeExpDistances(weight);
    });

    /* Normalize the weights */
    const T weigths_sum = std::reduce(std::execution::par_unseq, weights.begin(), weights.end());

    std::transform(std::execution::par_unseq, weights.begin(), weights.end(), weights.begin(), [&weigths_sum](const T& weight){
        return weight / weigths_sum;
    });

    /* Compute prediction */
    const T prediction = static_cast<T>(std::transform_reduce(std::execution::par_unseq, weights.begin(), weights.end(), points.cbegin(), static_cast<FloatType>(0), std::plus<FloatType>(),
                                               [](const FloatType& weight, const PointData<T, DIM>& point_data){
                                                return weight * static_cast<FloatType>(point_data.data);
                                               }));

    return prediction;
}

template<int32_t DIM>
requires Dimension<DIM>
std::array<FloatType, DIM> 
FillElementIDWCoordinates(t8_forest_t mesh, t8_locidx_t tree_idx, const t8_element_t* elem)
{
    std::array<double, 3> coords{};
    //t8_forest_element_coordinate (mesh, tree_idx, elem, 0, coords.data());
    t8_forest_element_linear_centroid (mesh, tree_idx, elem, coords.data());
    std::array<FloatType, DIM> coordinates;
    for_constexpr<0, DIM, 1>([&](auto IDX)
    {
        coordinates[IDX] = static_cast<FloatType>(coords[IDX]);
    });

    return coordinates;
}


#endif

}

#endif /* !CMC_AMR_LOSSY_PAR_MULTI_RES_EXTRACTION_UTIL_HXX */
