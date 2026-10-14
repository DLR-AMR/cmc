#ifndef UTILITIES_CMC_IDW_INTERPOLATION_UTIL_HXX
#define UTILITIES_CMC_IDW_INTERPOLATION_UTIL_HXX

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

namespace cmc::util::idw
{

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


}

#endif /* !UTILITIES_CMC_IDW_INTERPOLATION_UTIL_HXX */
