#ifndef CMC_CART_COORDINATE_HXX
#define CMC_CART_COORDINATE_HXX

#include "utilities/cmc_dimension_interval.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"

#include <array>
#include <vector>
#include <functional>

namespace cmc
{

struct CartesianCoordinate
{
    constexpr CartesianCoordinate(const int64_t x_coordinate = 0, const int64_t y_coordinate = 0, const int64_t z_coordinate = 0, const int64_t time_coordinate = 0)
    : coordinate{x_coordinate,y_coordinate,z_coordinate,time_coordinate}{};

    int64_t GetDimensionCoordinate(const Dimension dimension) const
    {
        return coordinate[dimension];
    };

    std::array<int64_t, Dimension::NumCoordinates> coordinate{0,0,0,0};
};

using TransformCartesianToIndexAccessorFn = std::function<DomainIndex(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)>;

TransformCartesianToIndexAccessorFn
GetCartesianCoordsToLinearIndexFunction (const DataLayout layout);

DomainIndex TransformCartesianCoordinateToLinearIndexLonLat(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

DomainIndex TransformCartesianCoordinateToLinearIndexLatLon(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

DomainIndex TransformCartesianCoordinateToLinearIndexLatLev(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

DomainIndex TransformCartesianCoordinateToLinearIndexLevLat(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

DomainIndex TransformCartesianCoordinateToLinearIndexLonLev(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

DomainIndex TransformCartesianCoordinateToLinearIndexLevLon(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

DomainIndex TransformCartesianCoordinateToLinearIndexLonLatLev(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

DomainIndex TransformCartesianCoordinateToLinearIndexLevLonLat(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

DomainIndex TransformCartesianCoordinateToLinearIndexLonLevLat(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

DomainIndex TransformCartesianCoordinateToLinearIndexLevLatLon(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

DomainIndex TransformCartesianCoordinateToLinearIndexLatLevLon(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

DomainIndex TransformCartesianCoordinateToLinearIndexLatLonLev(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain);

[[maybe_unused]] DomainIndex TransformCartesianCoordinateToLinearIndexError([[maybe_unused]] const std::vector<DomainIndex>& cartesian_coords, [[maybe_unused]] const GeoDomain& domain);

}

#endif /* !CMC_CART_COORDINATE_HXX */
