#include "utilities/cmc_cart_coordinate.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.h"

namespace cmc
{

DomainIndex TransformCartesianCoordinateToLinearIndexLonLat(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lat_length = domain.GetDimensionLength(Dimension::Lat);
    const DomainIndex lon_offset = cartesian_coords[0] * lat_length;
    const DomainIndex lat_offset = cartesian_coords[1];
    return lon_offset + lat_offset;
}

DomainIndex TransformCartesianCoordinateToLinearIndexLatLon(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lon_length = domain.GetDimensionLength(Dimension::Lon);
    const DomainIndex lon_offset = cartesian_coords[0];
    const DomainIndex lat_offset = cartesian_coords[1] * lon_length;
    return lon_offset + lat_offset;
}

DomainIndex TransformCartesianCoordinateToLinearIndexLatLev(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lev_length = domain.GetDimensionLength(Dimension::Lev);
    const DomainIndex lat_offset = cartesian_coords[0] * lev_length;
    const DomainIndex lev_offset = cartesian_coords[1];
    return lat_offset + lev_offset;
}

DomainIndex TransformCartesianCoordinateToLinearIndexLevLat(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lat_length = domain.GetDimensionLength(Dimension::Lat);
    const DomainIndex lat_offset = cartesian_coords[0];
    const DomainIndex lev_offset = cartesian_coords[1] * lat_length;
    return lat_offset + lev_offset;
}

DomainIndex TransformCartesianCoordinateToLinearIndexLonLev(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lev_length = domain.GetDimensionLength(Dimension::Lev);
    const DomainIndex lon_offset = cartesian_coords[0] * lev_length;
    const DomainIndex lev_offset = cartesian_coords[1];
    return lon_offset + lev_offset;
}

DomainIndex TransformCartesianCoordinateToLinearIndexLevLon(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lon_length = domain.GetDimensionLength(Dimension::Lon);
    const DomainIndex lon_offset = cartesian_coords[0];
    const DomainIndex lev_offset = cartesian_coords[1] * lon_length;
    return lon_offset + lev_offset;
}

DomainIndex TransformCartesianCoordinateToLinearIndexLonLatLev(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lat_length = domain.GetDimensionLength(Dimension::Lat);
    const DomainIndex lev_length = domain.GetDimensionLength(Dimension::Lev);
    const DomainIndex lon_offset = cartesian_coords[0] * lat_length * lev_length;
    const DomainIndex lat_offset = cartesian_coords[1] * lev_length;
    const DomainIndex lev_offset = cartesian_coords[2];
    return lon_offset + lat_offset + lev_offset;
}

DomainIndex TransformCartesianCoordinateToLinearIndexLevLonLat(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lon_length = domain.GetDimensionLength(Dimension::Lon);
    const DomainIndex lat_length = domain.GetDimensionLength(Dimension::Lat);
    const DomainIndex lon_offset = cartesian_coords[0] * lat_length;
    const DomainIndex lat_offset = cartesian_coords[1];
    const DomainIndex lev_offset = cartesian_coords[2] * lon_length * lat_length;
    return lon_offset + lat_offset + lev_offset;
}

DomainIndex TransformCartesianCoordinateToLinearIndexLonLevLat(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lat_length = domain.GetDimensionLength(Dimension::Lat);
    const DomainIndex lev_length = domain.GetDimensionLength(Dimension::Lev);
    const DomainIndex lon_offset = cartesian_coords[0] * lat_length * lev_length;
    const DomainIndex lat_offset = cartesian_coords[1];
    const DomainIndex lev_offset = cartesian_coords[2] * lat_length;
    return lon_offset + lat_offset + lev_offset;
}

DomainIndex TransformCartesianCoordinateToLinearIndexLevLatLon(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lon_length = domain.GetDimensionLength(Dimension::Lon);
    const DomainIndex lat_length = domain.GetDimensionLength(Dimension::Lat);
    const DomainIndex lon_offset = cartesian_coords[0];
    const DomainIndex lat_offset = cartesian_coords[1] * lon_length;
    const DomainIndex lev_offset = cartesian_coords[2] * lon_length * lat_length;
    return lon_offset + lat_offset + lev_offset;
}

DomainIndex TransformCartesianCoordinateToLinearIndexLatLevLon(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lon_length = domain.GetDimensionLength(Dimension::Lon);
    const DomainIndex lev_length = domain.GetDimensionLength(Dimension::Lev);
    const DomainIndex lon_offset = cartesian_coords[0];
    const DomainIndex lat_offset = cartesian_coords[1] * lon_length * lev_length;
    const DomainIndex lev_offset = cartesian_coords[2] * lon_length;
    return lon_offset + lat_offset + lev_offset;
}

DomainIndex TransformCartesianCoordinateToLinearIndexLatLonLev(const std::vector<DomainIndex>& cartesian_coords, const GeoDomain& domain)
{
    const DomainIndex lon_length = domain.GetDimensionLength(Dimension::Lon);
    const DomainIndex lev_length = domain.GetDimensionLength(Dimension::Lev);
    const DomainIndex lon_offset = cartesian_coords[0] * lev_length;
    const DomainIndex lat_offset = cartesian_coords[1] * lon_length * lev_length;
    const DomainIndex lev_offset = cartesian_coords[2];
    return lon_offset + lat_offset + lev_offset;
}

[[maybe_unused]] DomainIndex TransformCartesianCoordinateToLinearIndexError([[maybe_unused]] const std::vector<DomainIndex>& cartesian_coords, [[maybe_unused]] const GeoDomain& domain)
{
    return CMC_ERR;
}

TransformCartesianToIndexAccessorFn
GetCartesianCoordsToLinearIndexFunction(const DataLayout layout)
{
    switch (layout)
    {
        case DataLayout::Lon_Lat:
            return TransformCartesianCoordinateToLinearIndexLonLat;
        break;
        case DataLayout::Lat_Lon:
            return TransformCartesianCoordinateToLinearIndexLatLon;
        break;
        case DataLayout::Lat_Lev:
            return TransformCartesianCoordinateToLinearIndexLatLev;
        break;
        case DataLayout::Lev_Lat:
            return TransformCartesianCoordinateToLinearIndexLevLat;
        break;
        case DataLayout::Lon_Lev:
            return TransformCartesianCoordinateToLinearIndexLonLev;
        break;
        case DataLayout::Lev_Lon:
            return TransformCartesianCoordinateToLinearIndexLevLon;
        break;
        case DataLayout::Lon_Lat_Lev:
            return TransformCartesianCoordinateToLinearIndexLonLatLev;
        break;
        case DataLayout::Lev_Lon_Lat:
            return TransformCartesianCoordinateToLinearIndexLevLonLat;
        break;
        case DataLayout::Lon_Lev_Lat:
            return TransformCartesianCoordinateToLinearIndexLonLevLat;
        break;
        case DataLayout::Lev_Lat_Lon:
            return TransformCartesianCoordinateToLinearIndexLevLatLon;
        break;
        case DataLayout::Lat_Lev_Lon:
            return TransformCartesianCoordinateToLinearIndexLatLevLon;
        break;
        case DataLayout::Lat_Lon_Lev:
            return TransformCartesianCoordinateToLinearIndexLatLonLev;
        break;
        default :
            cmc_err_msg("The variable contains an undefined data layout.");
            return TransformCartesianCoordinateToLinearIndexError;
    }
}

}
