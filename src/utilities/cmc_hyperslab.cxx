#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_log_functions.h"

#include <numeric>

namespace cmc
{

inline void
Hyperslab::SetupDimension(const DimensionInterval& dimension)
{
    switch (dimension.dim)
    {
        case Dimension::Lon:
            start_indices_[Dimension::Lon] = dimension.start_index;
            count_indices_[Dimension::Lon] = dimension.end_index - dimension.start_index;
        break;
        case Dimension::Lat:
            start_indices_[Dimension::Lat] = dimension.start_index;
            count_indices_[Dimension::Lat] = dimension.end_index - dimension.start_index;
        break;
        case Dimension::Lev:
            start_indices_[Dimension::Lev] = dimension.start_index;
            count_indices_[Dimension::Lev] = dimension.end_index - dimension.start_index;
        break;
        case Dimension::Time:
            start_indices_[Dimension::Time] = dimension.start_index;
            count_indices_[Dimension::Time] = dimension.end_index - dimension.start_index;
        break;
        default:
            cmc_err_msg("The given dimension is not considered.");
    }
}

Hyperslab::Hyperslab(const DimensionInterval& dimension1)
{
    SetupDimension(dimension1);
};
Hyperslab::Hyperslab(const DimensionInterval& dimension1, const DimensionInterval& dimension2)
{
    SetupDimension(dimension1);
    SetupDimension(dimension2);
};
Hyperslab::Hyperslab(const DimensionInterval& dimension1, const DimensionInterval& dimension2,
                     const DimensionInterval& dimension3)
{
    SetupDimension(dimension1);
    SetupDimension(dimension2);
    SetupDimension(dimension3);
};
Hyperslab::Hyperslab(const DimensionInterval& dimension1, const DimensionInterval& dimension2,
                     const DimensionInterval& dimension3, const DimensionInterval& dimension4)
{
    SetupDimension(dimension1);
    SetupDimension(dimension2);
    SetupDimension(dimension3);
    SetupDimension(dimension4);
};

Hyperslab::Hyperslab(const HyperslabIndex lon_start, const HyperslabIndex lon_count,
                     const HyperslabIndex lat_start, const HyperslabIndex lat_count,
                     const HyperslabIndex lev_start, const HyperslabIndex lev_count,
                     const HyperslabIndex time_start, const HyperslabIndex time_count)
: start_indices_{lon_start, lat_start, lev_start, time_start}, count_indices_{lon_count, lat_count, lev_count, time_count}{};

HyperslabIndex
Hyperslab::GetNumberCoordinates() const
{
    const HyperslabIndex product =  std::accumulate(count_indices_.begin(), count_indices_.end(), 1, std::multiplies<HyperslabIndex>());
    return product;
};

HyperslabIndex
Hyperslab::GetNumberCoordinatesWithoutCertainDimension(const Dimension excluded_dimension) const
{
    HyperslabIndex product = 1;
    for (auto iter = count_indices_.begin(); iter != count_indices_.end(); ++iter)
    {
        product *= (std::distance(count_indices_.begin(), iter) != excluded_dimension ? *iter : 1);
    }
    return product;
}

HyperslabIndex
Hyperslab::GetDimensionStart(const Dimension dimension) const
{
    switch (dimension)
    {
        case Dimension::Lon:
            return start_indices_[Dimension::Lon];
        break;
        case Dimension::Lat:
            return start_indices_[Dimension::Lat];
        break;
        case Dimension::Lev:
            return start_indices_[Dimension::Lev];
        break;
        case Dimension::Time:
            return start_indices_[Dimension::Time];
        break;
        default:
            cmc_err_msg("The given dimension is not considered.");
    }
}

int
Hyperslab::GetDimensionality() const
{
    int dim = 0;
    for (auto count_iter = count_indices_.begin(); count_iter != count_indices_.end(); ++count_iter)
    {
        if (*count_iter > 1)
        {
            ++dim;
        }
    }
    return dim;
}

HyperslabIndex
Hyperslab::GetDimensionLength(const Dimension dimension) const
{
    switch (dimension)
    {
        case Dimension::Lon:
            return count_indices_[Dimension::Lon];
        break;
        case Dimension::Lat:
            return count_indices_[Dimension::Lat];
        break;
        case Dimension::Lev:
            return count_indices_[Dimension::Lev];
        break;
        case Dimension::Time:
            return count_indices_[Dimension::Time];
        break;
        default:
            cmc_err_msg("The given dimension is not considered.");
    }
}

/* (currently not!) Latitude Coordinates are flipped, since most applications label (or maybe this is problematic with the decompression) */
void UpdateHyperslabIndexLonLat(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lon_offset = index / lat_length;
    const HyperslabIndex lat_offset = (index - lon_offset * lat_length);
    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lon) + lon_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lat) + lat_offset;
}

void UpdateHyperslabIndexLatLon(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lat_offset = index / lon_length;
    const HyperslabIndex lon_offset = (index - lat_offset * lon_length);
    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lon) + lon_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lat) + lat_offset;
}

void UpdateHyperslabIndexLatLev(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);
    const HyperslabIndex lat_offset = index / lev_length;
    const HyperslabIndex lev_offset = (index - lat_offset * lev_length);
    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lat) + lat_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lev) + lev_offset;
}

void UpdateHyperslabIndexLevLat(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lev_offset = index / lat_length;
    const HyperslabIndex lat_offset = (index - lev_offset * lat_length);
    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lat) + lat_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lev) + lev_offset;
}

void UpdateHyperslabIndexLonLev(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);
    const HyperslabIndex lon_offset = index / lev_length;
    const HyperslabIndex lev_offset = (index - lon_offset * lev_length);
    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lon) + lon_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lev) + lev_offset;
}

void UpdateHyperslabIndexLevLon(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lev_offset = index / lon_length;
    const HyperslabIndex lon_offset = (index - lev_offset * lon_length);
    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lon) + lon_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lat) + lev_offset;
}

void UpdateHyperslabIndexLonLatLev(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);

    const HyperslabIndex lon_offset = index / (lat_length * lev_length);
    const HyperslabIndex lat_offset = (index - lon_offset * lat_length * lev_length) / lev_length; 
    const HyperslabIndex lev_offset = index - lon_offset * lat_length * lev_length - lat_offset * lev_length;

    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lon) + lon_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lat) + lat_offset;
    linear_indices[2] = hyperslab.GetDimensionStart(Dimension::Lev) + lev_offset;
}

void UpdateHyperslabIndexLevLonLat(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);

    const HyperslabIndex lev_offset = index / (lon_length * lat_length);
    const HyperslabIndex lon_offset = (index - lev_offset * lon_length * lat_length) / lat_length; 
    const HyperslabIndex lat_offset = index - lev_offset * lon_length * lat_length - lon_offset * lat_length;

    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lon) + lon_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lat) + lat_offset;
    linear_indices[2] = hyperslab.GetDimensionStart(Dimension::Lev) + lev_offset;
}

void UpdateHyperslabIndexLonLevLat(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);

    const HyperslabIndex lon_offset = index / (lev_length * lat_length);
    const HyperslabIndex lev_offset = (index - lon_offset * lev_length * lat_length) / lat_length; 
    const HyperslabIndex lat_offset = index - lon_offset * lev_length * lat_length - lev_offset * lat_length;

    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lon) + lon_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lat) + lat_offset;
    linear_indices[2] = hyperslab.GetDimensionStart(Dimension::Lev) + lev_offset;
}

void UpdateHyperslabIndexLevLatLon(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);

    const HyperslabIndex lev_offset = index / (lat_length * lon_length);
    const HyperslabIndex lat_offset = (index - lev_offset * lat_length * lon_length) / lon_length; 
    const HyperslabIndex lon_offset = index - lev_offset * lat_length * lon_length - lat_offset * lon_length;

    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lon) + lon_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lat) + lat_offset;
    linear_indices[2] = hyperslab.GetDimensionStart(Dimension::Lev) + lev_offset;
}

void UpdateHyperslabIndexLatLevLon(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);

    const HyperslabIndex lat_offset = index / (lev_length * lon_length);
    const HyperslabIndex lev_offset = (index - lat_offset * lev_length * lon_length) / lon_length; 
    const HyperslabIndex lon_offset = index - lat_offset * lev_length * lon_length - lev_offset * lon_length;

    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lon) + lon_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lat) + lat_offset;
    linear_indices[2] = hyperslab.GetDimensionStart(Dimension::Lev) + lev_offset;
}

void UpdateHyperslabIndexLatLonLev(const Hyperslab& hyperslab, std::vector<HyperslabIndex>& linear_indices, const HyperslabIndex index)
{
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);

    const HyperslabIndex lat_offset = index / (lon_length * lev_length);
    const HyperslabIndex lon_offset = (index - lat_offset * lon_length * lev_length) / lev_length; 
    const HyperslabIndex lev_offset = index - lat_offset * lon_length * lev_length - lon_offset * lev_length;

    linear_indices[0] = hyperslab.GetDimensionStart(Dimension::Lon) + lon_offset;
    linear_indices[1] = hyperslab.GetDimensionStart(Dimension::Lat) + lat_offset;
    linear_indices[2] = hyperslab.GetDimensionStart(Dimension::Lev) + lev_offset;
}

[[maybe_unused]] void UpdateHyperslabIndexError([[maybe_unused]] const Hyperslab& hyperslab, [[maybe_unused]] std::vector<HyperslabIndex>& linear_indices, [[maybe_unused]] const HyperslabIndex index){}

UpdateHyperslabCoordinateFn
GetHyperslabCoordinatesIterationFunction(const DataLayout layout)
{
    switch (layout)
    {
        case DataLayout::Lon_Lat:
            return UpdateHyperslabIndexLonLat;
        break;
        case DataLayout::Lat_Lon:
            return UpdateHyperslabIndexLatLon;
        break;
        case DataLayout::Lat_Lev:
            return UpdateHyperslabIndexLatLev;
        break;
        case DataLayout::Lev_Lat:
            return UpdateHyperslabIndexLevLat;
        break;
        case DataLayout::Lon_Lev:
            return UpdateHyperslabIndexLonLev;
        break;
        case DataLayout::Lev_Lon:
            return UpdateHyperslabIndexLevLon;
        break;
        case DataLayout::Lon_Lat_Lev:
            return UpdateHyperslabIndexLonLatLev;
        break;
        case DataLayout::Lev_Lon_Lat:
            return UpdateHyperslabIndexLevLonLat;
        break;
        case DataLayout::Lon_Lev_Lat:
            return UpdateHyperslabIndexLonLevLat;
        break;
        case DataLayout::Lev_Lat_Lon:
            return UpdateHyperslabIndexLevLatLon;
        break;
        case DataLayout::Lat_Lev_Lon:
            return UpdateHyperslabIndexLatLevLon;
        break;
        case DataLayout::Lat_Lon_Lev:
            return UpdateHyperslabIndexLatLonLev;
        break;
        default :
            cmc_err_msg("The variable contains an undefined data layout.");
            return UpdateHyperslabIndexError;
    }
}

HyperslabIndex ReceiveDimensionIndexFromReferenceCoordsLonLat(const std::vector<HyperslabIndex>& reference_coordinates, const Dimension dimension)
{
    switch(dimension)
    {
        case Dimension::Lon:
            return reference_coordinates[0];
        break;
        case Dimension::Lat:
            return reference_coordinates[1];
        break;
        default:
            cmc_err_msg("The supplied dimension does not correspond to the initially supplied data layout.");
            return CMC_ERR;
    }
}


HyperslabIndex ReceiveDimensionIndexFromReferenceCoordsLatLev(const std::vector<HyperslabIndex>& reference_coordinates, const Dimension dimension)
{
    switch(dimension)
    {
        case Dimension::Lat:
            return reference_coordinates[0];
        break;
        case Dimension::Lev:
            return reference_coordinates[1];
        break;
        default:
            cmc_err_msg("The supplied dimension does not correspond to the initially supplied data layout.");
            return CMC_ERR;
    }
}

HyperslabIndex ReceiveDimensionIndexFromReferenceCoordsLonLev(const std::vector<HyperslabIndex>& reference_coordinates, const Dimension dimension)
{
    switch(dimension)
    {
        case Dimension::Lon:
            return reference_coordinates[0];
        break;
        case Dimension::Lev:
            return reference_coordinates[1];
        break;
        default:
            cmc_err_msg("The supplied dimension does not correspond to the initially supplied data layout.");
            return CMC_ERR;
    }
}

HyperslabIndex ReceiveDimensionIndexFromReferenceCoordsLonLatLev(const std::vector<HyperslabIndex>& reference_coordinates, const Dimension dimension)
{
    switch(dimension)
    {
        case Dimension::Lon:
            return reference_coordinates[0];
        break;
        case Dimension::Lat:
            return reference_coordinates[1];
        break;
        case Dimension::Lev:
            return reference_coordinates[2];
        break;
        default:
            cmc_err_msg("The supplied dimension does not correspond to the initially supplied data layout.");
            return CMC_ERR;
    }
}

[[maybe_unused]] HyperslabIndex ReceiveDimensionIndexFromReferenceCoordsError([[maybe_unused]] const std::vector<HyperslabIndex>& reference_coordinates, [[maybe_unused]] const Dimension dimension){}

DimensionValueExtractionFn
GetDimensionValueFunctionForReferenceCoords(const DataLayout layout)
{
    switch (layout)
    {
        case DataLayout::Lon_Lat:
            [[fallthrough]];
        case DataLayout::Lat_Lon:
            return ReceiveDimensionIndexFromReferenceCoordsLonLat;
        break;
        case DataLayout::Lat_Lev:
            [[fallthrough]];
        case DataLayout::Lev_Lat:
            return ReceiveDimensionIndexFromReferenceCoordsLatLev;
        break;
        case DataLayout::Lon_Lev:
            [[fallthrough]];
        case DataLayout::Lev_Lon:
            return ReceiveDimensionIndexFromReferenceCoordsLonLev;
        break;
        case DataLayout::Lon_Lat_Lev:
            [[fallthrough]];
        case DataLayout::Lev_Lon_Lat:
            [[fallthrough]];
        case DataLayout::Lon_Lev_Lat:
            [[fallthrough]];
        case DataLayout::Lev_Lat_Lon:
            [[fallthrough]];
        case DataLayout::Lat_Lev_Lon:
            [[fallthrough]];
        case DataLayout::Lat_Lon_Lev:
            return ReceiveDimensionIndexFromReferenceCoordsLonLatLev;
        break;
        default :
            cmc_err_msg("The variable contains an undefined data layout.");
            return ReceiveDimensionIndexFromReferenceCoordsError;
    }
}

std::vector<HyperslabIndex> TrimRefCoordsToLonLat(const std::vector<HyperslabIndex>& reference_coordinates)
{
    return std::vector<HyperslabIndex>{reference_coordinates[0], reference_coordinates[1]};
}
std::vector<HyperslabIndex> TrimRefCoordsToLatLev(const std::vector<HyperslabIndex>& reference_coordinates)
{
    return std::vector<HyperslabIndex>{reference_coordinates[1], reference_coordinates[2]};
}
std::vector<HyperslabIndex> TrimRefCoordsToLonLev(const std::vector<HyperslabIndex>& reference_coordinates)
{
    return std::vector<HyperslabIndex>{reference_coordinates[0], reference_coordinates[2]};
}
[[maybe_unused]] std::vector<HyperslabIndex> TrimRefCoordsError([[maybe_unused]] const std::vector<HyperslabIndex>& reference_coordinates)
{
    return std::vector<HyperslabIndex>();
}

TrimRefCoordVectorFn
GetReferenceCoordTrimmingFunction(const DataLayout trimmed_data_layout)
{
    switch (trimmed_data_layout)
    {
        case DataLayout::Lon_Lat:
            [[fallthrough]];
        case DataLayout::Lat_Lon:
            return TrimRefCoordsToLonLat;
        break;
        case DataLayout::Lat_Lev:
            [[fallthrough]];
        case DataLayout::Lev_Lat:
            return TrimRefCoordsToLatLev;
        break;
        case DataLayout::Lon_Lev:
            [[fallthrough]];
        case DataLayout::Lev_Lon:
            return TrimRefCoordsToLonLev;
        break;
        default :
            cmc_err_msg("The trimmed data layout has to be 2D.");
            return TrimRefCoordsError;
    }
}


}
