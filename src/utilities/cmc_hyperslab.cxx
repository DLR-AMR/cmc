#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_log_functions.hxx"

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
}
Hyperslab::Hyperslab(const DimensionInterval& dimension1, const DimensionInterval& dimension2)
{
    SetupDimension(dimension1);
    SetupDimension(dimension2);
}
Hyperslab::Hyperslab(const DimensionInterval& dimension1, const DimensionInterval& dimension2,
                     const DimensionInterval& dimension3)
{
    SetupDimension(dimension1);
    SetupDimension(dimension2);
    SetupDimension(dimension3);
}
Hyperslab::Hyperslab(const DimensionInterval& dimension1, const DimensionInterval& dimension2,
                     const DimensionInterval& dimension3, const DimensionInterval& dimension4)
{
    SetupDimension(dimension1);
    SetupDimension(dimension2);
    SetupDimension(dimension3);
    SetupDimension(dimension4);
}

Hyperslab::Hyperslab(const HyperslabIndex lon_start, const HyperslabIndex lon_count,
                     const HyperslabIndex lat_start, const HyperslabIndex lat_count,
                     const HyperslabIndex lev_start, const HyperslabIndex lev_count,
                     const HyperslabIndex time_start, const HyperslabIndex time_count)
: start_indices_{lon_start, lat_start, lev_start, time_start}, count_indices_{lon_count, lat_count, lev_count, time_count}{}

HyperslabIndex
Hyperslab::GetNumberCoordinates() const
{
    const HyperslabIndex product =  std::accumulate(count_indices_.begin(), count_indices_.end(), 1, std::multiplies<HyperslabIndex>());
    return product;
}

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

[[maybe_unused]] HyperslabIndex ReceiveDimensionIndexFromReferenceCoordsError([[maybe_unused]] const std::vector<HyperslabIndex>& reference_coordinates, [[maybe_unused]] const Dimension dimension){return CMC_ERR;}

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


/** Extraction of indices for the whole subdomain **/
/** If indices are outside the global domain, they are flagged with kOutsideOfHyperslabDomain **/

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LonLat(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lat_start = global_hyperslab.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex global_lon_start = global_hyperslab.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex global_lat_length = global_hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex global_lon_length = global_hyperslab.GetDimensionLength(Dimension::Lon);

    const HyperslabIndex lat_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex lon_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex lat_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lon_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lon);

    cmc_assert(global_lon_start <= lon_start);
    cmc_assert(global_lat_start <= lat_start);

    HyperslabIndex offset = (lon_start - global_lon_start) * global_lat_length + (lat_start - global_lat_start);

    for (HyperslabIndex lon_iter = 0; lon_iter < lon_length; ++lon_iter)
    {
        for (HyperslabIndex lat_iter = 0; lat_iter < lat_length; ++lat_iter)
        {
            if ((lat_iter + lat_start < global_lat_start + global_lat_length) &&
                (lon_iter + lon_start < global_lon_start + global_lon_length))
            {
                /* We are in the domain of the global hyperslab */
                indices.push_back(offset + lat_iter);
            } else
            {
                /* Indicate that this is outside of the domain */
                indices.push_back(kOutsideOfHyperslabDomain);
            }
        }

        /* Update the offset */
        offset += global_lat_length;
    }

    return indices;
}

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LatLon(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lat_start = global_hyperslab.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex global_lon_start = global_hyperslab.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex global_lat_length = global_hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex global_lon_length = global_hyperslab.GetDimensionLength(Dimension::Lon);

    const HyperslabIndex lat_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex lon_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex lat_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lon_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lon);

    cmc_assert(global_lon_start <= lon_start);
    cmc_assert(global_lat_start <= lat_start);

    HyperslabIndex offset = (lat_start - global_lat_start) * global_lon_length + (lon_start - global_lon_start);

    for (HyperslabIndex lat_iter = 0; lat_iter < lat_length; ++lat_iter)
    {
        for (HyperslabIndex lon_iter = 0; lon_iter < lon_length; ++lon_iter)
        {
            if ((lat_iter + lat_start < global_lat_start + global_lat_length) &&
                (lon_iter + lon_start < global_lon_start + global_lon_length))
            {
                /* We are in the domain of the global hyperslab */
                indices.push_back(offset + lon_iter);
            } else
            {
                /* Indicate that this is outside of the domain */
                indices.push_back(kOutsideOfHyperslabDomain);
            }
        }

        /* Update the offset */
        offset += global_lon_length;
    }

    return indices;
}

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LatLev(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lat_start = global_hyperslab.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex global_lev_start = global_hyperslab.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex global_lat_length = global_hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex global_lev_length = global_hyperslab.GetDimensionLength(Dimension::Lev);

    const HyperslabIndex lat_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex lev_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex lat_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lev_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lev);

    cmc_assert(global_lev_start <= lev_start);
    cmc_assert(global_lat_start <= lat_start);

    HyperslabIndex offset = (lat_start - global_lat_start) * global_lev_length + (lev_start - global_lev_start);

    for (HyperslabIndex lat_iter = 0; lat_iter < lat_length; ++lat_iter)
    {
        for (HyperslabIndex lev_iter = 0; lev_iter < lev_length; ++lev_iter)
        {
            if ((lat_iter + lat_start < global_lat_start + global_lat_length) &&
                (lev_iter + lev_start < global_lev_start + global_lev_length))
            {
                /* We are in the domain of the global hyperslab */
                indices.push_back(offset + lev_iter);
            } else
            {
                /* Indicate that this is outside of the domain */
                indices.push_back(kOutsideOfHyperslabDomain);
            }
        }

        /* Update the offset */
        offset += global_lev_length;
    }

    return indices;
}

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LevLat(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lat_start = global_hyperslab.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex global_lev_start = global_hyperslab.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex global_lat_length = global_hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex global_lev_length = global_hyperslab.GetDimensionLength(Dimension::Lev);

    const HyperslabIndex lat_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex lev_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex lat_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lev_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lev);

    cmc_assert(global_lev_start <= lev_start);
    cmc_assert(global_lat_start <= lat_start);

    HyperslabIndex offset = (lev_start - global_lev_start) * global_lat_length + (lat_start - global_lat_start);

    for (HyperslabIndex lev_iter = 0; lev_iter < lev_length; ++lev_iter)
    {
        for (HyperslabIndex lat_iter = 0; lat_iter < lat_length; ++lat_iter)
        {
            if ((lat_iter + lat_start < global_lat_start + global_lat_length) &&
                (lev_iter + lev_start < global_lev_start + global_lev_length))
            {
                /* We are in the domain of the global hyperslab */
                indices.push_back(offset + lat_iter);
            } else
            {
                /* Indicate that this is outside of the domain */
                indices.push_back(kOutsideOfHyperslabDomain);
            }
        }

        /* Update the offset */
        offset += global_lat_length;
    }

    return indices;
}

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LonLev(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lon_start = global_hyperslab.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex global_lev_start = global_hyperslab.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex global_lon_length = global_hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex global_lev_length = global_hyperslab.GetDimensionLength(Dimension::Lev);

    const HyperslabIndex lon_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex lev_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex lon_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lev_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lev);

    cmc_assert(global_lev_start <= lev_start);
    cmc_assert(global_lon_start <= lon_start);

    HyperslabIndex offset = (lon_start - global_lon_start) * global_lev_length + (lev_start - global_lev_start);

    for (HyperslabIndex lon_iter = 0; lon_iter < lon_length; ++lon_iter)
    {
        for (HyperslabIndex lev_iter = 0; lev_iter < lev_length; ++lev_iter)
        {
            if ((lon_iter + lon_start < global_lon_start + global_lon_length) &&
                (lev_iter + lev_start < global_lev_start + global_lev_length))
            {
                /* We are in the domain of the global hyperslab */
                indices.push_back(offset + lev_iter);
            } else
            {
                /* Indicate that this is outside of the domain */
                indices.push_back(kOutsideOfHyperslabDomain);
            }
        }

        /* Update the offset */
        offset += global_lev_length;
    }

    return indices;
}

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LevLon(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lev_start = global_hyperslab.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex global_lon_start = global_hyperslab.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex global_lev_length = global_hyperslab.GetDimensionLength(Dimension::Lev);
    const HyperslabIndex global_lon_length = global_hyperslab.GetDimensionLength(Dimension::Lon);

    const HyperslabIndex lev_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex lon_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex lev_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lev);
    const HyperslabIndex lon_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lon);

    cmc_assert(global_lon_start <= lon_start);
    cmc_assert(global_lev_start <= lev_start);

    HyperslabIndex offset = (lev_start - global_lev_start) * global_lon_length + (lon_start - global_lon_start);

    for (HyperslabIndex lev_iter = 0; lev_iter < lev_length; ++lev_iter)
    {
        for (HyperslabIndex lon_iter = 0; lon_iter < lon_length; ++lon_iter)
        {
            if ((lev_iter + lev_start < global_lev_start + global_lev_length) &&
                (lon_iter + lon_start < global_lon_start + global_lon_length))
            {
                /* We are in the domain of the global hyperslab */
                indices.push_back(offset + lon_iter);
            } else
            {
                /* Indicate that this is outside of the domain */
                indices.push_back(kOutsideOfHyperslabDomain);
            }
        }

        /* Update the offset */
        offset += global_lon_length;
    }

    return indices;
}

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LonLatLev(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lon_start = global_hyperslab.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex global_lat_start = global_hyperslab.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex global_lev_start = global_hyperslab.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex global_lon_length = global_hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex global_lat_length = global_hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex global_lev_length = global_hyperslab.GetDimensionLength(Dimension::Lev);

    const HyperslabIndex lon_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex lat_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex lev_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex lon_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lat_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lev_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lev);

    cmc_assert(global_lev_start <= lev_start);
    cmc_assert(global_lat_start <= lat_start);
    cmc_assert(global_lon_start <= lon_start);

    HyperslabIndex offset = (lon_start - global_lon_start) * global_lat_length * global_lev_length + (lat_start - global_lat_start) * global_lev_length + (lev_start - global_lev_start);

    for (HyperslabIndex lon_iter = 0; lon_iter < lon_length; ++lon_iter)
    {
        for (HyperslabIndex lat_iter = 0; lat_iter < lat_length; ++lat_iter)
        {
            for (HyperslabIndex lev_iter = 0; lev_iter < lev_length; ++lev_iter)
            {
                if ((lat_iter + lat_start < global_lat_start + global_lat_length) &&
                    (lev_iter + lev_start < global_lev_start + global_lev_length) &&
                    (lon_iter + lon_start < global_lon_start + global_lon_length))
                {
                    /* We are in the domain of the global hyperslab */
                    indices.push_back(offset + lev_iter);
                } else
                {
                    /* Indicate that this is outside of the domain */
                    indices.push_back(kOutsideOfHyperslabDomain);
                }
            }

            /* Update the offset */
            offset += global_lev_length;
        }

        /* Update the offset */
        offset += (global_lat_length - lat_length) * global_lev_length;
    }
    return indices;
}

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LevLonLat(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lon_start = global_hyperslab.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex global_lat_start = global_hyperslab.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex global_lev_start = global_hyperslab.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex global_lon_length = global_hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex global_lat_length = global_hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex global_lev_length = global_hyperslab.GetDimensionLength(Dimension::Lev);

    const HyperslabIndex lon_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex lat_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex lev_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex lon_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lat_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lev_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lev);

    cmc_assert(global_lev_start <= lev_start);
    cmc_assert(global_lat_start <= lat_start);
    cmc_assert(global_lon_start <= lon_start);

    HyperslabIndex offset = (lev_start - global_lev_start) * global_lon_length * global_lat_length + (lon_start - global_lon_start) * global_lat_length + (lat_start - global_lat_start);

    for (HyperslabIndex lev_iter = 0; lev_iter < lev_length; ++lev_iter)
    {
        for (HyperslabIndex lon_iter = 0; lon_iter < lon_length; ++lon_iter)
        {
            for (HyperslabIndex lat_iter = 0; lat_iter < lat_length; ++lat_iter)
            {
                if ((lat_iter + lat_start < global_lat_start + global_lat_length) &&
                    (lev_iter + lev_start < global_lev_start + global_lev_length) &&
                    (lon_iter + lon_start < global_lon_start + global_lon_length))
                {
                    /* We are in the domain of the global hyperslab */
                    indices.push_back(offset + lat_iter);
                } else
                {
                    /* Indicate that this is outside of the domain */
                    indices.push_back(kOutsideOfHyperslabDomain);
                }
            }

            /* Update the offset */
            offset += global_lat_length;
        }

        /* Update the offset */
        offset += (global_lon_length - lon_length) * global_lat_length;
    }
    return indices;
}

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LonLevLat(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lon_start = global_hyperslab.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex global_lat_start = global_hyperslab.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex global_lev_start = global_hyperslab.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex global_lon_length = global_hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex global_lat_length = global_hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex global_lev_length = global_hyperslab.GetDimensionLength(Dimension::Lev);

    const HyperslabIndex lon_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex lat_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex lev_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex lon_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lat_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lev_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lev);

    cmc_assert(global_lev_start <= lev_start);
    cmc_assert(global_lat_start <= lat_start);
    cmc_assert(global_lon_start <= lon_start);

    HyperslabIndex offset = (lon_start - global_lon_start) * global_lev_length * global_lat_length + (lev_start - global_lev_start) * global_lat_length + (lat_start - global_lat_start);

    for (HyperslabIndex lon_iter = 0; lon_iter < lon_length; ++lon_iter)
    {
        for (HyperslabIndex lev_iter = 0; lev_iter < lev_length; ++lev_iter)
        {
            for (HyperslabIndex lat_iter = 0; lat_iter < lat_length; ++lat_iter)
            {
                if ((lat_iter + lat_start < global_lat_start + global_lat_length) &&
                    (lev_iter + lev_start < global_lev_start + global_lev_length) &&
                    (lon_iter + lon_start < global_lon_start + global_lon_length))
                {
                    /* We are in the domain of the global hyperslab */
                    indices.push_back(offset + lat_iter);
                } else
                {
                    /* Indicate that this is outside of the domain */
                    indices.push_back(kOutsideOfHyperslabDomain);
                }
            }

            /* Update the offset */
            offset += global_lat_length;
        }

        /* Update the offset */
        offset += (global_lev_length - lev_length) * global_lat_length;
    }
    return indices;
}

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LevLatLon(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lon_start = global_hyperslab.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex global_lat_start = global_hyperslab.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex global_lev_start = global_hyperslab.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex global_lon_length = global_hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex global_lat_length = global_hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex global_lev_length = global_hyperslab.GetDimensionLength(Dimension::Lev);

    const HyperslabIndex lon_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex lat_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex lev_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex lon_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lat_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lev_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lev);

    cmc_assert(global_lev_start <= lev_start);
    cmc_assert(global_lat_start <= lat_start);
    cmc_assert(global_lon_start <= lon_start);

    HyperslabIndex offset = (lev_start - global_lev_start) * global_lat_length * global_lon_length + (lat_start - global_lat_start) * global_lon_length + (lon_start - global_lon_start);

    for (HyperslabIndex lev_iter = 0; lev_iter < lev_length; ++lev_iter)
    {
        for (HyperslabIndex lat_iter = 0; lat_iter < lat_length; ++lat_iter)
        {
            for (HyperslabIndex lon_iter = 0; lon_iter < lon_length; ++lon_iter)
            {
                if ((lat_iter + lat_start < global_lat_start + global_lat_length) &&
                    (lev_iter + lev_start < global_lev_start + global_lev_length) &&
                    (lon_iter + lon_start < global_lon_start + global_lon_length))
                {
                    /* We are in the domain of the global hyperslab */
                    indices.push_back(offset + lon_iter);
                } else
                {
                    /* Indicate that this is outside of the domain */
                    indices.push_back(kOutsideOfHyperslabDomain);
                }
            }

            /* Update the offset */
            offset += global_lon_length;
        }

        /* Update the offset */
        offset += (global_lat_length - lat_length) * global_lon_length;
    } 
    return indices;
}

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LatLevLon(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lon_start = global_hyperslab.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex global_lat_start = global_hyperslab.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex global_lev_start = global_hyperslab.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex global_lon_length = global_hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex global_lat_length = global_hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex global_lev_length = global_hyperslab.GetDimensionLength(Dimension::Lev);

    const HyperslabIndex lon_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex lat_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex lev_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex lon_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lat_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lev_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lev);

    cmc_assert(global_lev_start <= lev_start);
    cmc_assert(global_lat_start <= lat_start);
    cmc_assert(global_lon_start <= lon_start);

    HyperslabIndex offset = (lat_start - global_lat_start) * global_lev_length * global_lon_length + (lev_start - global_lev_start) * global_lon_length + (lon_start - global_lon_start);

    for (HyperslabIndex lat_iter = 0; lat_iter < lat_length; ++lat_iter)
    {
        for (HyperslabIndex lev_iter = 0; lev_iter < lev_length; ++lev_iter)
        {
            for (HyperslabIndex lon_iter = 0; lon_iter < lon_length; ++lon_iter)
            {
                if ((lat_iter + lat_start < global_lat_start + global_lat_length) &&
                    (lev_iter + lev_start < global_lev_start + global_lev_length) &&
                    (lon_iter + lon_start < global_lon_start + global_lon_length))
                {
                    /* We are in the domain of the global hyperslab */
                    indices.push_back(offset + lon_iter);
                } else
                {
                    /* Indicate that this is outside of the domain */
                    indices.push_back(kOutsideOfHyperslabDomain);
                }
            }

            /* Update the offset */
            offset += global_lon_length;
        }

        /* Update the offset */
        offset += (global_lev_length - lev_length) * global_lon_length;
    }
    return indices;
}

std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_LatLonLev(const Hyperslab& global_hyperslab, const Hyperslab& domain_indices_to_extract)
{
    std::vector<HyperslabIndex> indices;
    indices.reserve(domain_indices_to_extract.GetNumberCoordinates());

    const HyperslabIndex global_lon_start = global_hyperslab.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex global_lat_start = global_hyperslab.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex global_lev_start = global_hyperslab.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex global_lon_length = global_hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex global_lat_length = global_hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex global_lev_length = global_hyperslab.GetDimensionLength(Dimension::Lev);

    const HyperslabIndex lon_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lon);
    const HyperslabIndex lat_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lat);
    const HyperslabIndex lev_start = domain_indices_to_extract.GetDimensionStart(Dimension::Lev);
    const HyperslabIndex lon_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lat_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lev_length = domain_indices_to_extract.GetDimensionLength(Dimension::Lev);

    cmc_assert(global_lev_start <= lev_start);
    cmc_assert(global_lat_start <= lat_start);
    cmc_assert(global_lon_start <= lon_start);

    HyperslabIndex offset = (lat_start - global_lat_start) * global_lon_length * global_lev_length + (lon_start - global_lon_start) * global_lev_length + (lev_start - global_lev_start);

    for (HyperslabIndex lat_iter = 0; lat_iter < lat_length; ++lat_iter)
    {
        for (HyperslabIndex lon_iter = 0; lon_iter < lon_length; ++lon_iter)
        {
            for (HyperslabIndex lev_iter = 0; lev_iter < lev_length; ++lev_iter)
            {
                if ((lat_iter + lat_start < global_lat_start + global_lat_length) &&
                    (lev_iter + lev_start < global_lev_start + global_lev_length) &&
                    (lon_iter + lon_start < global_lon_start + global_lon_length))
                {
                    /* We are in the domain of the global hyperslab */
                    indices.push_back(offset + lev_iter);
                } else
                {
                    /* Indicate that this is outside of the domain */
                    indices.push_back(kOutsideOfHyperslabDomain);
                }
            }

            /* Update the offset */
            offset += global_lev_length;
        }

        /* Update the offset */
        offset += (global_lon_length - lon_length) * global_lev_length;
    }
    return indices;
}

[[maybe_unused]] std::vector<HyperslabIndex> GetIndicesForHyperslabDataExtraction_Error([[maybe_unused]] const Hyperslab& global_hyperslab, [[maybe_unused]] const Hyperslab& domain_indices_to_extract){return std::vector<HyperslabIndex>();}

LinearIndicesExtractionFn
GetIndicesExtractionFunction(const DataLayout layout)
{
    switch (layout)
    {
        case DataLayout::Lon_Lat:
            return GetIndicesForHyperslabDataExtraction_LonLat;
        break;
        case DataLayout::Lat_Lon:
            return GetIndicesForHyperslabDataExtraction_LatLon;
        break;
        case DataLayout::Lat_Lev:
            return GetIndicesForHyperslabDataExtraction_LatLev;
        break;
        case DataLayout::Lev_Lat:
            return GetIndicesForHyperslabDataExtraction_LevLat;
        break;
        case DataLayout::Lon_Lev:
            return GetIndicesForHyperslabDataExtraction_LonLev;
        break;
        case DataLayout::Lev_Lon:
            return GetIndicesForHyperslabDataExtraction_LevLon;
        break;
        case DataLayout::Lon_Lat_Lev:
            return GetIndicesForHyperslabDataExtraction_LonLatLev;
        break;
        case DataLayout::Lev_Lon_Lat:
            return GetIndicesForHyperslabDataExtraction_LevLonLat;
        break;
        case DataLayout::Lon_Lev_Lat:
            return GetIndicesForHyperslabDataExtraction_LonLevLat;
        break;
        case DataLayout::Lev_Lat_Lon:
            return GetIndicesForHyperslabDataExtraction_LevLatLon;
        break;
        case DataLayout::Lat_Lev_Lon:
            return GetIndicesForHyperslabDataExtraction_LatLevLon;
        break;
        case DataLayout::Lat_Lon_Lev:
            return GetIndicesForHyperslabDataExtraction_LatLonLev;
        break;
        default :
            cmc_err_msg("The variable contains an undefined data layout.");
            return GetIndicesForHyperslabDataExtraction_Error;
    }
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLonLat(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);
    return hs_indices[0] * lat_length + hs_indices[1];
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLatLon(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);
    return hs_indices[0] + hs_indices[1] * lon_length;
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLatLev(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);
    return hs_indices[0] * lev_length + hs_indices[1];
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLevLat(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);
    return hs_indices[0] + hs_indices[1] * lat_length;
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLonLev(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);
    return hs_indices[0] * lev_length + hs_indices[1];
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLevLon(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);
    return hs_indices[0] + hs_indices[1] * lon_length;
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLonLatLev(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);

    return hs_indices[0] * lat_length * lev_length + hs_indices[1] * lev_length + hs_indices[2];
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLevLonLat(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);

    return hs_indices[0] * lat_length + hs_indices[1] + hs_indices[2] * lon_length * lat_length;
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLonLevLat(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);

    return hs_indices[0] * lev_length * lat_length + hs_indices[1] + hs_indices[2] * lat_length;
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLevLatLon(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lat_length = hyperslab.GetDimensionLength(Dimension::Lat);
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);

    return hs_indices[0] + hs_indices[1] * lon_length + hs_indices[2] * lat_length * lon_length;
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLatLevLon(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);

    return hs_indices[0] + hs_indices[1] * lev_length * lon_length + hs_indices[2] * lon_length;
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsLatLonLev(const Hyperslab& hyperslab, const std::vector<HyperslabIndex>& hs_indices)
{
    const HyperslabIndex lon_length = hyperslab.GetDimensionLength(Dimension::Lon);
    const HyperslabIndex lev_length = hyperslab.GetDimensionLength(Dimension::Lev);

    return hs_indices[0] * lev_length + hs_indices[1] * lev_length * lon_length + hs_indices[2];
}

HyperslabIndex GetLinearizedIndexFromHyperslabCoordsError([[maybe_unused]] const Hyperslab& hyperslab, [[maybe_unused]] const std::vector<HyperslabIndex>& hs_indices){return CMC_ERR;}

LinearizeHyperslabCoordiantesFn
GetLinearizedIndexFromHyperslabCoordsFunction(const DataLayout layout)
{
    switch (layout)
    {
        case DataLayout::Lon_Lat:
            return GetLinearizedIndexFromHyperslabCoordsLonLat;
        break;
        case DataLayout::Lat_Lon:
            return GetLinearizedIndexFromHyperslabCoordsLatLon;
        break;
        case DataLayout::Lat_Lev:
            return GetLinearizedIndexFromHyperslabCoordsLatLev;
        break;
        case DataLayout::Lev_Lat:
            return GetLinearizedIndexFromHyperslabCoordsLevLat;
        break;
        case DataLayout::Lon_Lev:
            return GetLinearizedIndexFromHyperslabCoordsLonLev;
        break;
        case DataLayout::Lev_Lon:
            return GetLinearizedIndexFromHyperslabCoordsLevLon;
        break;
        case DataLayout::Lon_Lat_Lev:
            return GetLinearizedIndexFromHyperslabCoordsLonLatLev;
        break;
        case DataLayout::Lev_Lon_Lat:
            return GetLinearizedIndexFromHyperslabCoordsLevLonLat;
        break;
        case DataLayout::Lon_Lev_Lat:
            return GetLinearizedIndexFromHyperslabCoordsLonLevLat;
        break;
        case DataLayout::Lev_Lat_Lon:
            return GetLinearizedIndexFromHyperslabCoordsLevLatLon;
        break;
        case DataLayout::Lat_Lev_Lon:
            return GetLinearizedIndexFromHyperslabCoordsLatLevLon;
        break;
        case DataLayout::Lat_Lon_Lev:
            return GetLinearizedIndexFromHyperslabCoordsLatLonLev;
        break;
        default :
            cmc_err_msg("The variable contains an undefined data layout.");
            return GetLinearizedIndexFromHyperslabCoordsError;
    }
}




}
