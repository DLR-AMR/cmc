#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_log_functions.hxx"

namespace cmc
{

void
GeoDomain::SetupDimension(const DimensionInterval& dimension)
{
    switch (dimension.dim)
    {
        case Dimension::Lon:
            start_indices_[Dimension::Lon] = dimension.start_index;
            end_indices_[Dimension::Lon] = dimension.end_index;
        break;
        case Dimension::Lat:
            start_indices_[Dimension::Lat] = dimension.start_index;
            end_indices_[Dimension::Lat] = dimension.end_index;
        break;
        case Dimension::Lev:
            start_indices_[Dimension::Lev] = dimension.start_index;
            end_indices_[Dimension::Lev] = dimension.end_index;
        break;
        case Dimension::Time:
            start_indices_[Dimension::Time] = dimension.start_index;
            end_indices_[Dimension::Time] = dimension.end_index;
        break;
        default:
            cmc_err_msg("The dimension cannot be considered.");
    }
}

GeoDomain::GeoDomain(const DimensionInterval& dimension1)
{
    SetupDimension(dimension1);
}
GeoDomain::GeoDomain(const DimensionInterval& dimension1, const DimensionInterval& dimension2)
{
    SetupDimension(dimension1);
    SetupDimension(dimension2);
}
GeoDomain::GeoDomain(const DimensionInterval& dimension1, const DimensionInterval& dimension2,
                     const DimensionInterval& dimension3)
{
    SetupDimension(dimension1);
    SetupDimension(dimension2);
    SetupDimension(dimension3);
}
GeoDomain::GeoDomain(const DimensionInterval& dimension1, const DimensionInterval& dimension2,
                     const DimensionInterval& dimension3, const DimensionInterval& dimension4)
{
    SetupDimension(dimension1);
    SetupDimension(dimension2);
    SetupDimension(dimension3);
    SetupDimension(dimension4);
}

GeoDomain::GeoDomain(const DomainIndex lon_dimension_start, const DomainIndex lon_dimension_end,
              const DomainIndex lat_dimension_start, const DomainIndex lat_dimension_end,
              const DomainIndex lev_dimension_start, const DomainIndex lev_dimension_end,
              const DomainIndex time_dimension_start, const DomainIndex time_dimension_end)
: start_indices_{lon_dimension_start, lat_dimension_start, lev_dimension_start, time_dimension_start},
  end_indices_{lon_dimension_end, lat_dimension_end, lev_dimension_end, time_dimension_end} {}

DomainIndex
GeoDomain::GetDimensionLength(const Dimension dimension) const
{
    switch (dimension)
    {
        case Dimension::Lon:
            return end_indices_[Dimension::Lon] - start_indices_[Dimension::Lon];
        break;
        case Dimension::Lat:
            return end_indices_[Dimension::Lat] - start_indices_[Dimension::Lat];
        break;
        case Dimension::Lev:
            return end_indices_[Dimension::Lev] - start_indices_[Dimension::Lev];
        break;
        case Dimension::Time:
            return end_indices_[Dimension::Time] - start_indices_[Dimension::Time];
        break;
        default:
            cmc_err_msg("The supplied dimension is not considered.");
    }
}

int
GeoDomain::GetDimensionality() const
{
    int dimensionality{0};
    for (int i = 0; i < Dimension::NumCoordinates; ++i)
    {
        if (GetDimensionLength(static_cast<Dimension>(i)) > 1)
        {
            ++dimensionality;
        }
    }
    return dimensionality;
}

DomainIndex
GeoDomain::GetLargestDimensionLength() const
{
    DomainIndex largest_dim_length{0};

    for (int i = 0; i < Dimension::NumCoordinates; ++i)
    {
        const DomainIndex dim_length = GetDimensionLength(static_cast<Dimension>(i));
        if (dim_length > largest_dim_length)
        {
            largest_dim_length = dim_length;
        }
    }
    return largest_dim_length;
}

DomainIndex 
GeoDomain::GetNumberReferenceCoordsCovered() const
{
    DomainIndex num_covered_coords = 1;

    for (int i = 0; i < Dimension::NumCoordinates; ++i)
    {
        const DomainIndex dim_length = GetDimensionLength(static_cast<Dimension>(i));
        if (dim_length >= 1)
        {
            num_covered_coords *= dim_length;
        }
    }

    return num_covered_coords;
}

void
GeoDomain::UpdateDimension(const DimensionInterval& dimension)
{
    SetupDimension(dimension);
}

void
GeoDomain::ClearDimension(const Dimension removed_dimension)
{
    start_indices_[removed_dimension] = 0;
    end_indices_[removed_dimension] = 0;
}

DomainIndex
GeoDomain::GetDimensionStartIndex(const Dimension dimension) const
{
    return start_indices_[dimension];
}

DomainIndex
GeoDomain::GetDimensionEndIndex(const Dimension dimension) const
{
    return end_indices_[dimension];
}

bool
CompareGeoDomains(const GeoDomain& domain1, const GeoDomain& domain2)
{
    if (std::equal(domain1.start_indices_.begin(), domain1.start_indices_.end(), domain2.start_indices_.begin()) &&
        std::equal(domain1.start_indices_.begin(), domain1.start_indices_.end(), domain2.start_indices_.begin()))
    {
        return true;
    } else
    {
        return false;
    }
}

GeoDomain 
ExtendGeoDomain(const GeoDomain& domain, const DimensionInterval& add_dimension)
{
    GeoDomain extended_domain = domain;
    extended_domain.UpdateDimension(add_dimension);
    return extended_domain;
}

bool
GeoDomain::IsValid() const
{
    for (int iter = 0; iter < Dimension::NumCoordinates; ++iter)
    {
        if (start_indices_[iter] > end_indices_[iter])
        {
            return false;
        }
    }

    return true;
}

GeoDomain
GeoDomain::GetZeroOffsetDomain() const
{
    return GeoDomain(0, GetDimensionLength(Dimension::Lon),
                     0, GetDimensionLength(Dimension::Lat),
                     0, GetDimensionLength(Dimension::Lev),
                     0, GetDimensionLength(Dimension::Time));
}


GeoDomain
ReconstructGeoDomainFromStartAndEndIndices(const std::array<DomainIndex, Dimension::NumCoordinates>& start_indices, const std::array<DomainIndex, Dimension::NumCoordinates>& end_indices)
{
    GeoDomain domain;
    domain.start_indices_ = start_indices;
    domain.end_indices_ = end_indices;

    return domain;
}

}
