#include "utilities/cmc_geo_utilities.hxx"

namespace cmc 
{

Hyperslab
TransformGeoDomainToHyperslab(const GeoDomain& domain)
{
    return Hyperslab(domain.GetDimensionStartIndex(Dimension::Lon), (domain.GetDimensionLength(Dimension::Lon) > 1 ? domain.GetDimensionLength(Dimension::Lon) : 1),
                     domain.GetDimensionStartIndex(Dimension::Lat), (domain.GetDimensionLength(Dimension::Lat) > 1 ? domain.GetDimensionLength(Dimension::Lat) : 1),
                     domain.GetDimensionStartIndex(Dimension::Lev), (domain.GetDimensionLength(Dimension::Lev) > 1 ? domain.GetDimensionLength(Dimension::Lev) : 1),
                     domain.GetDimensionStartIndex(Dimension::Time), (domain.GetDimensionLength(Dimension::Time) > 1 ? domain.GetDimensionLength(Dimension::Time) : 1));
}

GeoDomain
TransformHyperslabToGeoDomain(const Hyperslab& hyperslab)
{
    return GeoDomain(hyperslab.GetDimensionStart(Dimension::Lon),  hyperslab.GetDimensionStart(Dimension::Lon)  + hyperslab.GetDimensionLength(Dimension::Lon),
                     hyperslab.GetDimensionStart(Dimension::Lat),  hyperslab.GetDimensionStart(Dimension::Lat)  + hyperslab.GetDimensionLength(Dimension::Lat),
                     hyperslab.GetDimensionStart(Dimension::Lev),  hyperslab.GetDimensionStart(Dimension::Lev)  + hyperslab.GetDimensionLength(Dimension::Lev),
                     hyperslab.GetDimensionStart(Dimension::Time), hyperslab.GetDimensionStart(Dimension::Time) + hyperslab.GetDimensionLength(Dimension::Time));
}

Hyperslab
SubtractGeoDomainOffset(const Hyperslab& hyperslab, const GeoDomain& domain)
{
    return Hyperslab(hyperslab.GetDimensionStart(Dimension::Lon) - domain.GetDimensionStartIndex(Dimension::Lon),   hyperslab.GetDimensionLength(Dimension::Lon),
                     hyperslab.GetDimensionStart(Dimension::Lat) - domain.GetDimensionStartIndex(Dimension::Lat),   hyperslab.GetDimensionLength(Dimension::Lat),
                     hyperslab.GetDimensionStart(Dimension::Lev) - domain.GetDimensionStartIndex(Dimension::Lev),   hyperslab.GetDimensionLength(Dimension::Lev),
                     hyperslab.GetDimensionStart(Dimension::Time) - domain.GetDimensionStartIndex(Dimension::Time), hyperslab.GetDimensionLength(Dimension::Time));
}

}
