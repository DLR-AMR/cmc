#include "utilities/cmc_geo_utilities.hxx"

namespace cmc 
{

Hyperslab
TransformGeoDomainToHyperslab(const GeoDomain& domain)
{
    return Hyperslab(domain.GetDimensionStartIndex(Dimension::Lon), domain.GetDimensionLength(Dimension::Lon),
                     domain.GetDimensionStartIndex(Dimension::Lat), domain.GetDimensionLength(Dimension::Lat),
                     domain.GetDimensionStartIndex(Dimension::Lev), domain.GetDimensionLength(Dimension::Lev),
                     domain.GetDimensionStartIndex(Dimension::Time), domain.GetDimensionLength(Dimension::Time));
}

GeoDomain
TransformHyperslabToGeoDomain(const Hyperslab& hyperslab)
{
    return GeoDomain(hyperslab.GetDimensionStart(Dimension::Lon),  hyperslab.GetDimensionStart(Dimension::Lon)  + hyperslab.GetDimensionLength(Dimension::Lon),
                     hyperslab.GetDimensionStart(Dimension::Lat),  hyperslab.GetDimensionStart(Dimension::Lat)  + hyperslab.GetDimensionLength(Dimension::Lat),
                     hyperslab.GetDimensionStart(Dimension::Lev),  hyperslab.GetDimensionStart(Dimension::Lev)  + hyperslab.GetDimensionLength(Dimension::Lev),
                     hyperslab.GetDimensionStart(Dimension::Time), hyperslab.GetDimensionStart(Dimension::Time) + hyperslab.GetDimensionLength(Dimension::Time));
}

}
