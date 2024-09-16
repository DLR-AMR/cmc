#ifndef CMC_GEO_UTILITIES_HXX
#define CMC_GEO_UTILITIES_HXX

#include "utilities/cmc_dimension_interval.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_hyperslab.hxx"


namespace cmc
{

Hyperslab TransformGeoDomainToHyperslab(const GeoDomain& domain);

GeoDomain TransformHyperslabToGeoDomain(const Hyperslab& hyperslab);

Hyperslab SubtractGeoDomainOffset(const Hyperslab& hyperslab, const GeoDomain& domain);

}

#endif /* !CMC_GEO_UTILITIES_HXX */
