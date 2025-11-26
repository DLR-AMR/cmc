#ifndef CMC_GEO_UTILITIES_HXX
#define CMC_GEO_UTILITIES_HXX

#include "utilities/cmc_dimension_interval.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_hyperslab.hxx"

#include <vector>

namespace cmc
{

Hyperslab TransformGeoDomainToHyperslab(const GeoDomain& domain);

GeoDomain TransformHyperslabToGeoDomain(const Hyperslab& hyperslab);

Hyperslab SubtractGeoDomainOffset(const Hyperslab& hyperslab, const GeoDomain& domain);

GeoDomain GetDefaultDomain(const DataLayout layout, const std::vector<size_t>& dimension_lengths);

}

#endif /* !CMC_GEO_UTILITIES_HXX */
