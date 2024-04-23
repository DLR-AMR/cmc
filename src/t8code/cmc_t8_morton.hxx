#ifndef CMC_T8_MORTON_HXX
#define CMC_T8_MORTON_HXX
/**
 * @file cmc_t8_morton.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_hyperslab.hxx"

namespace cmc
{

MortonIndex
GetMortonIndex(const std::vector<DomainIndex>& coordinates, const int dimensionality);

std::vector<MortonIndex>
TransformHyperslabCoordinatesToMortonIndices(const std::vector<Hyperslab>& hyperslabs, const DataLayout layout, const GeoDomain& global_domain);

}

#endif /* !CMC_T8_MORTON_HXX */
