#ifndef CMC_DIMENSION_INTERVAL_HXX
#define CMC_DIMENSION_INTERVAL_HXX

#include "utilities/cmc_utilities.hxx"

namespace cmc
{

struct DimensionInterval
{
    constexpr DimensionInterval(Dimension dimension, DomainIndex start_idx, DomainIndex end_idx)
    : dim{dimension}, start_index{start_idx}, end_index{end_idx}{};

    Dimension dim{Dimension::DimensionUndefined};
    DomainIndex start_index{0};
    DomainIndex end_index{0};
};


}

#endif /* !CMC_DIMENSION_INTERVAL_HXX */
