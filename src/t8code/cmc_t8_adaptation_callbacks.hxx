#ifndef CMC_T8_ADAPTATION_CALLBACKS_HXX
#define CMC_T8_ADAPTATION_CALLBACKS_HXX

#include "utilities/cmc_utilities.hxx"

#if CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest_general.h>
#endif

namespace cmc::t8
{

typedef t8_forest_adapt_t AdaptationFn;

/* Helper functions for return values during the t8code adaptation call */
constexpr t8_locidx_t kCoarsenElements = -1;
constexpr t8_locidx_t kRefineElement = 1;
constexpr t8_locidx_t kLeaveElementUnchanged = 0;

}

#endif /* !CMC_T8_ADAPTATION_CALLBACKS_HXX */
