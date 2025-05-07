#ifndef CMC_T8_ADAPTATION_CALLBACKS_HXX
#define CMC_T8_ADAPTATION_CALLBACKS_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_vector_view.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
//#include <t8_cmesh.h>
//include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_forest/t8_forest_general.h>
//#include <t8_schemes/t8_scheme.hxx> 
//#include <t8_schemes/t8_default/t8_default.hxx>
//#include <t8_forest/t8_forest_iterate.h>
//#include <p4est.h>
//#include <p8est.h>
#endif

#include <vector>
#include <functional>

namespace cmc::t8
{

//using AdaptationFn = t8_forest_adapt_t;
typedef t8_forest_adapt_t AdaptationFn;

/* Helper functions for return values during the t8code adaptation call */
constexpr t8_locidx_t kCoarsenElements = -1;
constexpr t8_locidx_t kRefineElement = 1;
constexpr t8_locidx_t kLeaveElementUnchanged = 0;

}

#endif /* !CMC_T8_ADAPTATION_CALLBACKS_HXX */
