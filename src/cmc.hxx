#ifndef CMC_HXX
#define CMC_HXX

/* Include the configuration macros describing the linkage status to other libraries */
#include "cmc_config.h"

namespace cmc
{
/* The minimum initilization only initializes MPI if available, the default initialization additionally initializes SC and t8code if available.
 * The same holds for the finilizaiton, therefore, those funcitons should be called consistently. */
constexpr bool kMinimumInitialization = false;
constexpr bool kDefaultInitialization = true;
constexpr bool kMinimumFinalization = false;
constexpr bool kDefaultFinalization = true;

void CmcInitialize(const bool initilization = kDefaultInitialization);

void CmcFinalize(const bool finalization = kDefaultFinalization);

}

#endif /* CMC_HXX */
