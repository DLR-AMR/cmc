#ifndef CMC_HXX
#define CMC_HXX

/* Include the configuration macros describing the linkage status to other libraries */
#include "cmc_config.h"

namespace cmc
{
/* The minimum initilization only initializes MPI if available, the default initialization additionally initializes SC and t8code if available. */
constexpr bool kMinimumInitialization = false;
constexpr bool kDefaultInitialization = true;

void CmcInitialize(const bool initialization = kDefaultInitialization);

void CmcFinalize();

}

#endif /* CMC_HXX */
