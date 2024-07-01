#ifndef CMC_HXX
#define CMC_HXX

/* Include the configuration macros describing the linkage status to other libraries */
#include "cmc_config.h"
#include "mpi/cmc_mpi.hxx"

namespace cmc
{

void CmcInitialize();

void CmcFinalize();

}

#endif /* CMC_HXX */
