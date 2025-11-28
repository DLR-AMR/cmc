/** This is cmc - a compression tool for climate data sets */

#include "cmc.hxx"

#ifdef CMC_ENABLE_MPI
#include "mpi/cmc_mpi.hxx"
#endif

#ifdef CMC_WITH_T8CODE
#include "t8.h"
#endif

#include <iostream>

namespace cmc
{

static bool CmcIsInitialized = false;
static bool CmcIsFinalized = false;

static bool kInitializationScheme;

void
CmcInitialize(const bool initialization)
{
    if (CmcIsInitialized)
    {
        std::cout << "[cmc] WARNING: cmc has already been initialized." << std::endl;
        return;
    }
    #ifdef CMC_ENABLE_MPI
    MPIInitialize();

    int err, rank;
    /* Get the rank of the process */
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPICheckError(err);

    #ifdef CMC_WITH_T8CODE
    if (initialization)
    {
      sc_init (MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
      t8_init (SC_LP_DEBUG);
    }
    #endif

    /* Only print out the initilization info on rank zero */
    if (rank == 0)
    {
      std::cout << "[cmc] cmc has been initialized." << std::endl; 
    }
    #else
    std::cout << "[cmc] cmc has been initialized." << std::endl; 
    #endif

    CmcIsInitialized = true;
    kInitializationScheme = initialization;
}

void
CmcFinalize()
{
    if (not CmcIsInitialized || CmcIsFinalized)
    {
        std::cout << "[cmc] WARNING: cmc has not been initialized before or is already finalized." << std::endl; 
        return;
    }

    #ifdef CMC_ENABLE_MPI
    int err, rank;
    /* Get the rank of the process */
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPICheckError(err);
  
    #ifdef CMC_WITH_T8CODE
    if (kInitializationScheme)
    {
       sc_finalize();
    }
    #endif
  
    MPIFinalize();
  
    /* Only print out the initilization info on rank zero */
    if (rank == 0)
    {
      std::cout << "[cmc] cmc has been finalized." << std::endl;
    }
    #else
    std::cout << "[cmc] cmc has been finalized." << std::endl;
    #endif

    CmcIsFinalized = true;
    CmcIsInitialized = false;
}

}
