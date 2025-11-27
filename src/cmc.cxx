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

void
CmcInitialize(const bool initilization)
{
    #ifdef CMC_ENABLE_MPI
    MPIInitialize();

    int err, rank;
    /* Get the rank of the process */
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPICheckError(err);

    #ifdef CMC_WITH_T8CODE
    if (initilization)
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

}

void
CmcFinalize(const bool finalization)
{
    #ifdef CMC_ENABLE_MPI
    int err, rank;
    /* Get the rank of the process */
    err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPICheckError(err);
  
    #ifdef CMC_WITH_T8CODE
    if (finalization)
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
}

}
