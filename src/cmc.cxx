/** This is cmc - a compression tool for climate data sets */

#include "cmc.hxx"

#include <iostream>

namespace cmc
{

void
CmcInitialize()
{
  MPIInitialize();

  //TODO: Initialize t8code

  #ifdef CMC_ENABLE_MPI
  int err, rank;
  /* Get the rank of the process */
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPICheckError(err);
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
CmcFinalize()
{
  #ifdef CMC_ENABLE_MPI
  int err, rank;
  /* Get the rank of the process */
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPICheckError(err);
  #endif

  //TODO: Finalize t8code

  MPIFinalize();

  #ifdef CMC_ENABLE_MPI
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
