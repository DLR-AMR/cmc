/** This is cmc - a compression tool for climate data sets */

#include "cmc.h"
#include "mpi/cmc_mpi.h"
#include "t8code/cmc_t8code.h"

/** Initialize cmc and it's sub-modules */
void
cmc_initialize()
{
  /* Initialize MPI */
  cmc_mpi_initialize();
  /* Initialize t8code */
  cmc_t8code_initialize(MPI_COMM_WORLD);
  std::cout << "[CMC] cmc has been initialized." << std::endl; 
}

/** Initialize cmc and it's sub-modules with a specific MPI Communicator */
void
cmc_initialize_mpi_comm(MPI_Comm comm)
{
  /* Initialize MPI */
  cmc_mpi_initialize();
  /* Initialize t8code */
  cmc_t8code_initialize(comm);
  std::cout << "[CMC] cmc has been initialized." << std::endl; 
}

/** Finalize cmc and it's submodules */
void
cmc_finalize()
{
  /* Finalize t8code */
  cmc_t8code_finalize(0);
  /* Finalize MPI */
  cmc_mpi_finalize();
  std::cout << "[CMC] cmc has been finalized." << std::endl;
}

/** Finalize cmc and it's submodules except MPI */
void
cmc_finalize_without_mpi()
{
  /* Finalize t8code */
  cmc_t8code_finalize(0);
  std::cout << "[CMC] cmc has been finalized." << std::endl;
}
