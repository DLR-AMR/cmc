/** This is cmc - a compression tool for climate data sets */

#include "cmc.h"

/** Initialize cmc and it's sub-modules */
void
cmc_initialize()
{
  /* Initialize MPI */
  cmc_mpi_initialize();
  /* Initialize t8code */
  //TODO:
  //cmc_t8code_initialize(MPI_COMM_WORLD);
  #ifdef CMC_ENABLE_MPI
  int err, rank;
  /* Get the rank of the process */
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  cmc_mpi_check_err(err);
  /* Only print out the initilization info on rank zero */
  if (rank == 0)
  {
    std::cout << "[cmc] cmc has been initialized." << std::endl; 
  }
  #else
  std::cout << "[cmc] cmc has been initialized." << std::endl; 
  #endif
}

/** Initialize cmc and it's sub-modules with a specific MPI Communicator */
void
cmc_initialize_mpi_comm(MPI_Comm comm)
{
  /* Initialize MPI */
  cmc_mpi_initialize();
  /* Initialize t8code */
  //TODO:
  //ccmc_t8code_initialize(comm);
  #ifdef CMC_ENABLE_MPI
  int err, rank;
  /* Get the rank of the process */
  err = MPI_Comm_rank(comm, &rank);
  cmc_mpi_check_err(err);
  /* Only print out the initilization info on rank zero */
  if (rank == 0)
  {
    std::cout << "[cmc] cmc has been initialized." << std::endl;
  }
  #else
  std::cout << "[cmc] cmc has been initialized." << std::endl;
  #endif
}

/** Finalize cmc and it's submodules */
void
cmc_finalize()
{
  #ifdef CMC_ENABLE_MPI
  int err, rank;
  /* Get the rank of the process */
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  cmc_mpi_check_err(err);
  #endif

  /* Finalize t8code */
  //TODO:
  //ccmc_t8code_finalize(0);
  /* Finalize MPI */
  cmc_mpi_finalize();
  #ifdef CMC_ENABLE_MPI
  /* Only print out the initilization info on rank zero */
  if (rank == 0)
  {
    std::cout << "[cmc] cmc has been finalized." << std::endl;
  }
  #else
  /* Print out the finalizing info */
  std::cout << "[cmc] cmc has been finalized." << std::endl;
  #endif
}

/** Finalize cmc and it's submodules except MPI */
void
cmc_finalize_without_mpi()
{
  /* Finalize t8code */
  //TODO:
  //ccmc_t8code_finalize(0);
  #ifdef CMC_ENABLE_MPI
  int err, rank; 
  err = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  cmc_mpi_check_err(err);
  /* Only print out the finalizing info on rank zero of the MPI_WORLD_COMM */
  if (rank == 0)
  {
    std::cout << "[cmc] cmc has been finalized." << std::endl;
  }
  #else
  std::cout << "[cmc] cmc has been finalized." << std::endl;
  #endif
}
