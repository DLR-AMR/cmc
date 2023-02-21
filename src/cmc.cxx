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
  cmc_t8code_initialize();
}

/** Finalize cmc and it's submodules */
void
cmc_finalize()
{
  /* Finalize t8code */
  cmc_t8code_finalize(0);
  /* Finalize MPI */
  cmc_mpi_finalize();
}

/** Finalize cmc and it's submodules except MPI */
void
cmc_finalize_without_mpi()
{
  /* Finalize t8code */
  cmc_t8code_finalize(0);
}