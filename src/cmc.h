#ifndef CMC_H
#define CMC_H

/* Include the configuration macros describing the linkage status to other libraries */
#include "cmc_config.h"
#include "mpi/cmc_mpi.h"

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
/* In case the Fortran interface is built, these function needs to be callable from C-wrapper scripts.
 * In order to prevent name-mangling, we have to declare them as extern "C"
 */
extern "C" {
#endif
#endif

/**
 * @brief Initialize cmc by initializing all sub-modules (e.g. MPI, t8code, etc..)
 * @note If MPI is enabled, the MPI Communicator to use, will be MPI_COMM_WORLD
 */
void
cmc_initialize();

/**
 * @brief Initialize cmc by initializing all sub-modules (e.g. MPI, t8code, etc..) with a specific MPI Communicator
 * 
 * @param comm The communicator to be used by cmc, t8code, etc.
 */
void
cmc_initialize_mpi_comm(MPI_Comm comm);

/**
 * @brief Finalize cmc by finalizing all sub-modules (e.g. MPI, t8code, etc..)
 */
void
cmc_finalize();

/**
 * @brief Finalize cmc by finalizing all sub-modules (e.g. t8code, etc..) except for MPI.
 */
void
cmc_finalize_without_mpi();

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
}
#endif
#endif

#endif /* CMC_H */
