#ifndef CMC_H
#define CMC_H

/* Include the configuration macros describing the linkage status to other libraries */
#include "cmc_config.h"

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
/* In case the Fortran interface is built, these function needs to be callable from C-wrapper scripts.
 * In order to prevent name-mangling, we have to declare them as extern "C"
 */
extern "C" {
#endif
#endif

/** \fn void cmc_initialize()
  * Initialize cmc by initializing all sub-modules (e.g. MPI, t8code, etc..)
*/
void
cmc_initialize();

/** \fn void cmc_finalize()
  * Finalize cmc by finalizing all sub-modules (e.g. MPI, t8code, etc..)
*/
void
cmc_finalize();

/** \fn void cmc_finalize()
  * Finalize cmc by finalizing all sub-modules (e.g. t8code, etc..) except for MPI.
  * \note This function may be particularily important when cmc is initialized and finalized as a sub-module, for example within an external software package 
*/
void
cmc_finalize_without_mpi();

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
}
#endif
#endif

#endif /* CMC_H */