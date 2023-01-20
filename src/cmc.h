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

void
cmc_initialize();

void
cmc_finalize();

void
cmc_finalize_without_mpi();

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
}
#endif
#endif

#endif /* CMC_H */