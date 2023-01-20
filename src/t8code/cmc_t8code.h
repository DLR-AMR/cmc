#ifndef CMC_T8CODE_H
#define CMC_T8CODE_H

#include "cmc_config.h"

/* Include some t8code dependent header files */
#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest.h>
#endif


#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
/* In case the Fortran interface is built, these function needs to be callable from C-wrapper scripts.
 * In order to prevent name-mangling, we have to declare them as extern "C"
 */
extern "C" {
#endif
#endif

void
cmc_t8code_initialize();

void
cmc_t8code_finalize(const int flag_finalize_mpi);

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
}
#endif
#endif

#endif /* CMC_T8CODE_H */