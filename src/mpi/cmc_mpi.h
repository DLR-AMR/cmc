#ifndef CMC_MPI_H
#define CMC_MPI_H

#include "cmc_config.h"
#include "utilities/cmc_constants_definitions.h"

#ifdef CMC_ENABLE_MPI
#include "mpi.h"
#endif

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
/* In case the Fortran interface is built, these function needs to be callable from C-wrapper scripts.
 * In order to prevent name-mangling, we have to declare them as extern "C"
 */
extern "C" {
#endif
#endif

/**
 * @brief Initialize MPI execution environment 
 * 
 * @note This functions only initializes MPI if it has not been initialized before
 */
void
cmc_mpi_initialize();

/**
 * @brief Terminate the MPI execution environment
 */
void
cmc_mpi_finalize();

/**
 * @brief Abort MPI if an error occured (on communicator MPI_COMM_WORLD) 
 * 
 * @param _err_code Error Code to return to the invoking environment
 * @param _location Location of the error which occured (File name and Line number)
 */
void cmc_mpi_abort(const int _err_code, const char* _location);


#define cmc_mpi_check_err(err) ((err) == MPI_SUCCESS ? (void) 0 : cmc_mpi_abort(err, CMC_FILE_LOCATION))

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
}
#endif
#endif

#endif /* CMC_MPI_H */
