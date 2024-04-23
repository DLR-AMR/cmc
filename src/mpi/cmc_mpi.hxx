#ifndef CMC_MPI_HXX
#define CMC_MPI_HXX
/**
 * @file cmc_mpi.hxx
 */

#include "utilities/cmc_utilities.hxx"

#ifdef CMC_ENABLE_MPI
#include "mpi.h"
#endif

namespace cmc
{

#ifdef CMC_ENABLE_MPI
#define MPI_MORTON_INDEX_T MPI_UINT64_T
#endif

/**
 * @brief Initialize MPI execution environment 
 * 
 * @note This functions only initializes MPI if it has not been initialized before
 */
void
MPIInitialize();

/**
 * @brief Terminate the MPI execution environment
 */
void
MPIFinalize();

/**
 * @brief Abort MPI if an error occured (on communicator MPI_COMM_WORLD) 
 * 
 * @param _err_code Error Code to return to the invoking environment
 * @param _location Location of the error which occured (File name and Line number)
 */
[[noreturn]] void
MPIAbort(const int _err_code, const char* _location);


//#define cmc_mpi_check_err(err) ((err) == MPI_SUCCESS ? (void) 0 : cmc_mpi_abort(err, CMC_FILE_LOCATION))

inline
void
MPICheckError(const int return_value)
{
    if (return_value != MPI_SUCCESS)
    {
        MPIAbort(return_value, CMC_FILE_LOCATION);
    }
}

}

#endif /* !CMC_MPI_HXX */
