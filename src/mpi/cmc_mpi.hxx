#ifndef CMC_MPI_HXX
#define CMC_MPI_HXX
/**
 * @file cmc_mpi.hxx
 */

#include "utilities/cmc_utilities.hxx"

#ifdef CMC_ENABLE_MPI
/* We define a macro to skip the inclusion of MPI-CXX interface when OpenMPI is used (but buitl with CXX support) */
#ifndef OMPI_SKIP_MPICXX
#define OMPI_SKIP_MPICXX
#endif
/* Include the MPi-Interface*/
#include "mpi.h"
#endif

namespace cmc
{

#ifdef CMC_ENABLE_MPI
#define MPI_MORTON_INDEX_T MPI_INT64_T
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

/**
 * @brief An MPI error function checking the return value of MPI functions and 
 * issues an error if the MPI operation was not successfull.
 */
#define MPICheckError(err) ((err) == MPI_SUCCESS ? (void) 0 : MPIAbort(err, CMC_FILE_LOCATION))

}

#endif /* !CMC_MPI_HXX */
