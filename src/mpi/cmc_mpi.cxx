#include "mpi/cmc_mpi.hxx"

namespace cmc
{

/** Initialize MPI (if has not been initialized beforehand) */
void
MPIInitialize()
{
    #ifdef CMC_ENABLE_MPI
    int init_flag{0};

    /* Check if MPI already has been initialized */
    int err{MPI_Initialized(&init_flag)};
    MPICheckError(err);

    if (init_flag == 0)
    {
        err = MPI_Init(NULL, NULL);
        MPICheckError(err);
    }
    #endif
}

/** Finalize MPI */
void
MPIFinalize()
{
    #ifdef CMC_ENABLE_MPI
    int err{MPI_Finalize()};
    MPICheckError(err);
    #endif
}

/** Abort MPI */
[[noreturn]]
void
MPIAbort(const int _err_code, const char* _location)
{
    #ifdef CMC_ENABLE_MPI
    int rank{0};
    int err{MPI_Comm_rank(MPI_COMM_WORLD, &rank)};
    if (err == MPI_SUCCESS)
    {
        if (rank == 0)
        {
            std::cout << "CMC_MPI_ABORT is invoked..." << std::endl << "The error-code: " << _err_code << " occured in: " << _location;
        }
    } else
    {
        std::cout << "CMC_MPI_ABORT is invoked..." << std::endl << "The error-code: " << _err_code << " occured in: " << CMC_FILE_LOCATION;
    }
    err = MPI_Abort(MPI_COMM_WORLD, _err_code);
    if (err == MPI_SUCCESS)
    {
        std::exit(EXIT_FAILURE);
    } else {
        std::abort();
    }
    #endif
}


}











#include "cmc.h"
#include "cmc_mpi.h"


/** Initialize MPI (if has not been initialized beforehand) */
void
cmc_mpi_initialize()
{
    #ifdef CMC_ENABLE_MPI
    int init_flag{0};

    /* Check if MPI already has been initialized */
    int err{MPI_Initialized(&init_flag)};
    cmc_mpi_check_err(err);

    if (init_flag == 0)
    {
        err = MPI_Init(NULL, NULL);
        cmc_mpi_check_err(err);
    }
    #endif
}

/** Finalize MPI */
void
cmc_mpi_finalize()
{
    #ifdef CMC_ENABLE_MPI
    int err{MPI_Finalize()};
    cmc_mpi_check_err(err);
    #endif
}

/** Abort MPI */
void
cmc_mpi_abort(const int _err_code, const char* _location)
{
    #ifdef CMC_ENABLE_MPI
    int rank{0};
    int err{MPI_Comm_rank(MPI_COMM_WORLD, &rank)};
    if (err == MPI_SUCCESS)
    {
        if (rank == 0)
        {
            std::cout << "CMC_MPI_ABORT is invoked..." << std::endl << "An MPI-Error with Error-Code: " << _err_code << " occured in: " << _location;
        }
    } else
    {
        std::cout << "CMC_MPI_ABORT is invoked..." << std::endl << "An MPI-Error with Error-Code: " << _err_code << " occured in: " << _location;
    }
    err = MPI_Abort(MPI_COMM_WORLD, _err_code);
    if (err == MPI_SUCCESS)
    {
        std::exit(EXIT_FAILURE);
    } else {
        std::abort();
    }
    #endif
}
