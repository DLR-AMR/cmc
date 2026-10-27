#include "mpi/cmc_mpi.hxx"
#include "utilities/cmc_log_functions.hxx"

namespace cmc
{

/** Initialize MPI (if has not been initialized beforehand) */
void
MPIInitialize()
{
#if CMC_ENABLE_MPI
    int init_flag{0};

    /* Check if MPI already has been initialized */
    const int rv_is_init = MPI_Initialized(&init_flag);
    MPICheckError(rv_is_init);

    if (init_flag == 0)
    {
        const int rv_init = MPI_Init(NULL, NULL);
        MPICheckError(rv_init);
    }
#else
    cmc_warn_msg("MPIInitialize is called although CMC is not configured with MPI!");
#endif
}

/** Finalize MPI */
void
MPIFinalize()
{
#if CMC_ENABLE_MPI
    int finalized_flag{0};
    const int rv_is_fin = MPI_Finalized(&finalized_flag);
    MPICheckError(rv_is_fin);
    if (finalized_flag == 0)
    {
        const int rv_fin = MPI_Finalize();
        MPICheckError(rv_fin);
    }
#else
    cmc_warn_msg("MPIFinalize is called although CMC is not configured with MPI!");
#endif
}

/** Abort MPI */
[[noreturn]]
void
MPIAbort(const int _err_code, const char* _location)
{
#if CMC_ENABLE_MPI
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
        std::cout << "CMC_MPI_ABORT is invoked..." << std::endl << "The error-code: " << _err_code << " occured in: " << _location;
    }
    err = MPI_Abort(MPI_COMM_WORLD, _err_code);
    if (err == MPI_SUCCESS)
    {
        std::exit(EXIT_FAILURE);
    } else
#endif
    {
        std::cout << "CMC_MPI_ABORT is invoked..." << std::endl << "The error-code: " << _err_code << " occured in: " << _location << std::endl;
        std::abort();
    }
}


}
