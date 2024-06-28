#ifndef CMC_LOG_FUNCTIONS_H
#define CMC_LOG_FUNCTIONS_H
#include "cmc.hxx"
#include "mpi/cmc_mpi.hxx"

namespace cmc
{

template<typename T>
void
cmc_print_args(const T& msg)
{
    std::cout << msg << std::endl;
}

template<typename T, typename... Ts>
void
cmc_print_args(const T& msg, Ts... msgs)
{
    std::cout << msg;
    cmc_print_args(msgs...);
}

template<typename... Ts>
void
cmc_global_msg(MPI_Comm comm, Ts... msgs)
{
    /* check for MPI_Comm_rank */
    #ifdef CMC_ENABLE_MPI
        int rank{0};
        int err{MPI_Comm_rank(comm, &rank)};
        MPICheckError(err);
        if (rank == 0)
        {
            std::cout << "[cmc] ";
            cmc_print_args(msgs...);
        }
    #else
        std::cout << "[cmc] ";
        cmc_print_args(msgs...);
    #endif
}

template<typename... Ts>
void
cmc_global_msg(Ts... msgs)
{
    /* check for MPI_Comm_rank */
    #ifdef CMC_ENABLE_MPI
        int rank{0};
        int err{MPI_Comm_rank(MPI_COMM_WORLD, &rank)};
        MPICheckError(err);
        if (rank == 0)
        {
            std::cout << "[cmc] ";
            cmc_print_args(msgs...);
        }
    #else
        std::cout << "[cmc] ";
        cmc_print_args(msgs...);
    #endif
}

template<typename... Ts>
void
cmc_debug_msg(MPI_Comm comm, Ts... msgs)
{
     #ifdef CMC_ENABLE_MPI
        int rank{0}, size{0};
        int err{MPI_Comm_rank(comm, &rank)};
        MPICheckError(err);
        err = MPI_Comm_size(comm, &size);
        MPICheckError(err);
        if (size > 1)
        {
            std::cout << "[cmc] [rank " << rank << "] DEBUG: ";
        } else {
            std::cout << "[cmc] DEBUG: ";
        }
        cmc_print_args(msgs...);
    #else
        std::cout << "[cmc] DEBUG: ";
        cmc_print_args(msgs...);
    #endif
}

template<typename... Ts>
void
cmc_debug_msg(Ts... msgs)
{
     #ifdef CMC_ENABLE_MPI
        int rank{0}, size{0};
        int err{MPI_Comm_rank(MPI_COMM_WORLD, &rank)};
        MPICheckError(err);
        err = MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPICheckError(err);
        if (size > 1)
        {
            std::cout << "[cmc] [rank " << rank << "] DEBUG: ";
        } else {
            std::cout << "[cmc] DEBUG: ";
        }
        cmc_print_args(msgs...);
    #else
        std::cout << "[cmc] DEBUG: ";
        cmc_print_args(msgs...);
    #endif
}

template<typename... Ts>
void
cmc_msg(MPI_Comm comm, Ts... msgs)
{
    #ifdef CMC_ENABLE_MPI
        int rank{0}, size{0};
        int err{MPI_Comm_rank(comm, &rank)};
        MPICheckError(err);
        err = MPI_Comm_size(comm, &size);
        MPICheckError(err);
        if (size > 1)
        {
            std::cout << "[cmc] [rank " << rank << "] ";
        } else {
            std::cout << "[cmc] ";
        }
        cmc_print_args(msgs...);
    #else
        std::cout << "[cmc] ";
        cmc_print_args(msgs...);
    #endif
}

template<typename... Ts>
void
cmc_msg(Ts... msgs)
{
    #ifdef CMC_ENABLE_MPI
        int rank{0}, size{0};
        int err{MPI_Comm_rank(MPI_COMM_WORLD, &rank)};
        MPICheckError(err);
        err = MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPICheckError(err);
        if (size > 1)
        {
            std::cout << "[cmc] [rank " << rank << "] ";
        } else {
            std::cout << "[cmc] ";
        }
        cmc_print_args(msgs...);
    #else
        std::cout << "[cmc] ";
        cmc_print_args(msgs...);
    #endif
}

template<typename... Ts>
[[noreturn]] void
cmc_err_msg(MPI_Comm comm, Ts... msgs)
{
    #ifdef CMC_ENABLE_MPI
        int rank{0};
        int err{MPI_Comm_rank(comm, &rank)};
        MPICheckError(err);
        std::cout << "[cmc] [rank " << rank << "] ERROR: ";
        cmc_print_args(msgs...);
        err = MPI_Abort(comm, MPI_ERR_OTHER);
        MPICheckError(err);
        exit(EXIT_FAILURE);
    #else
        std::cout << "[cmc] ERROR: ";
        cmc_print_args(msgs...);
        std::exit(EXIT_FAILURE);
    #endif
}

template<typename... Ts>
[[noreturn]] void
cmc_err_msg(Ts... msgs)
{
    #ifdef CMC_ENABLE_MPI
        int rank{0};
        int err{MPI_Comm_rank(MPI_COMM_WORLD, &rank)};
        MPICheckError(err);
        std::cout << "[cmc] [rank " << rank << "] ERROR: ";
        cmc_print_args(msgs...);
        err = MPI_Abort(MPI_COMM_WORLD, MPI_ERR_OTHER);
        MPICheckError(err);
        exit(EXIT_FAILURE);
    #else
        std::cout << "[cmc] ERROR: ";
        cmc_print_args(msgs...);
        std::exit(EXIT_FAILURE);
    #endif
    
}


template<typename... Ts>
void
cmc_warn_msg(MPI_Comm comm, Ts... msgs)
{
    #ifdef CMC_ENABLE_MPI
        int rank{0}, size{0};
        int err{MPI_Comm_rank(comm, &rank)};
        MPICheckError(err);
        err = MPI_Comm_size(comm, &size);
        MPICheckError(err);
        if (size > 1)
        {
            std::cout << "[cmc] [rank " << rank << "] WARNING: ";
        } else {
            std::cout << "[cmc] WARNING:";
        }
        cmc_print_args(msgs...);
    #else
        std::cout << "[cmc] WARNING: ";
        cmc_print_args(msgs...);
    #endif
}

template<typename... Ts>
void
cmc_warn_msg(Ts... msgs)
{
    #ifdef CMC_ENABLE_MPI
        int rank{0}, size{0};
        int err{MPI_Comm_rank(MPI_COMM_WORLD, &rank)};
        MPICheckError(err);
        err = MPI_Comm_size(MPI_COMM_WORLD, &size);
        MPICheckError(err);
        if (size > 1)
        {
            std::cout << "[cmc] [rank " << rank << "] WARNING: ";
        } else {
            std::cout << "[cmc] WARNING: ";
        }
        cmc_print_args(msgs...);
    #else
        std::cout << "[cmc] WARNING: ";
        cmc_print_args(msgs...);
    #endif
}

}

#endif /* CMC_LOG_FUNCTIONS_H */
