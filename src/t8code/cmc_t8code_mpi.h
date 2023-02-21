#ifndef CMC_T8_MPI_H
#define CMC_T8_MPI_H

#include "mpi/cmc_mpi.h"

struct cmc_t8_mpi_data
{
    cmc_t8_mpi_data(){};
    cmc_t8_mpi_data(const int rank, const int num_vars)
    : partner_rank{rank}
    {
        data.reserve(num_vars);
    };
    ~cmc_t8_mpi_data(){};

    int partner_rank{-1};
    std::vector<var_dynamic_array_t*> data;
    std::vector<uint64_t> morton_indices; //Currently, it is only possible to perform parallel computations if all variables are defined on the same domain.
};

/* Typedefs for sending and receiving data */
typedef struct cmc_t8_mpi_data cmc_t8_mpi_send_data;
typedef struct cmc_t8_mpi_data cmc_t8_mpi_recv_data;


#endif /* CMC_T8_MPI_H */