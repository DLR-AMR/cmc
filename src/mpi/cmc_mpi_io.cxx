#include "mpi/cmc_mpi_io.hxx"
#include "utilities/cmc_log_functions.hxx"

namespace cmc
{

std::vector<Hyperslab>
DetermineSlicedDataDistribution(const Hyperslab& global_hyperslab, const MPI_Comm comm)
{
    #ifdef CMC_ENABLE_MPI
    int err, rank, size;

    /* Get the size of the communicator */
    err = MPI_Comm_size(comm, &size);
    MPICheckError(err);

    /* Get the local rank id */
    err = MPI_Comm_rank(comm, &rank);
    MPICheckError(err);

    /* Create a new hyperslab variable defining the local domain for each process */
    Hyperslab local_hyperslab;

    const double drank = static_cast<double>(rank);
    const double dsize = static_cast<double>(size);
    
    auto new_start_val_iter = local_hyperslab.StartIndicesBegin();
    auto new_count_val_iter = local_hyperslab.CountIndicesBegin();

    bool first_dimension_has_been_split = false;
    
    /* Try to find an equally blocked distribution of the hyperslabs */
    for (auto start_val_iter = global_hyperslab.StartIndicesBegin(),
         count_val_iter = global_hyperslab.CountIndicesBegin();
         start_val_iter != global_hyperslab.StartIndicesEnd();
         ++start_val_iter, ++count_val_iter, ++new_start_val_iter, ++new_count_val_iter)
    {
        /* If the dimension is considered, we split the first one between all processes */
        if (*count_val_iter > 1)
        {
            /* Split the first dimension equally between all processes */
            if (!first_dimension_has_been_split)
            {
                /* Calculate the start offset */
                const HyperslabIndex start_offset = *start_val_iter + static_cast<HyperslabIndex>((drank * static_cast<long double>(*count_val_iter)) / dsize);

                /* Calculate the succeeding start offset */
                const HyperslabIndex next_start_offset = *start_val_iter + static_cast<HyperslabIndex>(((drank + 1.0) * static_cast<long double>(*count_val_iter)) / dsize);

                /* Calculate the count value */
                const HyperslabIndex count_value = next_start_offset - start_offset;

                /* Set the new start and count value within the local hyperslab */
                *new_start_val_iter = start_offset;
                *new_count_val_iter = count_value;

                /* Set the flag */
                first_dimension_has_been_split = true; 
            } else
            { 
                /* Set the new start and count value within the local hyperslab */
                *new_start_val_iter = *start_val_iter;
                *new_count_val_iter = *count_val_iter;
            }
        }
    }

    std::vector<Hyperslab> local_hyperslabs;
    local_hyperslabs.emplace_back(std::move(local_hyperslab));

    return local_hyperslabs;
    #else
    return std::vector<Hyperslab>(1, global_hyperslab);
    #endif
}


}

