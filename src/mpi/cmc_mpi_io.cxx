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


#if 0

#include "utilities/cmc_util.h"
#include "utilities/cmc_log_functions.hxx"

/* Calculate a blocked distribution based on supplied numbers for the amount of processes per dimension */
static
cmc_mpi_nc_data
cmc_mpi_calculate_nc_reading_data_distribution_num_procs_per_dim(const std::vector<uint64_t>& start_values, const std::vector<uint64_t>& count_values, const MPI_Comm comm, const std::vector<int>& num_procs_per_dim)
{
    int err, rank, size;

    /* Get the size of the communicator */
    err = MPI_Comm_size(comm, &size);
    cmc_mpi_check_err(err);

    /* Check if the blocked distribution given by the parameter @var num_procs_per_dim equals or is less than the size of the communicator */ 
    cmc_assert(std::reduce(num_procs_per_dim.begin(), num_procs_per_dim.end(), 1, std::multiplies<int>()) <= size);

    /* Get the (local) rank id */
    err = MPI_Comm_rank(comm, &rank);
    cmc_mpi_check_err(err);

    const int num_reading_procs = std::reduce(num_procs_per_dim.begin(), num_procs_per_dim.end(), 1, std::multiplies<int>());

    if (num_reading_procs < size)
    {
        /* We need to update the size if not all processes of the communicator are going to read data */
        size = num_reading_procs;
    }

    std::vector<uint64_t> start_vals;
    start_vals.reserve(start_values.size());
    std::vector<uint64_t> count_vals;
    count_vals.reserve(start_values.size());

    if (rank < num_reading_procs)
    {
        int current_dim_id{0};

        /* Variables for calculating the distribution */
        uint64_t global_offset{0};
        uint64_t box_length{0};

        int vrank = rank;
        int divisor = 1;
        uint64_t div = 1;

        /* Iterate over all dimensions and determine a blocked distribution */
        for (auto iter{count_values.begin()}; iter != count_values.end(); ++iter, ++current_dim_id)
        {
            cmc_assert(num_procs_per_dim[current_dim_id] <= size);

            /* If a dimension is not considered in the data hyperslab, just skip it */
            if (*iter <= 1UL)
            {
                /* These values indicate that we are not reading data from this dimension */
                start_vals.push_back(start_values[current_dim_id]);
                count_vals.push_back(1);
                continue;
            } else
            {
                if (num_procs_per_dim[current_dim_id] > 1)
                {
                    if (div > 1)
                    {
                        vrank = rank / std::pow(2, div - 1);
                    } else
                    {
                        vrank = rank;
                    }
                    divisor = num_procs_per_dim[current_dim_id];
                    global_offset = static_cast<uint64_t>((((double) (vrank % divisor) * (long double) (*iter)) / (double) divisor));

                    box_length = (vrank % divisor) + 1 >= num_procs_per_dim[current_dim_id] ? (*iter) - global_offset : static_cast<uint64_t>((((double) ((vrank % divisor) + 1) * (long double) (*iter)) / (double) divisor)) - global_offset;

                    ++div;
                } else
                {
                    global_offset = 0;
                    box_length = *iter;
                }

                /* Save the process-local start and count values for the current dimension */
                start_vals.push_back(start_values[current_dim_id] + global_offset);
                count_vals.push_back(box_length);
            }
        }

    }
    else
    {
        for (auto iter{count_values.begin()}; iter != count_values.end(); ++iter)
        {
            start_vals.push_back(0);
            count_vals.push_back((*iter != 1 ? 0 : 1));
        }
    }

    /* Return a struct dscribing the data which should be read in on each rank */
    return cmc_mpi_nc_data(std::move(start_vals), std::move(count_vals), data_distribution_t::CMC_BLOCKED);
}

cmc_mpi_nc_data
cmc_mpi_calculate_nc_reading_data_distribution_blocked(const std::vector<uint64_t>& start_values, const std::vector<uint64_t>& count_values, const MPI_Comm comm, const std::vector<int>& num_procs_per_dim)
{
    return cmc_mpi_calculate_nc_reading_data_distribution_num_procs_per_dim(start_values, count_values, comm, num_procs_per_dim);
}

#if 0
//Not implemented yet
cmc_mpi_nc_data
cmc_mpi_calculate_nc_reading_data_distribution_striped(const std::vector<uint64_t>& start_values, const std::vector<uint64_t>& count_values, const MPI_Comm comm, const std::vector<int>& num_procs_per_dim)
{
    return cmc_mpi_calculate_nc_reading_data_distribution_num_procs_per_dim(start_values, count_values, comm, num_procs_per_dim);
}

//Not implemented yet
cmc_mpi_nc_data
cmc_mpi_calculate_nc_reading_data_distribution_zcurve(const std::vector<uint64_t>& start_values, const std::vector<uint64_t>& count_values, const MPI_Comm comm, const std::vector<int>& num_procs_per_dim)
{
    ...
}
#endif


#endif
