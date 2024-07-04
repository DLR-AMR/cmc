#include "mpi/cmc_mpi_io.hxx"

namespace cmc
{


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
