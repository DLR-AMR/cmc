#include "cmc_mpi_io.h"
#include "utilities/cmc_util.h"
#include "utilities/cmc_log_functions.h"

static
cmc_mpi_nc_data
cmc_mpi_calculate_nc_reading_data_distribution_blocked_default(const std::vector<uint64_t>& start_values, const std::vector<uint64_t>& count_values, const MPI_Comm comm)
{
    int err, rank, size;

    /* Get the size of the communicator */
    err = MPI_Comm_size(comm, &size);
    cmc_mpi_check_err(err);
    
    /* Get the (local) rank id */
    err = MPI_Comm_rank(comm, &rank);
    cmc_mpi_check_err(err);

    /* Variables for calculating the distribution */
    uint64_t leftover_dim_points{0}, evenly_split_dim_points{0};
    uint64_t global_offset{0};

    /* Flag indicating that one diemnsion was already not split (in order to cover the whole domain) */
    bool non_splitted_dim_applied{false};

    std::vector<uint64_t> start_vals;
    std::vector<uint64_t> count_vals;

    int current_dim_id{0};

    /* Iterate over all dimensions and determine a blocked distribution */
    for (auto iter{count_values.begin()}; iter != count_values.end(); ++iter, ++current_dim_id)
    {
        /* If a dimension is not considered in the data hyperslab, just skip it */
        if (*iter <= 1UL)
        {
            /* These values indicate that we are not reading data from this dimension */
            start_vals.push_back(0);
            count_vals.push_back(1);
            continue;
        } else
        {
            //Currently, there is only one option to calculate a blocked distribution, therefore there are some assertion 
            cmc_assert(*iter >= static_cast<size_t>(size));
            if (non_splitted_dim_applied)
            {
                /* If the dimension is considered, calculate an offset for this dimension */
                //The assertion above indicates that at least one dimension point is assigned to each rank
                /* The leftover dimension points will be distributed */
                evenly_split_dim_points = static_cast<uint64_t>((*iter) / size);
                /* Calculate the offset for the start vector */
                global_offset = rank * evenly_split_dim_points;
                /* Calculate the points which cannot be evenly split between all processses */
                leftover_dim_points = (*iter) % size;
                /* Adjust the offset for the leftover points */
                /* Rank 0 starts at position zero nevertheless, therefore, the offset cannot change */
                if (rank != 0)
                {
                    if (rank < static_cast<int>(leftover_dim_points))
                    {
                        /* The ranks (starting from the lowest to the highest id) obtain another additional dimension point when the diemnsion cannot be split evenly */
                        ++evenly_split_dim_points;
                        /* Adjust the global offset */
                        global_offset += leftover_dim_points - rank;
                    } else
                    {
                        /* Adjust the global offset */
                        global_offset += leftover_dim_points;
                    }
                }

                /* Save the start and count values for this dimension */
                start_vals.push_back(start_values[current_dim_id] + global_offset);
                count_vals.push_back(evenly_split_dim_points);

            }
            else
            {
                /* Just copy the whole dimension length */
                start_vals.push_back(start_values[current_dim_id]);
                count_vals.push_back(*iter);

                /* The whole dimensions are already saved by default in start and count pointers */
                non_splitted_dim_applied = true;
            }
        }
    }
    std::cout << "Calculated mpi offsets netcdf\n";
    for (size_t i{0}; i < count_vals.size(); ++i)
    {
        std::cout << "start_val (" << i << "): " << start_vals[i] << " und count val: " << count_vals[i] << std::endl;
    }
    /* Return a struct dscribing the data which should be read in on each rank */
    return cmc_mpi_nc_data(std::move(start_vals), std::move(count_vals), data_distribution_t::CMC_BLOCKED);
}

/* Calculate a blocked distribution based on supplied numbers for the amount of processes per dimension */
static
cmc_mpi_nc_data
cmc_mpi_calculate_nc_reading_data_distribution_num_procs_per_dim(const std::vector<uint64_t>& start_values, const std::vector<uint64_t>& count_values, const MPI_Comm comm, const std::vector<int>& num_procs_per_dim)
{
    int err, rank, size;

    /* Get the size of the communicator */
    err = MPI_Comm_size(comm, &size);
    cmc_mpi_check_err(err);

    /* Check if the blocked distribution given by the parameter @var num_procs_per_dim equals the size of the communicator */ 
    cmc_assert(std::reduce(num_procs_per_dim.begin(), num_procs_per_dim.end(), 1, std::multiplies<int>()) == size);

    /* Get the (local) rank id */
    err = MPI_Comm_rank(comm, &rank);
    cmc_mpi_check_err(err);

    std::vector<uint64_t> start_vals;
    std::vector<uint64_t> count_vals;

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


cmc_mpi_nc_data
cmc_mpi_calculate_nc_reading_data_distribution(const std::vector<uint64_t>& start_values, const std::vector<uint64_t>& count_values, const MPI_Comm comm, const enum data_distribution_t preferred_data_distribution)
{
    switch(preferred_data_distribution)
    {
        case CMC_BLOCKED:
            return cmc_mpi_calculate_nc_reading_data_distribution_blocked_default(start_values, count_values, comm);
        break;
        default:
            cmc_err_msg("An not yet implemented or unknown data distribution type was supplied.");
    }

    /* Empty return value in order to prevent issues during compilation */
    return cmc_mpi_nc_data();
}
