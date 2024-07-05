#ifndef CMC_MPI_IO_HXX
#define CMC_MPI_IO_HXX

#include "mpi/cmc_mpi.hxx"
#include "cmc.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_geo_domain.hxx"

#include <vector>

namespace cmc
{


std::vector<Hyperslab>
DetermineSlicedDataDistribution(const Hyperslab& global_hyperslab, const MPI_Comm comm);

}

#if 0

#include "cmc_config.h"
#include "utilities/cmc_constants_definitions.h"
#include "cmc_mpi.h"
#include "utilities/cmc_util.h"
#include "utilities/cmc_container.h"

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
/* In case the Fortran interface is built, these function needs to be callable from C-wrapper scripts.
 * In order to prevent name-mangling, we have to declare them as extern "C"
 */
extern "C" {
#endif
#endif

enum data_distribution_t {DISTRIBUTION_UNDEFINED = 0, CMC_SERIAL_DATA, CMC_LINEARIZED, CMC_BLOCKED, CMC_MESSY, CMC_ZCURVE};

struct cmc_mpi_nc_data
{
    cmc_mpi_nc_data(){};
    cmc_mpi_nc_data(std::vector<uint64_t>&& _start_values, std::vector<uint64_t>&& _count_values, const data_distribution_t _data_dist)
    : start_values{_start_values}, count_values{_count_values}, data_dist{_data_dist} {};

    int num_process_local_blocks{0};
    std::vector<uint64_t> start_values;
    std::vector<uint64_t> count_values;
    data_distribution_t data_dist{data_distribution_t::DISTRIBUTION_UNDEFINED};
};

cmc_mpi_nc_data
cmc_mpi_calculate_nc_reading_data_distribution(const std::vector<uint64_t>& start_values, const std::vector<uint64_t>& count_values, const MPI_Comm comm, const enum data_distribution_t preferred_data_distribution);

cmc_mpi_nc_data
cmc_mpi_calculate_nc_reading_data_distribution_blocked(const std::vector<uint64_t>& start_values, const std::vector<uint64_t>& count_values, const MPI_Comm comm, const std::vector<int>& num_procs_per_dim);

struct cmc_mpi_t8_data
{
    cmc_mpi_t8_data(){};
    cmc_mpi_t8_data(cmc_mpi_t8_data const&) = default;
    cmc_mpi_t8_data(cmc_mpi_t8_data&& other) noexcept :
        variable_indices(std::move(other.variable_indices)),
        data(std::move(other.data)), morton_indices(std::move(other.morton_indices)), received_morton_indices(std::move(other.received_morton_indices))
    {};
    cmc_mpi_t8_data& operator=(const cmc_mpi_t8_data&) = default;
    #if 0
    cmc_mpi_t8_data(const int rank, const int num_vars)
    : partner_rank{rank}
    {
        data.reserve(num_vars);
    };
    #endif
    ~cmc_mpi_t8_data(){};

    //Change some of the vectors (i.e. morton_indices) to C-style arrays since the vector is useless for receiving MPI objects
    std::vector<int> variable_indices;
    std::vector<var_dynamic_array_t*> data;
    std::vector<uint64_t> morton_indices; //Currently, it is only possible to perform parallel computations if all variables are defined on the same domain.
    int received_morton_indices{0}; //The number of received Morton indices
};

/* Typedefs for sending and receiving data */
typedef struct cmc_mpi_t8_data cmc_mpi_t8_send_data;
//typedef struct cmc_mpi_t8_data cmc_mpi_t8_recv_data;

struct cmc_mpi_t8_recv_data
{
    cmc_mpi_t8_recv_data(){};

    std::vector<std::vector<var_dynamic_array_t*>> data;
    std::vector<std::vector<int>> received_variable_ids;
    std::vector<uint64_t*> morton_indices;
    std::vector<int> received_morton_indices;

    std::vector<int> reference_group_id_sender;
};


#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
}
#endif
#endif


#endif

#endif /* CMC_MPI_IO_H */
