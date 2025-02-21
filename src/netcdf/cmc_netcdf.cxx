#include "netcdf/cmc_netcdf.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "mpi/cmc_mpi_io.hxx"

#include <utility>

namespace cmc::nc
{

[[noreturn]] void
NcExit(const int _err_code, const std::string file, const int line)
{
    std::cout << "CMC_NETCDF_EXIT is invoked..." << std::endl << "A netCDF-Error occured, Code: " << _err_code << std::endl << nc_strerror(_err_code) << std::endl << "In: " << file << ": " << line << std::endl;
    std::exit(EXIT_FAILURE);
}

int
OpenSerial(const char* path_to_file)
{
    #ifdef CMC_WITH_NETCDF
    int err;
    int ncid;

    /* Open the file without explicit parallel access */
    err = nc__open(path_to_file, NC_NOWRITE, NULL, &ncid);
    CheckError(err);

    return ncid;
    #else
    cmc_err_msg("CMC is not linked against netCDF. Please recompile the build with netCDF\n.");
    return CMC_ERR;
    #endif
}

int
OpenParallel(const char* path_to_file, MPI_Comm comm)
{
    #ifdef CMC_WITH_NETCDF_PAR
    int err;
    int ncid;

    int comm_size;
    err = MPI_Comm_size(comm, &comm_size);
    MPICheckError(err);

    MPI_Info info = MPI_INFO_NULL;
    err = nc_open_par(path_to_file, NC_NOWRITE, comm, info, &ncid);
    CheckError(err);


    return ncid;
    #else
    cmc_err_msg("CMC is not linked against netCDF (at leat no parallel access functions are available). Please recompile the build with netCDF\n.");
    return CMC_ERR;
    #endif
}

}
