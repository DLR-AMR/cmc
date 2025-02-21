#ifndef CMC_NETCDF_HXX
#define CMC_NETCDF_HXX
/**
 * @file cmc_netcdf.hxx
 * @brief Via the 'include' of @file cmc_netcdf.h, this file collects all functions used for accessing netCDF files and storing the data of the netCDF variables as well as the geo-spatial domain on which the variables are defined.
 * Additionally, this file supplies C++ only functions for inquiring variables from a netCDF file
 */

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_coordinate_array.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_input_variable.hxx"

#ifdef CMC_WITH_NETCDF
#include "netcdf.h"
#endif
#ifdef CMC_WITH_NETCDF_PAR
#include "netcdf_par.h"
#endif

#include <vector>
#include <string>

namespace cmc::nc
{

/* Enum describing an access mode for netCDF files */
enum OpeningMode {Serial = 0, Parallel};

[[noreturn]] void
NcExit(const int err, const std::string file, const int line);

inline void CheckError_ (int err, const std::string file, const int line)
{
    if (err != NC_NOERR)
    {
        NcExit(err, file, line);
    }
}

#define CheckError(err) CheckError_(err, __FILE__, __LINE__)

int
OpenSerial(const char* path_to_file);

int
OpenParallel(const char* path_to_file, MPI_Comm comm);

}

#endif /* CMC_NETCDF_HXX */
