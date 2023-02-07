#ifndef CMC_NETCDF_H
#define CMC_NETCDF_H

#include "cmc.h"
#include "mpi/cmc_mpi.h"

#ifdef CMC_WITH_NETCDF
#include "netcdf.h"
#endif
#ifdef CMC_WITH_NETCDF_PAR
#include "netcdf_par.h"
#endif

#define CMC_NC_MAX_VAR_DIMS NC_MAX_VAR_DIMS
#define CMC_NC_VAR_NOT_CONSIDERED -1

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
/* In case the Fortran interface is built, these function needs to be callable from C-wrapper scripts.
 * In order to prevent name-mangling, we have to declare them as extern "C"
 */
extern "C" {
#endif
#endif

/* Opaque pointer typedefs for netCDF data structs */
typedef struct cmc_nc_data* cmc_nc_data_t;

#ifdef __cplusplus
[[noreturn]]
#else
_Noreturn
#endif
void
cmc_netcdf_exit(const int _err_code, const char* _location);

#ifdef CMC_WITH_NETCDF
#define cmc_nc_check_err(err) ((err) == NC_NOERR ? (void) 0 : cmc_netcdf_exit(err, CMC_FILE_LOCATION))
#endif

cmc_nc_data_t
cmc_nc_create(const int _ncid);

void
cmc_nc_set_mpi_communicator(cmc_nc_data_t nc_data, MPI_Comm comm);

void
cmc_nc_destroy(const cmc_nc_data_t _nc_data);

#ifdef __cplusplus
int
cmc_nc_open(const char* path_to_file, MPI_Comm comm = MPI_COMM_WORLD);
#else
int
cmc_nc_open(const char* path_to_file, MPI_Comm comm);
#endif

void
cmc_nc_close(int ncid);

void
cmc_inquire_coordinates(cmc_nc_data_t nc_data);

#ifdef __cplusplus
void
cmc_nc_inquire_var_data(cmc_nc_data_t nc_data, const size_t* start_ptr = nullptr, const size_t* count_ptr = nullptr);
#else
void
cmc_nc_inquire_var_data(cmc_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr);
#endif

void
cmc_nc_add_variable(cmc_nc_data_t nc_data, const char* var_name);

void
cmcc_nc_inquire_vars(cmc_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr, const int var_count, ...);

void
cmcc_nc_inquire_added_vars(cmc_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr);

/* If non-standard names for coordinate dimensions/variables are present, the netCDF dimension id for the corresponding dimension can be set via these functions */
void
cmc_nc_set_hint_latitude_dim(cmc_nc_data_t nc_data, const int lat_dim_id);

void
cmc_nc_set_hint_longitude_dim(cmc_nc_data_t nc_data, const int lon_dim_id);

void
cmc_nc_set_hint_elevation_dim(cmc_nc_data_t nc_data, const int lev_dim_id);

void
cmc_nc_set_hint_time_dim(cmc_nc_data_t nc_data, const int time_dim_id);

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
}
#endif
#endif

/* In the following, only internally used functions are declared (currently, not considered for general usage) */
/** \note They reside here in order to retain an opaque pointer to the netCDF data class and to strictly sepearte C++ functions from a C compliant interface */
#ifdef __cplusplus
void cmc_nc_preallocate_vars_by_name(cmc_nc_data_t, std::string&&);
#ifdef CMC_WITH_T8CODE
#include "t8code/cmc_t8code_data.hxx"
void _cmc_transform_nc_data_to_t8code_data(cmc_nc_data_t, cmc_t8_data_t);
#endif
#endif

#endif /* CMC_NETCDF_H */