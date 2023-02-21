#ifndef CMC_NETCDF_H
#define CMC_NETCDF_H
/**
 * @file cmc_netcdf.h
 * @brief This file collects all functions used for accessing netCDF files and storing the data of the netCDF variables as well as the geo-spatial domain on which the variables are defined
 */

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

/**
 * @brief This function allocates a @struct cmc_nc_data and returns a pointer to it
 * 
 * @param _ncid The id returned by a call to 'cmc_nc_open(...)'
 * @return A pointer to the allocated struct
 */
cmc_nc_data_t
cmc_nc_create(const int _ncid);

/**
 * @brief This function sets an MPI communicator to use for parallel file access via netCDF
 * 
 * @param nc_data The pointer to the @struct cmc_nc_data which should use the communicator @var comm
 * @param comm The MPI communicator to use in a parallel environment
 */
void
cmc_nc_set_mpi_communicator(cmc_nc_data_t nc_data, MPI_Comm comm);

/**
 * @brief This function destroys/deallocates a previously allocated @struct cmc_nc_data 
 * 
 * @param _nc_data A pointer to a @struct cmc_nc_data which should be deallocated
 */
void
cmc_nc_destroy(const cmc_nc_data_t _nc_data);

/**
 * @brief This function (@fn int cmc_nc_open(const char* path_to_file, MPI_Comm comm)) opens a netCDF file
 * 
 * @param path_to_file The path to the netCDF file which will be opnened
 * @param comm The MPI communicator to use in a parallel environment
 * @return An integer resembling the unique id of the opened netCDF file
 *
 * @note If compiled with a C++ compiler, the the default value for the @var comm is 'MPI_COMM_WORLD'
 */
 #ifdef __cplusplus
int
cmc_nc_open(const char* path_to_file, MPI_Comm comm = MPI_COMM_WORLD);
#else
int
cmc_nc_open(const char* path_to_file, MPI_Comm comm);
#endif


/**
 * @brief This function closes a (previously opened) netCDF file
 * 
 * @param ncid The id of the netCDF file which will be closed
 */
void
cmc_nc_close(int ncid);


/**
 * @brief This function inquires information about the geo-spatial domain on which the variables may be defined.
 * Therefore, the coordinate dimensions and coordinate variables are inquired by this function and stored within the @struct cmc_nc_data to which the @var nc_data points to.
 * 
 * @param nc_data A pointer to a @struct cmc_nc_data (the data is read from the netCDF file corresponding to the id with which the @struct cmc_nc_data was created)
 */
void
cmc_inquire_coordinates(cmc_nc_data_t nc_data);

#ifdef __cplusplus
/**
 * @brief This function inquires information about the geo-spatial domain on which the variables may be defined.
 * Therefore, the coordinate dimensions and coordinate variables are inquired by this function and stored within the @struct cmc_nc_data to which the @var nc_data points to.
 * 
 * @param nc_data A pointer to a @struct cmc_nc_data (the data is read from the netCDF file corresponding to the id with which the @struct cmc_nc_data was created)
 * @param start_ptr A pointer to an array containing the start index of each dimension (The length of the array has to coincide with the number of the dimensions of the variables)
 * @param count_ptr A pointer to an array containing the length for each dimension (The length of the array has to coincide with the number of the dimensions of the variables) 
 *
 * @note If nullptr's for @var start_ptr and @var count_ptr are supplied, the whole variables will be read from the file. Otherwise only the supplied hyperslab (defined by the @var start_ptr and @var count_ptr) will be read for each variable.
 * @note If compiled with a C++ compiler, default values for @var start_ptr and @var count_ptr are 'nullptr's
 */
void
cmc_nc_inquire_var_data(cmc_nc_data_t nc_data, const size_t* start_ptr = nullptr, const size_t* count_ptr = nullptr);
#else
void
cmc_nc_inquire_var_data(cmc_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr);
#endif

/**
 * @brief This function adds a the variable (with the supplied name) to the list of variables which will be read from the netCDF file
 * 
 * @param nc_data A pointer to a @struct cmc_nc_data (the data is read from the netCDF file corresponding to the id with which the @struct cmc_nc_data was created)
 * @param var_name The name of the variable which should be read 
 *
 * @note Variables have to be either added only by this function or with a single call to @fn cmcc_nc_inquire_vars or with a single call to @fn cmc_nc_inquire_vars (@see @file cmc_netcdf.hxx)
 * @note After (several times of) calling this function, the 'added vars' has to be inquired by a call to @fn cmcc_nc_inquire_added_vars)
 */
void
cmc_nc_add_variable(cmc_nc_data_t nc_data, const char* var_name);

/**
 * @brief This function inquires variables corressponding to the names supplied
 * 
 * @param nc_data A pointer to a @struct cmc_nc_data (the data is read from the netCDF file corresponding to the id with which the @struct cmc_nc_data was created)
 * @param start_ptr A pointer to an array containing the start index of each dimension (The length of the array has to coincide with the number of the dimensions of the variables)
 * @param count_ptr A pointer to an array containing the length for each dimension (The length of the array has to coincide with the number of the dimensions of the variables) 
 * @param var_count The number of variables which are supplied
 *
 * @note The hyberslab (defined by @var start_ptr and @var count_ptr) which will be read from each variable's data
 * @note This function uses 'va_args', all variable names should be appeneded as the last parameters in this function call (the amount of variables have to coincide with the number @var var_count)
 * @note This is a C-equivalent of the C++ function (@see @fn cmc_nc_inquire_vars defined in @file cmc_netcdf.hxx) using variadic templates
 */
void
cmcc_nc_inquire_vars(cmc_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr, const int var_count, ...);

/**
 * @brief This function inquires the data of all variables which were previously added by calls to @fn cmc_nc_add_variable
 * 
 * @param nc_data A pointer to a @struct cmc_nc_data (the data is read from the netCDF file corresponding to the id with which the @struct cmc_nc_data was created)
 * @param start_ptr A pointer to an array containing the start index of each dimension (The length of the array has to coincide with the number of the dimensions of the variables)
 * @param count_ptr A pointer to an array containing the length for each dimension (The length of the array has to coincide with the number of the dimensions of the variables)

 * @note Variables have to be either added only by this function or with a single call to @fn cmcc_nc_inquire_vars or with a single call to @fn cmc_nc_inquire_vars (@see @file cmc_netcdf.hxx)
 * @note After (several times of) calling this function, the 'added vars' has to be inquired by a call to @fn cmcc_nc_inquire_added_vars)
 */
void
cmcc_nc_inquire_added_vars(cmc_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr);

/**
 * @brief This function gives a hint which netCDF dimension id corresponds to the latitude dimension
 * 
 * @param nc_data A pointer to a @struct cmc_nc_data (the data is read from the netCDF file corresponding to the id with which the @struct cmc_nc_data was created)
 * @param lat_dim_id The internal netCDF dimension id which will be used as the latitude dimension
 *
 * @note If these functions are not used, the file is searched for standard names of coordinate dimension and coordinate variable names
 * Therefore, if non-standard names for coordinate dimensions/variables are present, the netCDF dimension id for the corresponding dimension can be set via these functions
 */
void
cmc_nc_set_hint_latitude_dim(cmc_nc_data_t nc_data, const int lat_dim_id);

/**
 * @brief This function gives a hint which netCDF dimension id corresponds to the longitude dimension
 * 
 * @param nc_data A pointer to a @struct cmc_nc_data (the data is read from the netCDF file corresponding to the id with which the @struct cmc_nc_data was created)
 * @param lat_dim_id The internal netCDF dimension id which will be used as the longitude dimension
 *
 * @note If these functions are not used, the file is searched for standard names of coordinate dimension and coordinate variable names
 * Therefore, if non-standard names for coordinate dimensions/variables are present, the netCDF dimension id for the corresponding dimension can be set via these functions
 */
void
cmc_nc_set_hint_longitude_dim(cmc_nc_data_t nc_data, const int lon_dim_id);

/**
 * @brief This function gives a hint which netCDF dimension id corresponds to the elevation dimension
 * 
 * @param nc_data A pointer to a @struct cmc_nc_data (the data is read from the netCDF file corresponding to the id with which the @struct cmc_nc_data was created)
 * @param lat_dim_id The internal netCDF dimension id which will be used as the elevation dimension
 *
 * @note If these functions are not used, the file is searched for standard names of coordinate dimension and coordinate variable names
 * Therefore, if non-standard names for coordinate dimensions/variables are present, the netCDF dimension id for the corresponding dimension can be set via these functions
 */
void
cmc_nc_set_hint_elevation_dim(cmc_nc_data_t nc_data, const int lev_dim_id);

/**
 * @brief This function gives a hint which netCDF dimension id corresponds to the 'time' dimension
 *
 * @note Currently, the 'time'-coordinate is not considered
 *
 * @param nc_data A pointer to a @struct cmc_nc_data (the data is read from the netCDF file corresponding to the id with which the @struct cmc_nc_data was created)
 * @param lat_dim_id The internal netCDF dimension id which will be used as the 'time' dimension
 *
 * @note If these functions are not used, the file is searched for standard names of coordinate dimension and coordinate variable names
 * Therefore, if non-standard names for coordinate dimensions/variables are present, the netCDF dimension id for the corresponding dimension can be set via these functions
 */
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