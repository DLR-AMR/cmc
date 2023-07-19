#ifndef CMC_NETCDF_HXX
#define CMC_NETCDF_HXX
/**
 * @file cmc_netcdf.hxx
 * @brief Via the 'include' of @file cmc_netcdf.h, this file collects all functions used for accessing netCDF files and storing the data of the netCDF variables as well as the geo-spatial domain on which the variables are defined.
 * Additionally, this file supplies C++ only functions for inquiring variables from a netCDF file
 */

#include "cmc_netcdf.h"
#include "utilities/cmc_geo_util.h"

void _cmc_nc_push_back_var(cmc_nc_data_t nc_data, std::string&& var_name);
void _cmc_nc_reserve_vars(cmc_nc_data_t nc_data, const size_t num_variables);

/**
 * @brief This function is internally used and allocates variables based on the name of the @var var_name
 * 
 * @param nc_data A pointer to a @struct cmc_nc_data (the data is read from the netCDF file corresponding to the id with which the @struct cmc_nc_data was created)
 * @param var_name The name of a variable
 * @tparam var_names The parameter pack containing the variable names
 */
template<typename... Ts>
void
cmc_nc_preallocate_vars_by_name(cmc_nc_data_t nc_data, std::string&& var_name, Ts&&... var_names)
{
    #ifdef CMC_WITH_NETCDF
    /* Pre-allocate variables by constructing new meta classes and save the given names of these variables */
    _cmc_nc_push_back_var(nc_data, std::move(var_name));
    cmc_nc_preallocate_vars_by_name(nc_data, std::forward<Ts>(var_names)...);
    #endif
};


/* Inquire information about the coordinate dimension and variables, as well as information about the actual supplied netCDF variables and store their data slice which is defined by the given hyperslab */
/**
 * @brief This function inquires information about the coordinate dimension and variables, as well as information about the actual supplied netCDF variables and store their data slice which is defined by the given hyperslab
 * 
 * @param nc_data A pointer to a @struct cmc_nc_data (the data is read from the netCDF file corresponding to the id with which the @struct cmc_nc_data was created)
 * @param start_ptr A pointer to an array containing the start index of each dimension (The length of the array has to coincide with the number of the dimensions of the variables)
 * @param count_ptr A pointer to an array containing the length for each dimension (The length of the array has to coincide with the number of the dimensions of the variables) 
 * @tparam var_names The parameter pack containing all the vairbale's names

 * @note There is a C-equivalent function @see @fn cmcc_nc_inquire_vars (in @file cmc_netcdf.h)
 */
template<typename... Ts>
void
cmc_nc_inquire_vars(cmc_nc_data_t nc_data, const size_t* start_ptr, const size_t* count_ptr, Ts&&... var_names)
{
    #ifdef CMC_WITH_NETCDF
    /* Inquire information about dimensions and read coordinate variables */
    cmc_inquire_coordinates(nc_data);

    /* Reserve memory for the variables */
    _cmc_nc_reserve_vars(nc_data, (sizeof...(Ts)));

    /* Allocate the variables based on the given data */
    cmc_nc_preallocate_vars_by_name(nc_data, std::forward<Ts>(var_names)...);

    /* Inquire the data of these variables */
    cmc_nc_inquire_var_data(nc_data, start_ptr, count_ptr);

    #else
    cmc_err_msg("CMC is not compiled with netCDF, please reconfigure with netCDF linkage in order to use this function.");
    #endif
};

/* Set a blocked reading if the file is processed in parallel */
void
cmc_nc_set_blocked_reading(cmc_nc_data_t nc_data, const std::vector<int> blocked_domain_num_processes_per_dimension);


#endif /* CMC_NETCDF_HXX */
