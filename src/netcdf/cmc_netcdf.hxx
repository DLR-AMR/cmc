#ifndef CMC_NETCDF_HXX
#define CMC_NETCDF_HXX

#include "cmc_netcdf.h"
#include "utilities/cmc_geo_util.h"

void _cmc_nc_push_back_var(cmc_nc_data_t nc_data, std::string&& var_name);
void _cmc_nc_reserve_vars(cmc_nc_data_t nc_data, const size_t num_variables);

/** Begin TEMPLATED Functions **/

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

/** End TEMPLATED Functions **/


#endif /* CMC_NETCDF_HXX */
