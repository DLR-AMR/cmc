#include "cmc_t8_nc.h"

void
cmc_t8_nc_setup_compression(cmc_nc_data_t nc_data, cmc_t8_data_t t8_data){
    #if CMC_WITH_T8CODE && CMC_WITH_NETCDF
    _cmc_transform_nc_data_to_t8code_data(nc_data, t8_data);
    #endif
};
