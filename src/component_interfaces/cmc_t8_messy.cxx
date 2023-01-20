#include "cmc_t8_messy.h"

#if 0
//TODO: update
void
cmc_t8_messy_setup_compression(CMC_MESSY_DATA_T messy_data,  CMC_T8_DATA* t8_data)
{
    #ifdef CMC_WITH_T8CODE
    #ifdef CMC_ENABLE_FORTRAN
    /* Call the intern setup compression function in 'cmc_messy.cxx' */
    _cmc_transform_messy_data_to_t8code_data(messy_data, t8_data);
    #endif
    #endif
}
#endif
