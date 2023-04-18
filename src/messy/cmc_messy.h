#ifndef CMC_MESSY_H
#define CMC_MESSY_H

//#include "fortran_interfaces/cmcc.h"
#include "cmc.h"
#include "utilities/cmc_constants_definitions.h"

#ifdef CMC_ENABLE_MPI
#include "mpi.h"
#endif

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
/* In case the Fortran interface is built, these function needs to be callable from C-wrapper scripts.
 * In order to prevent name-mangling, we have to declare them as extern "C"
 */
extern "C" {
#endif
#endif


/* Opaque pointer to data structure */
typedef struct cmc_messy_data* cmc_messy_data_t;

cmc_messy_data_t
cmc_setup_messy_data(MPI_Fint comm);

void
cmc_destroy_messy_data(cmc_messy_data_t messy_data);



#if 0
//TODO: update


/* Opaque pointer to data struct */
typedef struct CMC_MESSY_DATA* CMC_MESSY_DATA_T;

CMC_MESSY_DATA_T cmc_setup_messy(int* _dimension_sizes, const char* axis_representation, const int _num_variables, const double _missing_value, MPI_Fint comm_f);

void
cmc_messy_set_tracer_variable(CMC_MESSY_DATA_T messy_data, const char*  tracer_name, const int tracer_name_length, const double* data_array);

void
cmc_messy_set_reference_tracer(CMC_MESSY_DATA_T messy_data, const char*  tracer_name, const int tracer_name_length, const double* data_array);

void
cmc_messy_free_data(CMC_MESSY_DATA_T messy_data);

#endif

#ifdef __cplusplus
#ifdef CMC_ENABLE_FORTRAN
}
#endif
#endif


/* In the following, only internally used functions are declared (currently, not considered for general usage) */
#ifdef __cplusplus
#ifdef CMC_WITH_T8CODE
#include "t8code/cmc_t8code_data.h"
//TODO: update
//void _cmc_transform_messy_data_to_t8code_data(CMC_MESSY_DATA_T messy_data, CMC_T8_DATA* t8_data);
#endif
#endif

#endif /* CMC_MESSY_H */