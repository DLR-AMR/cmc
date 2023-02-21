#ifndef CMCC_WRAPPER_H
#define CMCC_WRAPPER_H

#include "fortran_interfaces/cmcc.h"
#include "netcdf/cmc_netcdf.h"
#include "lossy/cmc_amr_compressor.h"

#ifdef CMC_ENABLE_MPI
#include "mpi.h"
#endif

int
cmcc_nc_open(const char* path_to_file_f, const int length, MPI_Fint comm_f);

void
cmcc_nc_add_variable(cmc_nc_data_t nc_data, const char* var_name_f, const int length);

cmc_amr_data_t
cmcc_create_amr_compression_data(cmc_nc_data_t nc_data, const MPI_Fint comm_f);

void
cmcc_amr_write_vtk_file(cmc_amr_data_t amr_data, const char* file_prefix, const int file_prefix_length);

//TODO: update messy functions
//cmc_amr_data_t
//cmcc_create_amr_compression_data_messy(cmc_messy_data_t messy_data, const MPI_Fint comm_f);

#endif /* CMCC_WRAPPER_H */
