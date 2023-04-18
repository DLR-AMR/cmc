#ifndef CMCC_WRAPPER_H
#define CMCC_WRAPPER_H
/**
 * @file cmcc_wrapper.h
 * @brief This files provides wrapper functions for specific functions which may not be suitable for a direct Fortran interface
 */
#include "fortran_interfaces/cmcc.h"
#include "netcdf/cmc_netcdf.h"
#include "lossy/cmc_amr_compressor.h"

#ifdef CMC_ENABLE_MPI
#include "mpi.h"
#endif

/**
 * @brief This functions opens a netCDF-File
 * 
 * @param path_to_file_f The path to the netCDF-File which should be opened
 * @param length The length of the pathname
 * @param comm_f The Fortran MPI Communicator (if the file will be opened in parallel access)
 * @return int A unique id corresponding to this file 
 */
int
cmcc_nc_open(const char* path_to_file_f, const int length, MPI_Fint comm_f);

/**
 * @brief This functions adds a variable to the list of variables which will be read from the netCDF-File
 * 
 * @param nc_data The pointer to the netCDF data structure
 * @param var_name_f The name of the varibale
 * @param length The length of the variable's name
 */
void
cmcc_nc_add_variable(cmc_nc_data_t nc_data, const char* var_name_f, const int length);

/**
 * @brief This functions creates a @struct cmc_amr_data based on the netCDF input data supplied via @var nc_data and returns a pointer of the struct
 * 
 * @param nc_data The input data which was read from a netCDF-File
 * @param comm_f The Fotran MPI Communicator
 * @return cmc_amr_data_t 
 */
cmc_amr_data_t
cmcc_create_amr_compression_data(cmc_nc_data_t nc_data, const MPI_Fint comm_f);

/**
 * @brief This functions writes a vtk-File of the of the data residing in the @var amr_data 
 * 
 * @param amr_data The 'amr compression data'
 * @param file_prefix The desired output filename (without suffix/file extension)
 * @param file_prefix_length The length of the file's prefix name
 *
 * @note This function can be called before as well as after the compression of the data.
 *       In a 'One-For-One' compression has been used, each variable will be written into it's own file
 */
void
cmcc_amr_write_vtk_file(cmc_amr_data_t amr_data, const char* file_prefix, const int file_prefix_length);

/**
 * @brief Initialize cmc and it's submodules given a specific MPI Communicator
 * 
 * @param comm_f The MPI Communicator to use for cmc and it's submodules
 */
void
cmcc_initialize_mpi_comm(const MPI_Fint comm_f);

//TODO: update messy functions
//cmc_amr_data_t
//cmcc_create_amr_compression_data_messy(cmc_messy_data_t messy_data, const MPI_Fint comm_f);

#endif /* CMCC_WRAPPER_H */
