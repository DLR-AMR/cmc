#include "cmcc_wrapper.h"

int
cmcc_nc_open(const char* path_to_file_f, const int length, MPI_Fint comm_f)
{
    #ifdef CMC_WITH_NETCDF
    char path_to_file_c[length + 1];
    memcpy(path_to_file_c, path_to_file_f, length * sizeof(char));
    path_to_file_c[length] = '\0';
    return cmc_nc_open(path_to_file_c, MPI_Comm_f2c(comm_f));
    #endif
}

void
cmcc_nc_add_variable(cmc_nc_data_t nc_data, const char* var_name_f, const int length)
{
    #ifdef CMC_WITH_NETCDF
    char var_name_c[length + 1];
    memcpy(var_name_c, var_name_f, length * sizeof(char));
    var_name_c[length] = '\0';
    cmc_nc_add_variable(nc_data, var_name_c);
    #endif
}

cmc_amr_data_t
cmcc_create_amr_compression_data(cmc_nc_data_t nc_data, const MPI_Fint comm_f)
{
    #ifdef CMC_WITH_T8CODE
    #ifdef CMC_WITH_NETCDF
    cmc_create_amr_compression_data(nc_data, MPI_Comm_f2c(comm_f));
    #endif
    #endif
}

void
cmcc_amr_write_vtk_file(cmc_amr_data_t amr_data, const char* file_prefix, const int file_prefix_length)
{
    #ifdef CMC_WITH_T8CODE
    /* Allocate memory fot the file-prefix */
    char* fprefix = (char*) malloc(sizeof(char) * (file_prefix_length + 1));
    /* Copy the 'Fortran'-name */
    memcpy(fprefix, file_prefix, sizeof(char) * file_prefix_length);
    /* Add a trailing null-terminating character */
    fprefix[file_prefix_length] = '\0';

    /* Now we can call the C function with the C-compliant 'file_prefix' */
    cmc_amr_write_vtk_file(amr_data, fprefix);

    /* Deallocate the fpreix */
    free(fprefix);
    #endif
}

//TODO: update messy functions
#if 0
cmc_amr_data_t
cmcc_create_amr_compression_data_messy(cmc_messy_data_t messy_data, const MPI_Fint comm_f)
{
    #ifdef CMC_WITH_T8CODE
    #ifdef CMC_ENABLE_FORTRAN
    cmc_create_amr_compression_data_messy(messy_data, MPI_Comm_f2c(comm_f));
    #endif
    #endif
}
#endif
