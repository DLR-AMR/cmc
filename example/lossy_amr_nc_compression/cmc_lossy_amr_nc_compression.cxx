#include "cmc.h"
#include "mpi/cmc_mpi.h"
#include "utilities/cmc_log_functions.h"
#include "netcdf/cmc_netcdf.h"
#include "netcdf/cmc_netcdf.hxx"
#include "t8code/cmc_t8code.h"
#include "t8code/cmc_t8code_data.hxx"
#include "t8code/cmc_t8code_geo_mesh.h"
#include "t8code/cmc_t8code_geo_data.h"
#include "t8_forest/t8_forest_vtk.h"
#include "component_interfaces/cmc_t8_nc.h"
#include "lossy/cmc_amr_compressor.h"


int
main(int argc, char* argv[])
{
  /* Initialize cmc */
  cmc_initialize();

  {
    /* Open a netCDF file with geo-spatial data */
    //For opening out of the 'cmc_v03' directory
    //int ncid = cmc_nc_open("../data/ECMWF_ERA-40_subset.nc");
    /* Open a MESSy simulation file */
    int ncid = cmc_nc_open("../data/MESSy_DATA/MESSy2/raw/tracer/RC1-base-07_0028_restart_0001_tracer_gp.nc");

    /* Create a class holding the data from the netCDF-File */
    //CMC_NC_DATA_T nc_data{ncid};
    cmc_nc_data_t nc_data{cmc_nc_create(ncid)};

    /* Define a hyperslab of the data */
    /*This hyperslab will be read in and used for the compression */
    //const size_t start_ptr[3] = {0,0,0};  //Example netCDF File
    //const size_t count_ptr[3] = {1,73,144}; //Example netCDF File
    //const size_t start_ptr[4] = {0,0,0,0}; //Example MESSy netCDF File 
    //const size_t count_ptr[4] = {1,19,32,64}; // Example MESSy netCDF File
    
    const size_t start_ptr[3] = {0,0,0};  //MESSy Tracer Initialization File
    const size_t count_ptr[3] = {5,64,128}; //MESSy Tracer Initialization File

    /* Inquire given data/variables which are defined on a geo-spatial domain (latitude, longitude, height, (time)) */
    //cmc_nc_inquire_vars(nc_data, start_ptr, count_ptr, "tco3", "p2t");
    cmc_nc_inquire_vars(nc_data, start_ptr, count_ptr, "O3", "CH4");

    /* Define data classes holding the variable data and forests as well as additional information during the compression process */
    cmc_amr_data_t amr_data;

    /* Create/Allocate a new 'AMR_DATA' class for the compression of netCDF inquired data */
    amr_data = cmc_create_amr_compression_data(nc_data, MPI_COMM_WORLD);

    /* Split the 3D variable */
    //cmc_amr_pre_setup_split_3D_variable(amr_data, 0, DATA_LAYOUT::CMC_2D_LAT_LON);

    /* Set a compression criterium - e.g. error threshold woth a predefined tolerance */
    cmc_amr_pre_setup_set_compression_criterium_error_threshold(amr_data, 0.02);
    
    /* Write out a netCDF File containing the uncompressed data */
    //cmc_amr_write_netcdf_file(amr_data, "example_initial.nc");
    
    /* Setup the compression for a given 'compression mode' */
    cmc_amr_setup_compression(amr_data, CMC_T8_COMPRESSION_MODE::ONE_FOR_ONE_2D);

    /* Write out a netCDF File containing the uncompressed data */
    //cmc_amr_write_netcdf_file(amr_data, "example_uncompressed.nc");
    
    /* Execute the adaptation/compression */
    cmc_amr_compress(amr_data);

    //cmc_amr_write_netcdf_file(amr_data, "example_compressed.nc");

    /* Decompress the data */
    cmc_amr_decompress(amr_data);

  
  	/* Write a vtk file of the decompressed data */
    cmc_amr_write_vtk_file(amr_data, "cmc_decompressed_data");

    /* Write a netCDF file */
    /* Write out a netCDF File containing the uncompressed data */
    cmc_amr_write_netcdf_file(amr_data, "example_decompressed.nc", CMC_AMR_WRITE_ALL_VARS_TO_NETCDF);

    /* Close the netCDF file */
    cmc_nc_close(ncid);

    /* Deallocate the netCDF data */
    cmc_nc_destroy(nc_data);

    /* Deallocate the Lossy AMR Compression data */
    cmc_amr_destroy(amr_data);
  }

  /* Finalize cmc */
  cmc_finalize();

  return 0;
}
