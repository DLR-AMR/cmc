#include "cmc.h"
#include "mpi/cmc_mpi.h"
#include "utilities/cmc_log_functions.h"
#include "netcdf/cmc_netcdf.h"
#include "netcdf/cmc_netcdf.hxx"
#include "t8code/cmc_t8code.h"
#include "t8code/cmc_t8code_data.hxx"
#include "t8code/cmc_t8code_geo_mesh.h"
#include "t8code/cmc_t8code_geo_data.h"
#include "component_interfaces/cmc_t8_nc.h"
#include "lossy/cmc_amr_compressor.h"


int
main(int argc, char* argv[])
{
  /* Initialize cmc */
  cmc_initialize();

  {
  
  #if 1
    //New test for parallel netCDF input
    cmc_nc_data_t nc_data = cmc_nc_start("../../data/test_nc4_file.nc", cmc_nc_opening_mode::CMC_NC_PARALLEL, MPI_COMM_WORLD);
    
    const size_t start_ptr[3] = {0,0,0};  //Example netCDF File
    const size_t count_ptr[3] = {1,73,144}; //Example netCDF File

    //std::vector<int> p_dist{1,2,2};
    //cmc_nc_set_blocked_reading(nc_data, p_dist);

    cmc_nc_inquire_vars(nc_data, start_ptr, count_ptr, "p2t");

    /* Define data classes holding the variable data and forests as well as additional information during the compression process */
    cmc_amr_data_t amr_data;

    /* Create/Allocate a new 'AMR_DATA' class for the compression of netCDF inquired data */
    amr_data = cmc_create_amr_compression_data(nc_data, MPI_COMM_WORLD);

    /* Close the netCDF file and deallocate nc_data */
    cmc_nc_finish(nc_data);

    /* Set a compression criterium - e.g. error threshold with a predefined tolerance */
    cmc_amr_pre_setup_set_compression_criterium_error_threshold(amr_data, 0.03);

    /* Keep the initial data in order to check the actual introduced data inaccurcy after the decompression */
    cmc_amr_pre_setup_set_flag_in_order_to_keep_the_initial_data(amr_data, 1);

    /* Setup the compression for a given 'compression mode' */
    cmc_amr_setup_compression(amr_data, CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL);

    /* Execute the adaptation/compression */
    cmc_amr_compress(amr_data);

    cmc_amr_write_vtk_file(amr_data, "example_compressed");

    cmc_amr_decompress(amr_data);
    
    cmc_amr_write_vtk_file(amr_data, "example_decompressed");

    /* Deallocate the Lossy AMR Compression data */
    cmc_amr_destroy(amr_data);

#else
    /* Open a netCDF file with geo-spatial data */
    //For opening out of the 'cmc_v03' directory
    //int ncid = cmc_nc_open("../../data/ECMWF_ERA-40_subset.nc");
    //int ncid = cmc_nc_open("../../data/test_nc4_file.nc");
    /* Open a MESSy simulation file */
    int ncid = cmc_nc_open("../../data/MESSy_DATA/MESSy2/raw/tracer/RC1-base-07_0028_restart_0001_tracer_gp.nc");

    /* Create a class holding the data from the netCDF-File */
    //CMC_NC_DATA_T nc_data{ncid};
    cmc_nc_data_t nc_data{cmc_nc_create(ncid)};

    /* Define a hyperslab of the data */
    /* This hyperslab will be read in and used for the compression */
    //const size_t start_ptr[3] = {0,0,0};  //Example netCDF File
    //const size_t count_ptr[3] = {1,73,144}; //Example netCDF File
    //const size_t start_ptr[4] = {0,0,0,0}; //Example MESSy netCDF File 
    //const size_t count_ptr[4] = {1,19,32,64}; // Example MESSy netCDF File
    
    const size_t start_ptr[3] = {0,0,0};  //MESSy Tracer Initialization File
    const size_t count_ptr[3] = {90,64,128}; //MESSy Tracer Initialization File
    //cmc_nc_set_mpi_communicator(nc_data, MPI_COMM_WORLD);
    /* Inquire given data/variables which are defined on a geo-spatial domain (latitude, longitude, height, (time)) */
    //cmc_nc_inquire_vars(nc_data, start_ptr, count_ptr, "p2t");
    cmc_nc_inquire_vars(nc_data, start_ptr, count_ptr, "O3");

    //int err = MPI_Barrier(MPI_COMM_WORLD);
    //cmc_mpi_check_err(err);
    /* Close the netCDF file */
    cmc_nc_close(ncid);
#endif

#if 0
    /* Define data classes holding the variable data and forests as well as additional information during the compression process */
    cmc_amr_data_t amr_data;

    /* Create/Allocate a new 'AMR_DATA' class for the compression of netCDF inquired data */
    amr_data = cmc_create_amr_compression_data(nc_data, MPI_COMM_WORLD);

    /* Split the 3D variable */
    //cmc_amr_pre_setup_split_3D_variable(amr_data, 0, DATA_LAYOUT::CMC_2D_LAT_LON);

    //cmc_universal_type_t lon_min = static_cast<float>(275);
    //cmc_universal_type_t lon_max = static_cast<float>(330);
    //cmc_universal_type_t lat_min = static_cast<float>(22);
    //cmc_universal_type_t lat_max = static_cast<float>(-58);
    //cmc_amr_pre_setup_set_compression_criterium_exclude_area(amr_data, CMC_COORD_IDS::CMC_LON, lon_min, lon_max);
    //cmc_amr_pre_setup_set_compression_criterium_exclude_area(amr_data, CMC_COORD_IDS::CMC_LAT, lat_min, lat_max);

    //cmc_universal_type_t lev_min = static_cast<double>(20);
    //cmc_universal_type_t lev_max = static_cast<double>(30);
    //cmc_amr_pre_setup_set_compression_criterium_exclude_area(amr_data, CMC_COORD_IDS::CMC_LEV, lev_min, lev_max);

    /* Set a compression criterium - e.g. error threshold with a predefined tolerance */
    cmc_amr_pre_setup_set_compression_criterium_error_threshold(amr_data, 0.09);
    
    /* Write out a netCDF File containing the uncompressed data */
    //cmc_amr_write_netcdf_file(amr_data, "example_initial.nc", CMC_AMR_WRITE_ALL_VARS_TO_NETCDF);
    
    /* Setup the compression for a given 'compression mode' */
    cmc_amr_setup_compression(amr_data, CMC_T8_COMPRESSION_MODE::ONE_FOR_ONE);

    /* Write a vtk file of the decompressed data */
    cmc_amr_write_vtk_file(amr_data, "cmc_initial_data");

    /* Write out a netCDF File containing the uncompressed data */
    //cmc_amr_write_netcdf_file(amr_data, "example_uncompressed.nc");
    
    /* Execute the adaptation/compression */
    //cmc_amr_compress(amr_data);

    //cmc_amr_write_netcdf_file(amr_data, "example_compressed.nc");

    /* Decompress the data */
    //cmc_amr_decompress(amr_data);

  
  	/* Write a vtk file of the decompressed data */
    //cmc_amr_write_vtk_file(amr_data, "cmc_decompressed_data");

    /* Write a netCDF file */
    /* Write out a netCDF File containing the uncompressed data */
    //cmc_amr_write_netcdf_file(amr_data, "example_decompressed.nc", CMC_AMR_WRITE_ALL_VARS_TO_NETCDF);

    /* Deallocate the netCDF data */
    cmc_nc_destroy(nc_data);

    /* Deallocate the Lossy AMR Compression data */
    cmc_amr_destroy(amr_data);

  #endif
  }

  /* Finalize cmc */
  cmc_finalize();

  return 0;
}
