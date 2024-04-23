#include "cmc.h"
#include "mpi/cmc_mpi.h"
#include "utilities/cmc_log_functions.h"
//#include "netcdf/cmc_netcdf.h"
#include "netcdf/cmc_netcdf.hxx"
//#include "t8code/cmc_t8code.h"
//#include "t8code/cmc_t8code_data.hxx"
//#include "t8code/cmc_t8code_geo_mesh.h"
//#include "t8code/cmc_t8code_geo_data.h"
//#include "component_interfaces/cmc_t8_nc.h"
//#include "lossy/cmc_amr_compressor.h"

#include "utilities/cmc_utilities.hxx"
#include "t8code/cmc_t8_data_variables.hxx"
#include "utilities/cmc_span.hxx"
#include "t8code/cmc_t8_interpolation.hxx"

int
main(int argc, char* argv[])
{
  /* Initialize cmc */
  cmc_initialize();

  {

    #if 0
    cmc::Variable<int> var;
    var.push_back(3);
    std::cout << "var Value at zero is: " << var[0] << std::endl;
    cmc::VariableAttributes<int>& attr = var.GetVariableAttributes();
    attr.SetMissingValue(1);

    cmc::Variable<double> var2;
    std::cout << "Missing value is of var is: " << var.GetVariableAttributes().GetMissingValue() << std::endl;

    //std::vector<cmc::CmcVariable> varvec;
    //varvec.push_back(cmc::MakeUniversalVariable(std::move(var)));
    //std::cout << "Missing value is: " << attr.GetMissingValue() << std::endl;
    //varvec.push_back(cmc::MakeUniversalVariable(std::move(var2)));

    cmc::CmcVariable var4(cmc::Variable<int>());

    cmc::Var var6(cmc::CmcType::Int32_t);
    cmc::Var var7{cmc::Variable<int>()};

    int z = 10;

    //var6.push_back(z);

    //std::vector<cmc::Var> vec;
    //vec.push_back(var6);
    //vec.push_back(std::move(var7));
//
    std::vector<double> testvec;
    testvec.push_back(1.0);
    testvec.push_back(2.0);
    testvec.push_back(3.0);
    testvec.push_back(4.0);
    testvec.push_back(5.0);
    testvec.push_back(6.0);

    //cmc::span<double> tspan(testvec.begin(), std::size_t(4));
    cmc::span<double, 4> sp;

    constexpr size_t ff = 4;
    cmc::span<double,ff> sp2(testvec.data(), ff);

    for (auto iter1 = sp2.begin(); iter1 != sp2.end(); ++iter1)
    {
      cmc_debug_msg("Hat val: ", *iter1);
    }

    cmc::Interpolate<double,4> interpol(cmc::InterpoalteToArithmeticMean);

    double res = interpol(sp2, -10.0);

    cmc_debug_msg("Interpoaltion results in: ", res);

    cmc::VectorView vv(testvec.data() + 1, 4);

    for (auto iter1 = vv.begin(); iter1 != vv.end(); ++iter1)
    {
      cmc_debug_msg("Hat val: ", *iter1);
    }

    //var.SetInterpolation(cmc::InterpolateTest);

    //std::vector<double> ia = var6.InterTest();

    //cmc_debug_msg("size of ia: ", ia.size());

    var2.SetInaccuracyStorage(cmc::TrackingOption::TrackFullInaccuracy);


    
    //var2.SetInaccuracyStorage(cmc::TrackingOption::TrackMinimalWorkingInaccuracy);

    //var7.Interpolate();

    //cmc::span sp3(testvec.data(), ff);

    //cmc::Variable<int32_t>& var6c = var6.GetVar<int32_t>();

    //std::cout << "var6c Value at zero is: " <<  var6c[0] << std::endl;

    //cmc::CmcUniversalType uv = var6[0];


    //std::cout << "var7 Value at zero is: " <<  std::get<int32_t>(uv) << std::endl;

  #endif
  
    #if 0
    //Inquire netCDF data, will be updated later
    cmc_nc_data_t nc_data = cmc_nc_start("../../data/MESSy_DATA/MESSy2/raw/tracer/tracer_in_par.nc", cmc_nc_opening_mode::CMC_NC_SERIAL);
    const size_t start_ptr[3] = {0,0,0};  //MESSy Tracer Initialization File
    const size_t count_ptr[3] = {90,64,128}; //MESSy Tracer Initialization File
    cmc_nc_inquire_vars(nc_data, start_ptr, count_ptr, "O3");

    //Gain data from a netCDF file
    cmc::NcData nc_data();

    cmc::CompressionSettings settings;

    settings.SetAbsoluteErrorCriterion(1.0);

    cmc::CompressionData compression_data(nc_data);

    comrpession_data.Setup();

    compression_data.Compress();

    #endif





  #if 0


  #if 1
    //New test for parallel netCDF input
    //cmc_nc_data_t nc_data = cmc_nc_start("../../data/test_nc4_file.nc", cmc_nc_opening_mode::CMC_NC_PARALLEL, MPI_COMM_WORLD);
    //const size_t start_ptr[3] = {0,0,0};  //Example netCDF File
    //const size_t count_ptr[3] = {1,73,144}; //Example netCDF File
    //cmc_nc_data_t nc_data = cmc_nc_start("../../data/tas_decreg_europe_v20140120_20010101_20010131.nc", cmc_nc_opening_mode::CMC_NC_SERIAL);
    //const size_t start_ptr[3] = {0,0,0};  //Example netCDF File
    //const size_t count_ptr[3] = {1,1026,1056}; //Example netCDF File
    //
    //std::vector<int> p_dist{1,2,1};
    //cmc_nc_set_blocked_reading(nc_data, p_dist);
    //
    //cmc_nc_inquire_vars(nc_data, start_ptr, count_ptr, "tas");
    //cmc_nc_inquire_vars(nc_data, start_ptr, count_ptr, "p2t");
    
    //cmc_nc_data_t nc_data = cmc_nc_start("../../data/MESSy_DATA/MESSy2/raw/tracer/RC1-base-07_0028_restart_0001_tracer_gp.nc", cmc_nc_opening_mode::CMC_NC_PARALLEL);
    
    cmc_nc_data_t nc_data = cmc_nc_start("../../data/MESSy_DATA/MESSy2/raw/tracer/tracer_in_par.nc", cmc_nc_opening_mode::CMC_NC_SERIAL);
    const size_t start_ptr[3] = {0,0,0};  //MESSy Tracer Initialization File
    const size_t count_ptr[3] = {90,64,128}; //MESSy Tracer Initialization File
    //std::vector<int> p_dist{2,1,1};
    //cmc_nc_set_blocked_reading(nc_data, p_dist);
    cmc_nc_inquire_vars(nc_data, start_ptr, count_ptr, "O3");

    /* Define data classes holding the variable data and forests as well as additional information during the compression process */
    cmc_amr_data_t amr_data;

    /* Create/Allocate a new 'AMR_DATA' class for the compression of netCDF inquired data */
    amr_data = cmc_create_amr_compression_data(nc_data, MPI_COMM_WORLD);

    /* Close the netCDF file and deallocate nc_data */
    cmc_nc_finish(nc_data);

    cmc_amr_pre_setup_split_3D_variable(amr_data, 0, DATA_LAYOUT::CMC_2D_LAT_LON);

    /* Set a compression criterium - e.g. error threshold with a predefined tolerance */
    cmc_amr_pre_setup_set_compression_criterium_relative_error_threshold(amr_data, 0.25);
    //cmc_amr_pre_setup_set_compression_criterium_absolute_error_threshold(amr_data, 0.00000005);

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
    cmc_amr_pre_setup_set_compression_criterium_relative_error_threshold(amr_data, 0.09);
    
    /* Write out a netCDF File containing the uncompressed data */
    //cmc_amr_write_netcdf_file(amr_data, "example_initial.nc", CMC_AMR_WRITE_ALL_VARS_TO_NETCDF);
    
    /* Setup the compression for a given 'compression mode' */
    cmc_amr_setup_compression(amr_data, CMC_T8_COMPRESSION_MODE::ONE_FOR_ALL);

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

    
#if 0
/* Open the netCDF file either for serial or parallel reading */
cmc_nc_data_t nc_data = cmc_nc_start ("path/to/file.nc", cmc_nc_opening_mode::CMC_NC_PARALLEL);

//Data Ordering: Level, Latitude, Longitude
const size_t start_ptr[3] = {0,0,0};  
const size_t count_ptr[3] = {90,64,128}; 

/* Inquire the coordinates, attributes, and data of the given variables */
cmc_nc_inquire_vars (nc_data, start_ptr, count_ptr, "O3", "CO2", "CH4");

//Do something with nc_data...

/* Close the netCDF file and deallocate nc_data */
cmc_nc_finish (nc_data);
#endif



  #endif





  }

  /* Finalize cmc */
  cmc_finalize();

  return 0;
}
