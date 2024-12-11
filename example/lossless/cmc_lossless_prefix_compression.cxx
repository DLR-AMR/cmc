#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "netcdf/cmc_netcdf.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "lossless/cmc_prefix_compression.hxx"
#include "decompression/cmc_prefix_decompression.hxx"
#include "utilities/cmc_binary_reader.hxx"

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {
    
    #if 0
    const std::string file = "../../data/era5_reanalysis_t2m_tc03_13_12_23.nc";
    
    /* Create an object which interatcs with the file and opens it */
    cmc::NcData nc_data(file, cmc::NcOpeningMode::Serial);

    nc_data.SetHintHeightDimension(2); //Set time as height dimension

    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 512),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, 512),
                             cmc::DimensionInterval(cmc::Dimension::Lev, 0, 1)
                             );

    /* Inquire the hyperslab of data for the given variables */
    nc_data.InquireVariables(hyperslab, "t2m");

    /* Close the file, since we have gathered the data we wanted */
    nc_data.CloseFileHandle();

    #else

    #if 0
    const std::string file = "../../data/SDRBENCH-EXASKY-NYX-512x512x512/velocity_x.f32";

    cmc::CmcType type(cmc::CmcType::Float);
    const std::string name("velocity_x");
    const int id = 0;
    const size_t num_elements = 512 * 512 * 512;
    cmc::CmcUniversalType missing_value(static_cast<float>(31866790.000000));
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 512),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 512),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 256)
                          );

    cmc::bin_reader::Reader binary_reader(file);

    cmc::InputVar variable = binary_reader.CreateVariableFromBinaryData(type, name, id, num_elements, missing_value, layout, domain);

    std::vector<cmc::InputVar> vars;
    vars.emplace_back(std::move(variable));

    #else
    //const std::string file = "../../data/SDRBENCH-EXASKY-NYX-512x512x512/temperature.f32";
    const std::string file = "../../data/100x500x500/Pf48.bin.f32";
    cmc::CmcType type(cmc::CmcType::Float);
    const std::string name("qcloudf_log");
    const int id = 0;
    const size_t num_elements = 500 * 500 * 100;
    cmc::CmcUniversalType missing_value(static_cast<float>(13.5));
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 100)
                          );

    cmc::bin_reader::Reader binary_reader(file);

    cmc::InputVar variable = binary_reader.CreateVariableFromBinaryData(type, name, id, num_elements, missing_value, layout, domain);

    std::vector<cmc::InputVar> vars;
    vars.emplace_back(std::move(variable));

    #endif

    #endif


    {
    
    /* Create the compression data */
    //cmc::prefix::Compressor compression_data(nc_data.TransferData()); //netCDF data
    cmc::prefix::Compressor compression_data(std::move(vars)); //Binary data

    //compression_data.SetSplitVariable(cmc::SplitVariable(cmc::kSplitAllVariables, cmc::Dimension::Lev));

    /* Setup the example data for the compression */
    compression_data.Setup();

    cmc::cmc_debug_msg("\n\nSetup is finished\n\n");

    compression_data.Compress();

    cmc::cmc_debug_msg("\n\nCompression is finished\n\n");

    cmc::cmc_debug_msg("Write Compressed Data");
    compression_data.WriteCompressedData("CLOUDf48.bin.f32.new_scheme.cmc", 0);

    }

    {
    #if 0
    cmc::cmc_debug_msg("\n\nDecompression Start\n\n");

    /* Decompress */
    cmc::prefix::Decompressor decoder("nyx_velocity_x");

    decoder.Setup();
    decoder.DecompressVariable("velocity_x");
    
    #endif
    }

    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
