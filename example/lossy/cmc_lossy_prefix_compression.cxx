#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "netcdf/cmc_netcdf.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "lossy/cmc_prefix_compression.hxx"
#include "decompression/cmc_prefix_decompression.hxx"
#include "utilities/cmc_binary_reader.hxx"

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {
    //const std::string file = "../../data/SDRBENCH-EXASKY-NYX-512x512x512/temperature.f32";
    const std::string file = "../../data/100x500x500/Pf48.bin.f32";
    cmc::CmcType type(cmc::CmcType::Float);
    const std::string name("qcloudf_log");
    const int id = 0;
    const size_t num_elements = 500 * 500 * 100;
    cmc::CmcUniversalType missing_value(static_cast<float>(3224.5));
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 100)
                          );

    cmc::bin_reader::Reader binary_reader(file);

    cmc::InputVar variable = binary_reader.CreateVariableFromBinaryData(type, name, id, num_elements, missing_value, layout, domain);

    std::vector<cmc::InputVar> vars;
    vars.emplace_back(std::move(variable));

     /* Create compression settings */
    cmc::CompressionSettings settings;

    //const double rel_max_err = 0.01;
    //settings.SetAbsoluteErrorCriterion(rel_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    const double abs_max_err = 0.05;
    settings.SetRelativeErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);


    {
    
    /* Create the compression data */
    //cmc::prefix::Compressor compression_data(nc_data.TransferData()); //netCDF data
    cmc::prefix::lossy::Compressor compression_data(std::move(vars), std::move(settings)); //Binary data

    //compression_data.SetSplitVariable(cmc::SplitVariable(cmc::kSplitAllVariables, cmc::Dimension::Lev));

    /* Setup the example data for the compression */
    compression_data.Setup();

    cmc::cmc_debug_msg("\n\nSetup is finished\n\n");

    compression_data.Compress();

    cmc::cmc_debug_msg("\n\nCompression is finished\n\n");

    cmc::cmc_debug_msg("Write Compressed Data");
    compression_data.WriteCompressedData("Pf48.bin.f32.new_scheme2.cmc", 0);

    }

    {
    #if 0
    cmc::cmc_debug_msg("\n\nDecompression Start\n\n");

    /* Decompress */
    cmc::prefix::lossy::Decompressor decoder("nyx_velocity_x");

    decoder.Setup();
    decoder.DecompressVariable("velocity_x");
    
    #endif
    }

    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
