#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "netcdf/cmc_netcdf.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "lossy/cmc_full_amr_compression.hxx"
#include "decompression/cmc_prefix_decompression.hxx"
#include "utilities/cmc_binary_reader.hxx"

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {

    const std::string file = "../../data/100x500x500/Pf48.bin.f32";
    cmc::CmcType type(cmc::CmcType::Float);
    const std::string name("precipf48");
    const int id = 0;
    const size_t num_elements = 500 * 500 * 100;
    cmc::CmcUniversalType missing_value(static_cast<float>(3224.4));
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

    const double abs_max_err = 3.0;
    settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    //const double rel_max_err = 0.0005;
    //settings.SetRelativeErrorCriterion(rel_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    //cmc::SplitVariable split(cmc::kSplitAllVariables, cmc::Dimension::Lev);
    //settings.SplitVariableByDimension(split);


    {
    /* Create the compression data */
    cmc::full_amr::Compressor compression_data(std::move(vars), std::move(settings));

    /* Setup the example data for the compression */
    compression_data.Setup();

    cmc::cmc_debug_msg("\n\nSetup is finished\n\n");

    compression_data.Compress();

    cmc::cmc_debug_msg("\n\nCompression is finished\n\n");
    cmc::cmc_debug_msg("Write Compressed Data");
    compression_data.WriteCompressedData("AMR_AC_PREF_Pf48_AbsErr_3_0.cmc", 0);

    }

    {
        #if 0
    cmc::cmc_debug_msg("\n\nDecompression Start\n\n");

    /* Decompress */
    cmc::prefix::Decompressor decoder("newly_example_compr_pref_t2m");

    decoder.Setup();
    decoder.DecompressVariable("t2m");
        #endif
    }


    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
