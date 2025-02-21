#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "input/cmc_netcdf.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "lossy/cmc_amr_compression.hxx"
#include "decompression/cmc_prefix_decompression.hxx"
#include "utilities/cmc_binary_reader.hxx"

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {
    
    #if 0
    //const std::string file = "../../data/era5_reanalysis_t2m_tc03_13_12_23.nc";
    const std::string file = "../../data/mptrac_era5_2021_07_01_00.nc";

    /* Create an object which interatcs with the file and opens it */
    cmc::input::netcdf::Data nc_data(file, cmc::nc::OpeningMode::Serial);

    //nc_data.SetHintHeightDimension(2); //Set time as height dimension
    nc_data.SetHintHeightDimension(6); //Set plev as height dimension

    //cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1440),
    //                         cmc::DimensionInterval(cmc::Dimension::Lat, 0, 721),
    //                         cmc::DimensionInterval(cmc::Dimension::Lev, 0, 1)
    //                         );

    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1200),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, 601),
                             cmc::DimensionInterval(cmc::Dimension::Lev, 0, 137)
                             );

    /* Inquire the hyperslab of data for the given variables */
    //nc_data.InquireVariables(hyperslab, "tco3");
    nc_data.InquireVariables(hyperslab, "t");

    /* Close the file, since we have gathered the data we wanted */
    nc_data.CloseFileHandle();

    #else

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

    #endif

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
    #if 0
    cmc::amr::Compressor compression_data(nc_data.TransferData(), std::move(settings));
    #else
    cmc::amr::Compressor compression_data(std::move(vars), std::move(settings));
    #endif
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
