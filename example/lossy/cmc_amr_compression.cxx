#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "netcdf/cmc_netcdf.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "lossy/cmc_amr_compression.hxx"
#include "decompression/cmc_prefix_decompression.hxx"

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {
    
    //const std::string file = "../../data/era5_reanalysis_t2m_tc03_13_12_23.nc";
    const std::string file = "../../data/mptrac_era5_2021_07_01_00.nc";

    /* Create an object which interatcs with the file and opens it */
    cmc::NcData nc_data(file, cmc::NcOpeningMode::Serial);

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

    /* Create compression settings */
    cmc::CompressionSettings settings;

    const double abs_max_err = 0.03125;
    settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    //const double rel_max_err = 0.0005;
    //settings.SetRelativeErrorCriterion(rel_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    cmc::SplitVariable split(cmc::kSplitAllVariables, cmc::Dimension::Lev);
    settings.SplitVariableByDimension(split);


    {
    /* Create the compression data */
    cmc::amr::Compressor compression_data(nc_data.TransferData(), std::move(settings));

    /* Setup the example data for the compression */
    compression_data.Setup();

    cmc::cmc_debug_msg("\n\nSetup is finished\n\n");

    compression_data.Compress();

    cmc::cmc_debug_msg("\n\nCompression is finished\n\n");
    cmc::cmc_debug_msg("Write Compressed Data");
    compression_data.WriteCompressedData("newly_example_compr_pref_t2m", 0);

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
