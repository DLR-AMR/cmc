#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "netcdf/cmc_netcdf.hxx"
#include "lossy/cmc_amr_lossy_compression_settings.hxx"
#include "lossless/cmc_prefix_compression.hxx"
#include "decompression/cmc_prefix_decompression.hxx"

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {
    
    const std::string file = "../../data/era5_reanalysis_t2m_tc03_13_12_23.nc";
    
    /* Create an object which interatcs with the file and opens it */
    cmc::NcData nc_data(file, cmc::NcOpeningMode::Serial);

    nc_data.SetHintHeightDimension(2); //Set plev as height dimension

    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1440),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, 721),
                             cmc::DimensionInterval(cmc::Dimension::Lev, 0, 1)
                             );

    /* Inquire the hyperslab of data for the given variables */
    nc_data.InquireVariables(hyperslab, "t2m");

    /* Close the file, since we have gathered the data we wanted */
    nc_data.CloseFileHandle();

    {
    /* Create the compression data */
    cmc::prefix::Compressor compression_data(nc_data.TransferData());

    /* Setup the example data for the compression */
    compression_data.Setup();

    cmc::cmc_debug_msg("\n\nSetup is finished\n\n");

    compression_data.Compress();

    cmc::cmc_debug_msg("\n\nCompression is finished\n\n");

    cmc::cmc_debug_msg("Write Compressed Data");
    compression_data.WriteCompressedData("newly_example_compr_pref_t2m", 0);

    }

    {
    cmc::cmc_debug_msg("\n\nDecompression Start\n\n");

    /* Decompress */
    cmc::prefix::Decompressor decoder("newly_example_compr_pref_t2m");

    decoder.Setup();
    decoder.DecompressVariable("t2m");
    
    }

    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
