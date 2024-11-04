#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "netcdf/cmc_netcdf.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "lossy/cmc_amr_compression.hxx"
#include "decompression/cmc_prefix_decompression.hxx"

#include <cfenv>
int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {
    
    #if 1

    #if 1
    const std::string file = "../../data/era5_reanalysis_t2m_tc03_13_12_23.nc";
    //const std::string file = "../../data/tas_decreg_europe_v20140120_20010101_20010131.nc";
    //const std::string file = "../../data/era5_reanalysis_data_16_11_23.nc";
    //const std::string file = "../../data/era5_reanalysis_pressure_lvls_fixed_time.nc";
    //const std::string file = "../../data/mptrac_era5_2021_07_01_00.nc";
    
    //const std::string file = "../../data/mptrac_era5_2021_07_01_00.nc";
    
    /* Create an object which interatcs with the file and opens it */
    cmc::NcData nc_data(file, cmc::NcOpeningMode::Serial);

    //cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1440),
    //                         cmc::DimensionInterval(cmc::Dimension::Lat, 0, 721)
    //                         );
    ///* Inquire the hyperslab of data for the given variables */
    //nc_data.InquireVariables(hyperslab, "t2m");

    //cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1440),
    //                         cmc::DimensionInterval(cmc::Dimension::Lat, 0, 721),
    //                         cmc::DimensionInterval(cmc::Dimension::Lev, 0, 37)
    //                         );

    nc_data.SetHintHeightDimension(2); //Set plev as height dimension

    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1440),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, 721),
                             cmc::DimensionInterval(cmc::Dimension::Lev, 0, 1)
                             );


    //Test domain
    //cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 8),
    //                         cmc::DimensionInterval(cmc::Dimension::Lat, 0, 8)
    //                         );
    //For testing rn
    //cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 16),
    //                         cmc::DimensionInterval(cmc::Dimension::Lat, 0, 16)
    //                         );

    /* Inquire the hyperslab of data for the given variables */
    //nc_data.InquireVariables(hyperslab, "tco3");
    nc_data.InquireVariables(hyperslab, "t2m");

    /* Close the file, since we have gathered the data we wanted */
    nc_data.CloseFileHandle();

    /* Create compression settings */
    cmc::CompressionSettings settings;

    const double abs_max_err = 0.0;
    settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    //const double rel_max_err = 0.0005;
    //settings.SetRelativeErrorCriterion(rel_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    //cmc::SplitVariable split(cmc::kSplitAllVariables, cmc::Dimension::Lev);
    //settings.SplitVariableByDimension(split);

    #endif

    {
    /* Create the compression data */
    cmc::amr::Compressor compression_data(nc_data.TransferData(), std::move(settings));

    /* Setup the example data for the compression */
    const bool perform_default_lossy_compression_as_well = false;
    compression_data.Setup(perform_default_lossy_compression_as_well);

    cmc::cmc_debug_msg("\n\nSetup is finished\n\n");

    compression_data.Compress();

    cmc::cmc_debug_msg("\n\nCompression is finished\n\n");

    cmc::cmc_debug_msg("Write Compressed Data");
    compression_data.WriteCompressedData("newly_example_compr_pref_t2m", 0);

    }
    #endif
    #if 1

    {
    cmc::cmc_debug_msg("\n\nDecompression Start\n\n");

    /* Decompress */
    cmc::prefix::Decompressor decoder("newly_example_compr_pref_t2m");

    decoder.Setup();
    decoder.DecompressVariable("t2m");
    
    }
    #endif
    //compression_data.WriteCompressedDataEGU("prefix_compressed_data.nc");
    
    //compression_data.WriteCompressedData("prefix_compressed_data.nc");

    //compression_data.WriteVTKFile("example_compr_vtk");

    //cmc::AmrMesh mesh33;
    //cmc::Var hihi = cmc::Var(cmc::CmcType::Int32_t);
    //hihi = compression_data.GetVariable(0);
    //hihi.SetAmrMesh(mesh33);
    /* Comrpessed data shpould look like: 
        103, 113, 105, 115, 123, 133, 125, 135, 106.5, 116.5, -2, -2, 126.5, 136.5, -2, -2, 144, -2, 146.5, -2, -2, -2,
    */
    //std::cout << "Jetzt Write COMPRESSED" << std::endl;
    //compression_data.WriteCompressedData("era5_land_t2m_abs_compr.nc");

    //compression_data.SupplementarySZLikeCompression();

    #if 0
    01000011 00000000 00000000 00000000
    00000000 00011010 01000000 10100000

    01000011 10100000 01000000 00011010

    01000011 10100000 01000000 00011000
    01000011010000001000000000110100
    #endif
    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
