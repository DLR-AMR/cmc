#include "cmc.h"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "netcdf/cmc_netcdf.hxx"
#include "lossy/cmc_amr_lossy_compression_settings.hxx"
#include "lossy/cmc_amr_lossy_compression.hxx"

#include <cfenv>
int
main(void)
{
    /* Initialize cmc */
    cmc_initialize();
    
    {
    
    #if 0
    //const std::string file = "../../data/era5_reanalysis_t2m_tc03_13_12_23.nc";
    //const std::string file = "../../data/tas_decreg_europe_v20140120_20010101_20010131.nc";
    const std::string file = "../../data/era5_reanalysis_data_16_11_23.nc";
    //const std::string file = "../../data/era5_reanalysis_pressure_lvls_fixed_time.nc";

    /* Create an object which interatcs with the file and opens it */
    cmc::NcData nc_data(file, cmc::NcOpeningMode::Serial);

    //cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1440),
    //                         cmc::DimensionInterval(cmc::Dimension::Lat, 0, 721)
    //                         );
    ///* Inquire the hyperslab of data for the given variables */
    //nc_data.InquireVariables(hyperslab, "t2m");

    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1440),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, 721)
                             );
    /* Inquire the hyperslab of data for the given variables */
    nc_data.InquireVariables(hyperslab, "t2m");

    /* Close the file, since we have gathered the data we wanted */
    nc_data.CloseFileHandle();

    /* Create compression settings */
    cmc::CompressionSettings settings;

    //const double abs_max_err = 1.0;
    //settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    const double rel_max_err = 0.01;
    settings.SetRelativeErrorCriterion(rel_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    #endif


    #if 1
    /* Specify the path of the netCDF file */
    const std::string file = "../../data/era5_reanalysis_t2m_tc03_13_12_23.nc";
    //const std::string file = "../../data/era5_reanalysis_pressure_lvls_fixed_time.nc";

    /* Create an object which interatcs with the file and opens it */
    cmc::NcData nc_data(file, cmc::NcOpeningMode::Serial);

    /* Define a hyperlsab corresponding to the data we want to read */
    //cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1440),
    //                         cmc::DimensionInterval(cmc::Dimension::Lat, 0, 721),
    //                         cmc::DimensionInterval(cmc::Dimension::Lev, 36, 37)
    //                         );

    nc_data.SetHintHeightDimension(2); //Set time as height dimension

    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1440),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, 721),
                             cmc::DimensionInterval(cmc::Dimension::Lev, 0, 360)
                             );
    /* Inquire the hyperslab of data for the given variables */
    nc_data.InquireVariables(hyperslab, "tco3");

    /* Close the file, since we have gathered the data we wanted */
    nc_data.CloseFileHandle();

    /* Create compression settings */
    cmc::CompressionSettings settings;

    //const double abs_max_err = 2.5;
    //settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    const double rel_max_err = 0.01;
    settings.SetRelativeErrorCriterion(rel_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    cmc::SplitVariable split(cmc::kSplitAllVariables, cmc::Dimension::Lev);
    settings.SplitVariableByDimension(split);

//
    //cmc::SplitVariable split(cmc::kSplitAllVariables, cmc::Dimension::Lev);
    //settings.SplitVariableByDimension(split);
#endif

#if 0
    //Gutes rel error ergebnis

    /* Specify the path of the netCDF file */
    const std::string file = "../../data/mptrac_era5_2021_07_01_00.nc";
    //const std::string file = "../../data/era5_reanalysis_pressure_lvls_fixed_time.nc";

    /* Create an object which interatcs with the file and opens it */
    cmc::NcData nc_data(file, cmc::NcOpeningMode::Serial);

    nc_data.SetHintHeightDimension(6); //Set plev as height dimension

    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1200),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, 601),
                             cmc::DimensionInterval(cmc::Dimension::Lev, 0, 137)
                             );
    /* Inquire the hyperslab of data for the given variables */
    nc_data.InquireVariables(hyperslab, "t");

    /* Close the file, since we have gathered the data we wanted */
    nc_data.CloseFileHandle();

    
    /* Create compression settings */
    cmc::CompressionSettings settings;

    const double abs_max_err = 0.25;
    settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    //const double rel_max_err = 0.00125;
    //settings.SetRelativeErrorCriterion(rel_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    cmc::SplitVariable split(cmc::kSplitAllVariables, cmc::Dimension::Lev);
    settings.SplitVariableByDimension(split);
#endif
    //const double abs_specific_err = 0.01; 
    //const cmc::GeoDomain spec_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 400, 900),
    //                                 cmc::DimensionInterval(cmc::Dimension::Lat, 200, 500));
    //settings.SetCertainErrorForDomain(cmc::CompressionCriterion::RelativeErrorThreshold, abs_specific_err, spec_domain, cmc::kErrorCriterionHoldsForAllVariables);
    //
    //const double abs_specific_err2 = 0.5;  
    //const cmc::GeoDomain spec_domain2(cmc::DimensionInterval(cmc::Dimension::Lon, 600, 800),
    //                                 cmc::DimensionInterval(cmc::Dimension::Lat, 300, 400));
    //settings.SetCertainErrorForDomain(cmc::CompressionCriterion::AbsoluteErrorThreshold, abs_specific_err2, spec_domain2, cmc::kErrorCriterionHoldsForAllVariables);

    //cmc::SplitVariable split(cmc::kSplitAllVariables, cmc::Dimension::Lev);
    //settings.SplitVariableByDimension(split);

    //const double rel_max_err = 0.01;
    //settings.SetRelativeErrorCriterion(rel_max_err, cmc::kErrorCriterionHoldsForAllVariables);

    /* Create the compression data */
    cmc::CompressionData compression_data(nc_data.TransferData(), std::move(settings));

    /* Setup the example data for the compression */
    compression_data.Setup();

    compression_data.Compress(cmc::CompressionMode::OneForOne);

    //compression_data.WriteVTKFile("example_compr_vtk");

    //cmc::AmrMesh mesh33;
    //cmc::Var hihi = cmc::Var(cmc::CmcType::Int32_t);
    //hihi = compression_data.GetVariable(0);
    //hihi.SetAmrMesh(mesh33);
    /* Comrpessed data shpould look like: 
        103, 113, 105, 115, 123, 133, 125, 135, 106.5, 116.5, -2, -2, 126.5, 136.5, -2, -2, 144, -2, 146.5, -2, -2, -2,
    */
    std::cout << "Jetzt Write COMPRESSED" << std::endl;
    compression_data.WriteCompressedData("test_era5_tco3_compr.nc");

    //compression_data.SupplementarySZLikeCompression();


    }

    /* Finalize cmc */
    cmc_finalize();

    return 0;
}
