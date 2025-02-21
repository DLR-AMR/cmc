#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "input/cmc_netcdf.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "lossless/cmc_diff_compression.hxx"
#include "lossless/cmc_residual_compression.hxx"
#include "lossless/cmc_diff_decompression.hxx"
#include "utilities/cmc_binary_reader.hxx"
#include "comparison/cmc_test_comparison_compression.hxx"
#include "comparison/cmc_test_comparison_decompression.hxx"

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {
    
    #if 1
    //const std::string file = "../../data/era5_reanalysis_t2m_tc03_13_12_23.nc";
    const std::string file = "../../data/mptrac_era5_2021_07_01_00.nc";
    //const std::string file = "../../data/v8.0_FT2022_GHG_CO2_2022_IND_PROCESSES_emi.nc";
    /* Create an object which interatcs with the file and opens it */
    cmc::input::netcdf::Data nc_data(file, cmc::nc::OpeningMode::Serial);

    //nc_data.SetHintHeightDimension(2); //Set time as height dimension
    nc_data.SetHintHeightDimension(6); //Set time as height dimension

    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1200),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, 601),
                             cmc::DimensionInterval(cmc::Dimension::Lev, 0, 137)
                             );

    const std::string var_name = "t";
    /* Inquire the hyperslab of data for the given variables */
    nc_data.InquireVariables(hyperslab, std::string(var_name));

    /* Close the file, since we have gathered the data we wanted */
    nc_data.CloseFileHandle();

    #else

    #if 1
    const std::string file = "../../data/SDRBENCH-EXASKY-NYX-512x512x512/dark_matter_density.f32";
    //const std::string file = "../../data/SDRBENCH-CESM-ATM-cleared-1800x3600/CLDLOW_1_1800_3600.dat";
    cmc::CmcType type(cmc::CmcType::Float);
    const std::string var_name("precipf48log10");
    const int id = 0;
    
    const size_t num_elements = 512 * 512 * 256;
    //const size_t num_elements = 1800 * 3600;
    
    //cmc::CmcUniversalType missing_value(static_cast<float>(31866790.000000)); //velocity_x
    //cmc::CmcUniversalType missing_value(static_cast<float>(4782584.0)); //Temp
    //cmc::CmcUniversalType missing_value(static_cast<float>(31866787.0)); //velocity_x
    //cmc::CmcUniversalType missing_value(static_cast<float>(56505650.0)); //velocity_y
    //cmc::CmcUniversalType missing_value(static_cast<float>(33386285.0)); //velocity_z
    cmc::CmcUniversalType missing_value(static_cast<float>(13780.0)); //dark_matter_density.f32
    //cmc::CmcUniversalType missing_value(static_cast<float>(115864.0)); //baryon_density.f32


    //cmc::CmcUniversalType missing_value(static_cast<float>(0.5725));//AEROD_v_1_1800_3600.dat
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.925));//CLDHGH_1_1800_3600.dat
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.98));//CLDLOW_1_1800_3600.dat
    


    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    //cmc::DataLayout layout(cmc::DataLayout::Lat_Lon);
    
    
    #if 0
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 512),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 512)
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 512)
                          );

    cmc::bin_reader::Reader binary_reader(file);

    //cmc::InputVar variable = binary_reader.CreateVariableFromBinaryData(type, var_name, id, num_elements, missing_value, layout, domain);
    cmc::InputVar variable = binary_reader.CreateVariableFromBinaryData(type, var_name, id, num_elements, missing_value, layout, domain);
    std::vector<cmc::InputVar> vars;
    vars.emplace_back(std::move(variable));

    #else


    cmc::GeoDomain global_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 512),
                                 cmc::DimensionInterval(cmc::Dimension::Lat, 0, 512),
                                 cmc::DimensionInterval(cmc::Dimension::Lev, 0, 512)
                                 );

    #if 0
    cmc::GeoDomain sub_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 512),
                              cmc::DimensionInterval(cmc::Dimension::Lat, 0, 512),
                              cmc::DimensionInterval(cmc::Dimension::Lev, 0, 256)
                              );
    #else
    cmc::GeoDomain sub_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 512),
                              cmc::DimensionInterval(cmc::Dimension::Lat, 0, 512),
                              cmc::DimensionInterval(cmc::Dimension::Lev, 256, 512)
                              );
    #endif

    cmc::bin_reader::Reader binary_reader(file);

    //cmc::InputVar variable = binary_reader.CreateVariableFromBinaryData(type, var_name, id, num_elements, missing_value, layout, domain);
    cmc::InputVar variable = binary_reader.CreateSubDomainVariableFromBinaryData(type, var_name, id, missing_value, layout, global_domain, sub_domain);
    std::vector<cmc::InputVar> vars;
    vars.emplace_back(std::move(variable));

    #endif

    #else


    #if 1
    const std::string file = "../../data/SDRBENCH-EXASKY-NYX-512x512x512/temperature.f32";
    //const std::string file = "../../data/100x500x500/Wf48.bin.f32";
    cmc::CmcType type(cmc::CmcType::Float);
    const std::string name("precipf48log10");
    const int id = 0;
    const size_t num_elements = 500 * 500 * 100;

    //Missing Values Hurricane ISABEL Dataset

    //cmc::CmcUniversalType missing_value(static_cast<float>(3224.4)); //P
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.00755)); //PRECIP
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.0)); //PRECIPf48.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.007295)); //QGraup
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.0)); //QGraup.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.0));//QRAINf48.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.0065));//QRAINf48
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.0025));//CLOUD
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.5));//CLOUD-log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.00205));//QCLOUD
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.5));//QCLOUD.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.00085));//QICE
    //cmc::CmcUniversalType missing_value(static_cast<float>(-3.0));//QICE.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(-3.0)); //QSNOW.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(30.0)); //QVAPOR
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.000875)); //QSNOW
    //cmc::CmcUniversalType missing_value(static_cast<float>(29.65)); //TC
    //cmc::CmcUniversalType missing_value(static_cast<float>(40.0)); //UF
    //cmc::CmcUniversalType missing_value(static_cast<float>(13.4)); //WF
    //cmc::CmcUniversalType missing_value(static_cast<float>(48.25)); //VF


    //Missing Values NYX Dataset
    //cmc::CmcUniversalType missing_value(static_cast<float>(31866790.000000)); //velocity_x
    cmc::CmcUniversalType missing_value(static_cast<float>(4782584.0)); //Temp
    //cmc::CmcUniversalType missing_value(static_cast<float>(31866787.0)); //velocity_x
    //cmc::CmcUniversalType missing_value(static_cast<float>(56505650.0)); //velocity_y
    //cmc::CmcUniversalType missing_value(static_cast<float>(33386285.0)); //velocity_z
    //cmc::CmcUniversalType missing_value(static_cast<float>(13780.0)); //dark_matter_density.f32
    //cmc::CmcUniversalType missing_value(static_cast<float>(115864.0)); //baryon_density.f32



    #if 0
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 512),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 512),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 128)
                          );
    #else
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 100)
                          );
    #endif

    cmc::bin_reader::Reader binary_reader(file);

    cmc::InputVar variable = binary_reader.CreateVariableFromBinaryData(type, name, id, num_elements, missing_value, layout, domain);

    std::vector<cmc::InputVar> vars;
    vars.emplace_back(std::move(variable));

    #endif

    #endif

    #endif


#if 0
0000000000000 1100010110001010011 = 404563 mit num significant bits: 19
0000000000000 1100000010100001001 = 394505 mit num significant bits: 19
0000000000000 1110010010110111110 = 468414 mit num significant bits: 19
0000000000000 1110011110010001111 = 474255 mit num significant bits: 19
0000000000000 0011011110010010011 = 113811 mit num significant bits: 19
0000000000000 0011001001000100001 = 102945 mit num significant bits: 19
0000000000000 0101010100011110001 = 174321 mit num significant bits: 19
0000000000000 0101011110101000100 = 179524 mit num significant bits: 19

LZC: 5 + 2 + 6 + 0 + 5 + 1 + 6 = 25 (+7 one bits induced)

Ordered increasingly
0000000000000 0011001001000100001 = 102945 mit num significant bits: 19
0000000000000 0011011110010010011 = 113811 mit num significant bits: 19
0000000000000 0101010100011110001 = 174321 mit num significant bits: 19
0000000000000 0101011110101000100 = 179524 mit num significant bits: 19
0000000000000 1100000010100001001 = 394505 mit num significant bits: 19
0000000000000 1100010110001010011 = 404563 mit num significant bits: 19
0000000000000 1110010010110111110 = 468414 mit num significant bits: 19
0000000000000 1110011110010001111 = 474255 mit num significant bits: 19

with xor 
LZC: 5 + 1 + 6 + 0 + 5 + 7 + 6 = 30 (+7 one bits induced)

Ordered increasingly + diff
    10101001110010 14
  1110110001011110 16
     1010001010011 13
110100011111000101 18
    10011101001010 14
  1111100101101011 16
     1011011010001 13

LZC: 29 + (7 one bits induced)
#endif

#if 1
    {
    
    /* Create the compression data */
    cmc::diff::Compressor compression_data(nc_data.TransferData()); //netCDF data
    //cmc::diff::Compressor compression_data(std::move(vars)); //Binary data
    
    //cmc::test_comparison::light_amr_pcp::Compressor compression_data(std::move(vars)); //Binary data

    //compression_data.SetSplitVariable(cmc::SplitVariable(cmc::kSplitAllVariables, cmc::Dimension::Lev));

    /* Setup the example data for the compression */
    compression_data.Setup();

    cmc::cmc_debug_msg("\n\nSetup is finished\n\n");

    compression_data.Compress();

    cmc::cmc_debug_msg("\n\nCompression is finished\n\n");

    cmc::cmc_debug_msg("Write Compressed Data");
    compression_data.WriteCompressedData("PRECIPf48.log10.bin.f32.new_diff_scheme.cmc", 0);

    }
#endif

    {
    #if 1
    cmc::cmc_debug_msg("\n\nDecompression Start\n\n");

    /* Decompress */
    //cmc::diff::Decompressor decoder("PRECIPf48.log10.bin.f32.new_diff_scheme.cmc");
    cmc::test_comparison::light_amr_pcp::Decompressor decoder("PRECIPf48.log10.bin.f32.new_diff_scheme.cmc");
    decoder.Setup();
    decoder.DecompressVariable("precipf48log10");
    //decoder.DecompressVariable(var_name);
    #endif
    }

    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
