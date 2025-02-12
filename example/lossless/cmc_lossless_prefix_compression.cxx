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
    //const std::string file = "../../data/SDRBENCH-EXASKY-NYX-512x512x512/baryon_density.f32";
    const std::string file = "../../data/SDRBENCH-CESM-ATM-cleared-1800x3600/CLDLOW_1_1800_3600.dat";
    //const std::string file = "../../data/100x500x500/PRECIPf48.bin.f32";
    cmc::CmcType type(cmc::CmcType::Float);
    const std::string name("velo");
    const int id = 0;
    //const size_t num_elements = 512 * 512 * 256;
    const size_t num_elements = 3600 * 1800;

    //cmc::CmcUniversalType missing_value(static_cast<float>(31866790.000000)); //velocity_x
    //cmc::CmcUniversalType missing_value(static_cast<float>(4782584.0)); //Temp
    //cmc::CmcUniversalType missing_value(static_cast<float>(31866787.0)); //velocity_x
    //cmc::CmcUniversalType missing_value(static_cast<float>(56505650.0)); //velocity_y
    //cmc::CmcUniversalType missing_value(static_cast<float>(33386285.0)); //velocity_z
    //cmc::CmcUniversalType missing_value(static_cast<float>(13780.0)); //dark_matter_density.f32
    //cmc::CmcUniversalType missing_value(static_cast<float>(115864.0)); //baryon_density.f32

    //cmc::CmcUniversalType missing_value(static_cast<float>(0.5725));//AEROD_v_1_1800_3600.dat
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.925));//CLDHGH_1_1800_3600.dat
    cmc::CmcUniversalType missing_value(static_cast<float>(0.98));//CLDLOW_1_1800_3600.dat
    
    
    cmc::DataLayout layout(cmc::DataLayout::Lat_Lon);

    #if 1
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 3600),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 1800)
                          //cmc::DimensionInterval(cmc::Dimension::Lev, 0, 512)
                          );
    #endif

    cmc::GeoDomain global_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 512),
                                 cmc::DimensionInterval(cmc::Dimension::Lat, 0, 512),
                                 cmc::DimensionInterval(cmc::Dimension::Lev, 0, 512)
                                 );

    cmc::GeoDomain sub_domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 512),
                              cmc::DimensionInterval(cmc::Dimension::Lat, 0, 512),
                              cmc::DimensionInterval(cmc::Dimension::Lev, 256, 512)
                              );

    cmc::bin_reader::Reader binary_reader(file);

    cmc::InputVar variable = binary_reader.CreateVariableFromBinaryData(type, name, id, num_elements, missing_value, layout, domain);
    //cmc::InputVar variable = binary_reader.CreateSubDomainVariableFromBinaryData(type, name, id, missing_value, layout, global_domain, sub_domain);

    std::vector<cmc::InputVar> vars;
    vars.emplace_back(std::move(variable));

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




[cmc] DEBUG: Index 3000032 is 0010111110 0100101110111010101111 = 1240751
[cmc] DEBUG: Index 3000033 is 0010111110 1000011011100100010000 = 2210064
[cmc] DEBUG: Index 3000034 is 0010111110 1010110001111011000001 = 2825921
[cmc] DEBUG: Index 3000035 is 0010111110 1100000001100111010100 = 3152340
[cmc] DEBUG: Index 3000036 is 0010111110 0001100110110111100110 =  421350
[cmc] DEBUG: Index 3000037 is 0010111110 0011010101100100001000 =  874760
[cmc] DEBUG: Index 3000038 is 0010111110 0110110110110001111101 = 1797245
[cmc] DEBUG: Index 3000039 is 0010111110 0110101000110000000100 = 1739780

1111111111111111111111 = 4194303

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
