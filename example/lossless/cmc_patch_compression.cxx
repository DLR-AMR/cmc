#include "cmc.hxx"
#include "patch/lossless/cmc_multi_res_extraction.hxx"
#include "patch/lossless/cmc_multi_res_decompression.hxx"

#include "input/cmc_binary_file_reader.hxx"

#include <vector>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {
    
    #if 0
    constexpr int32_t DIM  = 2;

    std::vector<float> init_data{2.34,2.79,3.10,2.56,
                                 3.21,3.07,2.99,2.76,
                                 3.01,2.88,2.34,2.19,
                                 2.27,2.45,2.76,2.56};
    std::array<int, DIM> init_dims{4,4};
    
    cmc::serial::patch::lossless::multi_res::CompressionVariable<float, DIM> compression_variable(std::move(init_data), init_dims);
 
    compression_variable.Compress();
    compression_variable.WriteData("cmc_patch_example_data.cmc");

    cmc::serial::patch::lossless::multi_res::CompressionInfo cr = cmc::serial::patch::lossless::multi_res::ReadCompressionData("cmc_patch_example_data.cmc");
    cmc::cmc_global_msg("\n\n\n");
    cmc::serial::patch::lossless::multi_res::DecompressionVariable<float, DIM> decompression_variable("cmc_patch_example_data.cmc");
  
    decompression_variable.Decompress();

    const std::vector<float> decompressed_data = decompression_variable.GetDecompressedData();
    cmc::cmc_global_msg("Decompressed data size: ", decompressed_data.size());
    for(const auto& val : decompressed_data)
    {
        //cmc::cmc_global_msg("Decompressed Value: ", val);
    }


    #endif

    cmc::cmc_global_msg("\n\n\nNext Test\n");

    const std::string file_name("/home/niklas/software/data/100x500x500/CLOUDf48.bin.f32");
    constexpr int DataFileDimension = 3;
    const std::array<int, DataFileDimension> dim_lengths{100, 500, 500};
    const std::endian file_endianness = std::endian::little;

    cmc::input::binary_file::Reader<float, DataFileDimension> bin_reader(file_name, dim_lengths, file_endianness);
    const std::vector<float> init_file_data = bin_reader.ReadData();
    std::vector<float> copy_init_data = init_file_data;

    cmc::serial::patch::lossless::multi_res::CompressionVariable<float, DataFileDimension> compression_variable2(std::move(copy_init_data), dim_lengths);
 
    compression_variable2.Compress();
    compression_variable2.WriteData("cmc_big_patch_example_data.cmc");

    cmc::serial::patch::lossless::multi_res::DecompressionVariable<float, DataFileDimension> decompression_variable2("cmc_big_patch_example_data.cmc");
  
    decompression_variable2.Decompress();

    const std::vector<float> decompressed_data2 = decompression_variable2.GetDecompressedData();

    cmc::cmc_debug_msg("Size init data: ", init_file_data.size(), ", Size decompressed: ", decompressed_data2.size());

    for (size_t idx{0}; idx < init_file_data.size(); ++idx)
    {
        if (init_file_data[idx] != decompressed_data2[idx])
        {
            cmc::cmc_err_msg("Unequal at pos ", idx);
        }
    }

#if 0

    0000000000000000000000000000000000000000000000000000001111111000

    0000000000000000000000000000000000000000000000000000001100011000
    /* Read in data */
    #if 1
    const std::string file = "../data/100x500x500/TCf48.bin.f32";
    //const std::string file = "../data/era5_reanalysis_t2m_tc03_13_12_23.nc";
    #else
    const std::string file = "../data/SDRBENCH-EXASKY-NYX-512x512x512/baryon_density.f32";
    #endif

    cmc::CmcType type(cmc::CmcType::Float);
    const std::string name("compr_test_var");
    const int id = 1;

    //Missing Values Hurricane ISABEL Dataset
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.0025));//CLOUD
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.5));//CLOUD.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(3224.4)); //P
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.00755)); //PRECIP
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.0)); //PRECIPf48.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.00205));//QCLOUD
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.5));//QCLOUD.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.007295)); //QGraup
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.0)); //QGraup.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.00085));//QICE
    //cmc::CmcUniversalType missing_value(static_cast<float>(-3.0));//QICE.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.0065));//QRAINf48
    //cmc::CmcUniversalType missing_value(static_cast<float>(-2.0));//QRAINf48.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(0.000875)); //QSNOW
    //cmc::CmcUniversalType missing_value(static_cast<float>(-3.0)); //QSNOW.log10
    //cmc::CmcUniversalType missing_value(static_cast<float>(30.0)); //QVAPOR
    //cmc::CmcUniversalType missing_value(static_cast<float>(29.65)); //TC
    //cmc::CmcUniversalType missing_value(static_cast<float>(40.0)); //UF
    //cmc::CmcUniversalType missing_value(static_cast<float>(48.25)); //VF
    //cmc::CmcUniversalType missing_value(static_cast<float>(13.4)); //WF
    
    //Missing values NYX data
    //cmc::CmcUniversalType missing_value(static_cast<float>(31868.0)); //velocity x
    //cmc::CmcUniversalType missing_value(static_cast<float>(56507.0)); //velocity y
    //cmc::CmcUniversalType missing_value(static_cast<float>(33387.0)); //velocity z
    //cmc::CmcUniversalType missing_value(static_cast<float>(4784.0)); //temp
    //cmc::CmcUniversalType missing_value(static_cast<float>(13780.0)); //dark matter
    //cmc::CmcUniversalType missing_value(static_cast<float>(115864.0)); //baryon denisty


    const size_t num_elements = 500 * 500 * 100;
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 100)
                          );

    cmc::CmcUniversalType arbitrary_missing_value(static_cast<float>(-1000.0));

    /* Generate input variables from the binary file */
    cmc::input::binary::Reader binary_reader(file);
    cmc::input::Var variable = binary_reader.CreateVariableFromBinaryData(type, name, id, num_elements, arbitrary_missing_value, layout, domain);
    //variable.SetMPIComm(MPI_COMM_SELF);
    std::vector<cmc::input::Var> input_variables{std::move(variable)};


    {
    /* Setup an embedded PrefixAMR (with plain suffix encoding) compression variable from the input variables */           
    //cmc::patch::lossless::prefix::plain_suffix::PatchCompressionVariable<float, 3> var(input_variables.front());
    cmc::patch::lossless::multi_res::PatchCompressionVariable<float, 3> var(input_variables.front());

    /* Perform the compression */
    var.Compress();

#if 1

    #if 1
    {
        /* Write out the compressed data to disk */
        cmc::compression_io::serial::Writer writer("prefix_example_lossless_compression_output");
        writer.SetVariable(&var);
        writer.Write();
    }
    #endif
    }

    #if 1

    /* Create a reader for the compressed output that has been stored */
    cmc::compression_io::serial::Reader reader("prefix_example_lossless_compression_output");
    
    /* Create an embedded decompressor from the compressed data */
    std::unique_ptr<cmc::patch::IPatchDecompressionVariable<float>> decompression_var = reader.ReadPatchVariableForDecompression<float>(name);

    /* Decompress the encoded data */
    decompression_var->Decompress();

    #if 1
    /* De-Mortonize the data and obtain the initial ordering of the input data */
    const std::vector<float> decompressed_data = decompression_var->GetDecompressedData();
    
    /* Write this decompressed data out to disk, in order to be able to compare it to the intiial data */
    FILE* file_out = fopen("decompressed_data.cmc", "wb");
    fwrite(decompressed_data.data(), sizeof(float), decompressed_data.size(), file_out);
    fclose(file_out);

    cmc::cmc_debug_msg("Size of decompressed data: ", decompressed_data.size());
    #endif

    #endif

#endif



#endif



    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
