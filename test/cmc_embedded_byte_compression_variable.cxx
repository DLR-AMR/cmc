#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "input/cmc_netcdf.hxx"
#include "lossless/cmc_embedded_prefix_extraction_compression.hxx"
#include "lossless/cmc_embedded_prefix_extraction_decompression.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "input/cmc_binary_reader.hxx"
//#include "lossless/cmc_embedded_byte_decompression_variable.hxx"

#include "compression_io/cmc_compression_output.hxx"
#include "compression_io/cmc_decompression_input.hxx"

#include <numeric>
#include <algorithm>
#include <memory>
#include <vector>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {
    /* Read in data */
    #if 0
    /* Specify the path of the netCDF file */
    const std::string file = "../../data/era5_reanalysis_data_16_11_23.nc";

    /* Create an object which interatcs with the file and opens it */
    cmc::input::netcdf::Data nc_data(file, cmc::nc::OpeningMode::Serial, MPI_COMM_WORLD);

    /* Define a hyperlsab corresponding to the data we want to read */
    cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1440),
                             cmc::DimensionInterval(cmc::Dimension::Lat, 0, 721));

    /* Inquire the hyperslab of data for the given variables */
    nc_data.InquireVariables(hyperslab, "t2m");

    /* Get the variable */
    std::vector<cmc::input::Var> input_variables = nc_data.TransferData();

    /* Close the file, since we have gathered the data we wanted */
    nc_data.CloseFileHandle();

    #else 


    const std::string file = "../../data/100x500x500/Wf48.bin.f32";
    cmc::CmcType type(cmc::CmcType::Float);
    const std::string name("precipf48");
    const int id = 1;
    //const size_t num_elements = 500 * 500 * 100;

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
    cmc::CmcUniversalType missing_value(static_cast<float>(13.4)); //WF
    //cmc::CmcUniversalType missing_value(static_cast<float>(48.25)); //VF

    const size_t num_elements = 500 * 500 * 100;
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 100)
                          );

    //const size_t num_elements = 500 * 500;
    //cmc::DataLayout layout(cmc::DataLayout::Lat_Lon);
    //cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 500),
    //                      cmc::DimensionInterval(cmc::Dimension::Lat, 0, 500)
    //                      );
    //
    
    cmc::input::binary::Reader binary_reader(file);

    cmc::input::Var variable = binary_reader.CreateVariableFromBinaryData(type, name, id, num_elements, missing_value, layout, domain);
    variable.SetMPIComm(MPI_COMM_WORLD);
    std::vector<cmc::input::Var> input_variables{std::move(variable)};

    #endif






    cmc::lossless::embedded::prefix::EmbeddedCompressionVariable<float> var(input_variables.front());


    #if 1
    const std::vector<cmc::SerializedCompressionValue<sizeof(float)>> initial_vals = var.GetData();

    var.Compress();


    cmc::compression_io::Writer writer("multi_res_example_lossless_compression_output.cmc");
    writer.SetVariable(&var);
    writer.Write();

    #endif

    #if 0
    cmc::compression_io::Reader reader("multi_res_example_lossless_compression_output.cmc", MPI_COMM_WORLD);
    std::unique_ptr<cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<float>> decompression_var = reader.ReadEmbeddedVariableForDecompression<float>("t2m");

    decompression_var->Decompress();

    const std::vector<cmc::SerializedCompressionValue<sizeof(float)>> decompressed_byte_values = decompression_var->GetDecompressedData();

    for (size_t idx = 0; idx < decompressed_byte_values.size(); ++idx)
    {
        const float initial_val = initial_vals[idx].ReinterpretDataAs<float>();
        const float decompressed_val = decompressed_byte_values[idx].ReinterpretDataAs<float>();
        //cmc::cmc_debug_msg("Decompressed value is: ", val, " and initial value is: ", data[idx], " count of significant bits: ", cr_iter->GetCountOfSignificantBits());
        //cmc::ExpectTrue(cmc::ApproxCompare(val, data[idx]));
        if (not cmc::ApproxCompare(initial_val, decompressed_val))
        {
            cmc::cmc_debug_msg("unequal values for idx: ", idx, ", initial val: ", initial_val, ", decompressed val: ", decompressed_val);
        }
    }

    #endif








    #if 0
    cmc::lossless::multi_res::CompressionVariable<float> var("test_var", forest, data);

    var.Compress();

    cmc::compression_io::Writer writer("multi_res_example_lossless_compression_output.cmc");
    
    writer.SetVariable(&var);

    writer.Write();

    cmc::compression_io::Reader reader("multi_res_example_lossless_compression_output.cmc", MPI_COMM_WORLD);

    std::unique_ptr<cmc::decompression::AbstractByteDecompressionVariable<float>> decompression_var = reader.ReadVariableForDecompression<float>("test_var");

    decompression_var->Decompress();
    cmc::ExpectTrue(decompression_var->Size() == data.size());

    std::vector<cmc::SerializedCompressionValue<sizeof(float)>> decompressed_byte_values = decompression_var->GetDecompressedData();

    int idx = 0;
    for (auto cr_iter = decompressed_byte_values.begin(); cr_iter != decompressed_byte_values.end(); ++cr_iter, ++idx)
    {
        const float val = cr_iter->ReinterpretDataAs<float>();
        //cmc::cmc_debug_msg("Decompressed value is: ", val, " and initial value is: ", data[idx], " count of significant bits: ", cr_iter->GetCountOfSignificantBits());
        cmc::ExpectTrue(cmc::ApproxCompare(val, data[idx]));
    }

    #endif

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
