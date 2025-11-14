#include "cmc.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "input/cmc_netcdf.hxx"
#include "lossy/cmc_embedded_trimmed_multi_res_extraction_compression.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "input/cmc_binary_reader.hxx"

#include "compression_io/cmc_compression_output.hxx"
#include "compression_io/cmc_decompression_input.hxx"

#include <numeric>
#include <algorithm>
#include <memory>
#include <vector>

const double kRelMaxError = 0.000001;
const double kAbsMaxError = 0.00005;

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {

    /* Read in data */
    #if 1
    const std::string file = "../../programs/data/100x500x500/Uf48.bin.f32";
    #else
    const std::string file = "../../programs/data/SDRBENCH-EXASKY-NYX-512x512x512/baryon_density.f32";
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
    cmc::CmcUniversalType missing_value(static_cast<float>(40.0)); //UF
    //cmc::CmcUniversalType missing_value(static_cast<float>(48.25)); //VF
    //cmc::CmcUniversalType missing_value(static_cast<float>(13.4)); //WF

    //Missing values NYX data
    //cmc::CmcUniversalType missing_value(static_cast<float>(31868.0)); //velocity x
    //cmc::CmcUniversalType missing_value(static_cast<float>(56507.0)); //velocity y
    //cmc::CmcUniversalType missing_value(static_cast<float>(33387.0)); //velocity z
    //cmc::CmcUniversalType missing_value(static_cast<float>(4784.0)); //temp
    //cmc::CmcUniversalType missing_value(static_cast<float>(13780.0)); //dark matter
    //cmc::CmcUniversalType missing_value(static_cast<float>(115864.0)); //baryon denisty

    {

    #if 1
    const size_t num_elements = 500 * 500 * 100;
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 500),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 100)
                          );
    #else
    const size_t num_elements = 512 * 512 * 512;
    cmc::DataLayout layout(cmc::DataLayout::Lev_Lat_Lon);
    cmc::GeoDomain domain(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 512),
                          cmc::DimensionInterval(cmc::Dimension::Lat, 0, 512),
                          cmc::DimensionInterval(cmc::Dimension::Lev, 0, 512)
                          );
    #endif
    
    /* Generate input variables from the binary file */
    cmc::input::binary::Reader binary_reader(file);
    cmc::input::Var variable = binary_reader.CreateVariableFromBinaryData(type, name, id, num_elements, missing_value, layout, domain);
    variable.SetMPIComm(MPI_COMM_SELF);
    std::vector<cmc::input::Var> input_variables{std::move(variable)};

    /* Allocate the compression settings */
    cmc::CompressionSettings settings;
    /* Set a relative error criterion */
    //settings.SetGeneralErrorCriterion(cmc::CompressionCriterion::RelativeErrorThreshold, kRelMaxError);
    /* Set an absolute error criterion */
    settings.SetGeneralErrorCriterion(cmc::CompressionCriterion::AbsoluteErrorThreshold, kAbsMaxError);


    /* Setup an embedded MultiResAMR compression variable from the input variables */
    cmc::lossy::embedded::multi_res::EmbeddedCompressionVariable<float> var(settings, input_variables.front());
    
    /* Perform the compression */
    var.Compress();

    {
        /* Write out the compressed data to disk */
        cmc::compression_io::Writer writer("trimmed_multi_res_example_lossy_compression_output.cmc", MPI_COMM_SELF);
        writer.SetVariable(&var);
        writer.Write();
    }

    
    }

    #if 0
    /* Create a reader for the compressed output that has been stored */
    cmc::compression_io::Reader reader("multi_res_example_lossless_compression_output.cmc", MPI_COMM_SELF);

    /* Create an embedded decompressor from the compressed data */
    std::unique_ptr<cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<float>> decompression_var = reader.ReadEmbeddedVariableForDecompression<float>("compr_test_var");

    /* Decompress the encoded data */
    decompression_var->Decompress();

    /* De-Mortonize the data and obtain the initial ordering of the input data */
    const std::vector<float> decompressed_data = decompression_var->DeMortonizeData();

    /* Write this decompressed data out to disk, in order to be able to compare it to the intiial data */
    FILE* file_out = fopen("decompressed_data.cmc", "wb");
    fwrite(decompressed_data.data(), sizeof(float), decompressed_data.size(), file_out);
    fclose(file_out);

    cmc::cmc_debug_msg("Size of decompressed data: ", decompressed_data.size());

    #endif

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}





#if 0

#include "cmc.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "lossy/cmc_multi_res_trimmed_residuals.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_utilities.hxx"
#include "compression_io/cmc_compression_output.hxx"
#include "compression_io/cmc_decompression_input.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#endif

#include <numeric>
#include <algorithm>
#include <memory>

const double kRelMaxError = 0.01;
const double kAbsMaxError = 0.0;

static t8_locidx_t
TestAdapt ([[maybe_unused]] t8_forest_t forest,
           [[maybe_unused]] t8_forest_t forest_from,
           t8_locidx_t which_tree,
           [[maybe_unused]] const t8_eclass_t tree_class,
           t8_locidx_t lelement_id,
           [[maybe_unused]] const t8_scheme_c * ts,
           const int is_family,
           [[maybe_unused]] const int num_elements,
           [[maybe_unused]] t8_element_t * elements[])
{

    if (is_family == 1 && which_tree == 1 && lelement_id >= 8 && lelement_id < 16)
    {
        return cmc::t8::kCoarsenElements;
    } else
    {
        return cmc::t8::kLeaveElementUnchanged;
    }
}

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {

    /* Create a mesh */
    const sc_MPI_Comm comm = MPI_COMM_SELF;
    t8_cmesh_t cmesh;
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, comm, 0, 0, 0);
    const t8_scheme *scheme = t8_scheme_new_default ();
    const int initial_level = 3;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);

    forest = t8_forest_new_adapt(forest, TestAdapt, 0, 0, NULL);

    /* Create a variable for containing a value for each element */
    const t8_locidx_t num_elements = t8_forest_get_local_num_leaf_elements(forest);
    std::vector<float> data(num_elements);
    std::iota(data.begin(), data.end(), 0.0);

    /* Allocate the compression settings */
    cmc::CompressionSettings settings;
    /* Set a relative error criterion */
    //settings.SetGeneralErrorCriterion(cmc::CompressionCriterion::RelativeErrorThreshold, kRelMaxError);
    /* Set an absolute error criterion */
    settings.SetGeneralErrorCriterion(cmc::CompressionCriterion::AbsoluteErrorThreshold, kAbsMaxError);

    /* Setup an embedded MultiResAMR compression variable from the input variables */
    cmc::lossy::multi_res::CompressionVariable<float> var(settings, "test_var", forest, data); 

    var.Compress();

    cmc::compression_io::Writer writer("example_lossless_compression_output.cmc", MPI_COMM_SELF);
    
    writer.SetVariable(&var);

    writer.Write();

    //cmc::compression_io::Reader reader("example_lossless_compression_output.cmc", MPI_COMM_WORLD);
    //
    //std::unique_ptr<cmc::decompression::AbstractByteDecompressionVariable<float>> decompression_var = reader.ReadVariableForDecompression<float>("test_var");
    //
    //decompression_var->Decompress();
    //cmc::ExpectTrue(decompression_var->Size() == data.size());
    //
    //std::vector<cmc::SerializedCompressionValue<sizeof(float)>> decompressed_byte_values = decompression_var->GetDecompressedData();
    //
    //int idx = 0;
    //for (auto cr_iter = decompressed_byte_values.begin(); cr_iter != decompressed_byte_values.end(); ++cr_iter, ++idx)
    //{
    //    const float val = cr_iter->ReinterpretDataAs<float>();
    //    //cmc::cmc_debug_msg("Decompressed value is: ", val, " and initial value is: ", data[idx], " count of significant bits: ", cr_iter->GetCountOfSignificantBits());
    //    cmc::ExpectTrue(cmc::ApproxCompare(val, data[idx]));
    //}

    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}

#endif

