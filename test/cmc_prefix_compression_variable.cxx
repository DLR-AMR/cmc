#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "lossless/cmc_prefix_extraction_compression.hxx"
#include "lossless/cmc_prefix_extraction_decompression.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_utilities.hxx"
#include "compression_io/cmc_compression_output.hxx"
#include "compression_io/cmc_decompression_input.hxx"

#include <numeric>
#include <algorithm>
#include <memory>
static t8_locidx_t
TestAdapt ([[maybe_unused]] t8_forest_t forest,
           [[maybe_unused]] t8_forest_t forest_from,
           t8_locidx_t which_tree,
           t8_locidx_t lelement_id,
           [[maybe_unused]] t8_eclass_scheme_c * ts,
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
    const sc_MPI_Comm comm = sc_MPI_COMM_WORLD;
    t8_cmesh_t cmesh;
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, comm, 0, 0, 0);
    t8_scheme_cxx_t* scheme = t8_scheme_new_default_cxx ();
    const int initial_level = 3;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);

    forest = t8_forest_new_adapt(forest, TestAdapt, 0, 0, NULL);

    /* Create a variable for containing a value for each element */
    const t8_locidx_t num_elements = t8_forest_get_local_num_elements(forest);
    std::vector<float> data(num_elements);
    std::iota(data.begin(), data.end(), 0.0);

    cmc::lossless::prefix::CompressionVariable<float> var("test_var", forest, data);

    var.Compress();

    cmc::compression_io::Writer writer("example_lossless_compression_output.cmc");
    
    writer.SetVariable(&var);

    writer.Write();

    cmc::compression_io::Reader reader("example_lossless_compression_output.cmc", MPI_COMM_WORLD);

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

    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
