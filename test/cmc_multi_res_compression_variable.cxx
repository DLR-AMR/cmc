#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "lossless/cmc_multi_res_extraction_compression.hxx"
#include "lossless/cmc_multi_res_extraction_decompression.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"

#include <numeric>
#include <algorithm>

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
    const sc_MPI_Comm comm = sc_MPI_COMM_SELF;
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

    cmc::lossless::multi_res::CompressionVariable<float> var("test_var", forest, data);

    var.Compress();

    std::vector<std::vector<uint8_t>> encoding;
    var.MoveEncodedDataInto(encoding);

    std::vector<std::vector<uint8_t>> mesh_encoding;
    var.MoveEncodedMeshInto(mesh_encoding);

    /* Serialize the data, as it would have been written out */
    std::vector<uint8_t> serialized_encoded_data;

    for(auto lvl_iter = encoding.rbegin(); lvl_iter != encoding.rend(); ++lvl_iter)
    {
        std::copy_n(lvl_iter->begin(), lvl_iter->size(), std::back_inserter(serialized_encoded_data));
    }    
    
    /* Serialize the encoded mesh, as it would have been written out */
    std::vector<uint8_t> serialized_mesh_encoded_data;

    for(auto lvl_iter = mesh_encoding.rbegin(); lvl_iter != mesh_encoding.rend(); ++lvl_iter)
    {
        cmc::cmc_debug_msg("Encoding of level: ", lvl_iter->size());
        std::copy_n(lvl_iter->begin(), lvl_iter->size(), std::back_inserter(serialized_mesh_encoded_data));
    } 

    /* Define the decompressor for the data stream */
    cmc::lossless::multi_res::DecompressionVariable<float> decompression_var("test_var", std::move(serialized_encoded_data), std::move(serialized_mesh_encoded_data));

    /* Decompress the data */
    decompression_var.Decompress();

    //cmc::cmc_debug_msg("Size of decompressed data: ", decompression_var.Size());
    cmc::ExpectTrue(decompression_var.Size() == data.size());

    std::vector<cmc::SerializedCompressionValue<sizeof(float)>> decompressed_byte_values = decompression_var.GetDecompressedData();

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
