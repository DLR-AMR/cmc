#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "lossless/cmc_prefix_extraction_compression.hxx"
#include "lossless/cmc_prefix_extraction_decompression.hxx"

#include <numeric>
#include <algorithm>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {

    /* Create a mesh */
    const sc_MPI_Comm comm = sc_MPI_COMM_SELF;
    t8_cmesh_t cmesh;
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
    t8_scheme_cxx_t* scheme = t8_scheme_new_default_cxx ();
    const int initial_level = 3;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);
    t8_eclass_scheme_c* ts = t8_forest_get_eclass_scheme(forest, t8_forest_get_eclass(forest, 0));

    /* Create a variable for containing a value for each element */
    const t8_locidx_t num_elements = t8_forest_get_local_num_elements(forest);
    std::vector<float> data(num_elements);
    std::iota(data.begin(), data.end(), 0.0);

    cmc::lossless::prefix::CompressionVariable<float> var("test_var", forest, data);

    var.Compress();

    std::vector<std::vector<uint8_t>> encoding;
    var.MoveEncodedDataInto(encoding);

    std::vector<std::vector<uint8_t>> mesh_encoding;
    var.MoveEncodedMeshInto(mesh_encoding);

    cmc::cmc_debug_msg("Levelwise encoding: levels: ", encoding.size());

    /* Serialize the data, as it would have been written out */
    std::vector<uint8_t> serialized_encoded_data;

    for(auto lvl_iter = encoding.rbegin(); lvl_iter != encoding.rend(); ++lvl_iter)
    {
        cmc::cmc_debug_msg("Encoding of level: ", lvl_iter->size());
        std::copy_n(lvl_iter->begin(), lvl_iter->size(), std::back_inserter(serialized_encoded_data));
    }    
    
    cmc::cmc_debug_msg("Size of serialized encoded data: ", serialized_encoded_data.size());

    /* Serialize the encoded mesh, as it would have been written out */
    std::vector<uint8_t> serialized_mesh_encoded_data;

    for(auto lvl_iter = mesh_encoding.rbegin(); lvl_iter != mesh_encoding.rend(); ++lvl_iter)
    {
        cmc::cmc_debug_msg("Encoding of level: ", lvl_iter->size());
        std::copy_n(lvl_iter->begin(), lvl_iter->size(), std::back_inserter(serialized_mesh_encoded_data));
    } 

    cmc::cmc_debug_msg("Size of serialized mesh data: ", serialized_mesh_encoded_data.size());


    /* Define the decompressor for the data stream */
    cmc::lossless::prefix::DecompressionVariable<float> decompression_var("test_var", std::move(serialized_encoded_data), std::move(serialized_mesh_encoded_data));

    /* Decompress the data */
    decompression_var.Decompress();

    cmc::cmc_debug_msg("Size of decompressed data: ", decompression_var.Size());

    std::vector<cmc::SerializedCompressionValue<sizeof(float)>> decompressed_byte_values = decompression_var.GetDecompressedData();

    int idx = 0;
    for (auto cr_iter = decompressed_byte_values.begin(); cr_iter != decompressed_byte_values.end(); ++cr_iter, ++idx)
    {
        const float val = cr_iter->ReinterpretDataAs<float>();
        cmc::cmc_debug_msg("Decompressed value is: ", val, " and initial value is: ", data[idx], " count of signif bits: ", cr_iter->GetCountOfSignificantBits());
    }

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
