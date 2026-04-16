#include "cmc.hxx"
#include "test/cmc_test.hxx"

#include "mpi/cmc_mpi.hxx"
#include "amr/lossless/cmc_par_multi_res_extraction_util.hxx"
#include "utilities/cmc_bits_stream_decoder.hxx"
#include "utilities/cmc_huffman_coder.hxx"
#include "amr/lossless/cmc_par_multi_res_extraction.hxx"

#include "amr/lossless/cmc_par_multi_res_decompression.hxx"

#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_forest/t8_forest_general.h>

#include <array>
#include <vector>
#include <span>
#include <numeric>
#include <algorithm>
#include <memory>
#include <limits>
#include <cmath>

void
Test2DFloatData()
{
    constexpr int32_t N = 16;
    constexpr int DIM = 2;

    /* We define some example data on one element */
    std::array<float, N> init_data{2.34,2.79,3.10,2.56,
                                   3.21,3.07,2.99,2.76,
                                   3.01,2.88,2.34,2.19,
                                   2.27,2.45,2.76,2.56};

    /* Compute the element encoding */
    const cmc::par::lossless::multi_res::IntraElementCoding<float, DIM, N> elem_coding = cmc::par::lossless::multi_res::ComputeElementEncoding<float, DIM, N>(std::span<float>(init_data));
    
    /* Get the coarse level predictor */
    const float coarse_predictor = elem_coding.GetCoarseValuePredictor();

    /* Setup of the Huffman Coder */
    ///////////////////////////////
    constexpr int num_entropy_symbols = cmc::par::lossless::multi_res::GetNumEntropySymbols<float>();
    std::vector<uint64_t> entropy_symbol_frequencies(num_entropy_symbols, 0);
    
    /* Iterate over the entropy codes and collect the frequencies */
    for (const auto& symbol : elem_coding.entropy_codes)
    {
        /* Update the frequency */
        ++entropy_symbol_frequencies[cmc::par::lossless::multi_res::MapEntropySymbolToArrayIndex<float>(symbol)];
    }

    std::vector<cmc::entropy_coding::huffman::EntropySymbol<cmc::par::lossless::multi_res::SymbolType>> global_symbol_frequencies;
    global_symbol_frequencies.reserve(num_entropy_symbols);

    for (int idx{0}; idx < num_entropy_symbols; ++idx)
    {
        /* Convert the index back to the entropy symbol */
        const cmc::par::lossless::multi_res::SymbolType entropy_symbol = cmc::par::lossless::multi_res::MapArrayIndexToEntropySymbol<float>(idx);

        /* Store the symbol with the global frequency */
        global_symbol_frequencies.emplace_back(entropy_symbol, entropy_symbol_frequencies[idx]);
    }

    /* We initialize the Huffman coder */
    cmc::entropy_coding::huffman::HuffmanCoder<cmc::par::lossless::multi_res::SymbolType> huffman_coder(global_symbol_frequencies);
    
    /* We serialize the Huffman coder */
    const std::vector<uint8_t> serialized_huffman_codes = huffman_coder.SerializeHuffmanCodes();
    ////////////////////////////////////
    /* End of the Huffman Coder setup */

    /* We encode the data */
    cmc::bits::vector elem_encoding;
    cmc::par::lossless::multi_res::PerformElementEncoding<float, DIM, N>(elem_encoding, huffman_coder, elem_coding);

    /* Serialize the stream of the encoded element data */
    const std::vector<uint64_t> encoded_elem_stream = elem_encoding.GetSerializedByteStreamBE();

    /* We setup a stream decoder for the data */
    /* Construct a stream decoder for the given serialized Huffman codes */
    cmc::bits::StreamDecoder<cmc::par::lossless::multi_res::SymbolType> stream_decoder;
    stream_decoder.StartHuffmanCodesDecoding(serialized_huffman_codes.data());

    /* Set the start of the decoder to the encoded stream */
    cmc::bits::vector_view encoded_stream_view(encoded_elem_stream.data());
    stream_decoder.StartDecoding(encoded_stream_view);

    /* We decompress the encoded data */
    std::array<float, N> decompressed_data = cmc::par::lossless::multi_res::PerformElementDecoding<float, DIM, N>(stream_decoder, coarse_predictor);

    /* We compare the initial and decompressed data for equality */
    for (int idx{0}; idx < N; ++idx)
    {
        cmc::ExpectTrue(init_data[idx] == decompressed_data[idx]); 
    }
}

//Lobatto Puntke auf [-1, 1] sind approx -0.57735026919, 0.57735026919
struct Point
{
    Point() = default;
    constexpr Point(const float x_, const float y_)
    : x{x_}, y{y_} {};

    float x{};
    float y{};
};

constexpr int kNumPointsPerElem = 4;
constexpr float p1 = -0.57735026919;
constexpr float p2 = 0.57735026919;
constexpr float a = -1.0;
constexpr float b = 1.0;

constexpr std::array<Point, kNumPointsPerElem> eval_coords{Point(p1, p1), Point(p2, p1), Point(p1, p2), Point(p2, p2)};

std::array<float, kNumPointsPerElem>
EvalPoints(t8_forest_t mesh, const int tree_id, const t8_element_t* elem)
{
    constexpr int face_id = 0;
    const float elem_face_length = static_cast<float>(t8_forest_element_face_area(mesh, tree_id, elem, face_id));

    constexpr int vertex_id = 0;
    std::array<double, 3> elem_coords{};
    t8_forest_element_coordinate (mesh, tree_id, elem, vertex_id, elem_coords.data());

    /* Transform reference coords */
    std::array<Point, kNumPointsPerElem> transformed_coords;
    for (int idx{0}; idx < kNumPointsPerElem; ++idx)
    {
        transformed_coords[idx] = Point(elem_coords[0] + elem_face_length * (std::fabs(a - eval_coords[idx].x) / std::fabs(b - a)),
                                        elem_coords[1] + elem_face_length * (std::fabs(a - eval_coords[idx].y) / std::fabs(b - a)));
    }

    /* Evalute the example function and return the values */
    std::array<float, kNumPointsPerElem> elem_values{};
    for (int idx{0}; idx < kNumPointsPerElem; ++idx)
    {
        elem_values[idx] = std::cos(transformed_coords[idx].x - 0.5) * std::cos(transformed_coords[idx].y - 0.5);
    }
    
    return elem_values;
}

static t8_locidx_t
TestAdapt ([[maybe_unused]] t8_forest_t forest,
           t8_forest_t forest_from,
           t8_locidx_t which_tree,
           [[maybe_unused]] const t8_eclass_t tree_class,
           t8_locidx_t lelement_id,
           [[maybe_unused]] const t8_scheme_c * ts,
           const int is_family,
           [[maybe_unused]] const int num_elements,
           [[maybe_unused]] t8_element_t * elements[])
{

    std::vector<double> elem_midpoint(3);
    t8_forest_element_centroid (forest_from, which_tree, elements[0], elem_midpoint.data());

    const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest_from, which_tree);
    if (gtree_id == 3 && is_family) {return -1;}
    if (std::fabs(elem_midpoint[0] - 0.5) <= 0.2 && std::fabs(elem_midpoint[1] - 0.5) <= 0.2)
    {
        return 1;
    } else
    {
        return 0;
    }
}

t8_forest_t
CreateInitMesh()
{
    /* Create a mesh */
    const sc_MPI_Comm comm = MPI_COMM_WORLD;
    t8_cmesh_t cmesh;
    #if 0
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
    #else
    cmesh = t8_cmesh_new_brick_2d (2, 2, 0, 0, MPI_COMM_WORLD);
    #endif
    const t8_scheme *scheme = t8_scheme_new_default ();
    //const int initial_level = 5;
    const int initial_level = 2;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);
    
    /* Partition For coarsening */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    const int partition_for_coarsening = 1; //TODO: change to 'one' when partition for coarsening is in t8code
    t8_forest_set_partition(partitioned_forest, forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    partitioned_forest = t8_forest_new_adapt(partitioned_forest, TestAdapt, 0, 0, NULL);

    return partitioned_forest;
}

std::pair<t8_cmesh_t, const t8_scheme *>
GetInitCMeshAndScheme()
{
    #if 0
    return std::make_pair(t8_cmesh_new_hypercube (T8_ECLASS_QUAD, MPI_COMM_WORLD, 0, 0, 0), t8_scheme_new_default());
    #else
    t8_cmesh_t cmesh = t8_cmesh_new_brick_2d (2, 2, 0, 0, MPI_COMM_WORLD);
    return std::make_pair(cmesh, t8_scheme_new_default());
    #endif
}

std::vector<float>
CreateInitData(const t8_forest_t mesh)
{
    const t8_locidx_t num_elems = t8_forest_get_local_num_leaf_elements(mesh);

    std::vector<float> data;
    data.reserve(num_elems * kNumPointsPerElem);

    const t8_scheme *scheme = t8_forest_get_scheme (mesh);

    t8_locidx_t num_trees = t8_forest_get_num_local_trees (mesh);

    for (t8_locidx_t itree = 0, idata = 0; itree < num_trees; itree++)
    {
        const t8_eclass_t tree_class = t8_forest_get_tree_class (mesh, itree);
        t8_locidx_t num_elems_in_tree = t8_forest_get_tree_num_leaf_elements (mesh, itree);

        for (t8_locidx_t ielement = 0; ielement < num_elems_in_tree; ielement++)
        {
            const t8_element_t *element = t8_forest_get_leaf_element_in_tree (mesh, itree, ielement);

            /* Evaluate the data for this element */
            std::array<float, kNumPointsPerElem> elem_data = EvalPoints(mesh, itree, element);

            /* Store the element data */
            std::copy_n(elem_data.begin(), kNumPointsPerElem, std::back_inserter(data));
        }
    }

    return data;
}

std::vector<float>
CreateInitData2(const t8_forest_t mesh)
{
    const t8_locidx_t num_elems = t8_forest_get_local_num_leaf_elements(mesh);

    std::vector<float> data;
    data.reserve(num_elems);

    const t8_scheme *scheme = t8_forest_get_scheme (mesh);

    t8_locidx_t num_trees = t8_forest_get_num_local_trees (mesh);

    for (t8_locidx_t itree = 0, idata = 0; itree < num_trees; itree++)
    {
        const t8_eclass_t tree_class = t8_forest_get_tree_class (mesh, itree);
        t8_locidx_t num_elems_in_tree = t8_forest_get_tree_num_leaf_elements (mesh, itree);

        for (t8_locidx_t ielement = 0; ielement < num_elems_in_tree; ielement++)
        {
            const t8_element_t *element = t8_forest_get_leaf_element_in_tree (mesh, itree, ielement);

            std::array<double, 3> elem_coords{};
            t8_forest_element_coordinate (mesh, itree, element, 0, elem_coords.data());

            data.push_back(std::cos(elem_coords[0] - 0.5) * std::cos(elem_coords[1] - 0.5));
        }
    }

    return data;
}

void
WriteDataToVTK(t8_forest_t mesh, const std::vector<float>& data)
{
    std::vector<double> double_data1;
    std::vector<double> double_data2;
    std::vector<double> double_data3;
    std::vector<double> double_data4;
    for (int idx{0}; idx < t8_forest_get_local_num_leaf_elements(mesh); ++idx)
    {
        double_data1.push_back(data[idx * kNumPointsPerElem]);
        double_data2.push_back(data[idx * kNumPointsPerElem + 1]);
        double_data3.push_back(data[idx * kNumPointsPerElem + 2]);
        double_data4.push_back(data[idx * kNumPointsPerElem + 3]);
    }
    cmc::cmc_global_msg("Num elems mesh: ", t8_forest_get_local_num_leaf_elements(mesh), " Num Points per elem: ", kNumPointsPerElem, ", Num Data: ", data.size());

    t8_vtk_data_field_t vtk_data[4];
    snprintf (vtk_data[0].description, BUFSIZ, "ExampleData1");
    vtk_data[0].type = T8_VTK_SCALAR;
    vtk_data[0].data = double_data1.data();
    snprintf (vtk_data[1].description, BUFSIZ, "ExampleData2");
    vtk_data[1].type = T8_VTK_SCALAR;
    vtk_data[1].data = double_data2.data();
    snprintf (vtk_data[2].description, BUFSIZ, "ExampleData3");
    vtk_data[2].type = T8_VTK_SCALAR;
    vtk_data[2].data = double_data3.data();
    snprintf (vtk_data[3].description, BUFSIZ, "ExampleData4");
    vtk_data[3].type = T8_VTK_SCALAR;
    vtk_data[3].data = double_data4.data();

    t8_forest_write_vtk_ext (mesh, "cmc_new_new_test_data_vis", 1, 1, 1, 1, 1, 1, 1, 4, vtk_data);
}

void
WriteDataToVTK2(t8_forest_t mesh, const std::vector<float>& data)
{
    std::vector<double> double_data1;
    for (int idx{0}; idx < t8_forest_get_local_num_leaf_elements(mesh); ++idx)
    {
        double_data1.push_back(data[idx]);
    }

    t8_vtk_data_field_t vtk_data[1];
    snprintf (vtk_data[0].description, BUFSIZ, "ExampleData1");
    vtk_data[0].type = T8_VTK_SCALAR;
    vtk_data[0].data = double_data1.data();

    t8_forest_write_vtk_ext (mesh, "cmc_new_new_test_data_vis", 1, 1, 1, 1, 1, 1, 1, 1, vtk_data);
}


static t8_locidx_t
CoarsenFourthTree ([[maybe_unused]] t8_forest_t forest,
                    t8_forest_t forest_from,
                    t8_locidx_t which_tree,
                    [[maybe_unused]] const t8_eclass_t tree_class,
                    [[maybe_unused]] t8_locidx_t lelement_id,
                    [[maybe_unused]] const t8_scheme_c * ts,
                    [[maybe_unused]] const int is_family,
                    [[maybe_unused]] const int num_elements,
                    [[maybe_unused]] t8_element_t * elements[])
{
    const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest_from, which_tree);

    if (gtree_id == 3)
    {
        return -1;
    } else
    {
        return 0;
    }
}

static t8_locidx_t
CoarsenAllFamilies ([[maybe_unused]] t8_forest_t forest,
                    [[maybe_unused]]t8_forest_t forest_from,
                    [[maybe_unused]] t8_locidx_t which_tree,
                    [[maybe_unused]] const t8_eclass_t tree_class,
                    [[maybe_unused]] t8_locidx_t lelement_id,
                    [[maybe_unused]] const t8_scheme_c * ts,
                    const int is_family,
                    [[maybe_unused]] const int num_elements,
                    [[maybe_unused]] t8_element_t * elements[])
{
    if (is_family)
    {
        return -1;
    } else
    {
        return 0;
    }
}

void
TestPFC()
{
    const MPI_Comm comm = MPI_COMM_WORLD;
    t8_cmesh_t cmesh = t8_cmesh_new_brick_2d (2, 2, 0, 0, comm);
    const t8_scheme* scheme = t8_scheme_new_default ();

    #if 0
    const int initial_level = 2;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);

    t8_forest_t adapted_forest = t8_forest_new_adapt(forest, CoarsenFourthTree, 0, 0, NULL);

    t8_forest_write_vtk_ext (adapted_forest, "test_pfc_mesh_adapted", 1, 1, 1, 1, 1, 1, 1, 0, NULL);

    t8_forest_t pfc_forest; 
    t8_forest_init(&pfc_forest);

    /* Partition the forest */
    const int partition_for_coarsening = 1; //TODO: change to 'one' when partition for coarsening is in t8code
    t8_forest_set_partition(pfc_forest, adapted_forest, partition_for_coarsening);
    t8_forest_commit(pfc_forest);

    t8_forest_write_vtk_ext (pfc_forest, "test_pfc_mesh_partitioned", 1, 1, 1, 1, 1, 1, 1, 0, NULL);
    #else
    
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, 0, 0, comm);
    t8_forest_t pfc_forest; 
    t8_forest_init(&pfc_forest);

    /* Partition the forest */
    const int partition_for_coarsening = 1; //TODO: change to 'one' when partition for coarsening is in t8code
    t8_forest_set_partition(pfc_forest, forest, partition_for_coarsening);
    t8_forest_commit(pfc_forest);

    #endif

}

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {

        //TestPFC();
        //cmc::cmc_global_msg("ENDDDDDDDDD");
        //cmc::cmc_err_msg("End of test"); 
        /* Setup a compression variable */
        constexpr int32_t DIM = 2;
        constexpr int32_t N = kNumPointsPerElem;
        const std::string file_name("test_par_multi_res_serial.cmc");

        std::vector<float> decompressed_data;

#if 1
        /* Create a cmesh */
        t8_forest_t mesh = CreateInitMesh();

        /* Create some data */
        //std::vector<float> data = CreateInitData(mesh);
        std::vector<float> data = CreateInitData2(mesh);

        /* Create a view on the data */
        const std::span<float> init_data(data);

        //WriteDataToVTK(mesh, data);
        WriteDataToVTK2(mesh, data);
        //cmc::cmc_err_msg("Stop here");
        {
        //cmc::par::lossless::multi_res::CompressionVariableMultiData<float, DIM, N> variable("test_var", mesh, init_data);
        cmc::par::lossless::multi_res::CompressionVariable<float, DIM> variable("test_var", mesh, init_data);
        
        /* Compress */
        variable.Compress(); 

        /* Write the encoded data */
        variable.WriteCompressedData(file_name);
        cmc::cmc_global_msg("Num elems mesh: ", t8_forest_get_local_num_leaf_elements(mesh));
        t8_forest_unref(&mesh);
        }
#endif
#if 0
        {
        cmc::par::lossless::multi_res::ReadCompressionInfo(file_name);
        
        /* Decompress */
        const auto [cmesh, scheme] = GetInitCMeshAndScheme();

        //cmc::par::lossless::multi_res::DecompressionVariableMultiData<float, DIM, N> decompression_variable(file_name, MPI_COMM_WORLD, cmesh, scheme);
        cmc::par::lossless::multi_res::DecompressionVariable<float, DIM> decompression_variable(file_name, MPI_COMM_WORLD, cmesh, scheme);
        decompression_variable.Decompress();

        auto [decompr_forest, decompressed_data_vec] = decompression_variable.GetDecompressedData(); 
        decompressed_data = decompressed_data_vec;

        //Write the data 
        WriteDataToVTK2(decompr_forest, decompressed_data);


        t8_forest_unref(&decompr_forest);
        }

#endif
#if 0
        /* Compare the initial and decompressed data */
        for (int idx{0}; idx < init_data.size(); ++idx)
        {
            if (init_data[idx] != decompressed_data[idx])
            {
                cmc::cmc_global_msg("Inequality at ", idx, ", Init: ", init_data[idx], ", Decompr: ", decompressed_data[idx]);
            }
        }

#endif

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
