#include "cmc.hxx"
#include "test/cmc_test.hxx"

#include "mpi/cmc_mpi.hxx"
#include "amr/lossless/cmc_par_multi_res_extraction_util.hxx"
#include "utilities/cmc_bits_stream_decoder.hxx"
#include "utilities/cmc_huffman_coder.hxx"

#include "amr/lossy/trixi/cmc_par_multi_res_extraction_idw.hxx"
#include "amr/lossy/trixi/cmc_par_multi_res_decompression_idw.hxx"

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
#include <filesystem>

namespace cmc::test
{

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

    cmesh = t8_cmesh_new_brick_2d (2, 2, 0, 0, MPI_COMM_WORLD);

    const t8_scheme *scheme = t8_scheme_new_default ();

    const int initial_level = 7;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);
    
    /* Partition For coarsening */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    const int partition_for_coarsening = 1;
    t8_forest_set_partition(partitioned_forest, forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    partitioned_forest = t8_forest_new_adapt(partitioned_forest, TestAdapt, 0, 0, NULL);

    return partitioned_forest;
}

std::pair<t8_cmesh_t, const t8_scheme *>
GetInitCMeshAndScheme()
{
    t8_cmesh_t cmesh = t8_cmesh_new_brick_2d (2, 2, 0, 0, MPI_COMM_WORLD);
    return std::make_pair(cmesh, t8_scheme_new_default());
}


std::vector<float>
CreateInitData(const t8_forest_t mesh)
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

[[maybe_unused]] void
WriteDataToVTK(t8_forest_t mesh, const std::vector<float>& data, const std::string& file_name)
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

    t8_forest_write_vtk_ext (mesh, file_name.c_str(), 1, 1, 1, 1, 0, 0, 0, 1, vtk_data);
}

std::pair<t8_forest_t, std::vector<float>>
RepartitionMeshAndDataToInitialPartition(t8_forest_t decompr_forest, std::vector<float>& decompressed_data_vec, const uint64_t init_mesh_offset)
{
    const t8_locidx_t current_local_elems = t8_forest_get_local_num_leaf_elements(decompr_forest);

    /* Keep the not-partitioned forest */
    t8_forest_ref(decompr_forest);

    /* Allocate a new forest */
    t8_forest_t partitioned_mesh;
    t8_forest_init (&partitioned_mesh);

    /* Set the forest for partitioning */
    t8_forest_set_partition (partitioned_mesh, decompr_forest, 0);

    /* Set the partition bound explicitly */
    t8_forest_set_partition_offset (partitioned_mesh, static_cast<t8_gloidx_t>(init_mesh_offset));

    /* Commit the forest, this step will perform the partitioning and ghost layer creation. */
    t8_forest_commit (partitioned_mesh);

    /* Get the number of local elements */
    const t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements(partitioned_mesh);
    
    /* Create an sc_array_t wrapper of the variable's local element data */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(decompressed_data_vec.data()), sizeof(float), current_local_elems);

    /* Allocate an output vector for the partitioned data */
    std::vector<float> partitioned_data(num_local_elements);

    /* Create a wrapper for the freshly allocated partitioned data */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(partitioned_data.data()), sizeof(float), num_local_elements);

    /* Partition the variables local element data */
    t8_forest_partition_data(decompr_forest, partitioned_mesh, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Free the former forest */
    t8_forest_unref(&decompr_forest);

    return std::make_pair(partitioned_mesh, std::move(partitioned_data));
}

}

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {
        /* Get the MPI rank */
        int rank{0};
        constexpr int kTestRootRank = 0;

        /* Get the rank within the initial communicator */
        const int rv_rank = MPI_Comm_rank(MPI_COMM_WORLD, &rank);
        if (rv_rank != MPI_SUCCESS)
        {
            cmc::cmc_err_msg("An error occured during the MPI rank request in MPI_COMM_WORLD!");
        }
        cmc::ExpectTrue(rv_rank == MPI_SUCCESS);

        constexpr int32_t DIM = 2;
        [[maybe_unused]] constexpr int32_t N = 1;
        const std::string file_name("cmc_out_test_trixi_multi_res_lossy_single_data_serial.cmc");

        std::vector<float> decompressed_data;

        /* Create an error domain */
        /* Set an absolute error criterion */
        const float kAbsMaxError = 0.05;
        cmc::ErrorDomain global_abs_error_domain(cmc::PermittedError(cmc::CompressionCriterion::AbsoluteErrorThreshold, kAbsMaxError), cmc::error_domain_fn::GeneralErrorCriterion);
        const float kRelMaxError = 0.01;
        cmc::ErrorDomain global_rel_error_domain(cmc::PermittedError(cmc::CompressionCriterion::RelativeErrorThreshold, kRelMaxError), cmc::error_domain_fn::GeneralErrorCriterion);
        
        std::vector<cmc::ErrorDomain> error_domains;
        error_domains.push_back(global_abs_error_domain);
        error_domains.push_back(global_rel_error_domain);

        /* Create a cmesh */
        t8_forest_t mesh = cmc::test::CreateInitMesh();

        /* Store the global number of elements */
        const uint64_t num_global_init_elems = t8_forest_get_global_num_leaf_elements(mesh);

        /* Store the initial data offset */
        const uint64_t init_mesh_offset = t8_forest_get_first_local_leaf_element_id(mesh);

        /* Create some data */
        std::vector<float> data = cmc::test::CreateInitData(mesh);

        /* Create a view on the data */
        const std::span<float> init_data(data);

        //cmc::test::WriteDataToVTK(mesh, data, "cmc_test_input_data_lossy_compr");

        //Compression
        {
        /* Perform the lossy compression with multiple data */
        cmc::par::lossy::multi_res::idw::trixi::CompressionVariable<float, DIM> variable("test_var", mesh, init_data, error_domains);

        /* Compress */
        variable.Compress();

        /* Write the encoded data */
        variable.WriteCompressedData(file_name);

        /* Deallocate the initial mesh */
        t8_forest_unref(&mesh);
        }

        //Decompression
        {
        cmc::par::lossy::multi_res::idw::trixi::ReadCompressionInfo(file_name);
        
        /* Get the initial base mesh */
        const auto [cmesh, scheme] = cmc::test::GetInitCMeshAndScheme();

        /* Decompress the mesh and the data */
        cmc::par::lossy::multi_res::idw::trixi::DecompressionVariable<float, DIM> decompression_variable(file_name, MPI_COMM_WORLD, cmesh, scheme);
        decompression_variable.Decompress();

        /* Access the decompressed components */
        auto [decompressed_forest, decompressed_data_vec] = decompression_variable.GetDecompressedData();
        
        /* Get the number of elements */
        const uint64_t num_global_decompr_elems = t8_forest_get_global_num_leaf_elements(decompressed_forest);

        cmc::ExpectTrue(num_global_decompr_elems == num_global_init_elems);

        /* Repartition the data accordingly to the initial data */
        auto [partitioned_decompressed_forest, partitioned_decompressed_data] =  cmc::test::RepartitionMeshAndDataToInitialPartition(decompressed_forest, decompressed_data_vec, init_mesh_offset);
        
        /* Save the decompressed data */
        decompressed_data = std::move(partitioned_decompressed_data);

        //cmc::test::WriteDataToVTK(partitioned_decompressed_forest, decompressed_data, "cmc_test_output_data_lossy_compr");

        /* Deallocate the decompressed forest */
        t8_forest_unref(&partitioned_decompressed_forest);
        }

        //Compare the initial and the decompressed data 
        int count_error{0};

        /* Compare the initial and decompressed data */
        for (int idx{0}; idx < data.size(); ++idx)
        {
            /* Compute the permitted minumum absolute deviation from the error criteria */
            const float rel_abs_dev = std::abs(data[idx] * kRelMaxError);
            const float abs_dev = kAbsMaxError;
            const float residual = std::abs(data[idx] - decompressed_data[idx]);

            if (residual > abs_dev || residual > rel_abs_dev)
            {
                //cmc::cmc_debug_msg("Inequality at ", idx, ", Init Data: ", data[idx], ", Decompressed Data: ", decompressed_data[idx], ", Residual: ", residual, ", Permitted Absolute Deviation ", abs_dev, ", Permitted Relative Deviation (in absolute terms): ", rel_abs_dev);
                //cmc::cmc_debug_msg("Bitset Init Data: , ", std::bitset<32>(std::bit_cast<uint32_t>(data[idx])),  ", Bitset Decompressed Data: ", std::bitset<32>(std::bit_cast<uint32_t>(decompressed_data[idx])));
                ++count_error;
            }

            cmc::ExpectTrue(residual <= abs_dev && residual <= rel_abs_dev);
        }
        cmc::cmc_debug_msg("Number of error domain violations (|init_data - decompressed_data| > permitted_error): ", count_error);
        cmc::ExpectTrue(count_error == 0);
        
        //Clean-Up
        /* Delete the test's output file */
        if (rank == kTestRootRank)
        {
            const std::filesystem::path output_file_path(file_name);
            if (std::filesystem::exists(output_file_path))
            {
                std::remove(output_file_path.c_str());
            }
        }
    
    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
