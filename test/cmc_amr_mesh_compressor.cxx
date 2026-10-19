#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "cmc_amr_mesh_compression.hxx"

#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_partition.h>

#include <string>

constexpr bool kWriteVTKDataTestAMRMeshCompression = false;

namespace cmc::test::mesh_compressor
{

static t8_locidx_t
TestAdapt ([[maybe_unused]] t8_forest_t forest,
           t8_forest_t forest_from,
           t8_locidx_t which_tree,
           [[maybe_unused]] const t8_eclass_t tree_class,
           [[maybe_unused]] t8_locidx_t lelement_id,
           [[maybe_unused]] const t8_scheme_c * ts,
           const int is_family,
           [[maybe_unused]] const int num_elements,
           [[maybe_unused]] t8_element_t * elements[])
{

    std::vector<double> elem_midpoint(3);
    t8_forest_element_centroid (forest_from, which_tree, elements[0], elem_midpoint.data());

    const t8_gloidx_t gtree_id = t8_forest_global_tree_id (forest_from, which_tree);
    if ((gtree_id == 3 || gtree_id == 7 || gtree_id == 10) && is_family) {return -1;}
    if (gtree_id == 4) {return 1;}
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
    t8_cmesh_init(&cmesh);
    t8_cmesh_new_hypercube_hybrid (cmesh, MPI_COMM_WORLD, 0);

    const t8_scheme *scheme = t8_scheme_new_default ();

    const int initial_level = 4;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);
    
    /* Partition For coarsening */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    constexpr int partition_for_coarsening = 1;
    t8_forest_set_partition(partitioned_forest, forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    partitioned_forest = t8_forest_new_adapt(partitioned_forest, TestAdapt, 0, 0, NULL);

    t8_forest_t adapted_partitioned_forest;
    t8_forest_init(&adapted_partitioned_forest);

    /* Partition the forest */
    constexpr int partition_for_coarsening2 = 0;
    t8_forest_set_partition(adapted_partitioned_forest, partitioned_forest, partition_for_coarsening2);
    t8_forest_commit(adapted_partitioned_forest);

    return adapted_partitioned_forest;
}

std::pair<t8_cmesh_t, const t8_scheme *>
GetInitCMeshAndScheme()
{
    t8_cmesh_t cmesh;
    t8_cmesh_init(&cmesh);
    t8_cmesh_new_hypercube_hybrid (cmesh, MPI_COMM_WORLD, 0);
    return std::make_pair(cmesh, t8_scheme_new_default());
}

}

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {
        /**** SETUP ****/
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

        const std::string file_name("cmc_out_test_amr_mesh_compressor.cmc");

        /* Create an intiial mesh */
        t8_forest_t mesh = cmc::test::mesh_compressor::CreateInitMesh();

        if constexpr (kWriteVTKDataTestAMRMeshCompression)
        {
            t8_forest_write_vtk_ext (mesh, "test_init_mesh", 1, 1, 1, 1, 0, 0, 0, 0, NULL);
        }

        /* Store the global number of elements */
        const uint64_t num_global_init_elems = static_cast<uint64_t>(t8_forest_get_global_num_leaf_elements(mesh));

        /* Store the initial data offset */
        const uint64_t init_mesh_offset = static_cast<uint64_t>(t8_forest_get_first_local_leaf_element_id(mesh));

        /* Declare a forest pointer for the reconstructed mesh */
        t8_forest_t reconstructed_mesh;

        /**** COMPRESSION ****/
        {
        /* Define the mesh compressor */
        cmc::mesh::Compressor mesh_compressor(mesh);

        /* Compress */
        mesh_compressor.CompressMeshStructure();

        /* Write the encoded data */
        mesh_compressor.WriteCompressedMeshStructure(file_name);

        }

        cmc::mesh::PrintMeshCompressionInfo(file_name);

        {
        /**** DECOMPRESSION ****/

        /* Get the initial cmesh and scheme */
        auto [cmesh, scheme] = cmc::test::mesh_compressor::GetInitCMeshAndScheme();

        /* Define the mesh decompressor */
        cmc::mesh::Decompressor mesh_decompressor(file_name, MPI_COMM_WORLD, cmesh, scheme);

        /* Decompress the mesh */
        reconstructed_mesh = mesh_decompressor.ReconstructMesh();

        const uint64_t num_global_reconstructed_elements = static_cast<uint64_t>(t8_forest_get_global_num_leaf_elements(reconstructed_mesh));
        cmc::cmc_debug_msg("The reconstructed mesh has ", num_global_reconstructed_elements, " elements.");

        cmc::ExpectTrue(num_global_reconstructed_elements == num_global_init_elems);

        }
    
        /* Partition the reconstructed mesh equal to the initial mesh */
        t8_forest_t reconstructed_init_partitioned_mesh;
        t8_forest_init (&reconstructed_init_partitioned_mesh);

        /* Set the forest for partitioning */
        t8_forest_set_partition (reconstructed_init_partitioned_mesh, reconstructed_mesh, 0);

        /* Set the partition bound explicitly */
        t8_forest_set_partition_offset (reconstructed_init_partitioned_mesh, static_cast<t8_gloidx_t>(init_mesh_offset));

        /* Commit the forest, this step will perform the partitioning */
        t8_forest_commit (reconstructed_init_partitioned_mesh);

        if constexpr (kWriteVTKDataTestAMRMeshCompression)
        {
            t8_forest_write_vtk_ext (reconstructed_init_partitioned_mesh, "test_reconstructed_mesh", 1, 1, 1, 1, 0, 0, 0, 0, NULL);
        }

        /**** COMPARISON OF INIT AND DECOMPRESED FOREST ****/
        /* Comapre both forests */
        const t8_locidx_t num_trees_init_mesh = t8_forest_get_num_local_trees (mesh);
        const t8_locidx_t num_trees_decompr_mesh = t8_forest_get_num_local_trees (reconstructed_init_partitioned_mesh);
        cmc::ExpectTrue(num_trees_init_mesh == num_trees_decompr_mesh);

        const t8_scheme* init_scheme = t8_forest_get_scheme (mesh);
        const t8_scheme* decompr_scheme = t8_forest_get_scheme (reconstructed_init_partitioned_mesh);
        
        for (t8_locidx_t itree = 0; itree < num_trees_init_mesh; ++itree)
        {
            const t8_eclass_t tree_class_init_mesh = t8_forest_get_tree_class (mesh, itree);
            const t8_eclass_t tree_class_decompr_mesh = t8_forest_get_tree_class (reconstructed_init_partitioned_mesh, itree);
            cmc::ExpectTrue(tree_class_init_mesh == tree_class_decompr_mesh);

            const t8_locidx_t num_elems_in_tree_init_mesh = t8_forest_get_tree_num_leaf_elements (mesh, itree);
            const t8_locidx_t num_elems_in_tree_decompr_mesh = t8_forest_get_tree_num_leaf_elements (reconstructed_init_partitioned_mesh, itree);
            cmc::ExpectTrue(num_elems_in_tree_init_mesh == num_elems_in_tree_decompr_mesh);

            for (t8_locidx_t ielement = 0; ielement < num_elems_in_tree_init_mesh; ++ielement)
            {
                const t8_element_t* init_element = t8_forest_get_leaf_element_in_tree (mesh, itree, ielement);
                const t8_element_t* decompr_element = t8_forest_get_leaf_element_in_tree (reconstructed_init_partitioned_mesh, itree, ielement);

                /* Get the element level */
                const int32_t init_elem_level = init_scheme->element_get_level(tree_class_init_mesh, init_element);
                const int32_t decompr_elem_level = decompr_scheme->element_get_level(tree_class_decompr_mesh, decompr_element);
                cmc::ExpectTrue(init_elem_level == decompr_elem_level);

                /* Check the element for equality */
                const int are_elems_equal = init_scheme->element_is_equal(tree_class_init_mesh, init_element, decompr_element);
                cmc::ExpectTrue(are_elems_equal != 0);
            }
        }

        /* Destroy all forests */
        t8_forest_unref(&mesh);
        t8_forest_unref(&reconstructed_init_partitioned_mesh);
        
        /* Delete the compressed test file */
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

    return 0;
}
