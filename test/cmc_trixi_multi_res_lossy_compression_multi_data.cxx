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

namespace cmc::test::lossy::multi_res
{

struct Point
{
    Point() = default;
    constexpr Point(const float x_, const float y_)
    : x{x_}, y{y_} {};

    float x{};
    float y{};
};

/* Distribution of Gauss-Lobatto poits order 4 (-1, -sqrt(1/5), +sqrt(1/5), +1) */
constexpr int kNumPointsPerDimension = 4;
constexpr int kNumPointsPerElem = 16;
constexpr float x0 = -1.0;
constexpr float x1 = -sqrt(1.0/5.0);
constexpr float x2 = +sqrt(1.0/5.0);
constexpr float x3 = +1.0;
constexpr float y0 = x0;
constexpr float y1 = x1;
constexpr float y2 = x2;
constexpr float y3 = x3;

constexpr float a = -1.0;
constexpr float b = 1.0;

constexpr std::array<Point, kNumPointsPerElem> eval_coords{Point(x0, y0), Point(x1, y0), Point(x2, y0), Point(x3, y0),
                                                           Point(x0, y1), Point(x1, y1), Point(x2, y1), Point(x3, y1),
                                                           Point(x0, y2), Point(x1, y2), Point(x2, y2), Point(x3, y2),
                                                           Point(x0, y3), Point(x1, y3), Point(x2, y3), Point(x3, y3)};

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

    cmesh = t8_cmesh_new_brick_2d (2, 2, 0, 0, MPI_COMM_WORLD);

    const t8_scheme *scheme = t8_scheme_new_default ();

    const int initial_level = 5;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);
    
    /* Partition For coarsening */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    const int partition_for_coarsening = 1; //TODO: change to 'one' when partition for coarsening is in t8code
    t8_forest_set_partition(partitioned_forest, forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    /* Adapt the uniform forest a little */
    //partitioned_forest = t8_forest_new_adapt(partitioned_forest, TestAdapt, 0, 0, NULL);

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

[[maybe_unused]] void
WriteDataToVTK(t8_forest_t mesh, const std::vector<float>& data, const std::string& file_name)
{
    std::vector<double> double_data1;
    std::vector<double> double_data2;
    std::vector<double> double_data3;
    std::vector<double> double_data4;
    std::vector<double> double_data5;
    std::vector<double> double_data6;
    std::vector<double> double_data7;
    std::vector<double> double_data8;
    std::vector<double> double_data9;
    std::vector<double> double_data10;
    std::vector<double> double_data11;
    std::vector<double> double_data12;
    std::vector<double> double_data13;
    std::vector<double> double_data14;
    std::vector<double> double_data15;
    std::vector<double> double_data16;

    for (int idx{0}; idx < t8_forest_get_local_num_leaf_elements(mesh); ++idx)
    {
        double_data1.push_back(data[idx * kNumPointsPerElem]);
        double_data2.push_back(data[idx * kNumPointsPerElem + 1]);
        double_data3.push_back(data[idx * kNumPointsPerElem + 2]);
        double_data4.push_back(data[idx * kNumPointsPerElem + 3]);
        double_data5.push_back(data[idx * kNumPointsPerElem + 4]);
        double_data6.push_back(data[idx * kNumPointsPerElem + 5]);
        double_data7.push_back(data[idx * kNumPointsPerElem + 6]);
        double_data8.push_back(data[idx * kNumPointsPerElem + 7]);
        double_data9.push_back(data[idx * kNumPointsPerElem + 8]);
        double_data10.push_back(data[idx * kNumPointsPerElem + 9]);
        double_data11.push_back(data[idx * kNumPointsPerElem + 10]);
        double_data12.push_back(data[idx * kNumPointsPerElem + 11]);
        double_data13.push_back(data[idx * kNumPointsPerElem + 12]);
        double_data14.push_back(data[idx * kNumPointsPerElem + 13]);
        double_data15.push_back(data[idx * kNumPointsPerElem + 14]);
        double_data16.push_back(data[idx * kNumPointsPerElem + 15]);
    }

    t8_vtk_data_field_t vtk_data[16];
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
    snprintf (vtk_data[4].description, BUFSIZ, "ExampleData5");
    vtk_data[4].type = T8_VTK_SCALAR;
    vtk_data[4].data = double_data5.data();
    snprintf (vtk_data[5].description, BUFSIZ, "ExampleData6");
    vtk_data[5].type = T8_VTK_SCALAR;
    vtk_data[5].data = double_data6.data();
    snprintf (vtk_data[6].description, BUFSIZ, "ExampleData7");
    vtk_data[6].type = T8_VTK_SCALAR;
    vtk_data[6].data = double_data7.data();
    snprintf (vtk_data[7].description, BUFSIZ, "ExampleData8");
    vtk_data[7].type = T8_VTK_SCALAR;
    vtk_data[7].data = double_data8.data();
    snprintf (vtk_data[8].description, BUFSIZ, "ExampleData9");
    vtk_data[8].type = T8_VTK_SCALAR;
    vtk_data[8].data = double_data9.data();
    snprintf (vtk_data[9].description, BUFSIZ, "ExampleData10");
    vtk_data[9].type = T8_VTK_SCALAR;
    vtk_data[9].data = double_data10.data();
    snprintf (vtk_data[10].description, BUFSIZ, "ExampleData11");
    vtk_data[10].type = T8_VTK_SCALAR;
    vtk_data[10].data = double_data11.data();
    snprintf (vtk_data[11].description, BUFSIZ, "ExampleData12");
    vtk_data[11].type = T8_VTK_SCALAR;
    vtk_data[11].data = double_data12.data();
    snprintf (vtk_data[12].description, BUFSIZ, "ExampleData13");
    vtk_data[12].type = T8_VTK_SCALAR;
    vtk_data[12].data = double_data13.data();
    snprintf (vtk_data[13].description, BUFSIZ, "ExampleData14");
    vtk_data[13].type = T8_VTK_SCALAR;
    vtk_data[13].data = double_data14.data();
    snprintf (vtk_data[14].description, BUFSIZ, "ExampleData15");
    vtk_data[14].type = T8_VTK_SCALAR;
    vtk_data[14].data = double_data15.data();
    snprintf (vtk_data[15].description, BUFSIZ, "ExampleData16");
    vtk_data[15].type = T8_VTK_SCALAR;
    vtk_data[15].data = double_data16.data();
    t8_forest_write_vtk_ext (mesh, file_name.c_str(), 1, 1, 1, 1, 0, 0, 0, 16, vtk_data);
}

std::pair<t8_forest_t, std::vector<float>>
RepartitionMeshAndDataToInitialPartition(t8_forest_t decompr_forest, std::vector<float>& decompressed_data_vec, const uint64_t init_mesh_offset)
{
    const t8_locidx_t current_local_elems = t8_forest_get_local_num_leaf_elements(decompr_forest);

    cmc::ExpectTrue(current_local_elems * cmc::test::lossy::multi_res::kNumPointsPerElem == decompressed_data_vec.size());

    /* Keep the not-partitioned forest */
    t8_forest_ref(decompr_forest);

    /** Partition the forest correctly and build a halo layer **/
    t8_forest_t partitioned_mesh;
    t8_forest_init (&partitioned_mesh);

    /* Set the forest for partitioning */
    t8_forest_set_partition (partitioned_mesh, decompr_forest, 0);

    /* Set the partition bound explicitly */
    t8_forest_set_partition_offset (partitioned_mesh, static_cast<t8_gloidx_t>(init_mesh_offset));

    /* Commit the forest, this step will perform the partitioning and ghost layer creation. */
    t8_forest_commit (partitioned_mesh);

    /** Exchange the data corectly to set up the halo layer **/
    /* Get the number of local elements */
    const t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements(partitioned_mesh);
    
    std::vector<float> init_component(current_local_elems);

    /* Allocate an output vector for the partitioned data */
    std::vector<float> partitioned_component(num_local_elements);

    /* Allocate an output vector for the partitioned data */
    std::vector<float> output_data(num_local_elements * kNumPointsPerElem);

    /* Loop over all components */
    for (int idx{0}; idx < kNumPointsPerElem; ++idx)
    {
        /* Fill input vector */
        for (int elem_idx{0}; elem_idx < current_local_elems; ++elem_idx)
        {
            init_component[elem_idx] = decompressed_data_vec[elem_idx * kNumPointsPerElem + idx];
        }

        const int rv_wait2 = MPI_Barrier(MPI_COMM_WORLD);
        if (rv_wait2 != MPI_SUCCESS)
        {
            cmc_err_msg("Error in MPI_Barrier during partitioning of the decompressed components!");
        }

        /* Create an sc_array_t wrapper of the variable's local element data */
        sc_array_t* in_data = sc_array_new_data (static_cast<void*>(init_component.data()), sizeof(float), current_local_elems);

        /* Create a wrapper for the freshly allocated partitioned data */
        sc_array_t* out_data = sc_array_new_data (static_cast<void*>(partitioned_component.data()), sizeof(float), num_local_elements);

        /* Partition the variables local element data */
        t8_forest_partition_data(decompr_forest, partitioned_mesh, in_data, out_data);

        /* Assign the components correctly to the output vector */
        for (int elem_idx{0}; elem_idx < num_local_elements; ++elem_idx)
        {
            output_data[elem_idx * kNumPointsPerElem + idx] = partitioned_component[elem_idx];
        }

        const int rv_wait = MPI_Barrier(MPI_COMM_WORLD);
        if (rv_wait != MPI_SUCCESS)
        {
            cmc_err_msg("Error in MPI_Barrier during partitioning of the decompressed components!");
        }

        /* Destroy the array wrappers */
        sc_array_destroy(in_data);
        sc_array_destroy(out_data);
    }

    /* Free the former forest */
    t8_forest_unref(&decompr_forest);

    return std::make_pair(partitioned_mesh, std::move(output_data));
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
        constexpr int32_t N = cmc::test::lossy::multi_res::kNumPointsPerDimension;
        const std::string file_name("cmc_out_test_trixi_multi_res_lossy_multi_data_serial.cmc");

        std::vector<float> decompressed_data;

        /* Create an error domain */
        const float kAbsMaxError = 0.025;
        cmc::ErrorDomain global_abs_error_domain(cmc::PermittedError(cmc::CompressionCriterion::AbsoluteErrorThreshold, kAbsMaxError), cmc::error_domain_fn::GeneralErrorCriterion);
        const float kRelMaxError = 0.01;
        cmc::ErrorDomain global_rel_error_domain(cmc::PermittedError(cmc::CompressionCriterion::RelativeErrorThreshold, kRelMaxError), cmc::error_domain_fn::GeneralErrorCriterion);
        
        std::vector<cmc::ErrorDomain> error_domains;
        error_domains.push_back(global_abs_error_domain);
        error_domains.push_back(global_rel_error_domain);

        /* Create a cmesh */
        t8_forest_t mesh = cmc::test::lossy::multi_res::CreateInitMesh();

        /* Store the global number of elements */
        const uint64_t num_global_init_elems = t8_forest_get_global_num_leaf_elements(mesh);

        /* Store the initial data offset */
        const uint64_t init_mesh_offset = t8_forest_get_first_local_leaf_element_id(mesh);

        /* Create some data */
        std::vector<float> data = cmc::test::lossy::multi_res::CreateInitData(mesh);

        /* Create a view on the data */
        const std::span<float> init_data(data);

        /* Write the data to a vtk file */
        //cmc::test::lossy::multi_res::WriteDataToVTK(mesh, data, "cmc_test_lossy_par_multi_res_serial_multi_data_input");

        //Compression
        {
        /* Supply the reference coordinates of the element */
        std::vector<cmc::par::lossy::idw::util::Coordinate<DIM>> ref_coords{cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x0, cmc::test::lossy::multi_res::y0}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x1, cmc::test::lossy::multi_res::y0}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x2, cmc::test::lossy::multi_res::y0}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x3, cmc::test::lossy::multi_res::y0}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x0, cmc::test::lossy::multi_res::y1}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x1, cmc::test::lossy::multi_res::y1}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x2, cmc::test::lossy::multi_res::y1}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x3, cmc::test::lossy::multi_res::y1}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x0, cmc::test::lossy::multi_res::y2}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x1, cmc::test::lossy::multi_res::y2}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x2, cmc::test::lossy::multi_res::y2}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x3, cmc::test::lossy::multi_res::y2}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x0, cmc::test::lossy::multi_res::y3}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x1, cmc::test::lossy::multi_res::y3}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x2, cmc::test::lossy::multi_res::y3}),
                                                                            cmc::par::lossy::idw::util::Coordinate<DIM>(std::array<float, DIM>{cmc::test::lossy::multi_res::x3, cmc::test::lossy::multi_res::y3}),
                                                                        };

        /* Perform the lossy compression with multiple data */
        cmc::par::lossy::multi_res::idw::trixi::CompressionVariableMultiData<float, DIM, N> variable("test_var", mesh, init_data, error_domains, ref_coords);

        /* Compress */
        variable.Compress();

        /* Write the encoded data */
        variable.WriteCompressedData(file_name);

        }

        /* Deallocate the mesh */
        t8_forest_unref(&mesh);
        mesh = nullptr;

        const int rv_wait = MPI_Barrier(MPI_COMM_WORLD);
        if (rv_wait != MPI_SUCCESS)
        {
            cmc::cmc_err_msg("Error in MPI_Barrier during partitioning of the decompressed components!");
        }


        //Decompression
        {
        cmc::par::lossy::multi_res::idw::trixi::ReadCompressionInfo(file_name);
        
        /* Decompress */
        const auto [cmesh, scheme] = cmc::test::lossy::multi_res::GetInitCMeshAndScheme();

        cmc::par::lossy::multi_res::idw::trixi::DecompressionVariableMultiData<float, DIM, N> decompression_variable(file_name, MPI_COMM_WORLD, cmesh, scheme);
        decompression_variable.Decompress();

        auto [decompressed_forest, decompressed_data_vec] = decompression_variable.GetDecompressedData();

        cmc::cmc_debug_msg("Size of decompressed data vec: ", decompressed_data_vec.size());
        /* Get the number of elements */
        const uint64_t num_global_decompr_elems = t8_forest_get_global_num_leaf_elements(decompressed_forest);

        cmc::ExpectTrue(num_global_decompr_elems == num_global_init_elems);

        /* Repartition the data accordingly to the initial data */
        auto [partitioned_decompressed_forest, partitioned_decompressed_data] =  cmc::test::lossy::multi_res::RepartitionMeshAndDataToInitialPartition(decompressed_forest, decompressed_data_vec, init_mesh_offset);
        
        /* Save the decompressed data */
        decompressed_data = std::move(partitioned_decompressed_data);

        /* Write the data to a vtk file */
        //cmc::test::lossy::multi_res::WriteDataToVTK(decompressed_forest, decompressed_data_vec, "cmc_test_lossy_par_multi_res_serial_multi_data_output");

        /* Deallocate the decompressed mesh */
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
