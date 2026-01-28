#ifndef CMC_T8_CONSTRUCT_MESH_AND_DATA_FOR_COMPRESSION_HXX
#define CMC_T8_CONSTRUCT_MESH_AND_DATA_FOR_COMPRESSION_HXX

#include "cmc.hxx"

#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh_readmshfile.h>
#include "t8code/cmc_t8_mesh.hxx"
#include "mpi/cmc_mpi.hxx"

#include "amr/lossless/cmc_multi_res_par_extraction_compression.hxx"

#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "compression_io/cmc_compression_mpi_io_output.hxx"
#include "compression_io/cmc_decompression_mpi_io_input.hxx"

#include "amr/lossless/cmc_multi_res_extraction_compression.hxx"
#include "compression_io/cmc_compression_nc_output.hxx"

#include <vector>
#include <cmath>
#include <numeric>
#include <utility>
#include <cmath>
#include <fstream>
#include <cstdlib>


namespace cmc::example
{

float
ExampleEvalFn(const std::vector<double> &x)
{
    return static_cast<float>(std::sin(x[0]) * std::sin(x[1]));
}

std::vector<float>
CreateExampleData(t8_forest_t forest)
{
    t8_locidx_t num_local_elems = t8_forest_get_local_num_leaf_elements(forest);

    std::vector<float> values;
    values.reserve(num_local_elems);

    const t8_scheme *scheme = t8_forest_get_scheme (forest);

    t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest);

    for (t8_locidx_t itree = 0, idata = 0; itree < num_trees; itree++)
    {
        const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
        t8_locidx_t num_elems_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);

        for (t8_locidx_t ielement = 0; ielement < num_elems_in_tree; ielement++)
        {
            const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);

            const int level = scheme->element_get_level (tree_class, element);

            std::vector<double> midpoint(3);
                
            t8_forest_element_centroid (forest, itree, element, midpoint.data());

            values.push_back(ExampleEvalFn(midpoint));
        }
    }

    if (values.size() != num_local_elems)
    {
        cmc::cmc_err_msg("The amount of local data does not match the amount of local elements!");
    }

    return values;
}
 
std::pair<t8_forest_t, std::vector<float>>
BuildInitialForest()
{
    /* Construct an example mesh */
    const sc_MPI_Comm comm = MPI_COMM_WORLD;
    t8_cmesh_t cmesh;
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_TRIANGLE, comm, 0, 0, 0);
    const t8_scheme *scheme = t8_scheme_new_default ();
    const int initial_level = 9;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);

    /* Partition the forest */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    const int partition_for_coarsening = 0; //TODO: change to 'one' when partition for coarsening is in t8code
    t8_forest_set_partition(partitioned_forest, forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    /* Create some example data for tests */
    std::vector<float> example_data = CreateExampleData(partitioned_forest);

    cmc::cmc_global_msg("The mesh and the data have been constructed.");

    return std::make_pair(partitioned_forest, example_data);
}

#if 0
template <typename T>
std::vector<T>
ReadDataFromStream(const std::string& file_name, const size_t num_datums)
{
    std::vector<T> values;
    size_t count{0};

    std::ifstream istrm(file_name, std::ios::binary);
    if (!istrm.is_open())
        std::cout << "failed to open " << file_name << '\n';
    else
    {
        while ( !istrm.eof() && count < num_datums) {
            T value;
            istrm.read(reinterpret_cast<char*>(&value), sizeof(T));

            values.push_back(value);
            ++count;
        }
        istrm.close();
    }

    return values;
}

void
WriteData(t8_forest_t forest, std::vector<double>& data)
{
    std::vector<double> double_data;
    double_data.reserve(data.size());

    for (auto val_iter = data.begin(); val_iter != data.end(); ++val_iter)
    {
        double_data.push_back(static_cast<double>(*val_iter));
    }

    t8_vtk_data_field_t vtk_data[1];
    snprintf (vtk_data[0].description, BUFSIZ, "ExampleData");
    vtk_data[0].type = T8_VTK_SCALAR;
    vtk_data[0].data = double_data.data();

    t8_forest_write_vtk_ext (forest, "cmc_t8_hybrid_circlesquare_example_mesh_and_data", 0, 0, 0, 0, 0, 0, 0, 1, vtk_data);
}
#endif

}

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {

    {
    
    /* Construct the initial mesh and the data*/
    auto [mesh, data] = cmc::example::BuildInitialForest();

    // Test parallel multi res with huffman entropy coding
    #if 1
    /* Create a compression variabvle from the mesh and the data */
    cmc::lossless::par::multi_res::CompressionVariable<float> var("test_var", mesh, data);

    /* Perform the compression */
    var.Compress();

    /* Write out the compressed data to disk */
    cmc::compression_io::mpi::Writer<float> writer("cmc_mpi_multi_res_example_data.cmc", MPI_COMM_WORLD);
    writer.SetVariable(var);
    writer.Write();

    cmc::cmc_debug_msg("Reader Start....");

    cmc::compression_io::mpi::Reader reader("cmc_mpi_multi_res_example_data.cmc", MPI_COMM_WORLD);
    std::vector<cmc::compression_io::mpi::VariableHull> variable_hulls = reader.ReadVariableHulls();

    cmc::cmc_debug_msg("num variables: ", variable_hulls.size());

    cmc::cmc_debug_msg("Var ID: ", variable_hulls[0].id, ", Name: ", variable_hulls[0].name, ", Mesh count: ", variable_hulls[0].global_num_mesh_bytes, ", Data Count: ", variable_hulls[0].global_num_data_bytes, ", Global Mesh Offset: ", variable_hulls[0].global_file_offset_mesh_start, ", Global Data Offset: ", variable_hulls[0].global_file_offset_data_start);

    #else
    // Test serial default Multi res compresison with entropy coding 
    cmc::lossless::multi_res::CompressionVariable<float> var("test_var", mesh, data);

    /* Perform the compression */
    var.Compress();

    /* Write out the compressed data to disk */
    cmc::compression_io::nc::Writer writer("cmc_serial_multi_res_example_data.cmc", MPI_COMM_SELF);
    writer.SetVariable(&var);
    writer.Write();

    #endif
    }

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}

#endif /* !CMC_T8_CONSTRUCT_MESH_AND_DATA_FOR_COMPRESSION_HXX */
