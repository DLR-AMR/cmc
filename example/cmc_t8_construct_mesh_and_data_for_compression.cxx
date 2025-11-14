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

#if 0
#include "lossless/cmc_multi_res_extraction_compression.hxx"
#include "lossless/cmc_multi_res_extraction_decompression.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "compression_io/cmc_compression_output.hxx"
#include "compression_io/cmc_decompression_input.hxx"
#else

#include "lossless/cmc_multi_res_extraction_compression.hxx"
#include "lossless/cmc_multi_res_extraction_decompression.hxx"

#include "lossless/cmc_prefix_extraction_compression.hxx"
#include "lossless/cmc_prefix_extraction_compression_plain_suffixes.hxx"
#include "lossless/cmc_prefix_extraction_decompression.hxx"

#include "lossless/cmc_test_pcp4_compression.hxx"

#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_utilities.hxx"
#include "compression_io/cmc_compression_output.hxx"
#include "compression_io/cmc_decompression_input.hxx"
#endif

#include <vector>
#include <cmath>
#include <numeric>
#include <utility>
#include <cmath>
#include <fstream>
#include <cstdlib>


namespace cmc::example
{

//constexpr int max_elem_level = 11;
//constexpr int init_elem_level = 7;
//constexpr double bandwidth = 64.0;

constexpr int max_elem_level = 11;
constexpr int init_elem_level = 7;
constexpr double bandwidth = 26.5;

typedef struct
{
  std::vector<double> M{0,0,0};
  /**< midpoint */
  double radius{0.25};
  /**< radius */
} levelset_sphere_data_t;

inline double
dist (const std::vector<double> &point_x, const std::vector<double> &point_y)
{
    cmc_assert(point_x.size() == point_y.size());
    double dist = std::inner_product (point_x.begin (), point_x.end (), point_y.begin (), 0.0, std::plus<double> (),
                                    [] (double x, double y) { return (x - y) * (x - y); });
    return std::sqrt (dist);
}

double
levelset_sphere_fn (const std::vector<double>  &x, void *data, const int level)
{
    levelset_sphere_data_t *ls_data = (levelset_sphere_data_t *) data;

  T8_ASSERT (ls_data->radius > 0);

  if (dist (x, ls_data->M) >= 2.0)
  {
    return (8/(dist (x, ls_data->M) * dist (x, ls_data->M) * dist (x, ls_data->M))) * fabs(std::sin(1.57079632679 * x[1])) * (fabs(x[0]) + 0.25 * fabs(x[1])) * std::cos(1.57079632679 - dist (x, ls_data->M) + 0.25) + 0.01*std::sin(1000*3.14*fabs(x[0]));
  }
  //return (8) *std::fabs(cos(0.5 * x[0]))*(dist (x, ls_data->M) - ls_data->radius);
  //return x[0] * std::log(8*dist (x, ls_data->M)) - 2.5;
  //return std::fabs(cos(0.55 * x[0]))*(dist (x, ls_data->M));
  //return 2*fabs(std::cos(1.57079632679 * x[0] * x[1]));
  return fabs(std::sin(1.57079632679 * x[1])) * (fabs(x[0]) + 0.25 * fabs(x[1])) * std::cos(1.57079632679 - dist (x, ls_data->M) + 0.25) + 0.01*std::sin(1000*3.14*fabs(x[0]));
}

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

inline std::vector<double>
t8_advect_problem_init_elements (t8_forest_t forest)
{
    levelset_sphere_data_t ls_data;
    ls_data.M[0] = 0.0;
    ls_data.M[1] = 0.0;
    ls_data.M[2] = 0.0;
    ls_data.radius = 1.00;

    t8_locidx_t num_local_elems = t8_forest_get_local_num_leaf_elements(forest);

    std::vector<double> values;
    values.reserve(num_local_elems);

    const t8_scheme *scheme = t8_forest_get_scheme (forest);

    t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest);

    for (t8_locidx_t itree = 0, idata = 0; itree < num_trees; itree++)
    {
        const t8_eclass_t tree_class = t8_forest_get_tree_class (forest, itree);
        t8_locidx_t num_elems_in_tree = t8_forest_get_tree_num_elements (forest, itree);

        for (t8_locidx_t ielement = 0; ielement < num_elems_in_tree; ielement++)
        {
            const t8_element_t *element = t8_forest_get_element_in_tree (forest, itree, ielement);

            const int level = scheme->element_get_level (tree_class, element);

            std::vector<double> midpoint(3);
                
            t8_forest_element_centroid (forest, itree, element, midpoint.data());

            values.push_back(static_cast<double>(levelset_sphere_fn(midpoint, &ls_data,level)));
        }
    }
    return values;
}

struct adapt_data_t
{
    std::vector<double> data;
};

std::vector<double> CenterPoint{0,0,0};

static int
t8_advect_adapt (t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree_id, const t8_eclass_t tree_class,
                 t8_locidx_t lelement_id, const t8_scheme *scheme, const int is_family,
                 [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
    adapt_data_t* adapt_data = (adapt_data_t*) t8_forest_get_user_data (forest);

    const int level = scheme->element_get_level (tree_class, elements[0]);
    if (level == max_elem_level) {
        /* It is not possible to refine this level */
        return 0;
    }
    
    const t8_locidx_t offset = t8_forest_get_tree_element_offset (forest_from, ltree_id);
    const float val = adapt_data->data[offset + lelement_id];

    const double elem_diam = t8_forest_element_diam (forest_from, ltree_id, elements[0]);
  
    std::vector<double> midpoint(3);
                
    t8_forest_element_centroid (forest_from, ltree_id, elements[0], midpoint.data());

    #if 0
    if (fabs(midpoint[0]) <= 0.7 && fabs(midpoint[1]) <= 0.7)
    {
         /* refine if level is not too large */
        return level < max_elem_level;       
    } else if (fabs(midpoint[0]) >= 0.9 && fabs(midpoint[1]>= 0.9))
    {
        return -(is_family && level > init_elem_level);
    } else
    {
        return 0;
    }
    #endif

    #if 1
    if (fabs (val) > 2 * bandwidth * elem_diam) {
        /* coarsen if this is a family and level is not too small */
        return -(is_family && level > init_elem_level);
    }
    else if (fabs (val) < bandwidth * elem_diam) {
        /* refine if level is not too large */
        return level < max_elem_level;
    }
    #endif

    return 0;
}

inline 
std::pair<t8_forest_t, std::vector<double>>
BuildInitialForest()
{
    t8_cmesh_t cmesh = t8_cmesh_from_msh_file ("../programs/t8code/example/IO/cmesh/gmsh/circlesquare_hybrid_hole",0,MPI_COMM_SELF,2,0,0);
    const t8_scheme *default_scheme = t8_scheme_new_default ();
    t8_forest_t forest = t8_forest_new_uniform (cmesh, default_scheme, init_elem_level, 1, MPI_COMM_SELF);
    
    adapt_data_t adapt_data;
    adapt_data.data = cmc::example::t8_advect_problem_init_elements(forest);

    const int num_iterations = max_elem_level - init_elem_level - 1;

    t8_forest_t forest_adapt;
    for (int i = 0; i < num_iterations; ++i)
    {
        t8_forest_init(&forest_adapt);
        t8_forest_set_user_data (forest_adapt, &adapt_data);
        t8_forest_set_adapt (forest_adapt, forest, t8_advect_adapt, 0);
        t8_forest_set_balance (forest_adapt, NULL, 1);
        t8_forest_commit (forest_adapt);
        forest = forest_adapt;
        adapt_data.data = cmc::example::t8_advect_problem_init_elements(forest);
    }

    //for (auto iter = adapt_data.data.begin(); iter != adapt_data.data.end(); ++iter)
    //{
    //    *iter = (-1.0) * (*iter);
    //}

    for (auto iter = adapt_data.data.begin(); iter != adapt_data.data.end(); ++iter)
    {
        *iter = (2.25 - (*iter));
    }

    double max = -1000.0;
    double min = +1000.0;

    for (auto iter = adapt_data.data.begin(); iter != adapt_data.data.end(); ++iter)
    {
        if (max < *iter)
        {
            max = *iter;
        }

        if (min > *iter)
        {
            min = *iter;
        }
    }

    cmc_debug_msg("Data Values: Max: ", max, ", Min: ", min);

    return std::make_pair(forest, std::move(adapt_data.data));
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
}

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {

    std::vector<double> init_data;

    {
    
    /* Construct the initial mesh and the data*/
    auto [mesh, data] = cmc::example::BuildInitialForest();
    //cmc::example::WriteData(mesh, data);

    /* Write out the data in order for other compressors to compress it */
    init_data = data;
    FILE* file_out = fopen("hybrid_square_circle_init_values.bin", "wb");
    fwrite(init_data.data(), sizeof(double), init_data.size(), file_out);
    fclose(file_out);

    /* Set up either a PrefixAMR or MultiResAMR compression */
    #if 0
    cmc::lossless::prefix::plain_suffix::CompressionVariable<double> var("test_var", mesh, data);
    #else
    cmc::lossless::multi_res::CompressionVariable<double> var("test_var", mesh, data);
    //cmc::lossless::test_pcp4::CompressionVariable<double> var("test_var", mesh, data);
    #endif

    /* Perform the compression */
    var.Compress();

    /* Write out the compressed data to disk */
    cmc::compression_io::Writer writer("example_lossless_compression_output.cmc", MPI_COMM_SELF);
    writer.SetVariable(&var);
    writer.Write();

    }

    /* Create a reader for the compressed output that has been stored */
    cmc::compression_io::Reader reader("example_lossless_compression_output.cmc", MPI_COMM_SELF);

    /* Create an embedded decompressor from the compressed data */
    std::unique_ptr<cmc::decompression::AbstractByteDecompressionVariable<double>> decompression_var = reader.ReadVariableForDecompression<double>("test_var");

    /* Decompress the encoded data */
    decompression_var->Decompress();

    if (not decompression_var->Size() == init_data.size()) {cmc::cmc_debug_msg("Error: Initial data and Decompressed data are unequal in their size");}

    /* Access the decompressed data */
    std::vector<cmc::SerializedCompressionValue<sizeof(double)>> decompressed_byte_values = decompression_var->GetDecompressedData();

    /* Compare the decompressed data with the initial data */
    bool inid = true;
    int idx = 0;
    for (auto cr_iter = decompressed_byte_values.begin(); cr_iter != decompressed_byte_values.end(); ++cr_iter, ++idx)
    {
        const double val = cr_iter->ReinterpretDataAs<double>();
        if (not cmc::ApproxCompare(val, init_data[idx]))
        {
            cmc::cmc_debug_msg("Not equal at index: ", idx, " decompressed value: ", val, " and initial value: ", init_data[idx]);
            inid = false;
        }
    }
    cmc::cmc_debug_msg("Comparison of initial and decompressed data yields: ", inid);

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}

#endif /* !CMC_T8_CONSTRUCT_MESH_AND_DATA_FOR_COMPRESSION_HXX */
