

#ifndef CMC_T8_CONSTRUCT_AMR_MESH_EXAMPLE_HXX
#define CMC_T8_CONSTRUCT_AMR_MESH_EXAMPLE_HXX

#include "cmc.hxx"

#include <t8_schemes/t8_default/t8_default.hxx>
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest_geometrical.h>
#include <t8_cmesh/t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_cmesh/t8_cmesh_io/t8_cmesh_readmshfile.h>
#include "t8code/cmc_t8_mesh.hxx"
#include "mpi/cmc_mpi.hxx"

#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <vector>
#include <cmath>
#include <numeric>
#include <utility>
#include <cmath>
#include <fstream>
#include <cstdlib>

#include "amr/lossy/cmc_par_multi_res_extraction_idw.hxx"
#include "amr/lossy/cmc_par_multi_res_decompression_idw.hxx"

#include "amr/lossless/cmc_par_multi_res_extraction.hxx"
#include "amr/lossless/cmc_par_multi_res_decompression.hxx"

constexpr bool kWriteVTKData = false;

namespace cmc::example
{

constexpr int init_elem_level = 6;
constexpr int adapt_iterations = 2;
constexpr float kCoarsenThresh = 0.006;
constexpr float kRefineThresh = 0.08;

struct vortex
{
    float x;
    float y;
    float ampl;
    float lx;
    float ly;
};

/* Define severeal vortices */
const std::vector<vortex> vortices = {
    {1.2f, 0.0f, 0.17f, 0.35f, 0.22},
    {1.7f, 0.22f, -0.15f, 0.4f, 0.25f},
    {1.7f, -0.22f, 0.15f, 0.4f, 0.25f}
};

constexpr float flowsp = 1.0f;
constexpr float cr = 1.0f;

std::pair<float, float>
GenerateData(const float x, const float y)
{
    const float r = (x*x + y*y) * (x*x + y*y);

    /* Potential flow */
    float ux = flowsp * (1.0f - cr*cr*(x*x -y*y)/r);
    float uy = -2.0f * flowsp * cr*cr*x*y/r;

    /* Apply Vortices */
    for (auto viter = vortices.begin(); viter != vortices.end(); ++viter)
    {
        const vortex& v = *viter;

        const float diffx = x - v.x;
        const float diffy = y - v.y;

        const float g = std::exp(-(diffx*diffx) / (v.lx*v.lx) - (diffy*diffy) / (v.ly*v.ly));

        ux += -2.0f * v.ampl * diffy / (v.ly * v.ly) * g;
        uy += 2.0f * v.ampl * diffx / (v.lx * v.lx) * g;
    }

    return std::make_pair(ux, uy);
}

float
AdaptationCriterion(const float x, const float y, const float h)
{
    constexpr float delta = 1e-4;
    const auto [xp1, xp2] = GenerateData(x + delta, y);
    const auto [xm1, xm2] = GenerateData(x - delta, y);
    const auto [yp1, yp2] = GenerateData(x, y + delta);
    const auto [ym1, ym2] = GenerateData(x, y - delta);

    /* Approximate differentiation */
    const float dux_dx = (xp1 - xm1) / (2.0f * delta);
    const float duy_dx = (xp2 - xm2) / (2.0f * delta);
    const float dux_dy = (yp1 - ym1) / (2.0f * delta);
    const float duy_dy = (yp2 - ym2) / (2.0f * delta);

    const float gn = std::sqrt(dux_dx * dux_dx + duy_dx * duy_dx + dux_dy * dux_dy + duy_dy * duy_dy);

    return h*gn / flowsp;
}

static int
t8_example_adapt ([[maybe_unused]] t8_forest_t forest, t8_forest_t forest_from, t8_locidx_t ltree_id, [[maybe_unused]] const t8_eclass_t tree_class,
                  [[maybe_unused]] t8_locidx_t lelement_id, [[maybe_unused]] const t8_scheme *scheme, const int is_family,
                  [[maybe_unused]] const int num_elements, t8_element_t *elements[])
{
    /* Get the element size */
    const double elem_diam = t8_forest_element_diam (forest_from, ltree_id, elements[0]);
  
    /* Get the midpoint */
    std::vector<double> midpoint(3);
    t8_forest_element_centroid (forest_from, ltree_id, elements[0], midpoint.data());

    /* if close to the cylinder, refine */
    if (std::sqrt(midpoint[0]*midpoint[0] + midpoint[1]*midpoint[1]) < 1.1)
    {
        return 1;
    }

    /* Evaluate crtierion */
    const float eta = AdaptationCriterion(midpoint[0], midpoint[1], elem_diam);

    if (eta > kRefineThresh)
    {
        return 1;
    } else if (eta < kCoarsenThresh)
    {
        if (is_family)
        {
            return -1;
        } else
        {
            return 0;
        }
    }

    return 0;
}

inline 
std::pair<t8_cmesh_t, const t8_scheme*>
GetInitCMeshAndScheme()
{
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_from_msh_file (&cmesh, "../sim_setup/t8code_build/example/IO/cmesh/gmsh/circlesquare_hybrid_hole",0,MPI_COMM_SELF,2,0,0);
    const t8_scheme *default_scheme = t8_scheme_new_default ();
    return std::make_pair(cmesh, default_scheme);
}

inline 
std::pair<t8_forest_t, std::vector<float>>
BuildInitialMeshAndData()
{
    cmc::cmc_global_msg("The intial mesh and data will be constructed.");
    t8_cmesh_t cmesh;
    t8_cmesh_init (&cmesh);
    t8_cmesh_from_msh_file (&cmesh, "../sim_setup/t8code_build/example/IO/cmesh/gmsh/circlesquare_hybrid_hole",0,MPI_COMM_SELF,2,0,0);
    const t8_scheme *default_scheme = t8_scheme_new_default ();
    t8_forest_t forest = t8_forest_new_uniform (cmesh, default_scheme, init_elem_level, 0, MPI_COMM_SELF);

    t8_forest_t forest_adapt;
    for (int i = 0; i < adapt_iterations; ++i)
    {
        cmc::cmc_global_msg("Mesh adaptation step ", i + 1);
        t8_forest_init(&forest_adapt);
        t8_forest_set_adapt (forest_adapt, forest, t8_example_adapt, 0);
        t8_forest_set_balance (forest_adapt, forest, 0);
        t8_forest_commit (forest_adapt);
        forest = forest_adapt;
    }

    /* Num mesh elements */
    const uint64_t num_mesh_elems = t8_forest_get_global_num_leaf_elements(forest);
    cmc::cmc_global_msg("Number of global mesh elements: ", num_mesh_elems);

    std::vector<float> data;
    data.reserve(num_mesh_elems);

    /* Iterate over the mesh and generate the data */
    t8_locidx_t num_trees = t8_forest_get_num_local_trees (forest);

    for (t8_locidx_t itree = 0; itree < num_trees; itree++)
    {
        t8_locidx_t num_elems_in_tree = t8_forest_get_tree_num_leaf_elements (forest, itree);

        for (t8_locidx_t ielement = 0; ielement < num_elems_in_tree; ielement++)
        {
            const t8_element_t *element = t8_forest_get_leaf_element_in_tree (forest, itree, ielement);

            /* Evaluate the data for this element */
            std::vector<double> midpoint(3);
            t8_forest_element_centroid (forest, itree, element, midpoint.data());

            const auto [vx, vy] = GenerateData(midpoint[0], midpoint[1]);
            const float u = std::sqrt(vx*vx + vy*vy);
            data.push_back(u);
        }
    }

    return std::make_pair(forest, std::move(data));
}

}

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {
     
    const std::string file_name = "cmc_amr_example_compressed.cmc";
    std::vector<float> decompressed_data;
    
    /* Construct the initial mesh and the data*/
    auto [mesh, data] = cmc::example::BuildInitialMeshAndData();

    std::FILE* file_out1 = std::fopen("cmc_amr_example_init_data.cmc", "wb");
    std::fwrite(data.data(), sizeof(float), data.size(), file_out1);
    std::fclose(file_out1);

    if constexpr (kWriteVTKData)
    {
        std::vector<double> ddata;
        for (const auto& val : data)
        {
            ddata.push_back(static_cast<double>(val));
        }
        t8_vtk_data_field_t vtk_data[1];
        constexpr int kCMCT8BufferSize{1024};
        snprintf (vtk_data[0].description, kCMCT8BufferSize, "values");
        vtk_data[0].type = T8_VTK_SCALAR;
        vtk_data[0].data = ddata.data();

        t8_forest_write_vtk_ext (mesh, "cmc_test_amr_example_init_mesh_and_data", 1, 1, 1, 1, 0, 0, 0, 1, vtk_data);
    }

    /* Create an error domain */
    /* Set a relative error criterion */

    const float kRelMaxError = 0.01;
    cmc::ErrorDomain global_rel_error_domain(cmc::PermittedError(cmc::CompressionCriterion::RelativeErrorThreshold, kRelMaxError), cmc::error_domain_fn::GeneralErrorCriterion);
    
    [[maybe_unused]] std::vector<cmc::ErrorDomain> error_domains;
    error_domains.push_back(global_rel_error_domain);

    //Compression
    {
    /* Perform the lossy compression with multiple data */
    cmc::par::lossy::multi_res::idw::CompressionVariable<float, 2> variable("test_var", mesh, data, error_domains);
    //cmc::par::lossless::multi_res::CompressionVariable<float, 2> variable("test_var", mesh, data);

    /* Compress */
    variable.Compress();

    /* Write the encoded data */
    variable.WriteCompressedData(file_name);

    /* Deallocate the initial mesh */
    t8_forest_unref(&mesh);
    }

    
    //Decompression
    {
    
    /* Get the initial base mesh */
    const auto [cmesh, scheme] = cmc::example::GetInitCMeshAndScheme();

    /* Decompress the mesh and the data */
    cmc::par::lossy::multi_res::idw::DecompressionVariable<float, 2> decompression_variable(file_name, MPI_COMM_WORLD, cmesh, scheme);
    //cmc::par::lossless::multi_res::DecompressionVariable<float, 2> decompression_variable(file_name, MPI_COMM_WORLD, cmesh, scheme);
    
    decompression_variable.Decompress();

    /* Access the decompressed components */
    auto [decompressed_forest, decompressed_data_vec] = decompression_variable.GetDecompressedData();
    
    /* Save the decompressed data */
    decompressed_data = std::move(decompressed_data_vec);

    std::FILE* file_out = std::fopen("cmc_amr_example_decompressed_data.cmc", "wb");
    std::fwrite(decompressed_data.data(), sizeof(float), decompressed_data.size(), file_out);
    std::fclose(file_out);

    if constexpr (kWriteVTKData)
    {
        std::vector<double> ddata;
        for (const auto& val : decompressed_data)
        {
            ddata.push_back(static_cast<double>(val));
        }
        t8_vtk_data_field_t vtk_data[1];
        constexpr int kCMCT8BufferSize{1024};
        snprintf (vtk_data[0].description, kCMCT8BufferSize, "values");
        vtk_data[0].type = T8_VTK_SCALAR;
        vtk_data[0].data = ddata.data();
        t8_forest_write_vtk_ext (mesh, "cmc_amr_example_decompressed_mesh_and_data", 1, 1, 1, 1, 0, 0, 0, 1, vtk_data);
    }

    /* Deallocate the decompressed forest */
    t8_forest_unref(&decompressed_forest);
    }

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}

#endif /* !CMC_T8_CONSTRUCT_AMR_MESH_EXAMPLE_HXX */
