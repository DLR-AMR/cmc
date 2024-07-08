#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "lossy/cmc_amr_lossy_compression_settings.hxx"
#include "lossy/cmc_amr_lossy_compression.hxx"
#include "netcdf/cmc_netcdf.hxx"

#include "t8_cmesh.h"
#include "t8_forest/t8_forest.h"
#include "t8_geometry/t8_geometry.h"
#if CMC_WITH_T8CODE
#include <t8_eclass.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_common/t8_default_common_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_quad/t8_default_quad_cxx.hxx>
#include <t8_schemes/t8_default/t8_default_hex/t8_default_hex_cxx.hxx>
#include "t8_schemes/t8_default/t8_default_c_interface.h"
#include "t8_element_c_interface.h"
#include <p4est.h>
#include <p8est.h>
#include <t8_geometry/t8_geometry_implementations/t8_geometry_linear.h>
#endif

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {
        #if 0
        /* Create a mesh */
        t8_cmesh_t cmesh;
        t8_cmesh_init (&cmesh);
        t8_geometry_c      *linear_geom = t8_geometry_linear_new (2);

        double vertices_quad[] = {0,0,0,
                                  1,0,0,
                                  0,1,0,
                                  1,1,0
                                  };

        double vertices_triangle[] = {0,1,0,
                                      1,1,0,
                                      0.5,1.5,0
                                      };
        #if 1
        t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_QUAD);
        t8_cmesh_set_tree_class (cmesh, 1, T8_ECLASS_TRIANGLE);
        
        t8_cmesh_register_geometry (cmesh, &linear_geom);

        t8_cmesh_set_tree_vertices (cmesh, 0, vertices_quad, 4);
        t8_cmesh_set_tree_vertices (cmesh, 1, vertices_triangle, 3);

        t8_cmesh_set_join (cmesh, 0, 1, 3, 2, 0);
        #else
        t8_cmesh_register_geometry (cmesh, linear_geom);
        t8_cmesh_set_tree_class (cmesh, 0, T8_ECLASS_TRIANGLE);
        t8_cmesh_set_tree_vertices (cmesh, 0, vertices_triangle, 3);
        #endif

        t8_cmesh_commit (cmesh, MPI_COMM_WORLD);

        t8_scheme_cxx_t *scheme = t8_scheme_new_default_cxx ();

        const int initial_ref_lvl = 5;

        t8_forest_t forest = t8_forest_new_uniform(cmesh, scheme, initial_ref_lvl, 0, MPI_COMM_WORLD);

        cmc::AmrMesh initial_mesh{forest, initial_ref_lvl};
        initial_mesh.IndicateWhetherDummyElementsArePresent(false);

        t8_forest_write_vtk(forest, "example_forest_diff_mesh");
        #endif

        /* Create input variables */
        const std::string file = "../../data/era5_reanalysis_pressure_lvls_fixed_time.nc";
        cmc::NcData nc_data(file, cmc::NcOpeningMode::Serial);

        cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 16),
                                 cmc::DimensionInterval(cmc::Dimension::Lat, 0, 8)
                                 );

         /* Inquire the hyperslab of data for the given variables */
        nc_data.InquireVariables(hyperslab, "t");

        /* Close the file, since we have gathered the data we wanted */
        nc_data.CloseFileHandle();

        #if 1
        /* Create compression settings */
        cmc::CompressionSettings settings;

        const double abs_max_err = 1.0;
        settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);

        /* Create the compression data */
        cmc::CompressionData compression_data(nc_data.TransferData(), std::move(settings));

        /* Setup the example data for the compression */
        //compression_data.Setup(initial_mesh);
        compression_data.Setup();
        
        compression_data.Compress(cmc::CompressionMode::OneForOne);

        compression_data.WriteVTKFile("ExTempVar");

    
        #endif

    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}