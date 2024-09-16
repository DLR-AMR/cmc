#include "cmc.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "lossy/cmc_amr_lossy_compression_settings.hxx"
#include "lossy/cmc_amr_lossy_compression.hxx"
#include "lossy/cmc_prefix_lossy_compression.hxx"
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
    
    sc_init (MPI_COMM_WORLD, 1, 1, NULL, SC_LP_ESSENTIAL);
    t8_init (SC_LP_DEBUG);

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
        const std::string file = "../../data/mptrac_era5_2021_07_01_00.nc";
        cmc::NcData nc_data(file, cmc::NcOpeningMode::Serial);

        //cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 700, 716),
        //                         cmc::DimensionInterval(cmc::Dimension::Lat, 344, 360)
        //                         //cmc::DimensionInterval(cmc::Dimension::Lev, 34, 36)
        //                         );
        //nc_data.SetHintHeightDimension(2);
        nc_data.SetHintHeightDimension(6); //Set plev as height dimension
        cmc::Hyperslab hyperslab(cmc::DimensionInterval(cmc::Dimension::Lon, 0, 1200),
                                 cmc::DimensionInterval(cmc::Dimension::Lat, 0, 601),
                                 cmc::DimensionInterval(cmc::Dimension::Lev, 136, 137)
                                 );

        /* Inquire the hyperslab of data for the given variables */
        nc_data.InquireVariables(hyperslab, "t");

        /* Close the file, since we have gathered the data we wanted */
        nc_data.CloseFileHandle();

        #if 1
        /* Create compression settings */
        cmc::CompressionSettings settings;

        const double abs_max_err = 0.0;
        settings.SetAbsoluteErrorCriterion(abs_max_err, cmc::kErrorCriterionHoldsForAllVariables);
        //const double rel_max_err = 0.01;
        //settings.SetRelativeErrorCriterion(rel_max_err, cmc::kErrorCriterionHoldsForAllVariables);
        //cmc::SplitVariable split(cmc::kSplitAllVariables, cmc::Dimension::Lev);
        //settings.SplitVariableByDimension(split);

        #if 1
        /* Create the compression data */
        cmc::CompressionData compression_data(nc_data.TransferData(), std::move(settings));

        /* Setup the example data for the compression */
        //compression_data.Setup(initial_mesh);
        compression_data.Setup();
        compression_data.WriteVTKFile("mptrac_t_initial");

        compression_data.Compress(cmc::CompressionMode::OneForOne);
        compression_data.WriteVTKFile("mptrac_t_abs_2_0");
        //compression_data.WriteCompressedData("co2_emi_data");
        #else
        /* Create the compression data */
        cmc::PrefixCompressionData compression_data(nc_data.TransferData(), std::move(settings));

        /* Setup the example data for the compression */
        const bool perform_default_lossy_compression_as_well = false;
        compression_data.Setup(perform_default_lossy_compression_as_well);

        compression_data.Compress();
        compression_data.WriteCompressedDataEGU("co2_emi_data_pref.nc");

        #endif
        #if 0
        compression_data.Compress(cmc::CompressionMode::OneForOne);

        compression_data.WriteVTKFile("ExTempVar");

        /* Decompress the variable and receive the external wrapper to it */
        cmc::OutputVar decomrpessed_float_var = compression_data.DecompressVariable(5);
        /* Obtain the actual decompressed variable based on it's data type */
        cmc::OutputVariable<float> decompressed_float_variable = decomrpessed_float_var.SeizeOutputVariable<float>();
        /* Get it's data */
        std::vector<float> decompressed_float_data = decompressed_float_variable.GetData();
        cmc::cmc_debug_msg("Num decompressed float data: ", decompressed_float_data.size());
        
        #endif

        //for (auto iter = decompressed_float_data.begin(); iter != decompressed_float_data.end(); ++iter)
        //{
        //    cmc::cmc_debug_msg(*iter, ", ");
        //}
        #endif

        #if 0
        compression_data.Decompress();
        compression_data.WriteVTKFile("ExTempDecompressed");
        #endif
    }
    sc_finalize ();
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}