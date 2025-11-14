#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "lossy/cmc_ac_compression_variable.hxx"
#include "lossy/cmc_adaptive_coarsening.hxx"
#include "utilities/cmc_compression_settings.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8_cmesh/t8_cmesh_examples.h>
#include <t8_schemes/t8_default/t8_default.hxx>
#endif

#include <numeric>

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
    const t8_scheme *scheme = t8_scheme_new_default();
    const int initial_level = 3;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);

    /* Create a variable for containing a value for each element */
    const t8_locidx_t num_elements = t8_forest_get_local_num_leaf_elements(forest);
    std::vector<float> data(num_elements);
    std::iota(data.begin(), data.end(), 0.0);

    /* Allocate the compression settings */
    cmc::CompressionSettings settings;

    /* Set a relative error criterion */
    const double rel_max_err = 0.01;
    settings.SetGeneralErrorCriterion(cmc::CompressionCriterion::RelativeErrorThreshold, rel_max_err);

    /* Compress the data */
    cmc::lossy::CompressionVariable variable("test_var", forest, data);


    variable.Compress(settings);


    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
