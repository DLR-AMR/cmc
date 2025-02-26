#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "lossy/cmc_ac_compression_variable.hxx"
#include "lossy/cmc_adaptive_coarsening.hxx"
#include "utilities/cmc_compression_settings.hxx"


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
    t8_scheme_cxx_t* scheme = t8_scheme_new_default_cxx ();
    const int initial_level = 3;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);
    t8_eclass_scheme_c* ts = t8_forest_get_eclass_scheme(forest, t8_forest_get_eclass(forest, 0));

    const t8_locidx_t num_elements = t8_forest_get_local_num_elements(forest);

    std::vector<float> data(num_elements);

    std::iota(data.begin(), data.end(), 0.0);

    cmc::CompressionSettings settings;

    const double rel_max_err = 0.01;

    settings.SetRelativeErrorCriterion(rel_max_err, cmc::kErrorCriterionHoldsForAllVariables);


    cmc::lossy::CompressionVariable variable("test_var", forest, data);

    variable.Compress(settings);


    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
