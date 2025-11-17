#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "amr/lossy/cmc_ac_compression_variable.hxx"
#include "amr/lossy/cmc_adaptive_coarsening.hxx"
#include "utilities/cmc_compression_settings.hxx"


#include <numeric>

const double kRelMaxError = 0.01;
const double kAbsMaxError = 1.0;

const double kNewLowerRelMaxError = 0.005;
const double kNewLowerAbsMaxError = 0.5;

static bool
EvaluateSetErrorCriteria(const cmc::PermittedError error)
{
    if (error.criterion == cmc::CompressionCriterion::AbsoluteErrorThreshold)
    {
        /* Absolute error criterion */
        if (error.error == kAbsMaxError)
        {
            return true;
        } else
        {
            return false;
        }
    } else
    {
        /* Relative error criterion */
        if (error.error == kRelMaxError)
        {
            return true;
        } else
        {
            return false;
        }
    }
}


static bool
EvaluateNewLowerSetErrorCriteria(const cmc::PermittedError error)
{
    if (error.criterion == cmc::CompressionCriterion::AbsoluteErrorThreshold)
    {
        /* Absolute error criterion */
        if (error.error == kNewLowerAbsMaxError)
        {
            return true;
        } else
        {
            return false;
        }
    } else
    {
        /* Relative error criterion */
        if (error.error == kNewLowerRelMaxError)
        {
            return true;
        } else
        {
            return false;
        }
    }
}

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
    const t8_scheme_c* ts = t8_forest_get_scheme(forest);

    /* Create a variable for containing a value for each element */
    const t8_locidx_t num_elements = t8_forest_get_local_num_leaf_elements(forest);
    std::vector<float> data(num_elements);
    std::iota(data.begin(), data.end(), 0.0);

    /* Allocate the compression settings */
    cmc::CompressionSettings settings;

    /* Set a relative error criterion */
    settings.SetGeneralErrorCriterion(cmc::CompressionCriterion::RelativeErrorThreshold, kRelMaxError);

    /* Set an absolute error criterion */
    cmc::ErrorDomain error_domain(cmc::PermittedError(cmc::CompressionCriterion::AbsoluteErrorThreshold, kAbsMaxError), cmc::error_domain_fn::GeneralErrorCriterion);
    settings.SetErrorDomain(error_domain);
    
    /* Get the first element from the mesh */
    const t8_element_t* first_elem = t8_forest_get_leaf_element_in_tree(forest, 0, 0);

    /* Get the permitted error criteria for the element */
    std::vector<cmc::PermittedError> errors = settings.FindRestrictingErrors(forest, 0, 0, ts, 1, &first_elem);
    
    /* There should be two errors associated to the element */
    cmc::ExpectTrue(errors.size() == 2);

    cmc::PermittedError first_error = errors.front();
    cmc::PermittedError second_error = errors.back();

    /* Check if the error criteria are set correctly */
    cmc::ExpectTrue(EvaluateSetErrorCriteria(first_error));
    cmc::ExpectTrue(EvaluateSetErrorCriteria(second_error));

    /** Set new lower error criteria **/
    /* Set an absolute error criterion */
    cmc::ErrorDomain error_domain2(cmc::PermittedError(cmc::CompressionCriterion::AbsoluteErrorThreshold, kNewLowerAbsMaxError), cmc::error_domain_fn::GeneralErrorCriterion);
    settings.SetErrorDomain(error_domain2);
    /* Set an absolute error criterion */
    cmc::ErrorDomain error_domain3(cmc::PermittedError(cmc::CompressionCriterion::RelativeErrorThreshold, kNewLowerRelMaxError), cmc::error_domain_fn::GeneralErrorCriterion);
    settings.SetErrorDomain(error_domain3);

    /* Find the restriciting error criteria for the element again */
    std::vector<cmc::PermittedError> new_errors = settings.FindRestrictingErrors(forest, 0, 0, ts, 1, &first_elem);
    
    /* There should be two errors associated to the element */
    cmc::ExpectTrue(new_errors.size() == 2);

    cmc::PermittedError new_first_error = new_errors.front();
    cmc::PermittedError new_second_error = new_errors.back();

    /* Check if the new lower error criteria are set correctly */
    cmc::ExpectTrue(EvaluateNewLowerSetErrorCriteria(new_first_error));
    cmc::ExpectTrue(EvaluateNewLowerSetErrorCriteria(new_second_error));

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
