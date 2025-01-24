#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_span.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_cmesh.h>
#include <t8_cmesh/t8_cmesh_examples.h> 
#include <t8_forest/t8_forest_general.h>
#include <t8_forest/t8_forest.h>
#include <t8_schemes/t8_default/t8_default_cxx.hxx> 
#include <t8_forest/t8_forest_iterate.h> 
#endif

#include <vector>

static int
coarsen_all_elements(t8_forest_t forest, t8_forest_t forest_from, int which_tree, int lelement_id, t8_eclass_scheme_c * ts,
                     const int is_family, const int num_elements, t8_element_t * elements[])
{
    if (is_family)
    {
        return -1;
    } else
    {
        return 0;
    }
}

int main(void)
{
#ifdef CMC_WITH_T8CODE
    /* Initialize cmc */
    cmc::CmcInitialize();

    const sc_MPI_Comm comm = sc_MPI_COMM_SELF;

    /* Create a mesh */
    t8_cmesh_t cmesh;
    cmesh = t8_cmesh_new_hypercube (T8_ECLASS_QUAD, comm, 0, 0, 0);
    t8_scheme_cxx_t* scheme = t8_scheme_new_default_cxx ();
    const int initial_level = 5;
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, initial_level, 0, comm);
    t8_eclass_scheme_c* ts = t8_forest_get_eclass_scheme(forest, t8_forest_get_eclass(forest, 0));

    /* We want to keep the initial forest */
    t8_forest_ref(forest);

    /* Create an AMR mesh */
    const int dimensionality = 2;
    cmc::AmrMesh mesh(forest, initial_level, dimensionality);

    const t8_locidx_t num_initial_elements = t8_forest_get_local_num_elements (forest);

    cmc::ExpectTrue(num_initial_elements == 1024);

    /* Get a certain element of the mesh */
    const t8_locidx_t elem_id1 = 100;
    t8_element_t* elem1 = t8_forest_get_element (forest, elem_id1, NULL);

    /* Get the intial element coverage of this element */
    const std::vector<cmc::DomainIndex> elem1_initial_coverage = GetInitialElementCoverage(mesh, ts, elem1);
    
    cmc::cmc_debug_msg("1) Num initial elem ids: ", elem1_initial_coverage.size());
    for (auto iter = elem1_initial_coverage.begin(); iter != elem1_initial_coverage.end(); ++iter)
    {
        cmc::cmc_debug_msg("Intial Elem Index: ", *iter);
    }

    cmc::ExpectTrue(elem1_initial_coverage.size() == 1);

    cmc::ExpectTrue(elem1_initial_coverage.front() == static_cast<cmc::DomainIndex>(elem_id1));

    /* Create a coarser forest */
    t8_forest_t forest_adapt = t8_forest_new_adapt (forest, coarsen_all_elements, 0, 0, NULL);

    cmc::ExpectTrue(t8_forest_get_local_num_elements(forest_adapt) == 256);

    /* Get a certain element of the mesh */
    const t8_locidx_t elem_id2 = 64;
    t8_element_t* elem2 = t8_forest_get_element (forest_adapt, elem_id2, NULL);

    /* Get the intial element coverage of this element */
    const std::vector<cmc::DomainIndex> elem2_initial_coverage = GetInitialElementCoverage(mesh, ts, elem2);

    cmc::cmc_debug_msg("2) Num initial elem ids: ", elem2_initial_coverage.size());
    for (auto iter = elem2_initial_coverage.begin(); iter != elem2_initial_coverage.end(); ++iter)
    {
        cmc::cmc_debug_msg("Intial Elem Index: ", *iter);
    }

    cmc::ExpectTrue(elem2_initial_coverage.size() == 4);

    cmc::ExpectTrue(elem2_initial_coverage[0] == 256);
    cmc::ExpectTrue(elem2_initial_coverage[1] == 257);
    cmc::ExpectTrue(elem2_initial_coverage[2] == 258);
    cmc::ExpectTrue(elem2_initial_coverage[3] == 259);

    /* Create an even coarser forest */
    forest_adapt = t8_forest_new_adapt (forest_adapt, coarsen_all_elements, 0, 0, NULL);

    cmc::ExpectTrue(t8_forest_get_local_num_elements(forest_adapt) == 64);

    /* Get a certain element of the mesh */
    const t8_locidx_t elem_id3 = 8;
    t8_element_t* elem3 = t8_forest_get_element (forest_adapt, elem_id3, NULL);

    /* Get the intial element coverage of this element */
    const std::vector<cmc::DomainIndex> elem3_initial_coverage = GetInitialElementCoverage(mesh, ts, elem3);

    cmc::cmc_debug_msg("3) Num initial elem ids: ", elem3_initial_coverage.size());
    for (auto iter = elem3_initial_coverage.begin(); iter != elem3_initial_coverage.end(); ++iter)
    {
        cmc::cmc_debug_msg("Intial Elem Index: ", *iter);
    }

    cmc::ExpectTrue(elem3_initial_coverage.size() == 16);

    cmc::ExpectTrue(elem3_initial_coverage[0] == 128);
    cmc::ExpectTrue(elem3_initial_coverage[1] == 129);
    cmc::ExpectTrue(elem3_initial_coverage[2] == 130);
    cmc::ExpectTrue(elem3_initial_coverage[3] == 131);
    cmc::ExpectTrue(elem3_initial_coverage[4] == 132);
    cmc::ExpectTrue(elem3_initial_coverage[5] == 133);
    cmc::ExpectTrue(elem3_initial_coverage[6] == 134);
    cmc::ExpectTrue(elem3_initial_coverage[7] == 135);
    cmc::ExpectTrue(elem3_initial_coverage[8] == 136);
    cmc::ExpectTrue(elem3_initial_coverage[9] == 137);
    cmc::ExpectTrue(elem3_initial_coverage[10] == 138);
    cmc::ExpectTrue(elem3_initial_coverage[11] == 139);
    cmc::ExpectTrue(elem3_initial_coverage[12] == 140);
    cmc::ExpectTrue(elem3_initial_coverage[13] == 141);
    cmc::ExpectTrue(elem3_initial_coverage[14] == 142);
    cmc::ExpectTrue(elem3_initial_coverage[15] == 143);

    /* Deallocate forest */
    t8_forest_unref (&forest_adapt);
    /* Finalize cmc */
    cmc::CmcFinalize();
#endif

    return cmc::CMC_TEST_SUCCESS;
}
