#include "t8code/cmc_t8_mpi.hxx"
#include "t8code/cmc_t8_morton.hxx"
#include "utilities/cmc_log_functions.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_schemes/t8_scheme.hxx> 
#include <t8_schemes/t8_default/t8_default.hxx>
#endif

namespace cmc
{

DataOffsets
GatherGlobalDataOffsets(const AmrMesh& mesh, const MPI_Comm comm)
{
    #ifdef CMC_ENABLE_MPI
    int ret_val, size, rank;

    ret_val = MPI_Comm_size(comm, &size);
    MPICheckError(ret_val);

    ret_val = MPI_Comm_rank(comm, &rank);
    MPICheckError(ret_val);

    const t8_forest_t forest = mesh.GetMesh();

    const t8_locidx_t ltree_id = 0; //The mesh is always built based on single tree
    const t8_eclass_t eclass = t8_forest_get_tree_class (forest, ltree_id);

    const t8_scheme_c* ts = t8_forest_get_scheme(forest);

    /* Get each process-local offset */
    MortonIndex elem_offset{-1}; 
    if (rank != 0)
    {
        elem_offset = static_cast<MortonIndex>(ts->element_get_linear_id(eclass, t8_forest_get_leaf_element_in_tree(forest, ltree_id, 0), mesh.GetInitialRefinementLevel()));
    }

    /* Define a vector capable o holding all offsets */
    DataOffsets offsets(static_cast<size_t>(size + 1));

    /* Gather all offsets */
    ret_val = MPI_Allgather(&elem_offset, 1, MPI_MORTON_INDEX_T, offsets.data(), 1, MPI_MORTON_INDEX_T, comm);
    MPICheckError(ret_val);

    /* Add the number of all global elements as a last entry */
    /* Add the possible upper bound for Morton indices as the last entry to the offset array */
    const DomainIndex last_domain_integer_coord = std::pow(static_cast<DomainIndex>(2), static_cast<DomainIndex>(mesh.GetInitialRefinementLevel())) - 1;
    std::vector<DomainIndex> first_outside_of_domain_index;
    std::fill_n(std::back_inserter(first_outside_of_domain_index), mesh.GetDimensionality(), last_domain_integer_coord);

    /* Choose the index directly after the last entry entry of the squared domain */
    offsets[size] = GetMortonIndex(first_outside_of_domain_index, mesh.GetDimensionality()) + 1;

    /* Store the communicator */
    offsets.SetMPIComm(comm);
    
    return offsets;

    #else
    /* Calculate the possible upper bound for Morton indices for the last entry to the offset array */
    const DomainIndex last_domain_integer_coord = std::pow(static_cast<DomainIndex>(2), static_cast<DomainIndex>(mesh.GetInitialRefinementLevel())) - 1;
    std::vector<DomainIndex> first_outside_of_domain_index;
    std::fill_n(std::back_inserter(first_outside_of_domain_index), mesh.GetDimensionality(), last_domain_integer_coord);

    /* Create an offset array for one process spanning the whole squared uniform intial domain */
    DataOffsets offsets(2);
    offsets[0] = -1;
    offsets[1] = GetMortonIndex(first_outside_of_domain_index, mesh.GetDimensionality()) + 1;
    return offsets;
    #endif
}

}
