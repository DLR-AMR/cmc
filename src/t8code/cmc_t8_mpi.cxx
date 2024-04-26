#include "t8code/cmc_t8_mpi.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_schemes/t8_default/t8_default_c_interface.h>
#include <t8_element_c_interface.h>
#endif

namespace cmc
{

DataOffsets
GatherGlobalDataOffsets(const AmrMesh mesh, const MPI_Comm comm)
{
    #ifdef CMC_ENABLE_MPI

    int ret_val, size, rank;

    ret_val = MPI_Comm_size(comm, &size);
    MPICheckError(ret_val);

    ret_val = MPI_Comm_rank(comm, &rank);
    MPICheckError(ret_val);

    const t8_locidx_t ltree_id = 0; //The mesh is always built based on single tree
    const t8_forest_t forest = mesh.GetMesh();

    t8_eclass_scheme_c* ts = t8_forest_get_eclass_scheme(forest, t8_forest_get_eclass(forest, ltree_id));

    const MortonIndex elem_offset = static_cast<MortonIndex>(t8_element_get_linear_id (ts,
                                                                                       t8_forest_get_element_in_tree(forest, 0, 0),
                                                                                       mesh.GetInitialRefinementLevel()));

    DataOffsets offsets(static_cast<size_t>(size));

    ret_val = MPI_Allgather(&elem_offset, 1, MPI_MORTON_INDEX_T, offsets.data(), 1, MPI_MORTON_INDEX_T, comm);
    MPICheckError(ret_val);

    return offsets;

    #else
    DataOffsets offsets(1);
    offsets[0] = 0;
    return offsets;
    #endif
}

}
