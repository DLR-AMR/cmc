#ifndef CMC_T8_MPI_HXX
#define CMC_T8_MPI_HXX
/**
 * @file cmc_t8_morton.hxx
 */

#include "utilities/cmc_utilities.hxx"
#include "mpi/cmc_mpi.hxx"
#include "mpi/cmc_mpi_data.hxx"
#include "t8code/cmc_t8_mesh.hxx"

#include <iterator>

namespace cmc
{

DataOffsets
GatherGlobalDataOffsets(const AmrMesh& mesh, const MPI_Comm comm);

}

#endif /* !CMC_T8_MPI_HXX */
