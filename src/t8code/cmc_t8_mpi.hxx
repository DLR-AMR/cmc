#ifndef CMC_T8_MPI_HXX
#define CMC_T8_MPI_HXX


#include "utilities/cmc_utilities.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#ifdef CMC_ENABLE_MPI
#include "mpi/cmc_mpi.hxx"
#include "mpi/cmc_mpi_data.hxx"
#endif

#include <iterator>

namespace cmc
{

DataOffsets
GatherGlobalDataOffsets(const AmrMesh& mesh, const MPI_Comm comm);

}

#endif /* !CMC_T8_MPI_HXX */
