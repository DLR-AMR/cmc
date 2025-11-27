#ifndef CMC_MPI_IO_HXX
#define CMC_MPI_IO_HXX

#include "mpi/cmc_mpi.hxx"
#include "cmc.hxx"
#include "utilities/cmc_hyperslab.hxx"
#include "utilities/cmc_geo_domain.hxx"

#include <vector>

namespace cmc
{

std::vector<Hyperslab>
DetermineSlicedDataDistribution(const Hyperslab& global_hyperslab, const MPI_Comm comm);

}

#endif /* CMC_MPI_IO_H */
