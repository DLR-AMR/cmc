#include "lossy/cmc_amr_lossy_compression.hxx"

namespace cmc
{

inline void
CompressionData::CompressionHasBeenApplied()
{
    is_compression_applied = true;
}

inline void
CompressionData::DecompressionHasBeenApplied()
{
    is_compression_applied = false;
}

inline void
CompressionData::Setup()
{
    //TODO: The setup may be skipped, if a forest and corresponding  variables have been sent via the constructor
    //Or if the variables are already assigned to an forest
    compression_data->BuildInitialMesh();

    compression_data->DistributeDataOnInitialMesh();
}

inline void
CompressionData::Compress()
{
    compression_data->CompressByAdaptiveCoarsening();

    CompressionHasBeenApplied();
}

inline void
CompressionData::SetMPICommunicator(const MPI_Comm communicator)
{
    comm = communicator;
}

void
CompressionData::Decompress()
{

}

void
CompressionData::WriteVTKData()
{


}

}
