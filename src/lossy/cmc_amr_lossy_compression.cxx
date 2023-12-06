#include "lossy/cmc_amr_lossy_compression.hxx"

namespace cmc
{

inline void
CompressionData::CompressionHasBeenApplied()
{
    is_compression_applied_ = true;
}

inline void
CompressionData::DecompressionHasBeenApplied()
{
    is_compression_applied_ = false;
}

inline void
CompressionData::Setup()
{
    //TODO: The setup may be skipped, if a forest and corresponding  variables have been sent via the constructor
    //Or if the variables are already assigned to an forest
    compression_data_->BuildInitialMesh();

    compression_data_->DistributeDataOnInitialMesh();
}

inline void
CompressionData::Compress()
{
    compression_data_->CompressByAdaptiveCoarsening();

    CompressionHasBeenApplied();
}

inline void
CompressionData::SetMPICommunicator(const MPI_Comm communicator)
{
    comm_ = communicator;
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
