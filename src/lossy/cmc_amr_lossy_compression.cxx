#include "lossy/cmc_amr_lossy_compression.hxx"
#include "utilities/cmc_output_variable.hxx"

namespace cmc
{

void
CompressionData::CompressionHasBeenApplied()
{
    is_compression_applied_ = true;
}

void
CompressionData::DecompressionHasBeenApplied()
{
    is_compression_applied_ = false;
}


void
CompressionData::Setup()
{
    compression_data_->SplitVariables();

    compression_data_->CheckConsistencyOfInputVariables();

    compression_data_->BuildInitialMesh();

    compression_data_->DistributeDataOnInitialMesh();//TODO: implement//serial for now

    compression_data_->ApplyScalingAndOffset();

    compression_data_->SetupVariablesForCompression();
}

void
CompressionData::Compress(const CompressionMode compression_mode)
{
    compression_data_->CompressByAdaptiveCoarsening(compression_mode);

    CompressionHasBeenApplied();
}

void
CompressionData::WriteVTKFile(const std::string& file_name) const
{
    compression_data_->WriteVTKFilePerVariable(file_name);
}

void
CompressionData::SetMPICommunicator(const MPI_Comm communicator)
{
    comm_ = communicator;
}

void
CompressionData::SetCompressionSettings(const CompressionSettings& settings)
{
    compression_settings_ = settings;
}

void
CompressionData::SetCompressionSettings(CompressionSettings&& settings)
{
    compression_settings_ = std::move(settings);
}

void
CompressionData::Decompress()
{
    cmc_assert(is_compression_applied_);

    compression_data_->DecompressToInitialRefinementLevel();    

    DecompressionHasBeenApplied();

}

OutputVar
CompressionData::DecompressVariable(const int variable_id)
{
    return compression_data_->DecompressVariable(variable_id);
}

#if 0
OutputVar
CompressionData::SeizeDecompressedVariable(const int variable_id)
{

}
#endif
#if 0
void
CompressionData::WriteVTKData()
{


}
#endif

bool
CompressionData::IsValidForCompression() const
{
    if (is_compression_applied_)
        return false;
    if (!compression_settings_.AreTheSettingsValid())
        return false;
    if (compression_data_ == nullptr)
        return false;
    if (!compression_data_->IsValidForCompression())
        return false;
    
    return true;
}

size_t
CompressionData::GetNumberOfInputVariables() const
{
    if (compression_data_ == nullptr)
    {
        cmc_err_msg("No input variables have been supplied.");
        return CMC_ERR;
    }

    return compression_data_->GetNumberOfInputVariables();
}

size_t
CompressionData::GetNumberOfCompressionVariables() const
{
    if (compression_data_ == nullptr)
    {
        cmc_err_msg("No input variables that could have been transformed to compression variables were supplied.");
        return CMC_ERR;
    }

    return compression_data_->GetNumberOfCompressionVariables();
}


void
CompressionData::WriteCompressedData(const std::string& file_name) const
{
    if (!is_compression_applied_)
    {
        cmc_warn_msg("The data has not been compressed. Therefore, it cannot be written to disk.");
        return;
    }

    compression_data_->WriteCompressedData(file_name);
}


void
CompressionData::SupplementarySZLikeCompression()
{
    if (!is_compression_applied_)
    {
        cmc_err_msg("Supplementary compression with SZ could not be applied.");
        return;
    }

    sz_compression_data_ = std::make_unique<SZCompressor>(compression_data_->GetCompressionVariables());

    sz_compression_data_->Compress();
}

}
