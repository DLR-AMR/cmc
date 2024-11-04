#include "lossy/cmc_amr_lossy_compression.hxx"
#include "utilities/cmc_output_variable.hxx"

#include "netcdf/cmc_netcdf.hxx"

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

    compression_data_->DistributeDataOnInitialMesh();

    compression_data_->ApplyScalingAndOffset();

    compression_data_->SetupVariablesForCompression();
   //cmc_err_msg("Hier ist ende");
}

void
CompressionData::Setup(const AmrMesh& mesh)
{
    //Check consistency of input variables and mesh
    compression_data_->SetInitialMesh(mesh);

    compression_data_->ApplyScalingAndOffset();
    
    compression_data_->SetupVariablesForCompression();
}

void
CompressionData::Compress(const CompressionMode compression_mode)
{
    compression_data_->CompressByAdaptiveCoarsening(compression_mode);

    CompressionHasBeenApplied();

    #if 0
    //Test session for now
    OutputVar decompressed_var = compression_data_->DecompressVariable(5);
    OutputVariable<float> var = decompressed_var.SeizeOutputVariable<float>();
    std::vector<float> data = var.GetData();
    cmc_debug_msg("Size of data of decompressed variable: ", data.size());

    int err;
    int ncid;

    /* Open the file without explicit parallel access */
    err = nc__open("../../data/era5_reanalysis_pressure_lvls_fixed_time.nc", NC_NOWRITE, NULL, &ncid);
    NcCheckError(err);

    std::vector<short> init_data(data.size(), 0);

    err = nc_get_var_short(ncid, 5, init_data.data());
    NcCheckError(err);

    //short max_err = 0;

    std::vector<float> fdata;
    fdata.reserve(data.size());

    const float scale = 0.00200295932540368;
    const float offset = 250.108564255201;
    for (auto siter = init_data.begin(); siter != init_data.end(); ++siter)
    {
        fdata.push_back( scale * static_cast<float>(*siter) + offset);
    }

    float max_err = 0.0;
    for (size_t i = 0; i < data.size(); ++i)
    {
        float err = std::abs(data[i] - fdata[i]);
        if (err > 1.0)
        {
            cmc_debug_msg("Errro ggreater then limit: ", err);
        }
        if (err > max_err)
        {
            max_err = err;
        }
    }
    cmc_debug_msg("Max error after decompression: ", max_err);

    #endif

    #if 0
    //Test session for now Rel Error
    OutputVar decompressed_var = compression_data_->DecompressVariable(4);
    OutputVariable<float> var = decompressed_var.SeizeOutputVariable<float>();
    std::vector<float> data = var.GetData();
    cmc_debug_msg("Size of data of decompressed variable: ", data.size());

    int err;
    int ncid;

    /* Open the file without explicit parallel access */
    err = nc__open("../../data/era5_reanalysis_t2m_tc03_13_12_23.nc", NC_NOWRITE, NULL, &ncid);
    NcCheckError(err);

    std::vector<short> init_data(data.size(), 0);

    err = nc_get_var_short(ncid, 4, init_data.data());
    NcCheckError(err);

    std::vector<float> fdata;
    fdata.reserve(data.size());

    const float scale = 9.25404335362387E-08;
    const float offset = 0.00790779508439177;
    for (auto siter = init_data.begin(); siter != init_data.end(); ++siter)
    {
        fdata.push_back( scale * static_cast<float>(*siter) + offset);
    }

    float max_err = 0.0;
    for (size_t i = 0; i < data.size(); ++i)
    {
        float err = std::abs(data[i] - fdata[i]);

        float rel_err = std::abs(fdata[i] - data[i]) / std::abs(fdata[i]);
        if (rel_err > 1.0)
        {
            cmc_debug_msg("Errro ggreater then limit: ", rel_err);
        }
        if (rel_err > max_err)
        {
            max_err = rel_err;
        }
    }
    cmc_debug_msg("Max error after decompression: ", max_err);
    #endif
    //WriteCompressedData("test_new_rel_err_era5_tco3.nc");
    //Decompress();
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

}
