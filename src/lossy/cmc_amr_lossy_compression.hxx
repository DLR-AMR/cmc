#ifndef CMC_AMR_LOSSY_COMPRESSION_HXX
#define CMC_AMR_LOSSY_COMPRESSION_HXX

/** @file cmc_amr_lossy_compressor.hxx
 * All functions in order to perform the compression/decompression of the AMR lossy compressor lie within this file
 */
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "t8code/cmc_t8_data.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "t8code/cmc_t8_data_variables_forward.hxx"
#include "utilities/cmc_output_variable_forward.hxx"

#ifdef CMC_ENABLE_MPI
#include <mpi.h>
#endif

#include <memory>

namespace cmc {

class CompressionData
{
public:
    CompressionData() = default;

    CompressionData(const std::vector<InputVar>& variables_to_compress, const CompressionSettings& settings)
    : compression_settings_{settings}
    {
        compression_data_ = std::make_unique<AmrData>(variables_to_compress, compression_settings_);
    };
    CompressionData(std::vector<InputVar>&& variables_to_compress, const CompressionSettings& settings)
    : compression_settings_{settings}
    {
        compression_data_ = std::make_unique<AmrData>(std::move(variables_to_compress), compression_settings_);
    };
    CompressionData(const std::vector<InputVar>& variables_to_compress, CompressionSettings&& settings)
    : compression_settings_(std::move(settings))
    {
        compression_data_ = std::make_unique<AmrData>(variables_to_compress, compression_settings_);
    };
    CompressionData(std::vector<InputVar>&& variables_to_compress, CompressionSettings&& settings)
    : compression_settings_(std::move(settings))
    {
        compression_data_ = std::make_unique<AmrData>(std::move(variables_to_compress), compression_settings_);
    };

    CompressionData(const CompressionData& other) = default;
    CompressionData& operator=(const CompressionData& other) = default;
    CompressionData(CompressionData&& other) = default;
    CompressionData& operator=(CompressionData&& other) = default;

    ~CompressionData() = default;

    void Setup();
    void Setup(const AmrMesh& mesh);
    void Compress(const CompressionMode compression_mode = CompressionMode::OneForOne);
    void SupplementarySZLikeCompression();
    //void Compress(const char* path_to_file);
    void Decompress();
    OutputVar DecompressVariable(const int variable_id);
    //void Decompress(const char* path_to_file);
    
    void SetMPICommunicator(const MPI_Comm communicator);
    void SetCompressionSettings(const CompressionSettings& settings);
    void SetCompressionSettings(CompressionSettings&& settings);
    //void WriteVTKData();

    void WriteCompressedData(const std::string& file_name) const;
    void WriteVTKFile(const std::string& file_name) const;

    bool IsValidForCompression() const;
    size_t GetNumberOfCompressionVariables() const;
    size_t GetNumberOfInputVariables() const;

    #if 0
    Var GetVariable(const int id) const
    {
        cmc_assert(compression_data_ != nullptr);
        return compression_data_->GetVariable(id);
    };
    Var&& GetVariable(const int id)
    {
        cmc_assert(compression_data_ != nullptr);
        return std::move(compression_data_->GetVariable(id));
    };
    #endif

    OutputVar SeizeDecompressedVariable(const int variable_id);
    
private:
    //friend class TestLossyCompression;
    
    void CompressionHasBeenApplied();
    void DecompressionHasBeenApplied();


    //std::unique_ptr<AmrData> compression_data_;
    std::unique_ptr<AmrData> compression_data_{nullptr};
    CompressionSettings compression_settings_;
    CompressionMode compression_mode_{CompressionMode::CompressionModeUndefined};
    MPI_Comm comm_{MPI_COMM_WORLD};
    bool is_compression_applied_{false};

    #if 0
    MPI_Comm comm_{MPI_COMM_WORLD};
    std::unique_ptr<AmrData> compression_data_;
    std::shared_ptr<CompressionSettings> compression_settings_;
    CompressionMode compression_mode_{CompressionMode::CompressionModeUndefined};
    bool is_compression_applied_{false};
    #endif
};


}

#endif /* !CMC_AMR_LOSSY_COMPRESSOR_HXX */
