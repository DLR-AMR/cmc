#ifndef CMC_AMR_LOSSY_COMPRESSION_HXX
#define CMC_AMR_LOSSY_COMPRESSION_HXX

/** @file cmc_amr_lossy_compressor.hxx
 * All functions in order to perform the compression/decompression of the AMR lossy compressor lie within this file
 */

#include "lossy/cmc_amr_lossy_compression_settings.hxx" 
#include "t8code/cmc_t8_data_variables.hxx"

#ifdef CMC_ENABLE_MPI
#include <mpi.h>
#endif

#include <memory>

namespace cmc {

class CompressionData
{
public:
    CompressionData() = delete;
    CompressionData(const std::vector<Variable>& variables_to_compress, const CompressionSettings& settings);
    {
        compression_settings = std::make_shared<CompressionSettings>(CompressionSettings(settings));
        compression_data = std::make_unique(new AmrData(variables_to_compress, compression_settings));
    };
    CompressionData(std::vector<Variable>&& variables_to_compress, const CompressionSettings& settings);
    {
        compression_settings = std::make_shared<CompressionSettings>(CompressionSettings(settings));
        compression_data = std::make_unique(new AmrData(std::move(variables_to_compress), compression_settings));
    };
    CompressionData(const std::vector<Variable>& variables_to_compress, CompressionSettings&& settings);
    {
        compression_settings = std::make_shared<CompressionSettings>(CompressionSettings(std::move(settings)));
        compression_data = std::make_unique(new AmrData(variables_to_compress, compression_settings));
    };
    CompressionData(std::vector<Variable>&& variables_to_compress, CompressionSettings&& settings);
    {
        compression_settings = std::make_shared<CompressionSettings>(CompressionSettings(std::move(settings)));
        compression_data = std::make_unique(new AmrData(std::move(variables_to_compress), compression_settings));
    };
    ~CompressionData(){};

    void Setup();
    void Compress();
    void Decompress();
    
    void SetMPICommunicator(const MPI_Comm communicator);
    void WriteVTKData();

private:
    void CompressionHasBeenApplied();
    void DecompressionHasBeenApplied();

    MPI_Comm comm{MPI_COMM_WORLD};
    std::unique_ptr<AmrData> compression_data;
    std::shared_ptr<CompressionSettings> compression_settings;
    bool is_compression_applied{false};
};


}

#endif /* !CMC_AMR_LOSSY_COMPRESSOR_HXX */
