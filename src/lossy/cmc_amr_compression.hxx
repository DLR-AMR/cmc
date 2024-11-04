#ifndef CMC_AMR_COMPRESSION_HXX
#define CMC_AMR_COMPRESSION_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "t8code/cmc_t8_data.hxx"
#include "mpi/cmc_mpi.hxx"
#include "lossy/cmc_compression_class.hxx"

#include <memory>

namespace cmc
{

namespace amr
{

class Compressor
{
public:
    Compressor() = default;

    Compressor(const std::vector<InputVar>& variables_to_compress, const CompressionSettings& settings)
    : compression_settings_{settings}
    {
        compression_data_ = std::make_unique<AmrData>(variables_to_compress, compression_settings_);
    };
    Compressor(std::vector<InputVar>&& variables_to_compress, const CompressionSettings& settings)
    : compression_settings_{settings}
    {
        compression_data_ = std::make_unique<AmrData>(std::move(variables_to_compress), compression_settings_);
    };
    Compressor(const std::vector<InputVar>& variables_to_compress, CompressionSettings&& settings)
    : compression_settings_(std::move(settings))
    {
        compression_data_ = std::make_unique<AmrData>(variables_to_compress, compression_settings_);
    };
    Compressor(std::vector<InputVar>&& variables_to_compress, CompressionSettings&& settings)
    : compression_settings_(std::move(settings))
    {
        compression_data_ = std::make_unique<AmrData>(std::move(variables_to_compress), compression_settings_);
    };

    Compressor(const Compressor& other) = default;
    Compressor& operator=(const Compressor& other) = default;
    Compressor(Compressor&& other) = default;
    Compressor& operator=(Compressor&& other) = default;

    ~Compressor() = default;

    int GetMpiSize() const;
    
    void Setup(const bool with_default_lossy_amr_compression = false);
    void Compress();
    void WriteCompressedData(const std::string& file_name, const int time_step) const;

private:
    std::unique_ptr<AmrData> compression_data_{nullptr};
    
    CompressionSettings compression_settings_;

    std::vector<ByteVar> compression_variables_;

    MPI_Comm comm_{MPI_COMM_WORLD};
    bool is_compression_applied_{false};


    bool perform_default_lossy_compression_{false};
};

}

}

#endif /* !CMC_AMR_COMPRESSION_HXX */
