#ifndef CMC_FULL_AMR_COMPRESSION_HXX
#define CMC_FULL_AMR_COMPRESSION_HXX
//Write lossy compressor that does
//  - adaptive coaresening on real initial values (not estamation)
//  - LSB Bit Truncation/Toggeling with remaining error
//  - PrefixAMR DataDeduplication (with refinement indication bits due to adaptive coarsening)
//      -> Stores the prefix lengths with an arithmetic encoder
//  - Stores the suffix length as index of the last "one" (tail bit), which allows to discard this last "one" bit in the encoded stream

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "t8code/cmc_t8_data.hxx"
#include "mpi/cmc_mpi.hxx"
#include "lossy/cmc_compression_class.hxx"
#include  "utilities/cmc_bit_map.hxx"

#include <vector>

namespace cmc
{

namespace full_amr
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
    
    void Setup();
    void Compress();
    void WriteCompressedData(const std::string& file_name, const int time_step) const;

private:
    std::unique_ptr<AmrData> compression_data_{nullptr};
    
    CompressionSettings compression_settings_;

    std::vector<ByteVar> compression_variables_;

    std::vector<AdaptiveCoarseningIndications> ac_indications_;

    MPI_Comm comm_{MPI_COMM_WORLD};
    bool is_compression_applied_{false};
};


}

}


#endif /* !CMC_FULL_AMR_COMPRESSION_HXX */
