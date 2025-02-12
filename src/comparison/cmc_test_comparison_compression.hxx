#ifndef CMC_LOSSLESS_TEST_COMPARISON_COMPRESSION_HXX
#define CMC_LOSSLESS_TEST_COMPARISON_COMPRESSION_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "t8code/cmc_t8_data.hxx"
#include "mpi/cmc_mpi.hxx"
#include "lossy/cmc_compression_class.hxx"

#include <memory>

namespace cmc
{

namespace test_comparison
{

namespace light_amr_pcp
{


class Compressor
{
public:
    Compressor() = default;

    Compressor(const std::vector<InputVar>& variables_to_compress)
    {
        compression_data_ = std::make_unique<AmrData>(variables_to_compress);
    };
    Compressor(std::vector<InputVar>&& variables_to_compress)
    {
        compression_data_ = std::make_unique<AmrData>(std::move(variables_to_compress));
    };

    Compressor(const Compressor& other) = default;
    Compressor& operator=(const Compressor& other) = default;
    Compressor(Compressor&& other) = default;
    Compressor& operator=(Compressor&& other) = default;

    ~Compressor() = default;

    int GetMpiSize() const;
    
    void SetSplitVariable(const SplitVariable& split_var);
    void SetSplitVariable(SplitVariable&& split_var);

    void Setup();
    void Compress();
    void WriteCompressedData(const std::string& file_name, const int time_step) const;

private:
    std::unique_ptr<AmrData> compression_data_{nullptr};

    std::vector<ByteVar> compression_variables_;

    std::vector<SplitVariable> split_variables_;

    MPI_Comm comm_{MPI_COMM_WORLD};
    bool is_compression_applied_{false};
};

}

}

}

#endif /* !CMC_LOSSLESS_TEST_COMPARISON_COMPRESSION_HXX */
