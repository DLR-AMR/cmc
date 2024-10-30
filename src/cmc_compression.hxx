#ifndef CMC_COMPRESSION_HXX
#define CMC_COMPRESSION_HXX

#include "cmc.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "lossy/cmc_amr_lossy_compression_settings.hxx"
#include "lossy/cmc_compression_class.hxx"
#include "lossy/cmc_prefix_lossy_compression.hxx"

#include <vector>
#include <memory>

namespace cmc
{

/*
Usage: 
Compressor(...) initialization

Compressor.Start();

Compressor.Compress();
//das in Compress Methode
while(data_is_available)
{
    Compressor.FeedData();

    Compressor.TransitionToNextTimeStep();
}

Compressor.Write();

*/
#if 1
class Compressor
{
public:
    Compressor() = delete;
    Compressor(std::unique_ptr<lossy::Compressor>&& compressor)
    : compressor_{std::move(compressor)} {};
    ~Compressor() = default;

    void Start();
    void Compress();
    void WriteCompressedData(const std::string& file_name);

private:
    //FeedData();
    //TransitionToNextTimeStep();

    std::vector<InputVar> input_data_;
    //std::vector<Variable> previous_state_;
    int timestep_{0};
    std::unique_ptr<lossy::Compressor> compressor_;
};

class Decompressor
{
public:
    Decompressor() = delete;
    Decompressor(const std::string& file_name)
    : file_name_{file_name} {};
    ~Decompressor() = default;

    void Start();
    void Decompress();
private:
    const std::string file_name_;
    std::unique_ptr<lossy::Decompressor> decompressor_;
};

#endif

//enum CompressionTechnique {};

//void Compress(std::vector<InputVar>&& input_variables, const CompressionTechnique compression_technique, const CompressionSettings& settings);

//TODO: maybe define some decompression settings, but maybe not
std::vector<OutputVar> Decompress(const std::string& file_name);

}

#endif /* !CMC_COMPRESSION_HXX */
