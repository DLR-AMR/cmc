#ifndef CMC_COMPRESSION_HXX
#define CMC_COMPRESSION_HXX

#include "cmc.hxx"
#include "utilities/cmc_input_variable.hxx"
#include "utilities/cmc_output_variable.hxx"
#include "lossy/cmc_amr_lossy_compression_settings.hxx"

#include <vector>

namespace cmc
{

//enum CompressionTechnique {};

//void Compress(std::vector<InputVar>&& input_variables, const CompressionTechnique compression_technique, const CompressionSettings& settings);

//TODO: maybe define some decompression settings, but maybe not
std::vector<OutputVar> Decompress(const std::string& file_name);

}

#endif /* !CMC_COMPRESSION_HXX */
