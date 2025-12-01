#ifndef CMC_DECOMPRESSION_SERIAL_INPUT_HXX
#define CMC_DECOMPRESSION_SERIAL_INPUT_HXX

#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_geo_domain.hxx"
#include "utilities/cmc_compression_schema.hxx"
#include "compression_io/cmc_serial_io_util.hxx"

#include "patch/lossless/cmc_patch_prefix_extraction_plain_suffixes_decompression.hxx"
#include "patch/lossless/cmc_patch_multi_res_extraction_decompression.hxx"

#include <string>
#include <cstdio>
#include <cstdio>
#include <memory>
#include <vector>

namespace cmc::compression_io::serial
{

class Reader
{
public:
    Reader() = delete;
    ~Reader() = default;

    Reader(const Reader& other) = default;
    Reader& operator=(const Reader& other) = default;
    Reader(Reader&& other) = default;
    Reader& operator=(Reader&& other) = default;

    Reader(const std::string& file_prefix)
    : file_prefix_{file_prefix} {};

    template <typename T> std::unique_ptr<cmc::patch::IPatchDecompressionVariable<T>> ReadPatchVariableForDecompression(const std::string& var_name);

private:
    const std::string file_prefix_;
};

template <typename T>
std::unique_ptr<cmc::patch::IPatchDecompressionVariable<T>>
Reader::ReadPatchVariableForDecompression(const std::string& variable_name)
{
    /* Reconstruct the file name */
    const std::string file_name = CreateFileName(file_prefix_, variable_name);

    /* Read the encoded variable */
    std::FILE* file_in = std::fopen(file_name.c_str(), "rb");
    if (file_in == NULL)
    {
        cmc_err_msg("The file ", file_name, " does not exist and therfore cannot be read for decompression.");
    }
    /* Read the additional information first */
    std::vector<attr_type_t> info(kNumHeaderInfo);
    const size_t header_retrieved_bytes = std::fread(info.data(), sizeof(attr_type_t), kNumHeaderInfo, file_in);

    if (header_retrieved_bytes != kNumHeaderInfo)
    {
        cmc_err_msg("The number of inquired information from the compressed file is not as expecetd.");
    }

    /* Get the information from the header */
    const CmcType data_type = static_cast<CmcType>(GetValueFromByteStream<attr_type_t>(info.data()));
    if (ConvertToCmcType<T>() != data_type) {cmc_err_msg("The template parameter of the Reader functionality does not match the data_type of the compressed variable.");}

    const cmc::CompressionSchema compression_scheme = static_cast<cmc::CompressionSchema>(GetValueFromByteStream<attr_type_t>(info.data() + 1));
    const cmc::DataLayout initial_data_layout = static_cast<DataLayout>(GetValueFromByteStream<attr_type_t>(info.data() + 2));
    const attr_type_t dimensionality = GetValueFromByteStream<attr_type_t>(info.data() + 3);
    const attr_type_t num_compression_iterations = GetValueFromByteStream<attr_type_t>(info.data() + 4);
    const attr_type_t num_pyramidal_dim_length_lvls = GetValueFromByteStream<attr_type_t>(info.data() + 5);

    /* Allocat data for the pyramidal dimension lengths */
    std::vector<std::vector<size_t>> pyramidal_dim_lengths;
    std::vector<attr_type_t> dim_lengths(dimensionality);

    for (attr_type_t lvl_iter{0}; lvl_iter < num_pyramidal_dim_length_lvls; ++lvl_iter)
    {
        /* Allocate the next level */
        pyramidal_dim_lengths.emplace_back(dimensionality);

        /* Read the dimension lengths for the current level */
        const size_t num_dim_lengths_read = std::fread(dim_lengths.data(), sizeof(attr_type_t), static_cast<size_t>(dimensionality), file_in);
        
        if (static_cast<size_t>(dimensionality) != num_dim_lengths_read)
        {
            cmc_err_msg("The expeceted number of dimension lengths could not be read from the file.");
        }

        for (attr_type_t dim_iter{0}; dim_iter < dimensionality; ++dim_iter)
        {
            const attr_type_t dim_length = GetValueFromByteStream<attr_type_t>(dim_lengths.data() + dim_iter);
            pyramidal_dim_lengths.back()[dim_iter] = static_cast<size_t>(dim_length);
        }
    }

    /* Compute the bytes to read */
    const auto curr_pos = std::ftell(file_in);
    const int seek_end_ret_val = std::fseek(file_in, 0, SEEK_END);
    if (seek_end_ret_val != 0)
    {
        cmc_err_msg("Could not move to the end of the file.");
    }

    /* Compute the amount of bytes for the encoded stream */
    const auto num_bytes = std::ftell(file_in) - curr_pos;

    /* Move back to the current positions */
    const int seek_curr_ret_val = std::fseek(file_in, curr_pos, SEEK_SET);
    if (seek_curr_ret_val != 0)
    {
        cmc_err_msg("Could not move to the previous position of the file.");
    }

    /* Afterwards, we read the remaining bytes from the file */
    std::vector<uint8_t> encoded_data(num_bytes);

    auto read_bytes_encoded_stream = std::fread(encoded_data.data(), sizeof(uint8_t), num_bytes, file_in);

    if (read_bytes_encoded_stream != num_bytes)
    {
        cmc_err_msg("Could not read the encoded data from the compressed file.");
    }
    
    GeoDomain init_domain;
    const std::vector<Dimension> dims_vec = GetDimensionVectorFromLayout(static_cast<DataLayout>(initial_data_layout));
    cmc_assert(dims_vec.size() == pyramidal_dim_lengths.back().size());

    int dim_int_idx{0};
    for (auto dim_iter = dims_vec.begin(); dim_iter != dims_vec.end(); ++dim_iter, ++dim_int_idx)
    {
        DimensionInterval dim(*dim_iter, 0, pyramidal_dim_lengths.back()[dim_int_idx]);
        init_domain.UpdateDimension(dim);
    }

    /* Invoke the correct decompressor */
    switch (compression_scheme)
    {
        case cmc::CompressionSchema::PatchPrefixExtractionPlainSuffixes:
            switch (dimensionality)
            {
                case 1:
                    return std::make_unique<patch::decompression::prefix::plain_suffix::DecompressionVariable<T, 1>>(std::move(encoded_data), variable_name, std::move(pyramidal_dim_lengths), init_domain, initial_data_layout, num_pyramidal_dim_length_lvls);
                    break;
                case 2:
                    return std::make_unique<patch::decompression::prefix::plain_suffix::DecompressionVariable<T, 2>>(std::move(encoded_data), variable_name, std::move(pyramidal_dim_lengths), init_domain, initial_data_layout, num_pyramidal_dim_length_lvls);
                    break;
                case 3:
                    return std::make_unique<patch::decompression::prefix::plain_suffix::DecompressionVariable<T, 3>>(std::move(encoded_data), variable_name, std::move(pyramidal_dim_lengths), init_domain, initial_data_layout, num_pyramidal_dim_length_lvls);
                    break;
                default:
                    cmc_err_msg("The dimensionality of the compressed variable is not recognized for a patch variable.");
                    return std::make_unique<patch::decompression::prefix::plain_suffix::DecompressionVariable<T, 3>>(std::move(encoded_data), variable_name, std::move(pyramidal_dim_lengths), init_domain, initial_data_layout, num_pyramidal_dim_length_lvls);
            }
        break;
        case cmc::CompressionSchema::PatchMultiResExtraction:
            switch (dimensionality)
            {
                case 1:
                    return std::make_unique<patch::decompression::multi_res::DecompressionVariable<T, 1>>(std::move(encoded_data), variable_name, std::move(pyramidal_dim_lengths), init_domain, initial_data_layout, num_pyramidal_dim_length_lvls);
                    break;
                case 2:
                    return std::make_unique<patch::decompression::multi_res::DecompressionVariable<T, 2>>(std::move(encoded_data), variable_name, std::move(pyramidal_dim_lengths), init_domain, initial_data_layout, num_pyramidal_dim_length_lvls);
                    break;
                case 3:
                    return std::make_unique<patch::decompression::multi_res::DecompressionVariable<T, 3>>(std::move(encoded_data), variable_name, std::move(pyramidal_dim_lengths), init_domain, initial_data_layout, num_pyramidal_dim_length_lvls);
                    break;
                default:
                    cmc_err_msg("The dimensionality of the compressed variable is not recognized for a patch variable.");
                    return std::make_unique<patch::decompression::multi_res::DecompressionVariable<T, 3>>(std::move(encoded_data), variable_name, std::move(pyramidal_dim_lengths), init_domain, initial_data_layout, num_pyramidal_dim_length_lvls);
            }
        break;
        default:
            cmc_err_msg("The compression schema of the compressed variable is not recognized for a patch variable.");
            return std::make_unique<patch::decompression::prefix::plain_suffix::DecompressionVariable<T, 3>>(std::move(encoded_data), variable_name, std::move(pyramidal_dim_lengths), init_domain, initial_data_layout, num_pyramidal_dim_length_lvls);
    }
}

}

#endif /* !CMC_DECOMPRESSION_SERIAL_INPUT_HXX */
