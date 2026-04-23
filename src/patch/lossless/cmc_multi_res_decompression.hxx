#ifndef CMC_PATCH_LOSSLESS_CMC_MULTI_RES_DECOMPRESSION_HXX
#define CMC_PATCH_LOSSLESS_CMC_MULTI_RES_DECOMPRESSION_HXX

#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "patch/lossless/cmc_patch_multi_res_extraction_util.hxx"
#include "utilities/cmc_compression_schema.hxx"

#include <stdexcept>
#include <filesystem>


#include <bitset>

namespace cmc::serial::patch::lossless::multi_res
{

struct CompressionInfo
{
    SizeType global_byte_count;
    SizeType encoding_start_byte_count;
    SizeType data_type;
    SizeType dimensionality;
    SizeType compression_scheme;
    SizeType num_compression_levels;
    std::vector<SizeType> global_level_byte_count;  
    std::vector<std::vector<int>> level_dim_lengths;
    SizeType num_bytes_serialized_dictionary;
    std::vector<SizeType> serialized_entropy_dictionary;
    std::vector<SizeType> variable_encoding;
};

inline CompressionInfo
ReadCompressionData(const std::string& file_name)
{
    /* Check if the output file exists */
    const std::filesystem::path output_file_path(file_name);
    if (not std::filesystem::exists(output_file_path))
    {
        cmc::cmc_err_msg("The compression output file ", output_file_path, " does not exist.");
    }

    std::FILE* file_in = std::fopen(file_name.c_str(), "rb");

    /* Check the file size */
    const auto file_size = std::filesystem::file_size(output_file_path);

    /* Read the first two SizeTypes from the file */
    std::array<uint64_t, 2> num_bytes{};
    std::fread(num_bytes.data(), sizeof(uint64_t), 2, file_in);
    
    const SizeType global_byte_count = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(num_bytes[0]);
    const SizeType preamble_byte_count = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(num_bytes[1]);

    if (preamble_byte_count <= 2 || preamble_byte_count > 8192)  {cmc_err_msg("The file does not hold compression information!");}

    /* Read the preamble */
    const int remaining_preamble_size = preamble_byte_count - 2 * sizeof(uint64_t);

    cmc_assert(remaining_preamble_size % sizeof(uint64_t) == 0);

    if (remaining_preamble_size % sizeof(uint64_t) != 0) {cmc_err_msg("The data in the file has an unexpected offset and therefore cannot be read!");}

    CompressionInfo compression_info;
    compression_info.global_byte_count = global_byte_count;
    compression_info.encoding_start_byte_count = preamble_byte_count;

    /* Read the remaining preamble bytes */
    std::vector<SizeType> preamble(remaining_preamble_size / sizeof(uint64_t));
    std::fread(preamble.data(), sizeof(uint64_t), remaining_preamble_size / sizeof(uint64_t), file_in);

    int offset = 0;
    /* Read the succeeding information from the file */
    compression_info.data_type = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(preamble[offset]);
    ++offset;

    compression_info.dimensionality = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(preamble[offset]);
    ++offset;

    if (compression_info.dimensionality <= 0 || compression_info.dimensionality > 4){cmc_err_msg("The dimensionality is not supported!");}

    compression_info.compression_scheme = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(preamble[offset]);
    ++offset;

    if (CompressionSchema::PatchMultiResExtraction != static_cast<CompressionSchema>(compression_info.compression_scheme)){cmc_err_msg("The compression scheme does not match!");}

    compression_info.num_compression_levels = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(preamble[offset]);
    ++offset;

    const int num_encdoing_levels = compression_info.num_compression_levels + 1;

    compression_info.global_level_byte_count.reserve(compression_info.num_compression_levels);
    for (SizeType idx{0}; idx < compression_info.num_compression_levels; ++idx)
    {
        compression_info.global_level_byte_count.push_back(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(preamble[offset]));
        ++offset;
    }

    compression_info.level_dim_lengths.reserve(num_encdoing_levels);
    for (SizeType idx{0}; idx < compression_info.num_compression_levels; ++idx)
    {
        compression_info.level_dim_lengths.emplace_back(compression_info.dimensionality, 0);
        for (SizeType dim_idx{0}; dim_idx < compression_info.dimensionality; ++dim_idx)
        {
            compression_info.level_dim_lengths.back()[dim_idx] = static_cast<int>(cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(preamble[offset]));
            ++offset;
        }
    }

    compression_info.num_bytes_serialized_dictionary = cmc::bits::ConvertBigEndianToNativeEndianness<SizeType>(preamble[offset]);
    ++offset;

    cmc_assert(compression_info.num_bytes_serialized_dictionary % sizeof(uint64_t) == 0);
    const int num_vals_serialized_entropy_dictionary = compression_info.num_bytes_serialized_dictionary / sizeof(uint64_t);

    compression_info.serialized_entropy_dictionary.reserve(num_vals_serialized_entropy_dictionary);
    std::copy_n(&preamble[offset], num_vals_serialized_entropy_dictionary, std::back_inserter(compression_info.serialized_entropy_dictionary));
    offset += num_vals_serialized_entropy_dictionary;

    /* Get the encoded level data */
    cmc_assert(global_byte_count <= file_size && global_byte_count > preamble_byte_count);
    cmc_assert((global_byte_count - preamble_byte_count) % sizeof(uint64_t) == 0);
    const SizeType encoded_level_values = (global_byte_count - preamble_byte_count) / sizeof(uint64_t);

    /* Read the compressed data from */
    std::vector<SizeType> encoded_data(encoded_level_values);
    const std::size_t encoded_vals_read = std::fread(encoded_data.data(), sizeof(uint64_t), encoded_level_values, file_in);
    if (static_cast<SizeType>(encoded_vals_read) != encoded_level_values)
    {
        cmc_err_msg("The encoded data has not been read correctly from the file!");
    }
    compression_info.variable_encoding = std::move(encoded_data);

    /* Close the file */
    std::fclose(file_in);

    return compression_info;
}


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
class DecompressionVariable
{
public:
    DecompressionVariable() = delete;
    DecompressionVariable(const std::string compressed_file_name)
    {
        CompressionInfo cr_info = ReadCompressionData(compressed_file_name);

        if (static_cast<int32_t>(cr_info.dimensionality) != DIM) {cmc_err_msg("The dimensionality of the compressed data (", cr_info.dimensionality,
                                                                              ") does and the decompression variable (", DIM, ") does not match!");}
        this->bytes_per_level_ = std::move(cr_info.global_level_byte_count);
        this->dim_lengths_ = std::move(cr_info.level_dim_lengths);
        this->encoded_data_ = std::move(cr_info.variable_encoding);
        this->num_decompression_steps_ = cr_info.num_compression_levels;
        this->stream_decoder_.StartHuffmanCodesDecoding(cr_info.serialized_entropy_dictionary.data());
    }

    void Decompress();
    void WriteData(const std::string file_name);
    std::vector<T> GetDecompressedData() const;
private:

    std::vector<T> data_;

    std::vector<uint64_t> bytes_per_level_;
    std::vector<std::vector<int>> dim_lengths_;
    std::vector<uint64_t> encoded_data_;
    cmc::bits::vector_view level_encoding_;
    cmc::bits::StreamDecoder<SymbolType> stream_decoder_;

    int num_decompression_steps_{0};
    int decompression_step_{0};

    bool has_been_decompressed_{false};

    static struct kTag4D{} tag4D;
    static struct kTag3D{} tag3D;
    static struct kTag2D{} tag2D;
    static struct kTag1D{} tag1D;

    void PerformDecompression(kTag1D);
    void PerformDecompression(kTag2D);
    void PerformDecompression(kTag3D);
    void PerformDecompression(kTag4D);
};


template <FourByteArithmeticType T>
inline T
RefineValue(cmc::bits::StreamDecoder<SymbolType>& stream_decoder, const T coarse_predictor)
{
    /* Decode the next entropy symbol */
    const SymbolType entropy_symbol = stream_decoder.DecodeNextEntropySymbol();

    /* Determine the leading zero count */
    const int lzc = GetLZCFromEntropySymbol(entropy_symbol);

    if (lzc < static_cast<int>(sizeof(T) * cmc::bits::kCharBit) - 1) [[likely]]
    {
        /* Compute the length of the significant residual bits */
        const int residual_length = sizeof(T) * cmc::bits::kCharBit - 1 - lzc;

        /* We obtain the residual and add the implicit one bit  */
        const FourByteResidualType residual = stream_decoder.GetNextBitSequence<FourByteResidualType>(residual_length) | (FourByteResidualType{1} << residual_length);

        /* We create the residual applied value */
        if (IsApproximationGreater(entropy_symbol))
        {
            /* Subtract the residual from the prediction */
            const FourByteResidualType value = cmc::bits::IntegerSubtraction(coarse_predictor, residual);

            return std::bit_cast<T>(value);
        } else
        {
            /* Add the residual to the prediction */
            const FourByteResidualType value = cmc::bits::IntegerAddition(coarse_predictor, residual);

            return std::bit_cast<T>(value);
        }
    } else
    {
        if (lzc == sizeof(T) * cmc::bits::kCharBit) [[likely]]
        {
            /* It is a zero residual */
            return coarse_predictor;
        }

        /* Compute the residual in case we do not need to extract a bit-sequence */
        constexpr FourByteResidualType residual{0x00000001};

        /* We create the residual applied value */
        if (IsApproximationGreater(entropy_symbol))
        {
            /* Subtract the residual from the prediction */
            const FourByteResidualType value = cmc::bits::IntegerSubtraction(coarse_predictor, residual);

            return std::bit_cast<T>(value);
        } else
        {
            /* Add the residual to the prediction */
            const FourByteResidualType value = cmc::bits::IntegerAddition(coarse_predictor, residual);

            return std::bit_cast<T>(value);
        }
    }
}

template<OneByteArithmeticType T>
inline std::vector<T>
DecodeRootLevelValue(cmc::bits::vector_view root_level_view)
{
    /* Get the not-encoded root level value */
    const T root_value = std::bit_cast<T>(root_level_view.GetNextBitSequence<OneByteResidualType>(sizeof(T) * cmc::bits::kCharBit));

    std::vector<T> root_data;
    root_data.push_back(root_value);

    return root_data;
}

template<TwoByteArithmeticType T>
inline std::vector<T>
DecodeRootLevelValue(cmc::bits::vector_view root_level_view)
{
    /* Get the not-encoded root level value */
    const T root_value = std::bit_cast<T>(root_level_view.GetNextBitSequence<TwoByteResidualType>(sizeof(T) * cmc::bits::kCharBit));

    std::vector<T> root_data;
    root_data.push_back(root_value);

    return root_data;
}

template<FourByteArithmeticType T>
inline std::vector<T>
DecodeRootLevelValue(cmc::bits::vector_view root_level_view)
{
    /* Get the not-encoded root level value */
    const T root_value = std::bit_cast<T>(root_level_view.GetNextBitSequence<FourByteResidualType>(sizeof(T) * cmc::bits::kCharBit));

    std::vector<T> root_data;
    root_data.push_back(root_value);

    return root_data;
}

template<EightByteArithmeticType T>
inline std::vector<T>
DecodeRootLevelValue(cmc::bits::vector_view root_level_view)
{
    /* Get the not-encoded root level value */
    const T root_value = std::bit_cast<T>(root_level_view.GetNextBitSequence<EightByteResidualType>(sizeof(T) * cmc::bits::kCharBit));

    std::vector<T> root_data;
    root_data.push_back(root_value);

    return root_data;
}


/** 2D Decompression **/
template <typename T>
inline T
GetCoarseValue(const std::vector<T>& coarse_data, const std::vector<int>& coarse_lvl_dims, const int next_lvl_lat, const int next_lvl_lon)
{
    cmc_assert(coarse_lvl_dims.size() == static_cast<size_t>(2));

    constexpr int kLatID = 0;
    constexpr int kLonID = 1;

    const int coarse_lvl_lat = next_lvl_lat / kDimReductionFactor;
    const int coarse_lvl_lon = next_lvl_lon / kDimReductionFactor;

    cmc_assert(coarse_lvl_lat < coarse_lvl_dims[kLatID]);
    cmc_assert(coarse_lvl_lon < coarse_lvl_dims[kLonID]);
    cmc_assert(coarse_lvl_lat * coarse_lvl_dims[kLonID] + coarse_lvl_lon < coarse_data.size());

    return coarse_data[coarse_lvl_lat * coarse_lvl_dims[kLonID] + coarse_lvl_lon];
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::PerformDecompression(DecompressionVariable<T, DIM>::kTag2D)
{
    cmc_debug_msg("Decompression of the variable starts...");

    [[maybe_unused]] constexpr int kDim = 2;
    constexpr int kLatID = 0;
    constexpr int kLonID = 1;

    /* Get a view on the root level */
    cmc::bits::vector_view root_level_encoding(this->encoded_data_.data());

    /* Decode the root level value */
    this->data_ = DecodeRootLevelValue<T>(root_level_encoding);

    cmc_assert(data_.size() == static_cast<size_t>(1));

    cmc_assert(this->bytes_per_level_[this->decompression_step_] % sizeof(uint64_t) == 0);

    size_t level_offset{this->bytes_per_level_[this->decompression_step_] / sizeof(uint64_t)};
    ++(this->decompression_step_);

    /* Perform the iterative decompression */
    for (int decompression_step_idx{0}; decompression_step_idx < this->num_decompression_steps_ -1; ++decompression_step_idx)
    {
        /* Set a view to this level's encoding */
        cmc::bits::vector_view level_encoding(this->encoded_data_.data() + level_offset);

        /* Set the view accordingly in the stream decoder */
        this->stream_decoder_.StartDecoding(level_encoding);

        /* Get the dimensions of this and the next level */
        const std::vector<int>& curr_lvl_dims = this->dim_lengths_[decompression_step_idx];
        const std::vector<int>& next_lvl_dims = this->dim_lengths_[decompression_step_idx + 1];

        /* Allocate a refined data vector */
        std::vector<T> data_new;
        data_new.reserve(next_lvl_dims[kLatID] * next_lvl_dims[kLonID]);

        /* Iterate over the finer level */
        for (int lat = 0; lat < next_lvl_dims[kLatID]; ++lat)
        {
            for (int lon = 0; lon < next_lvl_dims[kLonID]; ++lon)
            {
                /* Get the corresponding coarse level value */
                const T coarse_val = GetCoarseValue<T>(this->data_, curr_lvl_dims, lat, lon);
                
                /* Refine the value */
                const T fine_val = RefineValue<T>(this->stream_decoder_, coarse_val);

                /* Store the refined value */
                data_new.push_back(fine_val);
            }
        }

        /* Switch to the new refined data for the next iteration */
        data_ = std::move(data_new);

        /* Update the offset in the encoded data */
        level_offset += this->bytes_per_level_[this->decompression_step_] / sizeof(uint64_t);

        /* Update the decompression step */
        ++(this->decompression_step_);
    }

    has_been_decompressed_ = true;
    cmc_debug_msg("Decompression of variable is completed.");
}



/** 3D Compression **/
template <typename T>
inline T
GetCoarseValue(const std::vector<T>& coarse_data, const std::vector<int>& coarse_lvl_dims, const int next_lvl_lev, const int next_lvl_lat, const int next_lvl_lon)
{
    cmc_assert(coarse_lvl_dims.size() == static_cast<size_t>(3));
    constexpr int kLevID = 0;
    constexpr int kLatID = 1;
    constexpr int kLonID = 2;

    const int coarse_lvl_lev = next_lvl_lev / kDimReductionFactor;
    const int coarse_lvl_lat = next_lvl_lat / kDimReductionFactor;
    const int coarse_lvl_lon = next_lvl_lon / kDimReductionFactor;

    cmc_assert(coarse_lvl_lev < static_cast<int>(coarse_lvl_dims[kLevID]));
    cmc_assert(coarse_lvl_lat < static_cast<int>(coarse_lvl_dims[kLatID]));
    cmc_assert(coarse_lvl_lon < static_cast<int>(coarse_lvl_dims[kLonID]));
    cmc_assert(coarse_lvl_lev * coarse_lvl_dims[kLatID] * coarse_lvl_dims[kLonID] + coarse_lvl_lat * coarse_lvl_dims[kLonID] + coarse_lvl_lon < static_cast<int>(coarse_data.size()));

    return coarse_data[coarse_lvl_lev * coarse_lvl_dims[kLatID] * coarse_lvl_dims[kLonID] + coarse_lvl_lat * coarse_lvl_dims[kLonID] + coarse_lvl_lon];
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::PerformDecompression(DecompressionVariable<T, DIM>::kTag3D)
{
    cmc_debug_msg("Decompression of the variable starts...");

    [[maybe_unused]] constexpr int kDim = 3;
    constexpr int kLevID = 0;
    constexpr int kLatID = 1;
    constexpr int kLonID = 2;

    /* Get a view on the root level */
    cmc::bits::vector_view root_level_encoding(this->encoded_data_.data());

    /* Decode the root level value */
    this->data_ = DecodeRootLevelValue<T>(root_level_encoding);

    cmc_assert(data_.size() == static_cast<size_t>(1));

    cmc_assert(this->bytes_per_level_[this->decompression_step_] % sizeof(uint64_t) == 0);

    size_t level_offset{this->bytes_per_level_[this->decompression_step_] / sizeof(uint64_t)};
    ++(this->decompression_step_);

    /* Perform the iterative decompression */
    for (int decompression_step_idx{0}; decompression_step_idx < this->num_decompression_steps_ - 1; ++decompression_step_idx)
    {
        /* Set a view to this level's encoding */
        cmc::bits::vector_view level_encoding(this->encoded_data_.data() + level_offset);

        /* Set the view accordingly in the stream decoder */
        this->stream_decoder_.StartDecoding(level_encoding);

        /* Get the dimensions of this and the next level */
        const std::vector<int>& curr_lvl_dims = this->dim_lengths_[decompression_step_idx];
        const std::vector<int>& next_lvl_dims = this->dim_lengths_[decompression_step_idx + 1];

        /* Allocate a refined data vector */
        std::vector<T> data_new;
        data_new.reserve(next_lvl_dims[kLevID] * next_lvl_dims[kLatID] * next_lvl_dims[kLonID]);

        /* Iterate over the finer level */
        for (int lev = 0; lev < next_lvl_dims[kLevID]; ++lev)
        {
            for (int lat = 0; lat < next_lvl_dims[kLatID]; ++lat)
            {
                for (int lon = 0; lon < next_lvl_dims[kLonID]; ++lon)
                {
                    /* Get the corresponding coarse level value */
                    const T coarse_val = GetCoarseValue<T>(this->data_, curr_lvl_dims, lev, lat, lon);

                    /* Refine the value */
                    const T fine_val = RefineValue<T>(this->stream_decoder_, coarse_val);

                    /* Store the refined value */
                    data_new.push_back(fine_val);
                }
            }
        }

        /* Switch to the new refined data for the next iteration */
        data_ = std::move(data_new);

        /* Update the offset in the encoded data */
        level_offset += this->bytes_per_level_[this->decompression_step_] / sizeof(uint64_t);

        /* Update the decompression step */
        ++(this->decompression_step_);
    }

    has_been_decompressed_ = true;
    cmc_debug_msg("Decompression of variable is completed.");
}


/** 4D Compression **/
template <typename T>
inline T
GetCoarseValue(const std::vector<T>& coarse_data, const std::vector<int>& coarse_lvl_dims, const int next_lvl_time, const int next_lvl_lev, const int next_lvl_lat, const int next_lvl_lon)
{
    cmc_assert(coarse_lvl_dims.size() == static_cast<size_t>(4));
    constexpr int kTimeID = 0;
    constexpr int kLevID = 1;
    constexpr int kLatID = 2;
    constexpr int kLonID = 3;

    const int coarse_lvl_time = next_lvl_time / kDimReductionFactor;
    const int coarse_lvl_lev = next_lvl_lev / kDimReductionFactor;
    const int coarse_lvl_lat = next_lvl_lat / kDimReductionFactor;
    const int coarse_lvl_lon = next_lvl_lon / kDimReductionFactor;

    cmc_assert(coarse_lvl_time < coarse_lvl_dims[kTimeID]);
    cmc_assert(coarse_lvl_lev < coarse_lvl_dims[kLevID]);
    cmc_assert(coarse_lvl_lat < coarse_lvl_dims[kLatID]);
    cmc_assert(coarse_lvl_lon < coarse_lvl_dims[kLonID]);
    cmc_assert(coarse_lvl_time * coarse_lvl_dims[kLevID] * coarse_lvl_dims[kLatID] * coarse_lvl_dims[kLonID] + coarse_lvl_lev * coarse_lvl_dims[kLatID] * coarse_lvl_dims[kLonID] + coarse_lvl_lat * coarse_lvl_dims[kLonID] + coarse_lvl_lon < coarse_data.size());

    return coarse_data[coarse_lvl_time * coarse_lvl_dims[kLevID] * coarse_lvl_dims[kLatID] * coarse_lvl_dims[kLonID] + coarse_lvl_lev * coarse_lvl_dims[kLatID] * coarse_lvl_dims[kLonID] + coarse_lvl_lat * coarse_lvl_dims[kLonID] + coarse_lvl_lon];
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::PerformDecompression(DecompressionVariable<T, DIM>::kTag4D)
{
    cmc_debug_msg("Decompression of the variable starts...");

    [[maybe_unused]] constexpr int kDim = 4;
    constexpr int kTimeID = 0;
    constexpr int kLevID = 1;
    constexpr int kLatID = 2;
    constexpr int kLonID = 3;

    /* Get a view on the root level */
    cmc::bits::vector_view root_level_encoding(this->encoded_data_.data());

    /* Decode the root level value */
    this->data_ = DecodeRootLevelValue<T>(root_level_encoding);

    cmc_assert(data_.size() == static_cast<size_t>(1));

    cmc_assert(this->bytes_per_level_[this->decompression_step_] % sizeof(uint64_t) == 0);

    size_t level_offset{this->bytes_per_level_[this->decompression_step_] / sizeof(uint64_t)};
    ++(this->decompression_step_);

    /* Perform the iterative decompression */
    for (int decompression_step_idx{0}; decompression_step_idx < this->num_decompression_steps_ - 1; ++decompression_step_idx)
    {
        /* Set a view to this level's encoding */
        cmc::bits::vector_view level_encoding(this->encoded_data_.data() + level_offset);

        /* Set the view accordingly in the stream decoder */
        this->stream_decoder_.StartDecoding(level_encoding);

        /* Get the dimensions of this and the next level */
        const std::vector<int>& curr_lvl_dims = this->dim_lengths_[decompression_step_idx];
        const std::vector<int>& next_lvl_dims = this->dim_lengths_[decompression_step_idx + 1];

        /* Allocate a refined data vector */
        std::vector<T> data_new;
        data_new.reserve(next_lvl_dims[kTimeID] * next_lvl_dims[kLevID] * next_lvl_dims[kLatID] * next_lvl_dims[kLonID]);

        /* Iterate over the finer level */
        for (int time = 0; time < next_lvl_dims[kTimeID]; ++time)
        {
            for (int lev = 0; lev < next_lvl_dims[kLevID]; ++lev)
            {
                for (int lat = 0; lat < next_lvl_dims[kLatID]; ++lat)
                {
                    for (int lon = 0; lon < next_lvl_dims[kLonID]; ++lon)
                    {
                        /* Get the corresponding coarse level value */
                        const T coarse_val = GetCoarseValue<T>(this->data_, curr_lvl_dims, time, lev, lat, lon);

                        /* Refine the value */
                        const T fine_val = RefineValue<T>(this->stream_decoder_, coarse_val);

                        /* Store the refined value */
                        data_new.push_back(fine_val);
                    }
                }
            }
        }

        /* Switch to the new refined data for the next iteration */
        data_ = std::move(data_new);

        /* Update the offset in the encoded data */
        level_offset += this->bytes_per_level_[this->decompression_step_] / sizeof(uint64_t);

        /* Update the decompression step */
        ++(this->decompression_step_);
    }

    has_been_decompressed_ = true;
    cmc_debug_msg("Decompression of variable is completed.");
}


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::Decompress()
{
    if constexpr (DIM == 1)
    {
        cmc_debug_msg("Lossless 1D Decompression");
        this->PerformDecompression(DecompressionVariable<T, DIM>::tag1D);
    } else if constexpr (DIM == 2)
    {
        cmc_debug_msg("Lossless 2D Decompression");
        this->PerformDecompression(DecompressionVariable<T, DIM>::tag2D);
    } else if constexpr (DIM == 3)
    {
        cmc_debug_msg("Lossless 3D Decompression");
        this->PerformDecompression(DecompressionVariable<T, DIM>::tag3D);
    } else if constexpr (DIM == 4)
    {
        cmc_debug_msg("Lossless 4D Decompression");
        this->PerformDecompression(DecompressionVariable<T, DIM>::tag4D);
    } else
    {
        cmc_err_msg("Unsupported variable's dimensionality (Dim = ", DIM, ").");
    }
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
std::vector<T>
DecompressionVariable<T, DIM>::GetDecompressedData() const
{
    if (not has_been_decompressed_) {cmc_err_msg("The data has not yet been decompressed");}

    return this->data_;
}


}


#endif /* !CMC_PATCH_LOSSLESS_CMC_MULTI_RES_DECOMPRESSION_HXX */
