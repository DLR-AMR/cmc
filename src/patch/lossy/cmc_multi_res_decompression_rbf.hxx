#ifndef CMC_PATCH_LOSSY_MULTI_RES_DECOMPRESSION_RBF_HXX
#define CMC_PATCH_LOSSY_MULTI_RES_DECOMPRESSION_RBF_HXX

#include "cmc.hxx"
#include "patch/lossy/cmc_multi_res_util.hxx"
#include "amr/lossy/cmc_par_multi_res_error_mesh.hxx"
#include "amr/lossy/cmc_par_multi_res_idw_interpolation_util.hxx"
#include "utilities/cmc_compression_schema.hxx"

#include <string>
#include <span>
#include <filesystem>
#include <array>
#include <vector>
#include <algorithm>
#include <execution>

namespace cmc::patch::lossy::multi_res::rbf
{

/* Forward declaration of the general compression variable */
template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
class DecompressionVariable;

struct CompressionInfoStruct
{
    const uint64_t*
    GetCompressionHuffmanCodes() const
    {
        return encoded_huffman_codes.data();
    }

    uint64_t global_byte_count{0};
    uint64_t offset_start_encoding{0};
    uint64_t data_type{0};
    uint64_t dimensionality{0};
    uint64_t compression_scheme{0};
    uint64_t permitted_abs_error{0};
    uint64_t num_decompression_levels{0};
    uint64_t pack_size{0};
    std::vector<uint64_t> global_level_bytes;
    std::vector<std::vector<uint64_t>> dim_lengths_pyramid;
    uint64_t num_bytes_compression_huffman_codes{0};

    std::vector<uint64_t> header_data;
    std::vector<uint64_t> encoded_huffman_codes;
};

inline CompressionInfoStruct
ConstructCompressionInfoStruct(const std::vector<uint64_t>& encoded_preamble)
{
    /* Define a start pointer to the data */
    const uint64_t* start_ptr = encoded_preamble.data();
    size_t offset{0};

    CompressionInfoStruct info;

    info.global_byte_count = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(start_ptr + offset));
    ++offset;

    info.offset_start_encoding = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(start_ptr + offset));
    ++offset;

    info.data_type = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(start_ptr + offset));
    ++offset;

    info.dimensionality = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(start_ptr + offset));
    ++offset;

    info.compression_scheme = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(start_ptr + offset));
    ++offset;

    info.permitted_abs_error = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(start_ptr + offset));
    ++offset;

    info.num_decompression_levels = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(start_ptr + offset));
    ++offset;

    const int num_levels = info.num_decompression_levels;
    info.global_level_bytes.reserve(num_levels);
    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        info.global_level_bytes.push_back(cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(start_ptr + offset)));
        ++offset;
    }

    info.dim_lengths_pyramid.reserve(num_levels);
    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        info.dim_lengths_pyramid.emplace_back();
        info.dim_lengths_pyramid.back().reserve(info.dimensionality);
        for (uint64_t dim_idx{0}; dim_idx < info.dimensionality; ++dim_idx)
        {
            info.dim_lengths_pyramid.back().push_back(cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(start_ptr + offset)));
            ++offset;
        }
    }

    info.num_bytes_compression_huffman_codes = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(start_ptr + offset));
    ++offset;

    return info;
}

inline void
DecodeCompressionInfoStruct(CompressionInfoStruct& encoded_info)
{
    size_t offset{0};

    encoded_info.global_byte_count = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_info.header_data.data() + offset));
    ++offset;

    encoded_info.offset_start_encoding = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_info.header_data.data() + offset));
    ++offset;

    encoded_info.data_type = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_info.header_data.data() + offset));
    ++offset;

    encoded_info.dimensionality = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_info.header_data.data() + offset));
    ++offset;

    encoded_info.compression_scheme = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_info.header_data.data() + offset));
    ++offset;

    encoded_info.permitted_abs_error = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_info.header_data.data() + offset));
    ++offset;

    encoded_info.num_decompression_levels = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_info.header_data.data() + offset));
    ++offset;

    const int num_levels = encoded_info.num_decompression_levels;
    encoded_info.global_level_bytes.reserve(num_levels);
    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        encoded_info.global_level_bytes.push_back(cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_info.header_data.data() + offset)));
        ++offset;
    }

    encoded_info.dim_lengths_pyramid.reserve(num_levels);
    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        encoded_info.dim_lengths_pyramid.emplace_back();
        encoded_info.dim_lengths_pyramid.back().reserve(encoded_info.dimensionality);
        for (uint64_t dim_idx{0}; dim_idx < encoded_info.dimensionality; ++dim_idx)
        {
            encoded_info.dim_lengths_pyramid.back().push_back(cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_info.header_data.data() + offset)));
            ++offset;
        }
    }

    encoded_info.num_bytes_compression_huffman_codes = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_info.header_data.data() + offset));
    ++offset;
}

inline
void
PrintCompressionInfo(const std::vector<uint64_t>& encoded_preamble, const std::string& file_name)
{
    /* Get the compression info from the encoding */
    const CompressionInfoStruct info = ConstructCompressionInfoStruct(encoded_preamble);

    /* Print the gathered information */
    cmc_global_msg("Compression Information retrieved from: ", file_name);
    cmc_global_msg("\t Overall Byte Count: ", info.global_byte_count, " bytes");
    cmc_global_msg("\t Compression Preamble Byte Count: ", info.offset_start_encoding, " bytes");
    cmc_global_msg("\t DataType: ", info.data_type);
    cmc_global_msg("\t Dimensionality: ", info.dimensionality, "D");
    cmc_global_msg("\t Compression Scheme: ", info.compression_scheme);
    cmc_global_msg("\t Number of Compression LVLs: ", info.num_decompression_levels);

    cmc_global_msg("\t Compression Level Byte Count:");
    for (uint64_t lvl_idx{0}; lvl_idx < info.num_decompression_levels; ++lvl_idx)
    {
        cmc_global_msg("\t\t Level ", lvl_idx, ": ", info.global_level_bytes[lvl_idx], " bytes");
    }
    
    cmc_global_msg("\t Dimensionality Lengths per Level:");
    for (uint64_t lvl_idx{0}; lvl_idx < info.num_decompression_levels; ++lvl_idx)
    {
        std::cout << "\t\t Level " << lvl_idx << ": ";
        for (uint64_t dim_idx{0}; dim_idx < info.dimensionality; ++dim_idx)
        {
            std::cout << info.dim_lengths_pyramid[lvl_idx][dim_idx] << " ";
        }
        std::cout << std::endl;
    }

    cmc_global_msg("\t Number of Bytes Serialized Compression Huffman Codes: ", info.num_bytes_compression_huffman_codes, " bytes");
}

template <typename Int>
requires std::is_integral_v<Int>
inline void 
CheckFileError(const Int ret_val)
{
    if (ret_val != 0) [[unlikely]]
    {
        cmc_err_msg("An error during a file io operation occured.");
    }
}

inline void 
CheckFileReadCorrectness(const size_t ret_val, const size_t expected_num_elems)
{
    if (ret_val != expected_num_elems) [[unlikely]]
    {
        cmc_err_msg("An unexpected number of elements have been read from the file.");
    }
}

inline
void
ReadCompressionInfo(const std::string& file_name)
{
    if (const std::filesystem::path input_file_path(file_name); not std::filesystem::exists(input_file_path))
    {
        cmc_err_msg("The compressed file does not exist!");
    }

    /* Open the file */
    std::FILE* fhandle = std::fopen(file_name.c_str(), "rb");
    if (fhandle == nullptr)
    {
        cmc_err_msg("The file ", file_name, " could not be opened!");
    }

    /* Read the first two uint64_ts from the file */
    std::array<uint64_t, 2> num_bytes{};

    /* Move to the start of the file */
    const int rv_seek_start = std::fseek(fhandle, 0, SEEK_SET);
    CheckFileError<int>(rv_seek_start);

    /* Read the first two uint64_ts */
    const size_t rv_read_start_bytes = std::fread(num_bytes.data(), sizeof(uint64_t), 2, fhandle);
    CheckFileReadCorrectness(rv_read_start_bytes, 2);
    
    /* Read the second value from the file which gives the number of preamble bytes */
    const uint64_t num_preamble_bytes = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(num_bytes[1]);

    cmc_assert(num_preamble_bytes % sizeof(uint64_t) == 0); //The stream length is a multiple of 64 bit
    const size_t stream_length = num_preamble_bytes / sizeof(uint64_t);

    /* Allocate memory for the compression preamble */
    std::vector<uint64_t> preamble(stream_length);

    /* Move to the start of the file */
    const int rv_seek_start2 = std::fseek(fhandle, 0, SEEK_SET);
    CheckFileError<int>(rv_seek_start2);

    /* Read the complete preamble from the root rank of the shared communicator */
    const size_t rv_read_preamble = std::fread(preamble.data(), sizeof(uint64_t), stream_length, fhandle);
    CheckFileReadCorrectness(rv_read_preamble, stream_length);

    PrintCompressionInfo(preamble, file_name);

    /* Close the file */
    std::fclose(fhandle);

    cmc_debug_msg("The file ", file_name, " has been closed.");
}

constexpr int kMaxDecompressionLevelUndefined = -1;

/* Actual class definition of the decompression variable */
template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
class DecompressionVariable
{
public:
    DecompressionVariable() = delete;
    DecompressionVariable(const std::string& file_name);

    void Decompress();

    std::vector<T> GetDecompressedData() const;
    void MoveDecompressedDataInto(std::vector<T>& output_data);

    static struct kTag4D{} tag4D;
    static struct kTag3D{} tag3D;
    static struct kTag2D{} tag2D;
    static struct kTag1D{} tag1D;
    
private:
    void InquireCompressionInfo();
    const uint64_t* OpenSharedLevelDataWindow(const int level);
    void CloseSharedLevelDataWindow();
    void DecodeRootLevelValue();
    void PerformPredictionDecompression(const std::array<T, k1DNumControlValues>& control_values, const float permitted_abs_error, const int32_t lon, const std::array<int32_t, DIM>& next_lvl_dim_lengths);
    void PerformPredictionDecompression(const std::array<T, k2DNumControlValues>& control_values, const float permitted_abs_error, const int32_t lon, const int32_t lat, const std::array<int32_t, DIM>& next_lvl_dim_lengths);
    void PerformPredictionDecompression(const std::array<T, k3DNumControlValues>& control_values, const float permitted_abs_error, const int32_t lon, const int32_t lat, const int32_t lev, const std::array<int32_t, DIM>& next_lvl_dim_lengths);
    void PerformPredictionDecompression(const std::array<T, k4DNumControlValues>& control_values, const float permitted_abs_error, const int32_t lon, const int32_t lat, const int32_t lev, const int32_t time, const std::array<int32_t, DIM>& next_lvl_dim_lengths);
    void Decompress(kTag1D);
    void Decompress(kTag2D);
    void Decompress(kTag3D);
    void Decompress(kTag4D);

    /* File containing the compressed data */
    const std::string file_name_;

    /* The handle to the compressed file */
    std::FILE* fhandle_;

    /* Infos extarcted from the preamble */
    CompressionInfoStruct compression_info_;

    std::vector<std::array<int32_t, DIM>> dim_length_pyramid_;

    /* The encoding of the current comrpessed level */
    std::vector<uint64_t> level_encoding_;

    /* Previous coarse data */
    std::vector<T> coarse_data_;

    /* The current data during the decompression */
    std::vector<T> data_;

    /* The permitted absolute error criterion */
    float abs_error_;

    /* Number of decompression steps to perform */
    int32_t num_decompression_levels_;

    /* Indicator of the maximum desired decompression level */
    int32_t max_decompression_level{kMaxDecompressionLevelUndefined};

    /* A decoder */
    cmc::bits::StreamDecoder<SymbolType> stream_decoder_;

    /* A step counter for the compression */
    int32_t decompression_step_idx_{0};

    /* A flag whether the decompression has been carried out and the file been closed */
    bool is_already_decompressed_{false};
};

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
DecompressionVariable<T, DIM>::DecompressionVariable(const std::string& file_name)
: file_name_(file_name)
{
    /* Check if the compressed file exists */
    if (const std::filesystem::path input_file_path(this->file_name_); not std::filesystem::exists(input_file_path))
    {
        throw std::invalid_argument("The compressed file does not exist!");
    }

    /* Open the file */
    this->fhandle_ = std::fopen(file_name.c_str(), "rb");
    if (this->fhandle_ == nullptr)
    {
        cmc_err_msg("The file ", file_name, " could not be opened!");
    }
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::InquireCompressionInfo()
{
    /* Read the first two uint64_ts from the file */
    std::array<uint64_t, 2> num_bytes{};

    /* Move to the start of the file */
    const int rv_seek_start = std::fseek(this->fhandle_, 0, SEEK_SET);
    CheckFileError<int>(rv_seek_start);

    /* Read the first two uint64_ts */
    const size_t rv_read_start_bytes = std::fread(num_bytes.data(), sizeof(uint64_t), 2, this->fhandle_);
    CheckFileReadCorrectness(rv_read_start_bytes, 2);

    /* Read the second value from the file which gives the number of preamble bytes */
    const uint64_t num_preamble_bytes = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(num_bytes[1]);

    cmc_assert(num_preamble_bytes % sizeof(uint64_t) == 0); //The stream length is a multiple of 64 bit
    const size_t stream_length = num_preamble_bytes / sizeof(uint64_t);

    /* Allocate memory for the compression preamble */
    std::vector<uint64_t> preamble(stream_length);

    /* Move to the start of the file */
    const int rv_seek_start2 = std::fseek(this->fhandle_, 0, SEEK_SET);
    CheckFileError<int>(rv_seek_start2);

    /* Read the complete preamble from the root rank of the shared communicator */
    const size_t rv_read_preamble = std::fread(preamble.data(), sizeof(uint64_t), stream_length, this->fhandle_);
    CheckFileReadCorrectness(rv_read_preamble, stream_length);

    /* Store the preamble in the info struct */
    this->compression_info_.header_data = std::move(preamble);

    /* Decode the compression preamble */
    DecodeCompressionInfoStruct(this->compression_info_);

    /* Read the Huffman codes */
    const size_t num_vals_huffman_codes = this->compression_info_.num_bytes_compression_huffman_codes / sizeof(uint64_t);
    std::vector<uint64_t> encoded_huffman_codes(num_vals_huffman_codes);

    const size_t rv_read_huff_codes = std::fread(encoded_huffman_codes.data(), sizeof(uint64_t), num_vals_huffman_codes, this->fhandle_);
    CheckFileReadCorrectness(rv_read_huff_codes, num_vals_huffman_codes);

    /* Store the encoded Huffman codes */
    this->compression_info_.encoded_huffman_codes = std::move(encoded_huffman_codes);

    /* Construct the dimension length pyramid */
    this->dim_length_pyramid_.reserve(this->compression_info_.num_decompression_levels);
    for (uint64_t lvl_idx{0}; lvl_idx < this->compression_info_.num_decompression_levels; ++lvl_idx)
    {
        this->dim_length_pyramid_.emplace_back();
        for (int dim_idx{0}; dim_idx < DIM; ++dim_idx)
        {
            this->dim_length_pyramid_[lvl_idx][dim_idx] = this->compression_info_.dim_lengths_pyramid[lvl_idx][dim_idx];
        }
    }

    /* Store the number of decompression steps */
    this->num_decompression_levels_ = static_cast<int>(this->compression_info_.num_decompression_levels); 

    /* Store the permitted absolute error criterion */
    static_assert(sizeof(float) == sizeof(uint32_t));
    this->abs_error_ = std::bit_cast<float>(static_cast<uint32_t>(this->compression_info_.permitted_abs_error));
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
const uint64_t*
DecompressionVariable<T, DIM>::OpenSharedLevelDataWindow(const int level)
{
    cmc_assert(this->compression_info_.global_level_bytes.size() > static_cast<size_t>(level) && level >= 0);

    /* Get the number of bytes for this encoding level */
    const uint64_t level_bytes = this->compression_info_.global_level_bytes[level];
    cmc_assert(level_bytes % sizeof(uint64_t) == 0);
    const uint64_t level_vals = level_bytes / sizeof(uint64_t);

    /* Move to the correct position in the file for the current level on the current process */
    const uint64_t var_level_offset = std::accumulate(this->compression_info_.global_level_bytes.begin(), std::next(this->compression_info_.global_level_bytes.begin(), level), 0);
    const uint64_t file_offset = this->compression_info_.offset_start_encoding + this->compression_info_.num_bytes_compression_huffman_codes + var_level_offset;
    
    const int rv_seek_start = std::fseek(this->fhandle_, file_offset, SEEK_SET);
    CheckFileError<int>(rv_seek_start);

    this->level_encoding_ = std::vector<uint64_t>(level_vals);
    const size_t rv_read_level = std::fread(this->level_encoding_.data(), sizeof(uint64_t), level_vals, this->fhandle_);
    CheckFileReadCorrectness(rv_read_level, level_vals);

    return this->level_encoding_.data();
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::CloseSharedLevelDataWindow()
{
    this->level_encoding_.clear();
}

template<OneByteArithmeticType T>
inline std::vector<T>
GetRootLevelValueFromView(cmc::bits::vector_view lvl_data_start_view)
{
    std::vector<T> data;
    const OneByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<OneByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    data.push_back(std::bit_cast<T>(uvalue));

    return data;
}

template<TwoByteArithmeticType T>
inline std::vector<T>
GetRootLevelValueFromView(cmc::bits::vector_view lvl_data_start_view)
{
    std::vector<T> data;
    const TwoByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<TwoByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    data.push_back(std::bit_cast<T>(uvalue));

    return data;
}

template<FourByteArithmeticType T>
inline std::vector<T>
GetRootLevelValueFromView(cmc::bits::vector_view lvl_data_start_view)
{
    std::vector<T> data;
    const FourByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<FourByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    data.push_back(std::bit_cast<T>(uvalue));

    return data;
}

template<EightByteArithmeticType T>
inline std::vector<T>
GetRootLevelValueFromView(cmc::bits::vector_view lvl_data_start_view)
{
    std::vector<T> data;
    const EightByteResidualType uvalue = lvl_data_start_view.GetNextBitSequence<EightByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    data.push_back(std::bit_cast<T>(uvalue));

    return data;
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::DecodeRootLevelValue()
{
    /* Create a shared window on the root level */
    const int root_level{0};
    const uint64_t* shared_lvl_start_ptr = this->OpenSharedLevelDataWindow(root_level);

    /* Create a view on the data and extarct it, since there is no special encoding applied to root elvel values */
    cmc::bits::vector_view data_view(shared_lvl_start_ptr);

    /* De-Serialize the values */
    this->data_ = GetRootLevelValueFromView<T>(data_view);

    /* Close the shared level window after all values have been extarcted */
    this->CloseSharedLevelDataWindow();

    /* Update the decompression count */
    ++(this->decompression_step_idx_);
}

template<OneByteArithmeticType T>
inline T 
GetUnpredictableValue(cmc::bits::StreamDecoder<SymbolType>& stream_decoder)
{
    const auto value = stream_decoder.GetNextBitSequence<OneByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    return std::bit_cast<T>(value);
}

template<TwoByteArithmeticType T>
inline T 
GetUnpredictableValue(cmc::bits::StreamDecoder<SymbolType>& stream_decoder)
{
    const auto value = stream_decoder.GetNextBitSequence<TwoByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    return std::bit_cast<T>(value);
}

template<FourByteArithmeticType T>
inline T 
GetUnpredictableValue(cmc::bits::StreamDecoder<SymbolType>& stream_decoder)
{
    const auto value = stream_decoder.GetNextBitSequence<FourByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    return std::bit_cast<T>(value);
}

template<EightByteArithmeticType T>
inline T 
GetUnpredictableValue(cmc::bits::StreamDecoder<SymbolType>& stream_decoder)
{
    const auto value = stream_decoder.GetNextBitSequence<EightByteResidualType>(sizeof(T) * cmc::bits::kCharBit);
    return std::bit_cast<T>(value);
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::PerformPredictionDecompression(const std::array<T, k4DNumControlValues>& control_values, const float permitted_abs_error, const int32_t lon, const int32_t lat, const int32_t lev, const int32_t time, const std::array<int32_t, DIM>& next_lvl_dim_lengths)
{
    static_assert(DIM == 4);
    constexpr int32_t kTimeID = 0;
    constexpr int32_t kLevID = 1;
    constexpr int32_t kLatID = 2;
    constexpr int32_t kLonID = 3;

    /* Perform the prediction */
    std::array<T, kPackSize<DIM>> predicted_data = Perform4DRBFPrediction<T, DIM>(control_values);

    int32_t pred_access_idx{-1};
    /* Gather the values for this initial patch */
    for (int32_t time_idx = 0; time_idx < kDIMReductionFactor; ++time_idx)
    {
        for (int32_t lev_idx = 0; lev_idx < kDIMReductionFactor; ++lev_idx)
        {
            for (int32_t lat_idx = 0; lat_idx < kDIMReductionFactor; ++lat_idx)
            {
                for (int32_t lon_idx = 0; lon_idx < kDIMReductionFactor; ++lon_idx)
                {
                    ++pred_access_idx;
                    if (kDIMReductionFactor * time + time_idx >= next_lvl_dim_lengths[kTimeID] || kDIMReductionFactor * lev + lev_idx >= next_lvl_dim_lengths[kLevID] || kDIMReductionFactor * lat + lat_idx >= next_lvl_dim_lengths[kLatID] || kDIMReductionFactor * lon + lon_idx >= next_lvl_dim_lengths[kLonID]) [[unlikely]]
                    {
                        continue;
                    } else
                    {
                        /* Get the next entropy symbol */
                        const SymbolType next_symbol = this->stream_decoder_.DecodeNextEntropySymbol();

                        /* Check, if the value was un-predictable */
                        if (next_symbol == kFlagUnpredictable) [[unlikely]]
                        {
                            /* Get the unpredicted value from the stream */
                            const T unpredicted_value = GetUnpredictableValue<T>(this->stream_decoder_);

                            /* We set the correctly predicted data in the next level data */
                            SetValue<T>(this->data_, unpredicted_value, kDIMReductionFactor * time + time_idx, kDIMReductionFactor * lev + lev_idx,
                                        kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID], next_lvl_dim_lengths[kTimeID]);
                        } else
                        {
                            /* Get the quantization bin from the entropy code */
                            const auto [was_prediction_greater, quant_bin] = GetQuantizationBinFromEntropySymbol(next_symbol);

                            /* Apply the de-quantization */
                            const T decompressed_value = DequantizeValue<T>(predicted_data[pred_access_idx], permitted_abs_error, was_prediction_greater, quant_bin);

                            /* We set the correctly predicted data in the next level data */
                            SetValue<T>(this->data_, decompressed_value, kDIMReductionFactor * time + time_idx, kDIMReductionFactor * lev + lev_idx,
                                        kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID], next_lvl_dim_lengths[kTimeID]);
                        }
                    }
                }
            }
        }
    }
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::Decompress(DecompressionVariable<T, DIM>::kTag4D)
{
    constexpr int32_t kTimeID = 0;
    constexpr int32_t kLevID = 1;
    constexpr int32_t kLatID = 2;
    constexpr int32_t kLonID = 3;

    cmc_debug_msg("Lossy patch-based decompression starts...");

    cmc_assert(this->dim_length_pyramid_.size() == static_cast<size_t>(this->num_decompression_levels_));

    cmc_debug_msg("Number of decompression iterations to be performed: ", this->num_decompression_levels_);
    
    /* Set up the Huffman decoder within the stream decoder */
    this->stream_decoder_.StartHuffmanCodesDecoding(this->compression_info_.GetCompressionHuffmanCodes());

    /* Decode the root level */
    this->DecodeRootLevelValue();

    /* Define the iterators to the dimension lengths */
    auto coarse_dim_lengths_iter = this->dim_length_pyramid_.begin();

    /* Perform the iterative compression steps up until the root level */
    for (int32_t lvl_idx{1}; lvl_idx < this->num_decompression_levels_; ++lvl_idx, ++coarse_dim_lengths_iter)
    {
        /* Set the current data as the caorse level data */
        this->coarse_data_ = std::move(this->data_);

        /* Define a reference for the ease of notation */
        const std::vector<T>& coarse_data = this->coarse_data_;
        const std::array<int32_t, DIM>& coarse_dim_lengths = *coarse_dim_lengths_iter;
        
        /* Get the next level dimension lengths */
        const std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::next(coarse_dim_lengths_iter);

        /* Compute the nubmer of values on the next level */
        const int32_t num_next_level_elems = std::reduce(next_lvl_dim_lengths.begin(), next_lvl_dim_lengths.end(), int32_t{1}, std::multiplies<int32_t>());

        /* Allocate the next level data */
        this->data_ = std::vector<T>(num_next_level_elems);

        /* Open the corresponding encoding level */
        const uint64_t* shared_lvl_start_ptr = this->OpenSharedLevelDataWindow(lvl_idx);

        /* Create a view on the data */
        cmc::bits::vector_view level_view(shared_lvl_start_ptr);

        /* Start the decoding process from the view */
        this->stream_decoder_.StartDecoding(level_view);

        cmc_debug_msg("A refinement iteration is initialized (step: ", lvl_idx, ").");
        cmc_debug_msg("Level's data dimensions are: ", coarse_dim_lengths[kTimeID], ", ", coarse_dim_lengths[kLevID], ", ", coarse_dim_lengths[kLatID], ", ", coarse_dim_lengths[kLonID]);

        /* Iterate over patches */
        for (int32_t time = 0; time < coarse_dim_lengths[kTimeID]; ++time)
        {
            for (int32_t lev = 0; lev < coarse_dim_lengths[kLevID]; ++lev)
            {
                for (int32_t lat = 0; lat < coarse_dim_lengths[kLatID]; ++lat)
                {
                    for (int32_t lon = 0; lon < coarse_dim_lengths[kLonID]; ++lon)
                    {
                        /** Now this element will be refined via prediction **/

                        /* Gather the face control values */
                        const std::array<T, k4DNumControlValues> control_values = GetFaceControlValues<T, DIM>(coarse_data, time, lev, lat, lon, coarse_dim_lengths);

                        /* Perform the refinement prediction */
                        this->PerformPredictionDecompression(control_values, this->abs_error_, lon, lat, lev, time, next_lvl_dim_lengths);
                    }
                }
            }
        }

        /* Close the data level */
        this->CloseSharedLevelDataWindow();

        cmc_debug_msg("The refinement prediction iteration is finished.");
    }

    cmc_debug_msg("The lossy compression of is finished.");
}


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::PerformPredictionDecompression(const std::array<T, k3DNumControlValues>& control_values, const float permitted_abs_error, const int32_t lon, const int32_t lat, const int32_t lev, const std::array<int32_t, DIM>& next_lvl_dim_lengths)
{
    static_assert(DIM == 3);
    constexpr int32_t kLevID = 0;
    constexpr int32_t kLatID = 1;
    constexpr int32_t kLonID = 2;

    /* Perform the prediction */
    std::array<T, kPackSize<DIM>> predicted_data = Perform3DRBFPrediction<T, DIM>(control_values);

    int32_t pred_access_idx{-1};
    /* Gather the values for this initial patch */
    for (int32_t lev_idx = 0; lev_idx < kDIMReductionFactor; ++lev_idx)
    {
        for (int32_t lat_idx = 0; lat_idx < kDIMReductionFactor; ++lat_idx)
        {
            for (int32_t lon_idx = 0; lon_idx < kDIMReductionFactor; ++lon_idx)
            {
                ++pred_access_idx;
                if (kDIMReductionFactor * lev + lev_idx >= next_lvl_dim_lengths[kLevID] || kDIMReductionFactor * lat + lat_idx >= next_lvl_dim_lengths[kLatID] || kDIMReductionFactor * lon + lon_idx >= next_lvl_dim_lengths[kLonID]) [[unlikely]]
                {
                    continue;
                } else
                {
                    /* Get the next entropy symbol */
                    const SymbolType next_symbol = this->stream_decoder_.DecodeNextEntropySymbol();

                    /* Check, if the value was un-predictable */
                    if (next_symbol == kFlagUnpredictable) [[unlikely]]
                    {
                        /* Get the unpredicted value from the stream */
                        const T unpredicted_value = GetUnpredictableValue<T>(this->stream_decoder_);

                        /* We set the correctly predicted data in the next level data */
                        SetValue<T>(this->data_, unpredicted_value, kDIMReductionFactor * lev + lev_idx,
                                    kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID]);
                    } else
                    {
                        /* Get the quantization bin from the entropy code */
                        const auto [was_prediction_greater, quant_bin] = GetQuantizationBinFromEntropySymbol(next_symbol);

                        /* Apply the de-quantization */
                        const T decompressed_value = DequantizeValue<T>(predicted_data[pred_access_idx], permitted_abs_error, was_prediction_greater, quant_bin);

                        /* We set the correctly predicted data in the next level data */
                        SetValue<T>(this->data_, decompressed_value, kDIMReductionFactor * lev + lev_idx,
                                    kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID]);
                    }
                }
            }
        }
    }
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::Decompress(DecompressionVariable<T, DIM>::kTag3D)
{
    constexpr int32_t kLevID = 0;
    constexpr int32_t kLatID = 1;
    constexpr int32_t kLonID = 2;

    cmc_debug_msg("Lossy patch-based decompression starts...");

    cmc_assert(this->dim_length_pyramid_.size() == static_cast<size_t>(this->num_decompression_levels_));

    cmc_debug_msg("Number of decompression iterations to be performed: ", this->num_decompression_levels_);
    
    /* Set up the Huffman decoder within the stream decoder */
    this->stream_decoder_.StartHuffmanCodesDecoding(this->compression_info_.GetCompressionHuffmanCodes());

    /* Decode the root level */
    this->DecodeRootLevelValue();

    /* Define the iterators to the dimension lengths */
    auto coarse_dim_lengths_iter = this->dim_length_pyramid_.begin();

    /* Perform the iterative compression steps up until the root level */
    for (int32_t lvl_idx{1}; lvl_idx < this->num_decompression_levels_; ++lvl_idx, ++coarse_dim_lengths_iter)
    {
        /* Set the current data as the caorse level data */
        this->coarse_data_ = std::move(this->data_);

        /* Define a reference for the ease of notation */
        const std::vector<T>& coarse_data = this->coarse_data_;
        const std::array<int32_t, DIM>& coarse_dim_lengths = *coarse_dim_lengths_iter;
        
        /* Get the next level dimension lengths */
        const std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::next(coarse_dim_lengths_iter);

        /* Compute the nubmer of values on the next level */
        const int32_t num_next_level_elems = std::reduce(next_lvl_dim_lengths.begin(), next_lvl_dim_lengths.end(), int32_t{1}, std::multiplies<int32_t>());

        /* Allocate the next level data */
        this->data_ = std::vector<T>(num_next_level_elems);

        /* Open the corresponding encoding level */
        const uint64_t* shared_lvl_start_ptr = this->OpenSharedLevelDataWindow(lvl_idx);

        /* Create a view on the data */
        cmc::bits::vector_view level_view(shared_lvl_start_ptr);

        /* Start the decoding process from the view */
        this->stream_decoder_.StartDecoding(level_view);

        cmc_debug_msg("A refinement iteration is initialized (step: ", lvl_idx, ").");
        cmc_debug_msg("Level's data dimensions are: ", coarse_dim_lengths[kLevID], ", ", coarse_dim_lengths[kLatID], ", ", coarse_dim_lengths[kLonID]);

        /* Iterate over patches */
        for (int32_t lev = 0; lev < coarse_dim_lengths[kLevID]; ++lev)
        {
            for (int32_t lat = 0; lat < coarse_dim_lengths[kLatID]; ++lat)
            {
                for (int32_t lon = 0; lon < coarse_dim_lengths[kLonID]; ++lon)
                {
                    /** Now this element will be refined via prediction **/

                    /* Gather the face control values */
                    const std::array<T, k3DNumControlValues> control_values = GetFaceControlValues<T, DIM>(coarse_data, lev, lat, lon, coarse_dim_lengths);

                    /* Perform the refinement prediction */
                    this->PerformPredictionDecompression(control_values, this->abs_error_, lon, lat, lev, next_lvl_dim_lengths);
                }
            }
        }

        /* Close the data level */
        this->CloseSharedLevelDataWindow();

        cmc_debug_msg("The refinement prediction iteration is finished.");
    }

    cmc_debug_msg("The lossy compression of is finished.");
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::PerformPredictionDecompression(const std::array<T, k2DNumControlValues>& control_values, const float permitted_abs_error, const int32_t lon, const int32_t lat, const std::array<int32_t, DIM>& next_lvl_dim_lengths)
{
    static_assert(DIM == 2);
    constexpr int32_t kLatID = 0;
    constexpr int32_t kLonID = 1;

    /* Perform the prediction */
    std::array<T, kPackSize<DIM>> predicted_data = Perform2DRBFPrediction<T, DIM>(control_values);

    int32_t pred_access_idx{-1};
    /* Gather the values for this initial patch */
    for (int32_t lat_idx = 0; lat_idx < kDIMReductionFactor; ++lat_idx)
    {
        for (int32_t lon_idx = 0; lon_idx < kDIMReductionFactor; ++lon_idx)
        {
            ++pred_access_idx;
            if (kDIMReductionFactor * lat + lat_idx >= next_lvl_dim_lengths[kLatID] || kDIMReductionFactor * lon + lon_idx >= next_lvl_dim_lengths[kLonID]) [[unlikely]]
            {
                continue;
            } else
            {
                /* Get the next entropy symbol */
                const SymbolType next_symbol = this->stream_decoder_.DecodeNextEntropySymbol();

                /* Check, if the value was un-predictable */
                if (next_symbol == kFlagUnpredictable) [[unlikely]]
                {
                    /* Get the unpredicted value from the stream */
                    const T unpredicted_value = GetUnpredictableValue<T>(this->stream_decoder_);

                    /* We set the correctly predicted data in the next level data */
                    SetValue<T>(this->data_, unpredicted_value,
                                kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID]);
                } else
                {
                    /* Get the quantization bin from the entropy code */
                    const auto [was_prediction_greater, quant_bin] = GetQuantizationBinFromEntropySymbol(next_symbol);

                    /* Apply the de-quantization */
                    const T decompressed_value = DequantizeValue<T>(predicted_data[pred_access_idx], permitted_abs_error, was_prediction_greater, quant_bin);

                    /* We set the correctly predicted data in the next level data */
                    SetValue<T>(this->data_, decompressed_value,
                                kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID]);
                }
            }
        }
    }
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::Decompress(DecompressionVariable<T, DIM>::kTag2D)
{
    constexpr int32_t kLatID = 0;
    constexpr int32_t kLonID = 1;

    cmc_debug_msg("Lossy patch-based decompression starts...");

    cmc_assert(this->dim_length_pyramid_.size() == static_cast<size_t>(this->num_decompression_levels_));

    cmc_debug_msg("Number of decompression iterations to be performed: ", this->num_decompression_levels_);
    
    /* Set up the Huffman decoder within the stream decoder */
    this->stream_decoder_.StartHuffmanCodesDecoding(this->compression_info_.GetCompressionHuffmanCodes());

    /* Decode the root level */
    this->DecodeRootLevelValue();

    /* Define the iterators to the dimension lengths */
    auto coarse_dim_lengths_iter = this->dim_length_pyramid_.begin();

    /* Perform the iterative compression steps up until the root level */
    for (int32_t lvl_idx{1}; lvl_idx < this->num_decompression_levels_; ++lvl_idx, ++coarse_dim_lengths_iter)
    {
        /* Set the current data as the caorse level data */
        this->coarse_data_ = std::move(this->data_);

        /* Define a reference for the ease of notation */
        const std::vector<T>& coarse_data = this->coarse_data_;
        const std::array<int32_t, DIM>& coarse_dim_lengths = *coarse_dim_lengths_iter;
        
        /* Get the next level dimension lengths */
        const std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::next(coarse_dim_lengths_iter);

        /* Compute the nubmer of values on the next level */
        const int32_t num_next_level_elems = std::reduce(next_lvl_dim_lengths.begin(), next_lvl_dim_lengths.end(), int32_t{1}, std::multiplies<int32_t>());

        /* Allocate the next level data */
        this->data_ = std::vector<T>(num_next_level_elems);

        /* Open the corresponding encoding level */
        const uint64_t* shared_lvl_start_ptr = this->OpenSharedLevelDataWindow(lvl_idx);

        /* Create a view on the data */
        cmc::bits::vector_view level_view(shared_lvl_start_ptr);

        /* Start the decoding process from the view */
        this->stream_decoder_.StartDecoding(level_view);

        cmc_debug_msg("A refinement iteration is initialized (step: ", lvl_idx, ").");
        cmc_debug_msg("Level's data dimensions are: ", coarse_dim_lengths[kLatID], ", ", coarse_dim_lengths[kLonID]);

        /* Iterate over patches */
        for (int32_t lat = 0; lat < coarse_dim_lengths[kLatID]; ++lat)
        {
            for (int32_t lon = 0; lon < coarse_dim_lengths[kLonID]; ++lon)
            {
                /** Now this element will be refined via prediction **/

                /* Gather the face control values */
                const std::array<T, k2DNumControlValues> control_values = GetFaceControlValues<T, DIM>(coarse_data, lat, lon, coarse_dim_lengths);

                /* Perform the refinement prediction */
                this->PerformPredictionDecompression(control_values, this->abs_error_, lon, lat, next_lvl_dim_lengths);
            }
        }

        /* Close the data level */
        this->CloseSharedLevelDataWindow();

        cmc_debug_msg("The refinement prediction iteration is finished.");
    }

    cmc_debug_msg("The lossy decompression of is finished.");
}


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::PerformPredictionDecompression(const std::array<T, k1DNumControlValues>& control_values, const float permitted_abs_error, const int32_t lon, const std::array<int32_t, DIM>& next_lvl_dim_lengths)
{
    static_assert(DIM == 1);
    constexpr int32_t kLonID = 0;

    /* Perform the prediction */
    std::array<T, kPackSize<DIM>> predicted_data = Perform1DRBFPrediction<T, DIM>(control_values);

    int32_t pred_access_idx{-1};
    /* Gather the values for this initial patch */
    for (int32_t lon_idx = 0; lon_idx < kDIMReductionFactor; ++lon_idx)
    {
        ++pred_access_idx;
        if (kDIMReductionFactor * lon + lon_idx >= next_lvl_dim_lengths[kLonID]) [[unlikely]]
        {
            continue;
        } else
        {
            /* Get the next entropy symbol */
            const SymbolType next_symbol = this->stream_decoder_.DecodeNextEntropySymbol();

            /* Check, if the value was un-predictable */
            if (next_symbol == kFlagUnpredictable) [[unlikely]]
            {
                /* Get the unpredicted value from the stream */
                const T unpredicted_value = GetUnpredictableValue<T>(this->stream_decoder_);

                /* We set the correctly predicted data in the next level data */
                SetValue<T>(this->data_, unpredicted_value,
                            kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID]);
            } else
            {
                /* Get the quantization bin from the entropy code */
                const auto [was_prediction_greater, quant_bin] = GetQuantizationBinFromEntropySymbol(next_symbol);

                /* Apply the de-quantization */
                const T decompressed_value = DequantizeValue<T>(predicted_data[pred_access_idx], permitted_abs_error, was_prediction_greater, quant_bin);

                /* We set the correctly predicted data in the next level data */
                SetValue<T>(this->data_, decompressed_value,
                            kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID]);
            }
        }
    }
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::Decompress(DecompressionVariable<T, DIM>::kTag1D)
{
    constexpr int32_t kLonID = 0;

    cmc_debug_msg("Lossy patch-based decompression starts...");

    cmc_assert(this->dim_length_pyramid_.size() == static_cast<size_t>(this->num_decompression_levels_));

    cmc_debug_msg("Number of decompression iterations to be performed: ", this->num_decompression_levels_);
    
    /* Set up the Huffman decoder within the stream decoder */
    this->stream_decoder_.StartHuffmanCodesDecoding(this->compression_info_.GetCompressionHuffmanCodes());

    /* Decode the root level */
    this->DecodeRootLevelValue();

    /* Define the iterators to the dimension lengths */
    auto coarse_dim_lengths_iter = this->dim_length_pyramid_.begin();

    /* Perform the iterative compression steps up until the root level */
    for (int32_t lvl_idx{1}; lvl_idx < this->num_decompression_levels_; ++lvl_idx, ++coarse_dim_lengths_iter)
    {
        /* Set the current data as the caorse level data */
        this->coarse_data_ = std::move(this->data_);

        /* Define a reference for the ease of notation */
        const std::vector<T>& coarse_data = this->coarse_data_;
        const std::array<int32_t, DIM>& coarse_dim_lengths = *coarse_dim_lengths_iter;
        
        /* Get the next level dimension lengths */
        const std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::next(coarse_dim_lengths_iter);

        /* Compute the nubmer of values on the next level */
        const int32_t num_next_level_elems = std::reduce(next_lvl_dim_lengths.begin(), next_lvl_dim_lengths.end(), int32_t{1}, std::multiplies<int32_t>());

        /* Allocate the next level data */
        this->data_ = std::vector<T>(num_next_level_elems);

        /* Open the corresponding encoding level */
        const uint64_t* shared_lvl_start_ptr = this->OpenSharedLevelDataWindow(lvl_idx);

        /* Create a view on the data */
        cmc::bits::vector_view level_view(shared_lvl_start_ptr);

        /* Start the decoding process from the view */
        this->stream_decoder_.StartDecoding(level_view);

        cmc_debug_msg("A refinement iteration is initialized (step: ", lvl_idx, ").");
        cmc_debug_msg("Level's data dimensions are: ", coarse_dim_lengths[kLonID]);

        /* Iterate over patches */
        for (int32_t lon = 0; lon < coarse_dim_lengths[kLonID]; ++lon)
        {
            /** Now this element will be refined via prediction **/

            /* Gather the face control values */
            const std::array<T, k1DNumControlValues> control_values = GetFaceControlValues<T, DIM>(coarse_data, lon, coarse_dim_lengths);

            /* Perform the refinement prediction */
            this->PerformPredictionDecompression(control_values, this->abs_error_, lon, next_lvl_dim_lengths);
        }

        /* Close the data level */
        this->CloseSharedLevelDataWindow();

        cmc_debug_msg("The refinement prediction iteration is finished.");
    }

    cmc_debug_msg("The lossy decompression of is finished.");
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::Decompress()
{
    if (this->is_already_decompressed_) [[unlikely]]
    {
        cmc_err_msg("The variable has already been decompressed!");
    }

    /* Inquire the basic information about the compression */
    this->InquireCompressionInfo();

    if constexpr (DIM == 1)
    {
        cmc_global_msg("Lossy 1D Decompression");
        this->Decompress(DecompressionVariable<T, DIM>::tag1D);
    } else if constexpr (DIM == 2)
    {
        cmc_global_msg("Lossy 2D Decompression");
        this->Decompress(DecompressionVariable<T, DIM>::tag2D);
    } else if constexpr (DIM == 3)
    {
        cmc_global_msg("Lossy 3D Decompression");
        this->Decompress(DecompressionVariable<T, DIM>::tag3D);
    } else if constexpr (DIM == 4)
    {
        cmc_global_msg("Lossy 4D Decompression");
        this->Decompress(DecompressionVariable<T, DIM>::tag4D);
    } else
    {
        cmc_err_msg("Unsupported variable's dimensionality (DIM = ", DIM, ").");
    }

    /* Close the compressed file */
    std::fclose(this->fhandle_);

    /* Set the flag that the data has already been decompressed */
    this->is_already_decompressed_ = true;

    cmc_debug_msg("Decompression has been completed.");
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
std::vector<T>
DecompressionVariable<T, DIM>::GetDecompressedData() const 
{
    return this->data_;
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
DecompressionVariable<T, DIM>::MoveDecompressedDataInto(std::vector<T>& output_data)  
{
    output_data = std::move(this->data_);
    this->data_ = std::vector<T>();
}

}
#endif /* !CMC_PATCH_LOSSY_MULTI_RES_DECOMPRESSION_RBF_HXX */
