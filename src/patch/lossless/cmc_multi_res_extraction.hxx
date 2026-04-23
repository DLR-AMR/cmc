#ifndef CMC_PATCH_LOSSLESS_CMC_MULTI_RES_EXTRACTION_HXX
#define CMC_PATCH_LOSSLESS_CMC_MULTI_RES_EXTRACTION_HXX

#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "patch/lossless/cmc_patch_multi_res_extraction_util.hxx"
#include "utilities/cmc_compression_schema.hxx"

#include <stdexcept>
#include <filesystem>
#include <cstdio>

#include <bitset>

namespace cmc::serial::patch::lossless::multi_res
{

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
class CompressionVariable
{
public:
    CompressionVariable() = delete;
    explicit CompressionVariable(std::vector<T>&& init_data, const std::array<int, DIM> dimension_lengths)
    : data_(std::move(init_data)), init_dimension_lengths_(dimension_lengths)
    {
        const bool are_dimension_lengths_valid = std::invoke([&dimension_lengths]() -> bool {
            for (const int& dim_length : dimension_lengths)
            {
                if (dim_length <= 0){return false;}
            }
            return true;
        });
        if (not are_dimension_lengths_valid) {throw std::invalid_argument("Dimension lengths have to be non-negative!");}

        /* Get the largest dimension length */
        const auto max_dim_length_iter = std::max_element(dimension_lengths.begin(), dimension_lengths.end());
        if (max_dim_length_iter == dimension_lengths.end()){throw std::invalid_argument("dimension_lengths");}

        const int max_dim_length = *max_dim_length_iter;

        /* Determine the number of compression iterations */
        this->num_compression_lvls_ = ComputeNumCompressionLevels<DIM>(max_dim_length);
        if (this->num_compression_lvls_ <= 0) {throw std::invalid_argument("dimension_lengths");}

        /* Store the initial dimension lengths in the dimension length pyramid */
        this->dim_lengths_.reserve(this->num_compression_lvls_ + 1);
        this->dim_lengths_.push_back(dimension_lengths);
    }

    void Compress();
    void WriteData(const std::string file_name);

private:
    PatchEncodingData<T, DIM> PerformExtraction(const std::array<T, kNumMaxChildrenElements<DIM>>& patch_values, const int num_elements);

    std::vector<T> data_;
    const std::array<int, DIM> init_dimension_lengths_;

    std::vector<std::vector<SymbolType>> entropy_codes_;
    std::vector<std::vector<T>> residuals_;
    std::vector<std::vector<uint64_t>> encoded_data_;
    std::vector<uint64_t> serialized_entropy_dictionary_;

    std::vector<std::array<int, DIM>> dim_lengths_;
    int num_compression_lvls_;

    static struct kTag4D{} tag4D;
    static struct kTag3D{} tag3D;
    static struct kTag2D{} tag2D;
    static struct kTag1D{} tag1D;

    void PerformCompression(kTag1D);
    void PerformCompression(kTag2D);
    void PerformCompression(kTag3D);
    void PerformCompression(kTag4D);
};

template <typename T>
inline T
GetValue(const std::vector<T>& data, const int time, const int lev, const int lat, const int lon, const int kLonLength, const int kLatLength, const int kLevLength,  [[maybe_unused]] const int kTimeLength)
{
    cmc_assert(time * (kLevLength * kLatLength * kLonLength) + lev * (kLatLength * kLonLength) + lat * kLonLength + lon < static_cast<int>(data.size()));
    return data[time * (kLevLength * kLatLength * kLonLength) + lev * (kLatLength * kLonLength) + lat * kLonLength + lon];
}

template <typename T>
inline void
SetValue(std::vector<T>& data, const T value, const int time, const int lev, const int lat, const int lon, const int kLonLength, const int kLatLength,  const int kLevLength, [[maybe_unused]] const int kTimeLength)
{
    cmc_assert(time * (kLevLength * kLatLength * kLonLength) + lev * (kLatLength * kLonLength) + lat * kLonLength + lon < static_cast<int>(data.size()));
    data[time * (kLevLength * kLatLength * kLonLength) + lev * (kLatLength * kLonLength) + lat * kLonLength + lon] = value;
}

template <typename T>
inline T
GetValue(const std::vector<T>& data, const int lev, const int lat, const int lon, const int kLonLength, const int kLatLength, [[maybe_unused]] const int kLevLength)
{
    cmc_assert(lev * (kLatLength * kLonLength) + lat * kLonLength + lon < static_cast<int>(data.size()));
    return data[lev * (kLatLength * kLonLength) + lat * kLonLength + lon];
}

template <typename T>
inline void
SetValue(std::vector<T>& data, const T value, const int lev, const int lat, const int lon, const int kLonLength, const int kLatLength,  [[maybe_unused]] const int kLevLength)
{
    cmc_assert(lev * (kLatLength * kLonLength) + lat * kLonLength + lon < static_cast<int>(data.size()));
    data[lev * (kLatLength * kLonLength) + lat * kLonLength + lon] = value;
}

template <typename T>
inline T
GetValue(const std::vector<T>& data, const int lat, const int lon, const int kLonLength, [[maybe_unused]] const int kLatLength)
{
    cmc_assert(lat * kLonLength + lon < static_cast<int>(data.size()));
    return data[lat * kLonLength + lon];
}

template <typename T>
inline void
SetValue(std::vector<T>& data, const T value, const int lat, const int lon, const int kLonLength, [[maybe_unused]] const int kLatLength)
{
    cmc_assert(lat * kLonLength + lon < static_cast<int>(data.size()));
    data[lat * kLonLength + lon] = value;
}


template<typename T>
std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>>
CollectEntropySymbols(const std::vector<std::vector<SymbolType>>& level_entropy_symbols)
{
    /* Get the number of all possible entropy symbols */
    constexpr int num_entropy_symbols = GetNumEntropySymbols<T>();

    /* Set the array and zero intialiaze the frequencies */
    std::array<uint64_t, num_entropy_symbols> entropy_symbol_frequencies{};

    /* Iterate through all entropy codes and accumulate their frequencies */
    for (size_t lvl_idx{0}; lvl_idx < level_entropy_symbols.size(); ++lvl_idx)
    {
        for (size_t elem_idx{0}; elem_idx < level_entropy_symbols[lvl_idx].size(); ++elem_idx)
        {
            /* Convert the symbol to the corresponding array index */
            const int array_idx = MapEntropySymbolToArrayIndex<T>(level_entropy_symbols[lvl_idx][elem_idx]);

            /* Update the frequency */
            ++entropy_symbol_frequencies[array_idx];
        }
    }

    std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>> global_symbol_frequencies;
    global_symbol_frequencies.reserve(num_entropy_symbols);

    /* Iterate through all codes and construct the symbol frequency table */
    for (int idx{0}; idx < num_entropy_symbols; ++idx)
    {
        /* Convert the index back to the entropy symbol */
        const SymbolType entropy_symbol = MapArrayIndexToEntropySymbol<T>(idx);

        /* Store the symbol with the global frequency */
        global_symbol_frequencies.emplace_back(entropy_symbol, entropy_symbol_frequencies[idx]);
    }

    return global_symbol_frequencies;
}

template <OneByteArithmeticType T, int32_t DIM>
inline std::vector<uint64_t>
EncodeRootLevel(const std::vector<T>& root_level_data)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(root_level_data.size() * sizeof(T));
    
    for (size_t idx{0}; idx < root_level_data.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<OneByteResidualType>(std::bit_cast<OneByteResidualType>(root_level_data[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}

template <TwoByteArithmeticType T, int32_t DIM>
inline std::vector<uint64_t>
EncodeRootLevel(const std::vector<T>& root_level_data)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(root_level_data.size() * sizeof(T));
    
    for (size_t idx{0}; idx < root_level_data.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<TwoByteResidualType>(std::bit_cast<TwoByteResidualType>(root_level_data[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}

template <FourByteArithmeticType T, int32_t DIM>
inline std::vector<uint64_t>
EncodeRootLevel(const std::vector<T>& root_level_data)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(root_level_data.size() * sizeof(T));
    
    for (size_t idx{0}; idx < root_level_data.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<FourByteResidualType>(std::bit_cast<FourByteResidualType>(root_level_data[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}

template <EightByteArithmeticType T, int32_t DIM>
inline std::vector<uint64_t>
EncodeRootLevel(const std::vector<T>& root_level_data)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(root_level_data.size() * sizeof(T));
    
    for (size_t idx{0}; idx < root_level_data.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<EightByteResidualType>(std::bit_cast<EightByteResidualType>(root_level_data[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}

template <typename T, int32_t DIM>
requires Dimension<DIM>
std::vector<uint64_t>
GenerateVariableHeader(const std::vector<std::vector<uint64_t>>& encoded_data_, const std::vector<std::array<int, DIM>>& dim_lengths_, const std::vector<uint64_t>& serialized_entropy_dictionary_)
{
    const uint64_t var_header_size = 7 + encoded_data_.size() + DIM * dim_lengths_.size() + serialized_entropy_dictionary_.size();

    std::vector<uint64_t> var_header;
    var_header.reserve(var_header_size);

    cmc_assert(encoded_data_.size() >= 1);

    const int num_compression_levels = encoded_data_.size();

    uint64_t global_byte_count = var_header_size;
    for (const auto& lvl_encoding : encoded_data_)
    {
        global_byte_count += lvl_encoding.size();
    }
    global_byte_count *= sizeof(uint64_t);


    /* Global Bytes compressed variable */
    var_header.push_back(cmc::bits::ConvertToBigEndian<SizeType>(global_byte_count));

    /* Store the offset from the stream start to the start of the root level encoding */
    var_header.push_back(cmc::bits::ConvertToBigEndian<SizeType>(sizeof(uint64_t) * var_header_size));

    /* Store the data type */
    var_header.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(ConvertToCmcType<T>())));

    /* Store the dimensionality */
    var_header.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(DIM)));

    /* Store the compression scheme */
    var_header.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(CompressionSchema::PatchMultiResExtraction)));

    /* Store the number of mesh compresion levels */
    var_header.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(num_compression_levels)));

    /* Append the global level bytes */
    for (auto lvl_iter = encoded_data_.begin(); lvl_iter != encoded_data_.end(); ++lvl_iter)
    {
        var_header.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(lvl_iter->size() * sizeof(uint64_t))));
    }

    /* Store the dimension lengths per level (in reverse) */
    for (auto dim_lengths_iter = dim_lengths_.rbegin(); dim_lengths_iter != dim_lengths_.rend(); ++dim_lengths_iter)
    {
        for (int dim_idx{0}; dim_idx < DIM; ++dim_idx)
        {
            var_header.push_back(cmc::bits::ConvertToBigEndian<SizeType>(dim_lengths_iter->operator[](dim_idx)));
        }
    }

    /* Store the number of bytes for the entropy dictionary */
    var_header.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(serialized_entropy_dictionary_.size() * sizeof(uint64_t))));
    std::copy_n(serialized_entropy_dictionary_.begin(), serialized_entropy_dictionary_.size(), std::back_inserter(var_header));

    return var_header;
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::WriteData(const std::string file_name)
{
    /* Create an additional variable header */
    const std::vector<uint64_t> var_header = GenerateVariableHeader<T, DIM>(this->encoded_data_, this->dim_lengths_, this->serialized_entropy_dictionary_);

    /* Open the file, but delete it first, if it already exists */
    const std::filesystem::path output_file_path(file_name);
    if (std::filesystem::exists(output_file_path))
    {
        std::remove(output_file_path.c_str());
    }
    std::FILE* file_out = std::fopen(file_name.c_str(), "wb");

    /* Write the variable header */
    const std::size_t num_written_vals_header = std::fwrite(var_header.data(), sizeof(uint64_t), var_header.size(), file_out);
    if (num_written_vals_header != var_header.size()) {cmc::cmc_err_msg("An unexpecetd number of compressed elements has been written!");}

    /* Write the encoed level data */
    for (const std::vector<uint64_t>& encoded_level : this->encoded_data_)
    {
        const std::size_t num_written_lvl_vals = std::fwrite(encoded_level.data(), sizeof(uint64_t), encoded_level.size(), file_out);
        if (num_written_lvl_vals != encoded_level.size()) {cmc::cmc_err_msg("An unexpecetd number of compressed elements has been written!");}
    }

    /* Close the file */
    std::fclose(file_out);
}

template <OneByteArithmeticType T>
inline void
AppendResidualBits(cmc::bits::vector& level_encoding, const T residual, const int start_pos, const int end_pos)
{
    level_encoding.AppendBits(std::bit_cast<OneByteResidualType>(residual), start_pos, end_pos);
}

template <TwoByteArithmeticType T>
inline void
AppendResidualBits(cmc::bits::vector& level_encoding, const T residual, const int start_pos, const int end_pos)
{
    level_encoding.AppendBits(std::bit_cast<TwoByteResidualType>(residual), start_pos, end_pos);
}

template <FourByteArithmeticType T>
inline void
AppendResidualBits(cmc::bits::vector& level_encoding, const T residual, const int start_pos, const int end_pos)
{
    level_encoding.AppendBits(std::bit_cast<FourByteResidualType>(residual), start_pos, end_pos);
}

template <EightByteArithmeticType T>
inline void
AppendResidualBits(cmc::bits::vector& level_encoding, const T residual, const int start_pos, const int end_pos)
{
    level_encoding.AppendBits(std::bit_cast<EightByteResidualType>(residual), start_pos, end_pos);
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::Compress()
{
    if constexpr (DIM == 1)
    {
        cmc_debug_msg("Lossless 1D Compression");
        this->PerformCompression(CompressionVariable<T, DIM>::tag1D);
    } else if constexpr (DIM == 2)
    {
        cmc_debug_msg("Lossless 2D Compression");
        this->PerformCompression(CompressionVariable<T, DIM>::tag2D);
    } else if constexpr (DIM == 3)
    {
        cmc_debug_msg("Lossless 3D Compression");
        this->PerformCompression(CompressionVariable<T, DIM>::tag3D);
    } else if constexpr (DIM == 4)
    {
        cmc_debug_msg("Lossless 4D Compression");
        this->PerformCompression(CompressionVariable<T, DIM>::tag4D);
    } else
    {
        cmc_err_msg("Unsupported variable's dimensionality (Dim = ", DIM, ").");
    }
}

template<ArithmeticType T>
inline std::vector<T>
CreatePredictors(const std::span<const T> values)
{
    /* Utilize each individual value as a predictor */
    std::vector<T> predictors(values.size() + 2);
    std::copy_n(values.begin(), values.size(), predictors.begin());

    /* Append the arithmetic mean as a predictor */
    predictors[values.size()] = ComputeArithmeticMean<T>(std::span<T>(predictors.data(), values.size()));

    /* Append the mid-range as a predictor */
    predictors[values.size() + 1] = ComputeMidRange<T>(std::span<T>(predictors.data(), values.size()));

    return predictors;
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
PatchEncodingData<T, DIM>
CompressionVariable<T, DIM>::PerformExtraction(const std::array<T, kNumMaxChildrenElements<DIM>>& patch_values, const int num_elements)
{
    cmc_assert(num_elements > 0);
    if (num_elements == 1) [[unlikely]]
    {
        PatchEncodingData<T, DIM> elem_coding;
        elem_coding.coarse_value = patch_values.front();
        return elem_coding;
    }

    cmc_assert(num_elements > 1 && num_elements <= kNumMaxChildrenElements<DIM>);
    [[assume(num_elements > 1 && num_elements <= kNumMaxChildrenElements<DIM>)]];

    /* Perform the multi-resolution extraction */
    const std::vector<T> predictors = CreatePredictors<T>(std::span<const T>(patch_values.data(), num_elements));

    int current_lzc{-1};
    T lzc_maximizing_predictor{};

    PatchEncodingData<T, DIM> elem_coding;

    /* For all predictors, we evaluate the one that gives us the overall maximum number of leading zeros */
    for (int pred_idx{0}; pred_idx < num_elements + 2; ++pred_idx)
    {
        std::array<SymbolType, kNumMaxChildrenElements<DIM>> current_entropy_codes{};
        Residuals<T, DIM> current_residuals{};
        
        /* Check the predictor for all element values */
        for (int elem_idx{0}; elem_idx < num_elements; ++elem_idx)
        {
            /* Compute the residual */
            const auto [is_approx_greater, residual] = cmc::bits::ComputeIntegerResidual<T>(predictors[pred_idx], predictors[elem_idx]);
            
            /* Store the entropy code */
            current_entropy_codes[elem_idx] = CreateEntropySymbol(is_approx_greater, residual);
            
            /* Store the residual */
            current_residuals.residuals[elem_idx] = residual;
        }

        /* Compute the cumulative leading zero count */
        const int cumulative_lzc = std::transform_reduce(std::execution::par_unseq, current_residuals.residuals.cbegin(), current_residuals.residuals.cend(), static_cast<int>(0),
                                                         std::plus<>{}, [](auto res){return cmc::bits::GetLZC(res);});

        /* If the predictor maximizes the LZC, we store it */
        if (current_lzc < cumulative_lzc)
        {
            current_lzc = cumulative_lzc;
            lzc_maximizing_predictor = predictors[pred_idx];
            std::copy_n(current_residuals.residuals.cbegin(), kNumMaxChildrenElements<DIM>, elem_coding.residuals.begin());
            std::copy_n(current_entropy_codes.cbegin(), kNumMaxChildrenElements<DIM>, elem_coding.entropy_symbols.begin());
        }
    }

    /* We store the lzc_maximizing_predictor for the next coarse level */
    elem_coding.coarse_value = lzc_maximizing_predictor;


    return elem_coding;
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::PerformCompression(CompressionVariable<T, DIM>::kTag4D)
{
    [[maybe_unused]] constexpr int kDim = 4;
    constexpr int kTimeID = 0;
    constexpr int kLevID = 1;
    constexpr int kLatID = 2;
    constexpr int kLonID = 3;

    /* Perform the iterative compression steps up until the root level */
    for (int lvl_idx{0}; lvl_idx < num_compression_lvls_; ++lvl_idx)
    {
        /* Get the current dimension lengths */
        const int time_length = this->dim_lengths_.back()[kTimeID];
        const int lev_length = this->dim_lengths_.back()[kLevID];
        const int lat_length = this->dim_lengths_.back()[kLatID];
        const int lon_length = this->dim_lengths_.back()[kLonID];

        /* Allocate the entropy symbols on this level */
        this->entropy_codes_.emplace_back(time_length * lev_length * lat_length * lon_length, 0);

        /* Allocate the coarse level */
        std::vector<T> coarse_data;
        coarse_data.reserve(((time_length / kDimReductionFactor) + 1) * ((lev_length / kDimReductionFactor) + 1) * ((lat_length / kDimReductionFactor) + 1) * ((lon_length / kDimReductionFactor) + 1));

        std::array<T, kNumMaxChildrenElements<DIM>> patch_values{};

        /* Iterate over patches */
        for (int time = 0; time < time_length; time += kDimReductionFactor)
        {
            for (int lev = 0; lev < lev_length; lev += kDimReductionFactor)
            {
                for (int lat = 0; lat < lat_length; lat += kDimReductionFactor)
                {
                    for (int lon = 0; lon < lon_length; lon += kDimReductionFactor)
                    {
                        int elem_idx{0};

                        /* Gather the values for this patch */
                        for (int time_idx = 0; time_idx < kDimReductionFactor; ++time_idx)
                        {
                            for (int lev_idx = 0; lev_idx < kDimReductionFactor; ++lev_idx)
                            {
                                for (int lat_idx = 0; lat_idx < kDimReductionFactor; ++lat_idx)
                                {
                                    for (int lon_idx = 0; lon_idx < kDimReductionFactor; ++lon_idx)
                                    {
                                        if (time + time_idx >= time_length || lev + lev_idx >= lev_length || lat + lat_idx >= lat_length || lon + lon_idx >= lon_length) [[unlikely]]
                                        {
                                            continue;
                                        } else
                                        {
                                            patch_values[elem_idx] = GetValue<T>(this->data_, time + time_idx, lev + lev_idx, lat + lat_idx, lon + lon_idx, lon_length, lat_length, lev_length, time_length);
                                            ++elem_idx;
                                        }
                                    }
                                }
                            }
                        }

                        /* Perform an extraction operation on this patch of values */
                        const PatchEncodingData<T, DIM> extracted_values = this->PerformExtraction(patch_values, elem_idx);

                        /* Re-Assign the fine values */
                        int extracted_val_idx{0};
                        for (int time_idx = 0; time_idx < kDimReductionFactor; ++time_idx)
                        {
                            for (int lev_idx = 0; lev_idx < kDimReductionFactor; ++lev_idx)
                            {
                                for (int lat_idx = 0; lat_idx < kDimReductionFactor; ++lat_idx)
                                {
                                    for (int lon_idx = 0; lon_idx < kDimReductionFactor; ++lon_idx)
                                    {
                                        if (time + time_idx >= time_length || lev + lev_idx >= lev_length || lat + lat_idx >= lat_length || lon + lon_idx >= lon_length) [[unlikely]]
                                        {
                                            continue;
                                        } else
                                        {
                                            /* Update the fine value accordingly */
                                            SetValue<T>(this->data_, std::bit_cast<T>(extracted_values.residuals[extracted_val_idx]), time + time_idx, lev + lev_idx, lat + lat_idx, lon + lon_idx, lon_length, lat_length, lev_length, time_length);
                                            SetValue<SymbolType>(this->entropy_codes_.back(), extracted_values.entropy_symbols[extracted_val_idx], time + time_idx, lev + lev_idx, lat + lat_idx, lon + lon_idx, lon_length, lat_length, lev_length, time_length);
                                            ++extracted_val_idx;
                                        }
                                    }
                                }
                            }
                        }
                        /* Store the extracted information */
                        coarse_data.push_back(extracted_values.coarse_value);
                    }
                }
            }
        }

        /* Store the data */
        this->residuals_.push_back(std::move(this->data_));
        this->data_ = std::move(coarse_data);

        /* Update the dimension sizes */
        const int new_time_length = (time_length % kDimReductionFactor == 0 ? (time_length / kDimReductionFactor) : (time_length / kDimReductionFactor) + 1);
        const int new_lev_length = (lev_length % kDimReductionFactor == 0 ? (lev_length / kDimReductionFactor) : (lev_length / kDimReductionFactor) + 1);
        const int new_lat_length = (lat_length % kDimReductionFactor == 0 ? (lat_length / kDimReductionFactor) : (lat_length / kDimReductionFactor) + 1);
        const int new_lon_length = (lon_length % kDimReductionFactor == 0 ? (lon_length / kDimReductionFactor) : (lon_length / kDimReductionFactor) + 1);
    
        /* Store the dimensions for the next iteration */
        dim_lengths_.push_back(std::array<int, DIM>{new_time_length, new_lev_length, new_lat_length, new_lon_length});
    }

    /* Determine Entropy Codes */
    const std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>> entropy_symbols = CollectEntropySymbols<T>(this->entropy_codes_);

    /* Build a Huffman coder */
    cmc::entropy_coding::huffman::HuffmanCoder<SymbolType> entropy_coder(entropy_symbols);

    /* Store the serialized Huffman coder of the intra element compression */
    this->serialized_entropy_dictionary_ = entropy_coder.SerializeHuffmanCodesBEPadded();

    encoded_data_.reserve(num_compression_lvls_ + 1);

    /* Encode Root Level */
    encoded_data_.push_back(EncodeRootLevel<T, DIM>(this->data_));

    /* Encode data from higer levels in reverse */
    auto lvl_entropy_symbols_iter = this->entropy_codes_.rbegin();
    auto lvl_residual_iter = this->residuals_.rbegin();

    for (int lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx)
    {
        /* Define some references for the ease of notation */
        const std::vector<SymbolType>& entropy_symbols = *lvl_entropy_symbols_iter;
        const auto& residuals = *lvl_residual_iter;

        /* Allocate a cmc_bits_vector */
        cmc::bits::vector lvl_encoding;
        lvl_encoding.Reserve(sizeof(T) * cmc::bits::kCharBit * (residuals.size() / 2));

        cmc_assert(entropy_symbols.size() == residuals.size());
        const int num_elems = entropy_symbols.size();

        /* Encode the symbols and the entropy codes in an interleved fashion */
        for (int elem_idx{0}; elem_idx < num_elems; ++elem_idx)
        {
            /* Encode the entropy symbol */
            const cmc::entropy_coding::huffman::HuffmanCode code = entropy_coder.EncodeSymbol(entropy_symbols[elem_idx]);

            /* Retrieve the LZC from the entropy symbol */
            const int lzc = GetLZCFromEntropySymbol(entropy_symbols[elem_idx]);

            /* Serialize the encoded entropy symbol */
            lvl_encoding.AppendBits(code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - code.code_length), 0);
            
            /* We do not need to encode the implicit given one-bit following the LZC */
            if (lzc + 1 < static_cast<int>(sizeof(T) * cmc::bits::kCharBit)) [[likely]]
            {
                /* Append the significant reisdual bits */
                AppendResidualBits<T>(lvl_encoding, residuals[elem_idx], lzc + 1, 0);
            }
        }

        /* Store the level-wise encoding */
        encoded_data_.push_back(lvl_encoding.GetSerializedByteStreamBE());

        /* Update the iterators */
        ++lvl_entropy_symbols_iter;
        ++lvl_residual_iter;
    }
}


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::PerformCompression(CompressionVariable<T, DIM>::kTag3D)
{
    [[maybe_unused]] constexpr int kDim = 3;
    constexpr int kLevID = 0;
    constexpr int kLatID = 1;
    constexpr int kLonID = 2;

    /* Perform the iterative compression steps up until the root level */
    for (int lvl_idx{0}; lvl_idx < num_compression_lvls_; ++lvl_idx)
    {
        /* Get the current dimension lengths */
        const int lev_length = this->dim_lengths_.back()[kLevID];
        const int lat_length = this->dim_lengths_.back()[kLatID];
        const int lon_length = this->dim_lengths_.back()[kLonID];

        /* Allocate the entropy symbols on this level */
        this->entropy_codes_.emplace_back(lev_length * lat_length * lon_length, 0);

        /* Allocate the coarse level */
        std::vector<T> coarse_data;
        coarse_data.reserve(((lev_length / kDimReductionFactor) + 1) * ((lat_length / kDimReductionFactor) + 1) * ((lon_length / kDimReductionFactor) + 1));

        std::array<T, kNumMaxChildrenElements<DIM>> patch_values{};

        /* Iterate over patches */
        for (int lev = 0; lev < lev_length; lev += kDimReductionFactor)
        {
            for (int lat = 0; lat < lat_length; lat += kDimReductionFactor)
            {
                for (int lon = 0; lon < lon_length; lon += kDimReductionFactor)
                {
                    int elem_idx{0};

                    /* Gather the values for this patch */
                    for (int lev_idx = 0; lev_idx < kDimReductionFactor; ++lev_idx)
                    {
                        for (int lat_idx = 0; lat_idx < kDimReductionFactor; ++lat_idx)
                        {
                            for (int lon_idx = 0; lon_idx < kDimReductionFactor; ++lon_idx)
                            {
                                if (lev + lev_idx >= lev_length || lat + lat_idx >= lat_length || lon + lon_idx >= lon_length) [[unlikely]]
                                {
                                    continue;
                                } else
                                {
                                    patch_values[elem_idx] = GetValue<T>(this->data_, lev + lev_idx, lat + lat_idx, lon + lon_idx, lon_length, lat_length, lev_length);
                                    ++elem_idx;
                                }
                            }
                        }
                    }

                    /* Perform an extraction operation on this patch of values */
                    const PatchEncodingData<T, DIM> extracted_values = this->PerformExtraction(patch_values, elem_idx);

                    /* Re-Assign the fine values */
                    int extracted_val_idx{0};
                    for (int lev_idx = 0; lev_idx < kDimReductionFactor; ++lev_idx)
                    {
                        for (int lat_idx = 0; lat_idx < kDimReductionFactor; ++lat_idx)
                        {
                            for (int lon_idx = 0; lon_idx < kDimReductionFactor; ++lon_idx)
                            {
                                if (lev + lev_idx >= lev_length || lat + lat_idx >= lat_length || lon + lon_idx >= lon_length) [[unlikely]]
                                {
                                    continue;
                                } else
                                {
                                    /* Update the fine value accordingly */
                                    SetValue<T>(this->data_, std::bit_cast<T>(extracted_values.residuals[extracted_val_idx]), lev + lev_idx, lat + lat_idx, lon + lon_idx, lon_length, lat_length, lev_length);
                                    SetValue<SymbolType>(this->entropy_codes_.back(), extracted_values.entropy_symbols[extracted_val_idx], lev + lev_idx, lat + lat_idx, lon + lon_idx, lon_length, lat_length, lev_length);
                                    ++extracted_val_idx;
                                }
                            }
                        }
                    }

                    /* Store the extracted information */
                    coarse_data.push_back(extracted_values.coarse_value);
                }
            }
        }

        /* Store the data */
        this->residuals_.push_back(std::move(this->data_));
        this->data_ = std::move(coarse_data);

        /* Update the dimension sizes */
        const int new_lev_length = (lev_length % kDimReductionFactor == 0 ? (lev_length / kDimReductionFactor) : (lev_length / kDimReductionFactor) + 1);
        const int new_lat_length = (lat_length % kDimReductionFactor == 0 ? (lat_length / kDimReductionFactor) : (lat_length / kDimReductionFactor) + 1);
        const int new_lon_length = (lon_length % kDimReductionFactor == 0 ? (lon_length / kDimReductionFactor) : (lon_length / kDimReductionFactor) + 1);
    
        /* Store the dimensions for the next iteration */
        dim_lengths_.push_back(std::array<int, DIM>{new_lev_length, new_lat_length, new_lon_length});
    }

    /* Determine Entropy Codes */
    const std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>> entropy_symbols = CollectEntropySymbols<T>(this->entropy_codes_);

    /* Build a Huffman coder */
    cmc::entropy_coding::huffman::HuffmanCoder<SymbolType> entropy_coder(entropy_symbols);

    /* Store the serialized Huffman coder of the intra element compression */
    this->serialized_entropy_dictionary_ = entropy_coder.SerializeHuffmanCodesBEPadded();

    encoded_data_.reserve(num_compression_lvls_ + 1);

    /* Encode Root Level */
    encoded_data_.push_back(EncodeRootLevel<T, DIM>(this->data_));

    /* Encode data from higer levels in reverse */
    auto lvl_entropy_symbols_iter = this->entropy_codes_.rbegin();
    auto lvl_residual_iter = this->residuals_.rbegin();

    for (int lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx)
    {
        /* Define some references for the ease of notation */
        const std::vector<SymbolType>& entropy_symbols = *lvl_entropy_symbols_iter;
        const auto& residuals = *lvl_residual_iter;

        /* Allocate a cmc_bits_vector */
        cmc::bits::vector lvl_encoding;
        lvl_encoding.Reserve(sizeof(T) * cmc::bits::kCharBit * (residuals.size() / 2));

        cmc_assert(entropy_symbols.size() == residuals.size());
        const int num_elems = entropy_symbols.size();

        /* Encode the symbols and the entropy codes in an interleved fashion */
        for (int elem_idx{0}; elem_idx < num_elems; ++elem_idx)
        {
            /* Encode the entropy symbol */
            const cmc::entropy_coding::huffman::HuffmanCode code = entropy_coder.EncodeSymbol(entropy_symbols[elem_idx]);

            /* Retrieve the LZC from the entropy symbol */
            const int lzc = GetLZCFromEntropySymbol(entropy_symbols[elem_idx]);

            /* Serialize the encoded entropy symbol */
            lvl_encoding.AppendBits(code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - code.code_length), 0);
            
            /* We do not need to encode the implicit given one-bit following the LZC */
            if (lzc + 1 < static_cast<int>(sizeof(T) * cmc::bits::kCharBit)) [[likely]]
            {
                /* Append the significant reisdual bits */
                AppendResidualBits<T>(lvl_encoding, residuals[elem_idx], lzc + 1, 0);
            }
        }

        /* Store the level-wise encoding */
        encoded_data_.push_back(lvl_encoding.GetSerializedByteStreamBE());

        /* Update the iterators */
        ++lvl_entropy_symbols_iter;
        ++lvl_residual_iter;
    }
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::PerformCompression(CompressionVariable<T, DIM>::kTag2D)
{
    [[maybe_unused]] constexpr int kDim = 2;
    constexpr int kLatID = 0;
    constexpr int kLonID = 1;

    /* Perform the iterative compression steps up until the root level */
    for (int lvl_idx{0}; lvl_idx < num_compression_lvls_; ++lvl_idx)
    {
        /* Get the current dimension lengths */
        const int lat_length = this->dim_lengths_.back()[kLatID];
        const int lon_length = this->dim_lengths_.back()[kLonID];

        /* Allocate the entropy symbols on this level */
        this->entropy_codes_.emplace_back(lat_length * lon_length, 0);

        /* Allocate the coarse level */
        std::vector<T> coarse_data;
        coarse_data.reserve(((lat_length / kDimReductionFactor) + 1) * ((lon_length / kDimReductionFactor) + 1));

        std::array<T, kNumMaxChildrenElements<DIM>> patch_values{};

        /* Iterate over patches */
        for (int lat = 0; lat < lat_length; lat += kDimReductionFactor)
        {
            for (int lon = 0; lon < lon_length; lon += kDimReductionFactor)
            {
                int elem_idx{0};

                /* Gather the values for this patch */
                for (int lat_idx = 0; lat_idx < kDimReductionFactor; ++lat_idx)
                {
                    for (int lon_idx = 0; lon_idx < kDimReductionFactor; ++lon_idx)
                    {
                        if (lat + lat_idx >= lat_length || lon + lon_idx >= lon_length) [[unlikely]]
                        {
                            continue;
                        } else
                        {
                            patch_values[elem_idx] = GetValue<T>(data_, lat + lat_idx, lon + lon_idx, lon_length, lat_length);
                            ++elem_idx;
                        }
                    }
                }

                /* Perform an extraction operation on this patch of values */
                const PatchEncodingData<T, DIM> extracted_values = this->PerformExtraction(patch_values, elem_idx);

                /* Re-Assign the fine values */
                int extracted_val_idx{0};
                for (int lat_idx = 0; lat_idx < kDimReductionFactor; ++lat_idx)
                {
                    for (int lon_idx = 0; lon_idx < kDimReductionFactor; ++lon_idx)
                    {
                        if (lat + lat_idx >= lat_length || lon + lon_idx >= lon_length) [[unlikely]]
                        {
                            continue;
                        } else
                        {
                            /* Update the fine value accordingly */
                            SetValue<T>(this->data_, std::bit_cast<T>(extracted_values.residuals[extracted_val_idx]), lat + lat_idx, lon + lon_idx, lon_length, lat_length);
                            SetValue<SymbolType>(this->entropy_codes_.back(), extracted_values.entropy_symbols[extracted_val_idx], lat + lat_idx, lon + lon_idx, lon_length, lat_length);
                            ++extracted_val_idx;
                        }
                    }
                }
                /* Store the extracted information */
                coarse_data.push_back(extracted_values.coarse_value);
            }
        }

        /* Store the data */
        this->residuals_.push_back(std::move(this->data_));
        this->data_ = std::move(coarse_data);

        /* Update the dimension sizes */
        const int new_lat_length = (lat_length % kDimReductionFactor == 0 ? (lat_length / kDimReductionFactor) : (lat_length / kDimReductionFactor) + 1);
        const int new_lon_length = (lon_length % kDimReductionFactor == 0 ? (lon_length / kDimReductionFactor) : (lon_length / kDimReductionFactor) + 1);
    
        /* Store the dimensions for the next iteration */
        dim_lengths_.push_back(std::array<int, DIM>{new_lat_length, new_lon_length});
    }

    /* Determine Entropy Codes */
    const std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>> entropy_symbols = CollectEntropySymbols<T>(this->entropy_codes_);

    /* Build a Huffman coder */
    cmc::entropy_coding::huffman::HuffmanCoder<SymbolType> entropy_coder(entropy_symbols);

    /* Store the serialized Huffman coder of the intra element compression */
    this->serialized_entropy_dictionary_ = entropy_coder.SerializeHuffmanCodesBEPadded();

    encoded_data_.reserve(num_compression_lvls_ + 1);

    /* Encode Root Level */
    encoded_data_.push_back(EncodeRootLevel<T, DIM>(this->data_));

    /* Encode data from higer levels in reverse */
    auto lvl_entropy_symbols_iter = this->entropy_codes_.rbegin();
    auto lvl_residual_iter = this->residuals_.rbegin();

    cmc_assert(this->entropy_codes_.size() == this->residuals_.size());

    for (int lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx)
    {
        /* Define some references for the ease of notation */
        const std::vector<SymbolType>& entropy_symbols = *lvl_entropy_symbols_iter;
        const auto& residuals = *lvl_residual_iter;

        /* Allocate a cmc_bits_vector */
        cmc::bits::vector lvl_encoding;
        lvl_encoding.Reserve(sizeof(T) * cmc::bits::kCharBit * (residuals.size() / 2));

        cmc_assert(entropy_symbols.size() == residuals.size());
        const int num_elems = entropy_symbols.size();

        /* Encode the symbols and the entropy codes in an interleved fashion */
        for (int elem_idx{0}; elem_idx < num_elems; ++elem_idx)
        {
            /* Encode the entropy symbol */
            const cmc::entropy_coding::huffman::HuffmanCode code = entropy_coder.EncodeSymbol(entropy_symbols[elem_idx]);

            /* Retrieve the LZC from the entropy symbol */
            const int lzc = GetLZCFromEntropySymbol(entropy_symbols[elem_idx]);

            /* Serialize the encoded entropy symbol */
            lvl_encoding.AppendBits(code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - code.code_length), 0);
            
            /* We do not need to encode the implicit given one-bit following the LZC */
            if (lzc + 1 < static_cast<int>(sizeof(T) * cmc::bits::kCharBit)) [[likely]]
            {
                /* Append the significant reisdual bits */
                AppendResidualBits<T>(lvl_encoding, residuals[elem_idx], lzc + 1, 0);
            }
        }

        /* Store the level-wise encoding */
        encoded_data_.push_back(lvl_encoding.GetSerializedByteStreamBE());

        /* Update the iterators */
        ++lvl_entropy_symbols_iter;
        ++lvl_residual_iter;
    }
}

}

#endif /* !CMC_PATCH_LOSSLESS_CMC_MULTI_RES_EXTRACTION_HXX */
