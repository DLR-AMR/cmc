#ifndef CMC_PATCH_LOSSLESS_FAST_MULTI_RES_EXTRACTION_WITH_CORRECTION_HXX
#define CMC_PATCH_LOSSLESS_FAST_MULTI_RES_EXTRACTION_WITH_CORRECTION_HXX

#include "cmc.hxx"
#include "patch/lossless/cmc_multi_res_util.hxx"
#include "utilities/cmc_compression_schema.hxx"

#include <string>
#include <span>
#include <filesystem>
#include <array>
#include <vector>
#include <algorithm>
#include <execution>

namespace cmc::patch::lossless::multi_res::fast
{

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
class CompressionVariable
{
public:
    CompressionVariable() = delete;
    CompressionVariable(const std::span<const T> data, const std::array<int32_t, DIM>& dimension_lengths)
    : init_data_(data), dim_lengths_(dimension_lengths),
      max_dimension_length_{GetMaximumDimensionLength<T, DIM>(dimension_lengths)}, num_compression_lvls_{GetNumCompressionIterations<T, DIM>(dimension_lengths)}
      {
        /* Compute the expected amount of data */
        const int num_data = std::reduce(dimension_lengths.begin(), dimension_lengths.end(), 1, std::multiplies<int>());

        if (static_cast<size_t>(num_data) != data.size())
        {
            cmc_err_msg("The expected amount of data (", num_data, ") given the supplied dimension-length vector does not coincide with the initial data size (", data.size(), ")!");
        }

        /* Allocate and initialize the vectors */
        this->dim_length_pyramid_.reserve(this->num_compression_lvls_ + 1);
        this->dim_length_pyramid_.push_back(dimension_lengths);
        this->data_.reserve(num_data);
        std::copy_n(data.begin(), num_data, std::back_inserter(this->data_));
        this->levelwise_patch_encodings_.reserve(this->num_compression_lvls_);
    }

    void Compress();
    void WriteCompressedData(const std::string& file_name);

    static struct kTag4D{} tag4D;
    static struct kTag3D{} tag3D;
    static struct kTag2D{} tag2D;
    static struct kTag1D{} tag1D;

private:
    void Compress(kTag1D);
    void Compress(kTag2D);
    void Compress(kTag3D);
    void Compress(kTag4D);
    std::vector<uint64_t> EncodeRootLevelData() const;
    void EncodeData();

    const std::span<const T> init_data_;
    const std::array<int32_t, DIM> dim_lengths_;

    const int32_t max_dimension_length_{0};
    const int32_t num_compression_lvls_{0};

    std::vector<T> data_;
    std::vector<std::array<int32_t, DIM>> dim_length_pyramid_;

    std::vector<std::vector<PatchEncoding<T, DIM>>> levelwise_patch_encodings_;

    std::vector<uint64_t> serialized_entropy_dictionary_;
    std::vector<std::vector<uint64_t>> levelwise_encoded_data_;
};


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::Compress()
{
    if constexpr (DIM == 1)
    {
        cmc_global_msg("Lossless 1D Compression");
        this->Compress(CompressionVariable<T, DIM>::tag1D);
    } else if constexpr (DIM == 2)
    {
        cmc_global_msg("Lossless 2D Compression");
        this->Compress(CompressionVariable<T, DIM>::tag2D);
    } else if constexpr (DIM == 3)
    {
        cmc_global_msg("Lossless 3D Compression");
        this->Compress(CompressionVariable<T, DIM>::tag3D);
    } else if constexpr (DIM == 4)
    {
        cmc_global_msg("Lossless 4D Compression");
        this->Compress(CompressionVariable<T, DIM>::tag4D);
    } else
    {
        cmc_err_msg("Unsupported variable's dimensionality (DIM = ", DIM, ").");
    }
}


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
std::pair<T, PatchEncoding<T, DIM>>
PerformExtraction(const std::array<T, kPackSize<DIM>> patch_data, const int num_values, const bool is_init_extraction)
{
    cmc_assert(num_values >= 1 && num_values <= kPackSize<DIM>);
    [[assume(num_values >= 1 && num_values <= kPackSize<DIM>)]];

    if (num_values == 1) [[unlikely]]
    {
        PatchEncoding<T, DIM> zero_coding;
        zero_coding.entropy_symbols[0] = CreateEntropySymbol(false, TransformToUInteger<T>(0));
        zero_coding.residuals[0] = 0;
        zero_coding.num_elements = 1;
        return std::make_pair(patch_data[0], std::move(zero_coding));
    }

    /* Perform the multi-resolution extraction */
    const T mean_predictor = (is_init_extraction ? CreateMatchingMeanPredictor<T, kPackSize<DIM>>(patch_data, num_values) : CreateArithmeticMeanPredictor<T, kPackSize<DIM>>(patch_data, num_values));

    PatchEncoding<T, DIM> patch_coding;
    patch_coding.num_elements = num_values;

    /* All intermediate levels are encoded without tolerances */
    const int num_vals_to_encode = (is_init_extraction ? num_values - 1 : num_values);

    /* We store the residuals from all but the last value (which is then implicitly given via the mean value) */
    for (int elem_idx{0}; elem_idx < num_vals_to_encode; ++elem_idx)
    {
        /* Compute the residual */
        const auto [is_approx_greater, residual] = cmc::bits::ComputeIntegerResidual<T>(mean_predictor, patch_data[elem_idx]);
        
        /* Store the entropy code */
        patch_coding.entropy_symbols[elem_idx] = CreateEntropySymbol(is_approx_greater, residual);
        
        /* Store the residual */
        patch_coding.residuals[elem_idx] = residual;
    }

    return std::make_pair(mean_predictor, patch_coding);
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::Compress(CompressionVariable<T, DIM>::kTag4D)
{
    static_assert(DIM == 4);
    constexpr int32_t kTimeID = 0;
    constexpr int32_t kLevID = 1;
    constexpr int32_t kLatID = 2;
    constexpr int32_t kLonID = 3;

    int32_t paccess_idx{0};
    std::array<T, kPackSize<DIM>> patch_data{};

    for (int32_t lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx)
    {
        /* Allocate a new dimension array */
        this->dim_length_pyramid_.emplace_back();

        /* Get the current dimension lengths */
        const std::array<int32_t, DIM>& dim_lengths = *std::prev(this->dim_length_pyramid_.end(), 2);

        /* Get the next level dimension array */
        std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::prev(this->dim_length_pyramid_.end(), 1);
        next_lvl_dim_lengths[kTimeID] = (dim_lengths[kTimeID] % kDIMReductionFactor == 0 ? (dim_lengths[kTimeID] / kDIMReductionFactor) : (dim_lengths[kTimeID] / kDIMReductionFactor) + 1);
        next_lvl_dim_lengths[kLevID] = (dim_lengths[kLevID] % kDIMReductionFactor == 0 ? (dim_lengths[kLevID] / kDIMReductionFactor) : (dim_lengths[kLevID] / kDIMReductionFactor) + 1);
        next_lvl_dim_lengths[kLatID] = (dim_lengths[kLatID] % kDIMReductionFactor == 0 ? (dim_lengths[kLatID] / kDIMReductionFactor) : (dim_lengths[kLatID] / kDIMReductionFactor) + 1);
        next_lvl_dim_lengths[kLonID] = (dim_lengths[kLonID] % kDIMReductionFactor == 0 ? (dim_lengths[kLonID] / kDIMReductionFactor) : (dim_lengths[kLonID] / kDIMReductionFactor) + 1);

        /* Compute the number of next level values */
        const int32_t num_next_level_vals = next_lvl_dim_lengths[kTimeID] * next_lvl_dim_lengths[kLevID] * next_lvl_dim_lengths[kLatID] * next_lvl_dim_lengths[kLonID];
        /* Allocate a next level vector for the coarse data */
        std::vector<T> coarse_values;
        coarse_values.reserve(num_next_level_vals);

        /* Allocate a level-wise patch encdoing */
        this->levelwise_patch_encodings_.emplace_back();
        this->levelwise_patch_encodings_.back().reserve(num_next_level_vals);

        /* Only the init iteration is encoded with tolerances */
        const bool is_init_iteration = (kTryToCorrectMeanFastCompression ? (lvl_idx == 0) : false);

        /* Extract the coarse level values */
        for (int32_t time = 0; time < dim_lengths[kTimeID]; time += kDIMReductionFactor)
        {
            for (int32_t lev = 0; lev < dim_lengths[kLevID]; lev += kDIMReductionFactor)
            {
                for (int32_t lat = 0; lat < dim_lengths[kLatID]; lat += kDIMReductionFactor)
                {
                    for (int32_t lon = 0; lon < dim_lengths[kLonID]; lon += kDIMReductionFactor)
                    {
                        /* Reset the access idx for each new patch */
                        paccess_idx = 0;

                        /* Gather the values for this initial patch */
                        for (int32_t time_idx = 0; time_idx < kDIMReductionFactor; ++time_idx)
                        {
                            for (int32_t lev_idx = 0; lev_idx < kDIMReductionFactor; ++lev_idx)
                            {
                                for (int32_t lat_idx = 0; lat_idx < kDIMReductionFactor; ++lat_idx)
                                {
                                    for (int32_t lon_idx = 0; lon_idx < kDIMReductionFactor; ++lon_idx)
                                    {
                                        if (time + time_idx >= dim_lengths[kTimeID] || lev + lev_idx >= dim_lengths[kLevID] || lat + lat_idx >= dim_lengths[kLatID] || lon + lon_idx >= dim_lengths[kLonID]) [[unlikely]]
                                        {
                                            continue;
                                        } else
                                        {
                                            /* Retrieve the corresponding value from this patch */
                                            patch_data[paccess_idx] = GetValue<T>(this->data_, time + time_idx, lev + lev_idx, lat + lat_idx, lon + lon_idx,
                                                                                  dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID], dim_lengths[kTimeID]);
                                        }
                                        ++paccess_idx;
                                    }
                                }
                            }
                        }

                        /* Compute the coarse data and the residual encoding */
                        auto [coarse_value, patch_encoding] = PerformExtraction<T, DIM>(patch_data, paccess_idx, is_init_iteration);

                        /* Store the coarse value */
                        coarse_values.push_back(coarse_value);

                        /* Store the computed patch encoding */
                        this->levelwise_patch_encodings_.back().push_back(std::move(patch_encoding));
                    }
                }
            }
        }

        /* Switch to the coarse data for the next level extraction */
        this->data_ = std::move(coarse_values);

        cmc_debug_msg("The compression iteration step ", lvl_idx, " is finished.");
    }

    /* Encode the gathered data */
    this->EncodeData();

    cmc_debug_msg("The lossless compression of is finished.");
}


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::Compress(CompressionVariable<T, DIM>::kTag3D)
{
    static_assert(DIM == 3);
    constexpr int32_t kLevID = 0;
    constexpr int32_t kLatID = 1;
    constexpr int32_t kLonID = 2;

    int32_t paccess_idx{0};
    std::array<T, kPackSize<DIM>> patch_data{};

    for (int32_t lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx)
    {
        /* Allocate a new dimension array */
        this->dim_length_pyramid_.emplace_back();

        /* Get the current dimension lengths */
        const std::array<int32_t, DIM>& dim_lengths = *std::prev(this->dim_length_pyramid_.end(), 2);

        /* Get the next level dimension array */
        std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::prev(this->dim_length_pyramid_.end(), 1);
        next_lvl_dim_lengths[kLevID] = (dim_lengths[kLevID] % kDIMReductionFactor == 0 ? (dim_lengths[kLevID] / kDIMReductionFactor) : (dim_lengths[kLevID] / kDIMReductionFactor) + 1);
        next_lvl_dim_lengths[kLatID] = (dim_lengths[kLatID] % kDIMReductionFactor == 0 ? (dim_lengths[kLatID] / kDIMReductionFactor) : (dim_lengths[kLatID] / kDIMReductionFactor) + 1);
        next_lvl_dim_lengths[kLonID] = (dim_lengths[kLonID] % kDIMReductionFactor == 0 ? (dim_lengths[kLonID] / kDIMReductionFactor) : (dim_lengths[kLonID] / kDIMReductionFactor) + 1);

        /* Compute the number of next level values */
        const int32_t num_next_level_vals = next_lvl_dim_lengths[kLevID] * next_lvl_dim_lengths[kLatID] * next_lvl_dim_lengths[kLonID];
        /* Allocate a next level vector for the coarse data */
        std::vector<T> coarse_values;
        coarse_values.reserve(num_next_level_vals);

        /* Allocate a level-wise patch encdoing */
        this->levelwise_patch_encodings_.emplace_back();
        this->levelwise_patch_encodings_.back().reserve(num_next_level_vals);

        /* Only the init iteration is encoded with tolerances */
        const bool is_init_iteration = (kTryToCorrectMeanFastCompression ? (lvl_idx == 0) : false);

        /* Extract the coarse level values */
        for (int32_t lev = 0; lev < dim_lengths[kLevID]; lev += kDIMReductionFactor)
        {
            for (int32_t lat = 0; lat < dim_lengths[kLatID]; lat += kDIMReductionFactor)
            {
                for (int32_t lon = 0; lon < dim_lengths[kLonID]; lon += kDIMReductionFactor)
                {
                    /* Reset the access idx for each new patch */
                    paccess_idx = 0;

                    /* Gather the values for this initial patch */
                    for (int32_t lev_idx = 0; lev_idx < kDIMReductionFactor; ++lev_idx)
                    {
                        for (int32_t lat_idx = 0; lat_idx < kDIMReductionFactor; ++lat_idx)
                        {
                            for (int32_t lon_idx = 0; lon_idx < kDIMReductionFactor; ++lon_idx)
                            {
                                if (lev + lev_idx >= dim_lengths[kLevID] || lat + lat_idx >= dim_lengths[kLatID] || lon + lon_idx >= dim_lengths[kLonID]) [[unlikely]]
                                {
                                    continue;
                                } else
                                {
                                    /* Retrieve the corresponding value from this patch */
                                    patch_data[paccess_idx] = GetValue<T>(this->data_, lev + lev_idx, lat + lat_idx, lon + lon_idx,
                                                                          dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID]);
                                }
                                ++paccess_idx;
                            }
                        }
                    }

                    /* Compute the coarse data and the residual encoding */
                    auto [coarse_value, patch_encoding] = PerformExtraction<T, DIM>(patch_data, paccess_idx, is_init_iteration);

                    /* Store the coarse value */
                    coarse_values.push_back(coarse_value);

                    /* Store the computed patch encoding */
                    this->levelwise_patch_encodings_.back().push_back(std::move(patch_encoding));
                }
            }
        }

        /* Switch to the coarse data for the next level extraction */
        this->data_ = std::move(coarse_values);

        cmc_debug_msg("The compression iteration step ", lvl_idx, " is finished.");
    }

    /* Encode the gathered data */
    this->EncodeData();

    cmc_debug_msg("The lossless compression of is finished.");
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::Compress(CompressionVariable<T, DIM>::kTag2D)
{
    static_assert(DIM == 2);
    constexpr int32_t kLatID = 0;
    constexpr int32_t kLonID = 1;

    int32_t paccess_idx{0};
    std::array<T, kPackSize<DIM>> patch_data{};

    for (int32_t lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx)
    {
        /* Allocate a new dimension array */
        this->dim_length_pyramid_.emplace_back();

        /* Get the current dimension lengths */
        const std::array<int32_t, DIM>& dim_lengths = *std::prev(this->dim_length_pyramid_.end(), 2);

        /* Get the next level dimension array */
        std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::prev(this->dim_length_pyramid_.end(), 1);
        next_lvl_dim_lengths[kLatID] = (dim_lengths[kLatID] % kDIMReductionFactor == 0 ? (dim_lengths[kLatID] / kDIMReductionFactor) : (dim_lengths[kLatID] / kDIMReductionFactor) + 1);
        next_lvl_dim_lengths[kLonID] = (dim_lengths[kLonID] % kDIMReductionFactor == 0 ? (dim_lengths[kLonID] / kDIMReductionFactor) : (dim_lengths[kLonID] / kDIMReductionFactor) + 1);

        /* Compute the number of next level values */
        const int32_t num_next_level_vals = next_lvl_dim_lengths[kLatID] * next_lvl_dim_lengths[kLonID];
        /* Allocate a next level vector for the coarse data */
        std::vector<T> coarse_values;
        coarse_values.reserve(num_next_level_vals);

        /* Allocate a level-wise patch encdoing */
        this->levelwise_patch_encodings_.emplace_back();
        this->levelwise_patch_encodings_.back().reserve(num_next_level_vals);

        /* Only the init iteration is encoded with tolerances */
        const bool is_init_iteration = (kTryToCorrectMeanFastCompression ? (lvl_idx == 0) : false);

        /* Extract the coarse level values */
        for (int32_t lat = 0; lat < dim_lengths[kLatID]; lat += kDIMReductionFactor)
        {
            for (int32_t lon = 0; lon < dim_lengths[kLonID]; lon += kDIMReductionFactor)
            {
                /* Reset the access idx for each new patch */
                paccess_idx = 0;

                /* Gather the values for this initial patch */
                for (int32_t lat_idx = 0; lat_idx < kDIMReductionFactor; ++lat_idx)
                {
                    for (int32_t lon_idx = 0; lon_idx < kDIMReductionFactor; ++lon_idx)
                    {
                        if (lat + lat_idx >= dim_lengths[kLatID] || lon + lon_idx >= dim_lengths[kLonID]) [[unlikely]]
                        {
                            continue;
                        } else
                        {
                            /* Retrieve the corresponding value from this patch */
                            patch_data[paccess_idx] = GetValue<T>(this->data_, lat + lat_idx, lon + lon_idx,
                                                                  dim_lengths[kLonID], dim_lengths[kLatID]);
                        }
                        ++paccess_idx;
                    }
                }

                /* Compute the coarse data and the residual encoding */
                auto [coarse_value, patch_encoding] = PerformExtraction<T, DIM>(patch_data, paccess_idx, is_init_iteration);

                /* Store the coarse value */
                coarse_values.push_back(coarse_value);

                /* Store the computed patch encoding */
                this->levelwise_patch_encodings_.back().push_back(std::move(patch_encoding));
            }
        }

        /* Switch to the coarse data for the next level extraction */
        this->data_ = std::move(coarse_values);

        cmc_debug_msg("The compression iteration step ", lvl_idx, " is finished.");
    }

    /* Encode the gathered data */
    this->EncodeData();

    cmc_debug_msg("The lossless compression of is finished.");
}


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::Compress(CompressionVariable<T, DIM>::kTag1D)
{
    static_assert(DIM == 1);
    constexpr int32_t kLonID = 0;

    int32_t paccess_idx{0};
    std::array<T, kPackSize<DIM>> patch_data{};

    for (int32_t lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx)
    {
        /* Allocate a new dimension array */
        this->dim_length_pyramid_.emplace_back();

        /* Get the current dimension lengths */
        const std::array<int32_t, DIM>& dim_lengths = *std::prev(this->dim_length_pyramid_.end(), 2);

        /* Get the next level dimension array */
        std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::prev(this->dim_length_pyramid_.end(), 1);
        next_lvl_dim_lengths[kLonID] = (dim_lengths[kLonID] % kDIMReductionFactor == 0 ? (dim_lengths[kLonID] / kDIMReductionFactor) : (dim_lengths[kLonID] / kDIMReductionFactor) + 1);

        /* Compute the number of next level values */
        const int32_t num_next_level_vals = next_lvl_dim_lengths[kLonID];
        /* Allocate a next level vector for the coarse data */
        std::vector<T> coarse_values;
        coarse_values.reserve(num_next_level_vals);

        /* Allocate a level-wise patch encdoing */
        this->levelwise_patch_encodings_.emplace_back();
        this->levelwise_patch_encodings_.back().reserve(num_next_level_vals);

        /* Only the init iteration is encoded with tolerances */
        const bool is_init_iteration = (kTryToCorrectMeanFastCompression ? (lvl_idx == 0) : false);

        /* Extract the coarse level values */
        for (int32_t lon = 0; lon < dim_lengths[kLonID]; lon += kDIMReductionFactor)
        {
            /* Reset the access idx for each new patch */
            paccess_idx = 0;

            /* Gather the values for this initial patch */
            for (int32_t lon_idx = 0; lon_idx < kDIMReductionFactor; ++lon_idx)
            {
                if (lon + lon_idx >= dim_lengths[kLonID]) [[unlikely]]
                {
                    continue;
                } else
                {
                    /* Retrieve the corresponding value from this patch */
                    patch_data[paccess_idx] = GetValue<T>(this->data_, lon + lon_idx,
                                                          dim_lengths[kLonID]);
                }
                ++paccess_idx;
            }

            /* Compute the coarse data and the residual encoding */
            auto [coarse_value, patch_encoding] = PerformExtraction<T, DIM>(patch_data, paccess_idx, is_init_iteration);

            /* Store the coarse value */
            coarse_values.push_back(coarse_value);

            /* Store the computed patch encoding */
            this->levelwise_patch_encodings_.back().push_back(std::move(patch_encoding));
        }

        /* Switch to the coarse data for the next level extraction */
        this->data_ = std::move(coarse_values);

        cmc_debug_msg("The compression iteration step ", lvl_idx, " is finished.");
    }

    /* Encode the gathered data */
    this->EncodeData();

    cmc_debug_msg("The lossless compression of is finished.");
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>>
CollectEntropySymbols(const std::vector<std::vector<PatchEncoding<T, DIM>>>& levelwise_encodings)
{
    /* Get the number of all possible entropy symbols */
    constexpr int num_entropy_symbols = GetNumEntropySymbols<T>();

    /* Set the array and zero intialiaze the frequencies */
    std::array<uint64_t, num_entropy_symbols> entropy_symbol_frequencies{};

    /* Iterate through all entropy codes and accumulate their frequencies */
    for (size_t lvl_idx{0}; lvl_idx < levelwise_encodings.size(); ++lvl_idx)
    {
        /* Iterate through all coarsening data on this level */
        for (size_t coarsening_idx{0}; coarsening_idx < levelwise_encodings[lvl_idx].size(); ++coarsening_idx)
        {
            const uint32_t num_elems = ((kTryToCorrectMeanFastCompression && lvl_idx == 0) ? (levelwise_encodings[lvl_idx][coarsening_idx].num_elements > 1 ? levelwise_encodings[lvl_idx][coarsening_idx].num_elements - 1 : 1) : levelwise_encodings[lvl_idx][coarsening_idx].num_elements);

            /* Iterate over all entropy codes from this coarsening data */
            for (uint32_t entropy_sym_idx{0}; entropy_sym_idx < num_elems; ++entropy_sym_idx)
            {
                /* Convert the symbol to the corresponding array index */
                const int array_idx = MapEntropySymbolToArrayIndex<T>(levelwise_encodings[lvl_idx][coarsening_idx].entropy_symbols[entropy_sym_idx]);
                /* Update the frequency */
                ++entropy_symbol_frequencies[array_idx];
            }
        }
    }

    /* Add the process end symbol */
    AddProcessEndSymbol<T>(entropy_symbol_frequencies, 1);

    /* Replicate the global entropy frequencies over all levels */
    std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>> global_symbol_frequencies;
    global_symbol_frequencies.reserve(num_entropy_symbols);

    for (int idx{0}; idx < num_entropy_symbols; ++idx)
    {
        /* Convert the index back to the entropy symbol */
        const SymbolType entropy_symbol = MapArrayIndexToEntropySymbol<T>(idx);

        /* Store the symbol with the global frequency */
        global_symbol_frequencies.emplace_back(entropy_symbol, entropy_symbol_frequencies[idx]);
    }

    return global_symbol_frequencies;
}

template<OneByteArithmeticType T>
inline std::vector<uint64_t>
PerformRootLevelEncoding(const std::vector<T>& data_)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(data_.size() * sizeof(T));
    
    for (size_t idx{0}; idx < data_.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<OneByteResidualType>(std::bit_cast<OneByteResidualType>(data_[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}

template<TwoByteArithmeticType T>
inline std::vector<uint64_t>
PerformRootLevelEncoding(const std::vector<T>& data_)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(data_.size() * sizeof(T));
    
    for (size_t idx{0}; idx < data_.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<TwoByteResidualType>(std::bit_cast<TwoByteResidualType>(data_[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}

template<FourByteArithmeticType T>
inline std::vector<uint64_t>
PerformRootLevelEncoding(const std::vector<T>& data_)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(data_.size() * sizeof(T));
    
    for (size_t idx{0}; idx < data_.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<FourByteResidualType>(std::bit_cast<FourByteResidualType>(data_[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}

template<EightByteArithmeticType T>
inline std::vector<uint64_t>
PerformRootLevelEncoding(const std::vector<T>& data_)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(data_.size() * sizeof(T));
    
    for (size_t idx{0}; idx < data_.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<EightByteResidualType>(std::bit_cast<EightByteResidualType>(data_[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
std::vector<uint64_t>
CompressionVariable<T, DIM>::EncodeRootLevelData() const
{   
    /* The lastly added vector to the data pyramid resembles the root level */
    return PerformRootLevelEncoding<T>(this->data_);
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::EncodeData()
{
    /* Collect and exchange all entropy symbols */
    const std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>> entropy_symbols = CollectEntropySymbols<T>(this->levelwise_patch_encodings_);
    
    /* Create a Huffman encoder */
    cmc::entropy_coding::huffman::HuffmanCoder<SymbolType> entropy_coder(entropy_symbols);

    /* Store the serialized Huffman coder */
    this->serialized_entropy_dictionary_ = entropy_coder.SerializeHuffmanCodesBEPadded();

    /* Number of overall encoding steps */
    const int num_encoding_steps = this->num_compression_lvls_ + 1;

    /* Allocate an output vector for the levelwise encoding */
    this->levelwise_encoded_data_.reserve(num_encoding_steps);

    /* We need to encode the data resididng on the root level */
    this->levelwise_encoded_data_.push_back(this->EncodeRootLevelData());

    /* We encode the data from the root level to the leaf level */
    auto enc_iter = this->levelwise_patch_encodings_.rbegin();

    /* Iterate over all compression levels (Without the root level and the intra element level) */
    for (int32_t step_idx{1}; step_idx < num_encoding_steps; ++step_idx, ++enc_iter)
    {
        /* Allocate a bits::vector to store this level's encoded data */
        cmc::bits::vector lvl_data;
        lvl_data.Reserve((enc_iter->size() * sizeof(PatchEncoding<T, DIM>) * cmc::bits::kCharBit) / 2);

        /* Define a reference on the coarsening data for the ease of notation */
        const std::vector<PatchEncoding<T, DIM>>& lvl_encoding_data = *enc_iter;

        /* Number of data on this level */
        const int32_t num_elems = lvl_encoding_data.size();

        /* Iterate over this level's refinement indications */
        for (int32_t elem_idx{0}; elem_idx < num_elems; ++elem_idx)
        {
            /* Define a reference for the ease of notation */
            const PatchEncoding<T, DIM>& enc_data = lvl_encoding_data[elem_idx];

            /* Encode the quantization bins and potentially interleave the unpredicted values */
            const uint32_t num_elems = ((kTryToCorrectMeanFastCompression && step_idx == num_encoding_steps - 1) ? (enc_data.num_elements > 1 ? enc_data.num_elements - 1 : 1) : enc_data.num_elements);

            for (int val_idx{0}; val_idx < static_cast<int>(num_elems); ++val_idx)
            {
                /* Encode the entropy symbol */
                const cmc::entropy_coding::huffman::HuffmanCode code = entropy_coder.EncodeSymbol(enc_data.entropy_symbols[val_idx]);

                /* Retrieve the LZC from the entropy symbol */
                const int lzc = GetLZCFromEntropySymbol(enc_data.entropy_symbols[val_idx]);

                /* Serialize the encoded entropy symbol */
                lvl_data.AppendBits(code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - code.code_length), 0);
            
                /* We do not need to encode the implicit given one-bit following the LZC */
                if (lzc + 1 < static_cast<int>(sizeof(T) * cmc::bits::kCharBit)) [[likely]]
                {
                    /* Append the significant reisdual bits */
                    lvl_data.AppendBits(enc_data.residuals[val_idx], lzc + 1, 0);
                }

            }
        }

        /* At the end of the local encoding of the level, we append the process-end symbol */
        const cmc::entropy_coding::huffman::HuffmanCode process_lvl_end_code = entropy_coder.EncodeSymbol(kProcessEndSymbol<T>);

        /* Serialize the encoded process end symbol */
        lvl_data.AppendBits(process_lvl_end_code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - process_lvl_end_code.code_length), 0);

        /* We store this level's encoding in the variable's buffer */
        this->levelwise_encoded_data_.push_back(lvl_data.GetSerializedByteStreamBE());
    }
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::WriteCompressedData(const std::string& file_name)
{
    /* Check if the output file exists, if so, we delete it */
    const std::filesystem::path output_file_path(file_name);
    if (std::filesystem::exists(output_file_path))
    {
        std::filesystem::remove(output_file_path);
    }

    const uint64_t num_header_vals = 7 + this->levelwise_encoded_data_.size() + this->levelwise_encoded_data_.size() * DIM;
    const uint64_t header_size = sizeof(uint64_t) * num_header_vals;

    /* We construct the output stream */
    std::vector<uint64_t> header;
    header.reserve(num_header_vals);

    /* Global Bytes compressed variable */
    uint64_t level_bytes{0};
    for (size_t lvl_idx{0}; lvl_idx < this->levelwise_encoded_data_.size(); ++lvl_idx)
    {
        level_bytes += this->levelwise_encoded_data_[lvl_idx].size() * sizeof(uint64_t);
    }
    const uint64_t global_byte_count = header_size + this->serialized_entropy_dictionary_.size() * sizeof(uint64_t) + level_bytes;
    header.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(global_byte_count));

    /* Store the offset from the beginning to the start of the root level encoding */
    const uint64_t level_encoding_byte_offset = header_size;
    header.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(level_encoding_byte_offset));

    /* Store the data type */
    header.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(ConvertToCmcType<T>())));

    /* Store the dimensionality */
    header.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(DIM)));

    /* Store the compression scheme */
    header.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(CompressionSchema::ParallelMultiResExtraction)));

    /* Store the number of compression levels */
    header.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(this->levelwise_encoded_data_.size())));

    /* Append the global level bytes */
    const int32_t num_compr_level = static_cast<int32_t>(this->levelwise_encoded_data_.size());
    for (int32_t lvl_idx{0}; lvl_idx < num_compr_level; ++lvl_idx)
    {
        header.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(this->levelwise_encoded_data_[lvl_idx].size() * sizeof(uint64_t))));
    }

    /* Append the dimensionality pyramid (it needs to be traversed in reverse order) */
    const int32_t dim_pyra_size = static_cast<int32_t>(this->dim_length_pyramid_.size());
    for (size_t lvl_idx{0}; lvl_idx < this->dim_length_pyramid_.size(); ++lvl_idx)
    {
        for (int32_t dim_idx{0}; dim_idx < DIM; ++dim_idx)
        {
            header.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(this->dim_length_pyramid_[dim_pyra_size - 1 - lvl_idx][dim_idx])));
        }
    }

    /* Store the Huffman codes */
    const uint64_t num_bytes_huffman_codes = this->serialized_entropy_dictionary_.size() * sizeof(uint64_t);
    header.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(num_bytes_huffman_codes));

    /* Open the output file */
    std::FILE* file_out = std::fopen(file_name.c_str(), "wb");

    /* Write the header */
    std::fwrite(header.data(), sizeof(uint64_t), header.size(), file_out);

    /* Write the Huffman Codes */
    std::fwrite(this->serialized_entropy_dictionary_.data(), sizeof(uint64_t), this->serialized_entropy_dictionary_.size(), file_out);

    /* Write the levelwise encoding streams */
    for (size_t lvl_idx{0}; lvl_idx < this->levelwise_encoded_data_.size(); ++lvl_idx)
    {
        std::fwrite(this->levelwise_encoded_data_[lvl_idx].data(), sizeof(uint64_t), this->levelwise_encoded_data_[lvl_idx].size(), file_out);
    }

    /* And, finally, close the file */
    std::fclose(file_out);
}

}

#endif /* !CMC_PATCH_LOSSLESS_FAST_MULTI_RES_EXTRACTION_WITH_CORRECTION_HXX */
