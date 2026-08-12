#ifndef CMC_PATCH_LOSSY_MULTI_RES_RBF_HXX
#define CMC_PATCH_LOSSY_MULTI_RES_RBF_HXX

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

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
class CompressionVariable
{
public:

    CompressionVariable() = delete;
    CompressionVariable(const std::span<const T> data, const std::array<int32_t, DIM>& dimension_lengths, const float global_abs_permitted_error)
    : init_data_(data), dim_lengths_(dimension_lengths), abs_error_{global_abs_permitted_error},
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
        this->data_pyramid_.reserve(this->num_compression_lvls_ + 1);
        this->data_pyramid_.emplace_back();
        this->data_pyramid_.back().reserve(num_data);
        std::copy_n(data.begin(), num_data, std::back_inserter(this->data_pyramid_.back()));
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
    void CollectCoarseDataPyramid(kTag1D);
    void CollectCoarseDataPyramid(kTag2D);
    void CollectCoarseDataPyramid(kTag3D);
    void CollectCoarseDataPyramid(kTag4D);
    PatchEncoding<T, DIM> PerformPrediction(const std::array<T, k1DNumControlValues> control_values, const std::array<T, kPackSize<DIM>>& init_data, const float permitted_abs_error, const int32_t lon, std::vector<T>& next_level_data, const std::array<int32_t, DIM>& next_lvl_dim_lengths);
    PatchEncoding<T, DIM> PerformPrediction(const std::array<T, k2DNumControlValues> control_values, const std::array<T, kPackSize<DIM>>& init_data, const float permitted_abs_error, const int32_t lon, const int32_t lat, std::vector<T>& next_level_data, const std::array<int32_t, DIM>& next_lvl_dim_lengths);
    PatchEncoding<T, DIM> PerformPrediction(const std::array<T, k3DNumControlValues> control_values, const std::array<T, kPackSize<DIM>>& init_data, const float permitted_abs_error, const int32_t lon, const int32_t lat, const int32_t lev, std::vector<T>& next_level_data, const std::array<int32_t, DIM>& next_lvl_dim_lengths);
    PatchEncoding<T, DIM> PerformPrediction(const std::array<T, k4DNumControlValues> control_values, const std::array<T, kPackSize<DIM>>& init_data, const float permitted_abs_error, const int32_t lon, const int32_t lat, const int32_t lev, const int32_t time, std::vector<T>& next_level_data, const std::array<int32_t, DIM>& next_lvl_dim_lengths);
    std::vector<uint64_t> EncodeRootLevelData() const;
    void EncodeData();

    const std::span<const T> init_data_;
    const std::array<int32_t, DIM> dim_lengths_;
    const float abs_error_;

    const int32_t max_dimension_length_{0};
    const int32_t num_compression_lvls_{0};

    std::vector<T> data_;
    std::vector<std::array<int32_t, DIM>> dim_length_pyramid_;
    std::vector<std::vector<T>> data_pyramid_;

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
        cmc_global_msg("Lossy 1D Compression");
        this->Compress(CompressionVariable<T, DIM>::tag1D);
    } else if constexpr (DIM == 2)
    {
        cmc_global_msg("Lossy 2D Compression");
        this->Compress(CompressionVariable<T, DIM>::tag2D);
    } else if constexpr (DIM == 3)
    {
        cmc_global_msg("Lossy 3D Compression");
        this->Compress(CompressionVariable<T, DIM>::tag3D);
    } else if constexpr (DIM == 4)
    {
        cmc_global_msg("Lossy 4D Compression");
        this->Compress(CompressionVariable<T, DIM>::tag4D);
    } else
    {
        cmc_err_msg("Unsupported variable's dimensionality (DIM = ", DIM, ").");
    }
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::CollectCoarseDataPyramid(CompressionVariable<T, DIM>::kTag4D)
{
    static_assert(DIM == 4);
    constexpr int32_t kTimeID = 0;
    constexpr int32_t kLevID = 1;
    constexpr int32_t kLatID = 2;
    constexpr int32_t kLonID = 3;

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

        /* Allocate a next level vector for the coarse data */
        std::vector<T> coarse_values(next_lvl_dim_lengths[kTimeID] * next_lvl_dim_lengths[kLevID] * next_lvl_dim_lengths[kLatID] * next_lvl_dim_lengths[kLonID]);

        /* Extract the coarse level values */
        for (int32_t time = 0; time < dim_lengths[kTimeID]; time += kDIMReductionFactor)
        {
            for (int32_t lev = 0; lev < dim_lengths[kLevID]; lev += kDIMReductionFactor)
            {
                for (int32_t lat = 0; lat < dim_lengths[kLatID]; lat += kDIMReductionFactor)
                {
                    for (int32_t lon = 0; lon < dim_lengths[kLonID]; lon += kDIMReductionFactor)
                    {
                        /* Get the patch's first value */
                        const T coarse_patch_value = GetValue<T>(this->data_pyramid_.back(), time, lev, lat, lon, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID], dim_lengths[kTimeID]);

                        /* Store the coarse value */
                        SetValue<T>(coarse_values, coarse_patch_value, time / 2, lev / 2, lat / 2, lon / 2, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID], next_lvl_dim_lengths[kTimeID]);
                    }
                }
            }
        }

        /* Store this level's coarse values */
        this->data_pyramid_.push_back(std::move(coarse_values));
    }
    
    /* Now, the data pyramid and the dimension length pyramid have been filled accordingly */
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::CollectCoarseDataPyramid(CompressionVariable<T, DIM>::kTag3D)
{
    static_assert(DIM == 3);
    constexpr int32_t kLevID = 0;
    constexpr int32_t kLatID = 1;
    constexpr int32_t kLonID = 2;

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

        /* Allocate a next level vector for the coarse data */
        std::vector<T> coarse_values(next_lvl_dim_lengths[kLevID] * next_lvl_dim_lengths[kLatID] * next_lvl_dim_lengths[kLonID]);

        /* Extract the coarse level values */
        for (int32_t lev = 0; lev < dim_lengths[kLevID]; lev += kDIMReductionFactor)
        {
            for (int32_t lat = 0; lat < dim_lengths[kLatID]; lat += kDIMReductionFactor)
            {
                for (int32_t lon = 0; lon < dim_lengths[kLonID]; lon += kDIMReductionFactor)
                {
                    /* Get the patch's first value */
                    const T coarse_patch_value = GetValue<T>(this->data_pyramid_.back(), lev, lat, lon, dim_lengths[kLonID], dim_lengths[kLatID], dim_lengths[kLevID]);

                    /* Store the coarse value */
                    SetValue<T>(coarse_values, coarse_patch_value, lev / 2, lat / 2, lon / 2, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID]);
                }
            }
        }

        /* Store this level's coarse values */
        this->data_pyramid_.push_back(std::move(coarse_values));
    }
    
    /* Now, the data pyramid and the dimension length pyramid have been filled accordingly */
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::CollectCoarseDataPyramid(CompressionVariable<T, DIM>::kTag2D)
{
    static_assert(DIM == 2);
    constexpr int32_t kLatID = 0;
    constexpr int32_t kLonID = 1;

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

        /* Allocate a next level vector for the coarse data */
        std::vector<T> coarse_values(next_lvl_dim_lengths[kLatID] * next_lvl_dim_lengths[kLonID]);

        /* Extract the coarse level values */
        for (int32_t lat = 0; lat < dim_lengths[kLatID]; lat += kDIMReductionFactor)
        {
            for (int32_t lon = 0; lon < dim_lengths[kLonID]; lon += kDIMReductionFactor)
            {
                /* Get the patch's first value */
                const T coarse_patch_value = GetValue<T>(this->data_pyramid_.back(), lat, lon, dim_lengths[kLonID], dim_lengths[kLatID]);

                /* Store the coarse value */
                SetValue<T>(coarse_values, coarse_patch_value, lat / kDIMReductionFactor, lon / kDIMReductionFactor, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID]);
            }
        }

        /* Store this level's coarse values */
        this->data_pyramid_.push_back(std::move(coarse_values));
    }
    
    /* Now, the data pyramid and the dimension length pyramid have been filled accordingly */
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::CollectCoarseDataPyramid(CompressionVariable<T, DIM>::kTag1D)
{
    static_assert(DIM == 1);
    constexpr int32_t kLonID = 0;

    for (int32_t lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx)
    {
        /* Allocate a new dimension array */
        this->dim_length_pyramid_.emplace_back();

        /* Get the current dimension lengths */
        const std::array<int32_t, DIM>& dim_lengths = *std::prev(this->dim_length_pyramid_.end(), 2);

        /* Get the next level dimension array */
        std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::prev(this->dim_length_pyramid_.end(), 1);
        next_lvl_dim_lengths[kLonID] = (dim_lengths[kLonID] % kDIMReductionFactor == 0 ? (next_lvl_dim_lengths[kLonID] / kDIMReductionFactor) : (next_lvl_dim_lengths[kLonID] / kDIMReductionFactor) + 1);

        /* Allocate a next level vector for the coarse data */
        std::vector<T> coarse_values(next_lvl_dim_lengths[kLonID]);

        /* Extract the coarse level values */
        for (int32_t lon = 0; lon < dim_lengths[kLonID]; lon += kDIMReductionFactor)
        {
            /* Get the patch's first value */
            const T coarse_patch_value = GetValue<T>(this->data_pyramid_.back(), lon, dim_lengths[kLonID]);

            /* Store the coarse value */
            SetValue<T>(coarse_values, coarse_patch_value, lon / 2, next_lvl_dim_lengths[kLonID]);
        }

        /* Store this level's coarse values */
        this->data_pyramid_.push_back(std::move(coarse_values));
    }
    
    /* Now, the data pyramid and the dimension length pyramid have been filled accordingly */
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
PatchEncoding<T, DIM>
CompressionVariable<T, DIM>::PerformPrediction(const std::array<T, k4DNumControlValues> control_values, const std::array<T, kPackSize<DIM>>& init_data, const float permitted_abs_error, const int32_t lon, const int32_t lat, const int32_t lev, const int32_t time, std::vector<T>& next_level_data, const std::array<int32_t, DIM>& next_lvl_dim_lengths)
{
    static_assert(DIM == 4);
    constexpr int32_t kTimeID = 0;
    constexpr int32_t kLevID = 1;
    constexpr int32_t kLatID = 2;
    constexpr int32_t kLonID = 3;

    /* Perform the prediction */
    std::array<T, kPackSize<DIM>> predicted_data = Perform4DRBFPrediction<T, DIM>(control_values);

    /* Allocate the encoding struct */
    PatchEncoding<T, DIM> encoding;

    int32_t init_access_idx{0};
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
                       /* Compute the residual between prediction and actual value */
                       const T residual = GetAbsResidual<T>(init_data[init_access_idx], predicted_data[pred_access_idx]);

                        /* Determine whether the prediction is in line with the permitted error */
                        if (static_cast<float>(residual) <= permitted_abs_error)
                        {
                            /* If the residual is within the permitted error bound */
                            encoding.quantization_bins[init_access_idx] = kPredictionWithinBound;

                            /* We store the prediction */
                            SetValue<T>(next_level_data, predicted_data[pred_access_idx], kDIMReductionFactor * time + time_idx, kDIMReductionFactor * lev + lev_idx,
                                        kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID], next_lvl_dim_lengths[kTimeID]);

                        } else if (static_cast<float>(residual) <= kCheckResidualMaxDeviationFactor * permitted_abs_error)
                        {
                            /* If the residual is within the permitted error interval such that the error can be met by quantization */
                            const SymbolType bin = ComputeBin(static_cast<float>(residual), permitted_abs_error);

                            /* Check in which direction the quantization goes */
                            const bool is_prediction_greater = (predicted_data[pred_access_idx] >= init_data[init_access_idx]);

                            /* Create the entropy symbol from the information above */
                            const SymbolType entropy_symbol = CreateEntropySymbolFromQuantizationBin(bin, is_prediction_greater); 

                            /* In this case, we do not need to store anything apart the quantization bin */
                            encoding.quantization_bins[init_access_idx] = entropy_symbol;

                            /* De-Quantize the value */
                            const T decompressed_value = DequantizeValue<T>(predicted_data[pred_access_idx], permitted_abs_error, is_prediction_greater, bin);

                            /* We store the de-quantized value */
                            SetValue<T>(next_level_data, decompressed_value, kDIMReductionFactor * time + time_idx, kDIMReductionFactor * lev + lev_idx,
                                        kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID], next_lvl_dim_lengths[kTimeID]);

                        } else
                        {
                            /* In this case, the value is unpredictable and we store it, as it is */
                            /* We flag the value as unpredictable */
                            encoding.quantization_bins[init_access_idx] = kFlagUnpredictable;
                            
                            /* And we store the actual value */
                            encoding.unpredictable_values[init_access_idx] = init_data[init_access_idx];

                            SetValue<T>(next_level_data, init_data[init_access_idx], kDIMReductionFactor * time + time_idx, kDIMReductionFactor * lev + lev_idx,
                                        kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID], next_lvl_dim_lengths[kTimeID]);
                        }

                        
                        ++init_access_idx;
                    }
                }
            }
        }
    }

    /* Store the number of elements for this encoding strtuct */
    encoding.num_elements = init_access_idx;

    return encoding;
}


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
PatchEncoding<T, DIM>
CompressionVariable<T, DIM>::PerformPrediction(const std::array<T, k3DNumControlValues> control_values, const std::array<T, kPackSize<DIM>>& init_data, const float permitted_abs_error, const int32_t lon, const int32_t lat, const int32_t lev, std::vector<T>& next_level_data, const std::array<int32_t, DIM>& next_lvl_dim_lengths)
{
    static_assert(DIM == 3);
    constexpr int32_t kLevID = 0;
    constexpr int32_t kLatID = 1;
    constexpr int32_t kLonID = 2;

    /* Perform the prediction */
    std::array<T, kPackSize<DIM>> predicted_data = Perform3DRBFPrediction<T, DIM>(control_values);

    /* Allocate the encoding struct */
    PatchEncoding<T, DIM> encoding;

    int32_t init_access_idx{0};
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
                   /* Compute the residual between prediction and actual value */
                   const T residual = GetAbsResidual<T>(init_data[init_access_idx], predicted_data[pred_access_idx]);

                    /* Determine whether the prediction is in line with the permitted error */
                    if (static_cast<float>(residual) <= permitted_abs_error)
                    {
                        /* If the residual is within the permitted error bound */
                        encoding.quantization_bins[init_access_idx] = kPredictionWithinBound;

                        /* We store the prediction */
                        SetValue<T>(next_level_data, predicted_data[pred_access_idx], kDIMReductionFactor * lev + lev_idx,
                                    kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID]);

                    } else if (static_cast<float>(residual) <= kCheckResidualMaxDeviationFactor * permitted_abs_error)
                    {
                        /* If the residual is within the permitted error interval such that the error can be met by quantization */
                        const SymbolType bin = ComputeBin(static_cast<float>(residual), permitted_abs_error);

                        /* Check in which direction the quantization goes */
                        const bool is_prediction_greater = (predicted_data[pred_access_idx] >= init_data[init_access_idx]);

                        /* Create the entropy symbol from the information above */
                        const SymbolType entropy_symbol = CreateEntropySymbolFromQuantizationBin(bin, is_prediction_greater); 

                        /* In this case, we do not need to store anything apart the quantization bin */
                        encoding.quantization_bins[init_access_idx] = entropy_symbol;

                        /* De-Quantize the value */
                        const T decompressed_value = DequantizeValue<T>(predicted_data[pred_access_idx], permitted_abs_error, is_prediction_greater, bin);

                        /* We store the de-quantized value */
                        SetValue<T>(next_level_data, decompressed_value, kDIMReductionFactor * lev + lev_idx,
                                    kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID]);

                    } else
                    {
                        /* In this case, the value is unpredictable and we store it, as it is */
                        /* We flag the value as unpredictable */
                        encoding.quantization_bins[init_access_idx] = kFlagUnpredictable;
                        
                        /* And we store the actual value */
                        encoding.unpredictable_values[init_access_idx] = init_data[init_access_idx];

                        SetValue<T>(next_level_data, init_data[init_access_idx], kDIMReductionFactor * lev + lev_idx,
                                    kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID]);
                    }

                    ++init_access_idx;
                }
            }
        }
    }

    /* Store the number of elements for this encoding strtuct */
    encoding.num_elements = init_access_idx;

    return encoding;
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
PatchEncoding<T, DIM>
CompressionVariable<T, DIM>::PerformPrediction(const std::array<T, k2DNumControlValues> control_values, const std::array<T, kPackSize<DIM>>& init_data, const float permitted_abs_error, const int32_t lon, const int32_t lat, std::vector<T>& next_level_data, const std::array<int32_t, DIM>& next_lvl_dim_lengths)
{
    static_assert(DIM == 2);
    constexpr int32_t kLatID = 0;
    constexpr int32_t kLonID = 1;

    /* Perform the prediction */
    std::array<T, kPackSize<DIM>> predicted_data = Perform2DRBFPrediction<T, DIM>(control_values);

    /* Allocate the encoding struct */
    PatchEncoding<T, DIM> encoding;

    int32_t init_access_idx{0};
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
               /* Compute the residual between prediction and actual value */
               const T residual = GetAbsResidual<T>(init_data[init_access_idx], predicted_data[pred_access_idx]);

                /* Determine whether the prediction is in line with the permitted error */
                if (static_cast<float>(residual) <= permitted_abs_error)
                {
                    /* If the residual is within the permitted error bound */
                    encoding.quantization_bins[init_access_idx] = kPredictionWithinBound;

                    /* We store the prediction */
                    SetValue<T>(next_level_data, predicted_data[pred_access_idx],
                                kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID]);

                } else if (static_cast<float>(residual) <= kCheckResidualMaxDeviationFactor * permitted_abs_error)
                {
                    /* If the residual is within the permitted error interval such that the error can be met by quantization */
                    const SymbolType bin = ComputeBin(static_cast<float>(residual), permitted_abs_error);

                    /* Check in which direction the quantization goes */
                    const bool is_prediction_greater = (predicted_data[pred_access_idx] >= init_data[init_access_idx]);

                    /* Create the entropy symbol from the information above */
                    const SymbolType entropy_symbol = CreateEntropySymbolFromQuantizationBin(bin, is_prediction_greater); 

                    /* In this case, we do not need to store anything apart the quantization bin */
                    encoding.quantization_bins[init_access_idx] = entropy_symbol;

                    /* De-Quantize the value */
                    const T decompressed_value = DequantizeValue<T>(predicted_data[pred_access_idx], permitted_abs_error, is_prediction_greater, bin);

                    /* We store the de-quantized value */
                    SetValue<T>(next_level_data, decompressed_value,
                                kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID]);

                } else
                {
                    /* In this case, the value is unpredictable and we store it, as it is */
                    /* We flag the value as unpredictable */
                    encoding.quantization_bins[init_access_idx] = kFlagUnpredictable;
                    
                    /* And we store the actual value */
                    encoding.unpredictable_values[init_access_idx] = init_data[init_access_idx];

                    SetValue<T>(next_level_data, init_data[init_access_idx],
                                kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID]);
                }

                ++init_access_idx;
            }
        }
    }

    /* Store the number of elements for this encoding strtuct */
    encoding.num_elements = init_access_idx;

    return encoding;
}


template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
PatchEncoding<T, DIM>
CompressionVariable<T, DIM>::PerformPrediction(const std::array<T, k1DNumControlValues> control_values, const std::array<T, kPackSize<DIM>>& init_data, const float permitted_abs_error, const int32_t lon, std::vector<T>& next_level_data, const std::array<int32_t, DIM>& next_lvl_dim_lengths)
{
    static_assert(DIM == 1);
    constexpr int32_t kLonID = 0;

    /* Perform the prediction */
    std::array<T, kPackSize<DIM>> predicted_data = Perform1DRBFPrediction<T, DIM>(control_values);

    /* Allocate the encoding struct */
    PatchEncoding<T, DIM> encoding;

    int32_t init_access_idx{0};
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
           /* Compute the residual between prediction and actual value */
           const T residual = GetAbsResidual<T>(init_data[init_access_idx], predicted_data[pred_access_idx]);

            /* Determine whether the prediction is in line with the permitted error */
            if (static_cast<float>(residual) <= permitted_abs_error)
            {
                /* If the residual is within the permitted error bound */
                encoding.quantization_bins[init_access_idx] = kPredictionWithinBound;

                /* We store the prediction */
                SetValue<T>(next_level_data, predicted_data[pred_access_idx],
                            kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID]);

            } else if (static_cast<float>(residual) <= kCheckResidualMaxDeviationFactor * permitted_abs_error)
            {
                /* If the residual is within the permitted error interval such that the error can be met by quantization */
                const SymbolType bin = ComputeBin(static_cast<float>(residual), permitted_abs_error);

                /* Check in which direction the quantization goes */
                const bool is_prediction_greater = (predicted_data[pred_access_idx] >= init_data[init_access_idx]);

                /* Create the entropy symbol from the information above */
                const SymbolType entropy_symbol = CreateEntropySymbolFromQuantizationBin(bin, is_prediction_greater); 

                /* In this case, we do not need to store anything apart the quantization bin */
                encoding.quantization_bins[init_access_idx] = entropy_symbol;

                /* De-Quantize the value */
                const T decompressed_value = DequantizeValue<T>(predicted_data[pred_access_idx], permitted_abs_error, is_prediction_greater, bin);

                /* We store the de-quantized value */
                SetValue<T>(next_level_data, decompressed_value,
                            kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID]);

            } else
            {
                /* In this case, the value is unpredictable and we store it, as it is */
                /* We flag the value as unpredictable */
                encoding.quantization_bins[init_access_idx] = kFlagUnpredictable;
                
                /* And we store the actual value */
                encoding.unpredictable_values[init_access_idx] = init_data[init_access_idx];

                SetValue<T>(next_level_data, init_data[init_access_idx],
                            kDIMReductionFactor * lon + lon_idx, next_lvl_dim_lengths[kLonID]);
            }

            ++init_access_idx;
        }
    }

    /* Store the number of elements for this encoding strtuct */
    encoding.num_elements = init_access_idx;

    return encoding;
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>>
CollectEntropySymbols(const std::vector<std::vector<PatchEncoding<T, DIM>>>& levelwise_encodings)
{
    /* Get the number of all possible entropy symbols */
    constexpr int num_entropy_symbols = GetNumEntropySymbols();

    /* Set the array and zero intialiaze the frequencies */
    std::array<uint64_t, num_entropy_symbols> entropy_symbol_frequencies{};

    /* Iterate through all entropy codes and accumulate their frequencies */
    for (size_t lvl_idx{0}; lvl_idx < levelwise_encodings.size(); ++lvl_idx)
    {
        /* Iterate through all coarsening data on this level */
        for (size_t coarsening_idx{0}; coarsening_idx < levelwise_encodings[lvl_idx].size(); ++coarsening_idx)
        {
            /* Iterate over all entropy codes from this coarsening data */
            for (int32_t entropy_sym_idx{0}; entropy_sym_idx < levelwise_encodings[lvl_idx][coarsening_idx].num_elements; ++entropy_sym_idx)
            {
                /* Convert the symbol to the corresponding array index */
                const int array_idx = MapEntropySymbolToArrayIndex<T>(levelwise_encodings[lvl_idx][coarsening_idx].quantization_bins[entropy_sym_idx]);
                /* Update the frequency */
                ++entropy_symbol_frequencies[array_idx];
            }
        }
    }

    /* Add the process end symbol */
    AddProcessEndSymbol(entropy_symbol_frequencies, 1);

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
    return PerformRootLevelEncoding<T>(this->data_pyramid_.back());
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
    auto enc_iter = this->levelwise_patch_encodings_.begin();

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
            for (int child_elem_idx{0}; child_elem_idx < static_cast<int>(enc_data.num_elements); ++child_elem_idx)
            {
                /* Encode the entropy symbol */
                const cmc::entropy_coding::huffman::HuffmanCode code = entropy_coder.EncodeSymbol(enc_data.quantization_bins[child_elem_idx]);

                /* Serialize the encoded entropy symbol */
                lvl_data.AppendBits(code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - code.code_length), 0);
            
                /* If the data is not predicatble, we need to store the actual value */
                if (enc_data.quantization_bins[child_elem_idx] == kFlagUnpredictable) [[unlikely]]
                {
                    /* Append the not-predictable value fully */
                    lvl_data.AppendBits(TransformToUInteger<T>(enc_data.unpredictable_values[child_elem_idx]), 0, 0);
                }
            }
        }

        /* At the end of the local encoding of the level, we append the process-end symbol */
        const cmc::entropy_coding::huffman::HuffmanCode process_lvl_end_code = entropy_coder.EncodeSymbol(kProcessEndSymbol);

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

    const uint64_t num_header_vals = 8 + this->levelwise_encoded_data_.size() + this->levelwise_encoded_data_.size() * DIM;
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

    /* Store the permitted absolute error */
    const uint64_t abs_permitted_error = static_cast<uint64_t>(TransformToUInteger<float>(this->abs_error_));
    header.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(abs_permitted_error));

    /* Store the number of compression levels */
    header.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(this->levelwise_encoded_data_.size())));

    /* Append the global level bytes */
    for (size_t lvl_idx{0}; lvl_idx < this->levelwise_encoded_data_.size(); ++lvl_idx)
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

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::Compress(CompressionVariable<T, DIM>::kTag4D)
{
    constexpr int32_t kTimeID = 0;
    constexpr int32_t kLevID = 1;
    constexpr int32_t kLatID = 2;
    constexpr int32_t kLonID = 3;

    cmc_debug_msg("Lossy patch-based compression starts...");

    cmc_assert(this->dim_length_pyramid_.size() == 1ULL);

    /* Build up the coarse prediction pyramid */
    this->CollectCoarseDataPyramid(CompressionVariable<T, DIM>::tag4D);

    cmc_debug_msg("Number of compression iterations to be performed: ", this->num_compression_lvls_);
    
    int32_t paccess_idx{0};
    std::array<T, kPackSize<DIM>> init_patch_data{};

    /* Define the iterators to the dimension lengths and the coarse data */
    auto coarse_predictor_iter = this->data_pyramid_.rbegin();
    auto coarse_dim_lengths_iter = this->dim_length_pyramid_.rbegin();

    /* Set the start data */
    std::vector<T> data = *coarse_predictor_iter;

    /* Perform the iterative compression steps up until the root level */
    for (int32_t lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx, ++coarse_dim_lengths_iter, ++coarse_predictor_iter)
    {
        /* Define a reference for the ease of notation */
        const std::vector<T>& coarse_data = data;
        const std::array<int32_t, DIM>& coarse_dim_lengths = *coarse_dim_lengths_iter;
        
        const std::vector<T>& next_lvl_init_data = *std::next(coarse_predictor_iter);
        const std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::next(coarse_dim_lengths_iter);

        cmc_debug_msg("A coarsening iteration is initialized (step: ", lvl_idx, ").");
        cmc_debug_msg("Level's data dimensions are: ", coarse_dim_lengths[kTimeID], ", ", coarse_dim_lengths[kLevID], ", ", coarse_dim_lengths[kLatID], ", ", coarse_dim_lengths[kLonID]);
        
        /* Allocate the next level's data */
        std::vector<T> next_level_data(next_lvl_init_data.size());

        /* Set up a vector which will hold the encoding for this level */
        this->levelwise_patch_encodings_.emplace_back();
        this->levelwise_patch_encodings_.back().reserve(coarse_data.size());

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
                        const std::array<T, k4DNumControlValues> control_points = GetFaceControlValues<T, DIM>(coarse_data, time, lev, lat, lon, coarse_dim_lengths);

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
                                        if (kDIMReductionFactor * time + time_idx >= next_lvl_dim_lengths[kTimeID] || kDIMReductionFactor * lev + lev_idx >= next_lvl_dim_lengths[kLevID] || kDIMReductionFactor * lat + lat_idx >= next_lvl_dim_lengths[kLatID] || kDIMReductionFactor * lon + lon_idx >= next_lvl_dim_lengths[kLonID]) [[unlikely]]
                                        {
                                            continue;
                                        } else
                                        {
                                            /* Retrieve the corresponding value from this patch */
                                            init_patch_data[paccess_idx] = GetValue<T>(next_lvl_init_data, kDIMReductionFactor * time + time_idx, kDIMReductionFactor * lev + lev_idx, kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx,
                                                                                       next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID], next_lvl_dim_lengths[kTimeID]);
                                        }
                                        ++paccess_idx;
                                    }
                                }
                            }
                        }

                        /* Perform the prediction onto the next level */
                        PatchEncoding<T, DIM> elem_encoding = this->PerformPrediction(control_points, init_patch_data, this->abs_error_, lon, lat, lev, time, next_level_data, next_lvl_dim_lengths);

                        /* Store the encoding */
                        this->levelwise_patch_encodings_.back().push_back(std::move(elem_encoding));
                    }
                }
            }
        }

        /* Store the data for the next level */
        data = std::move(next_level_data);

        cmc_debug_msg("The refinement prediction iteration is finished.");
    }

    /* Encode the gathered data */
    this->EncodeData();

    cmc_debug_msg("The lossy compression of is finished.");
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

    cmc_debug_msg("Lossy patch-based compression starts...");

    cmc_assert(this->dim_length_pyramid_.size() == 1ULL);

    /* Build up the coarse prediction pyramid */
    this->CollectCoarseDataPyramid(CompressionVariable<T, DIM>::tag3D);

    cmc_debug_msg("Number of compression iterations to be performed: ", this->num_compression_lvls_);
    
    int32_t paccess_idx{0};
    std::array<T, kPackSize<DIM>> init_patch_data{};

    /* Define the iterators to the dimension lengths and the coarse data */
    auto coarse_predictor_iter = this->data_pyramid_.rbegin();
    auto coarse_dim_lengths_iter = this->dim_length_pyramid_.rbegin();

    /* Set the start data */
    std::vector<T> data = *coarse_predictor_iter;

    /* Perform the iterative compression steps up until the root level */
    for (int32_t lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx, ++coarse_dim_lengths_iter, ++coarse_predictor_iter)
    {
        /* Define a reference for the ease of notation */
        const std::vector<T>& coarse_data = data;
        const std::array<int32_t, DIM>& coarse_dim_lengths = *coarse_dim_lengths_iter;
        
        const std::vector<T>& next_lvl_init_data = *std::next(coarse_predictor_iter);
        const std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::next(coarse_dim_lengths_iter);

        cmc_debug_msg("A coarsening iteration is initialized (step: ", lvl_idx, ").");
        cmc_debug_msg("Level's data dimensions are: ", coarse_dim_lengths[kLevID], ", ", coarse_dim_lengths[kLatID], ", ", coarse_dim_lengths[kLonID]);
        
        /* Allocate the next level's data */
        std::vector<T> next_level_data(next_lvl_init_data.size());

        /* Set up a vector which will hold the encoding for this level */
        this->levelwise_patch_encodings_.emplace_back();
        this->levelwise_patch_encodings_.back().reserve(coarse_data.size());

        /* Iterate over patches */
        for (int32_t lev = 0; lev < coarse_dim_lengths[kLevID]; ++lev)
        {
            for (int32_t lat = 0; lat < coarse_dim_lengths[kLatID]; ++lat)
            {
                for (int32_t lon = 0; lon < coarse_dim_lengths[kLonID]; ++lon)
                {
                    /** Now this element will be refined via prediction **/

                    /* Gather the face control values */
                    const std::array<T, k3DNumControlValues> control_points = GetFaceControlValues<T, DIM>(coarse_data, lev, lat, lon, coarse_dim_lengths);

                    /* Reset the access idx for each new patch */
                    paccess_idx = 0;

                    /* Gather the values for this initial patch */
                    for (int32_t lev_idx = 0; lev_idx < kDIMReductionFactor; ++lev_idx)
                    {
                        for (int32_t lat_idx = 0; lat_idx < kDIMReductionFactor; ++lat_idx)
                        {
                            for (int32_t lon_idx = 0; lon_idx < kDIMReductionFactor; ++lon_idx)
                            {
                                if (kDIMReductionFactor * lev + lev_idx >= next_lvl_dim_lengths[kLevID] || kDIMReductionFactor * lat + lat_idx >= next_lvl_dim_lengths[kLatID] || kDIMReductionFactor * lon + lon_idx >= next_lvl_dim_lengths[kLonID]) [[unlikely]]
                                {
                                    continue;
                                } else
                                {
                                    /* Retrieve the corresponding value from this patch */
                                    init_patch_data[paccess_idx] = GetValue<T>(next_lvl_init_data, kDIMReductionFactor * lev + lev_idx, kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx,
                                                                               next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID], next_lvl_dim_lengths[kLevID]);
                                }
                                ++paccess_idx;
                            }
                        }
                    }

                    /* Perform the prediction onto the next level */
                    PatchEncoding<T, DIM> elem_encoding = this->PerformPrediction(control_points, init_patch_data, this->abs_error_, lon, lat, lev, next_level_data, next_lvl_dim_lengths);

                    /* Store the encoding */
                    this->levelwise_patch_encodings_.back().push_back(std::move(elem_encoding));
                }
            }
        }

        /* Store the data for the next level */
        data = std::move(next_level_data);

        cmc_debug_msg("The refinement prediction iteration is finished.");
    }

    /* Encode the gathered data */
    this->EncodeData();

    cmc_debug_msg("The lossy compression of is finished.");
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::Compress(CompressionVariable<T, DIM>::kTag2D)
{
    static_assert(DIM == 2);
    constexpr int32_t kLatID = 0;
    constexpr int32_t kLonID = 1;

    cmc_debug_msg("Lossy patch-based compression starts...");

    cmc_assert(this->dim_length_pyramid_.size() == 1ULL);

    /* Build up the coarse prediction pyramid */
    this->CollectCoarseDataPyramid(CompressionVariable<T, DIM>::tag2D);

    cmc_debug_msg("Number of compression iterations to be performed: ", this->num_compression_lvls_);
    
    int32_t paccess_idx{0};
    std::array<T, kPackSize<DIM>> init_patch_data{};

    /* Define the iterators to the dimension lengths and the coarse data */
    auto coarse_predictor_iter = this->data_pyramid_.rbegin();
    auto coarse_dim_lengths_iter = this->dim_length_pyramid_.rbegin();

    /* Set the start data */
    std::vector<T> data = *coarse_predictor_iter;

    /* Perform the iterative compression steps up until the root level */
    for (int32_t lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx, ++coarse_dim_lengths_iter, ++coarse_predictor_iter)
    {
        /* Define a reference for the ease of notation */
        const std::vector<T>& coarse_data = data;
        const std::array<int32_t, DIM>& coarse_dim_lengths = *coarse_dim_lengths_iter;
        
        const std::vector<T>& next_lvl_init_data = *std::next(coarse_predictor_iter);
        const std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::next(coarse_dim_lengths_iter);

        cmc_global_msg("A coarsening iteration is initialized (step: ", lvl_idx, ").");
        cmc_global_msg("Level's data dimensions are: ", coarse_dim_lengths[kLatID], ", ", coarse_dim_lengths[kLonID]);
        
        /* Allocate the next level's data */
        std::vector<T> next_level_data(next_lvl_init_data.size());

        /* Set up a vector which will hold the encoding for this level */
        this->levelwise_patch_encodings_.emplace_back();
        this->levelwise_patch_encodings_.back().reserve(coarse_data.size());

        /* Iterate over patches */
        for (int32_t lat = 0; lat < coarse_dim_lengths[kLatID]; ++lat)
        {
            for (int32_t lon = 0; lon < coarse_dim_lengths[kLonID]; ++lon)
            {
                /** Now this element will be refined via prediction **/

                /* Gather the face control values */
                const std::array<T, k2DNumControlValues> control_points = GetFaceControlValues<T, DIM>(coarse_data, lat, lon, coarse_dim_lengths);

                /* Reset the access idx for each new patch */
                paccess_idx = 0;

                /* Gather the values for this initial patch */
                for (int32_t lat_idx = 0; lat_idx < kDIMReductionFactor; ++lat_idx)
                {
                    for (int32_t lon_idx = 0; lon_idx < kDIMReductionFactor; ++lon_idx)
                    {
                        if (kDIMReductionFactor * lat + lat_idx >= next_lvl_dim_lengths[kLatID] || kDIMReductionFactor * lon + lon_idx >= next_lvl_dim_lengths[kLonID]) [[unlikely]]
                        {
                            continue;
                        } else
                        {
                            /* Retrieve the corresponding value from this patch */
                            init_patch_data[paccess_idx] = GetValue<T>(next_lvl_init_data, kDIMReductionFactor * lat + lat_idx, kDIMReductionFactor * lon + lon_idx,
                                                                       next_lvl_dim_lengths[kLonID], next_lvl_dim_lengths[kLatID]);
                        }
                        ++paccess_idx;
                    }
                }

                /* Perform the prediction onto the next level */
                PatchEncoding<T, DIM> elem_encoding = this->PerformPrediction(control_points, init_patch_data,  this->abs_error_, lon, lat, next_level_data, next_lvl_dim_lengths);

                /* Store the encoding */
                this->levelwise_patch_encodings_.back().push_back(std::move(elem_encoding));
            }
        }

        /* Store the data for the next level */
        data = std::move(next_level_data);

        cmc_debug_msg("The refinement prediction iteration is finished.");
    }

    /* Encode the gathered data */
    this->EncodeData();

    cmc_debug_msg("The lossy compression of is finished.");
}

template <ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionVariable<T, DIM>::Compress(CompressionVariable<T, DIM>::kTag1D)
{
    static_assert(DIM == 1);
    constexpr int32_t kLonID = 0;

    cmc_debug_msg("Lossy patch-based compression starts...");

    cmc_assert(this->dim_length_pyramid_.size() == 1ULL);

    /* Build up the coarse prediction pyramid */
    this->CollectCoarseDataPyramid(CompressionVariable<T, DIM>::tag1D);

    cmc_debug_msg("Number of compression iterations to be performed: ", this->num_compression_lvls_);
    
    int32_t paccess_idx{0};
    std::array<T, kPackSize<DIM>> init_patch_data{};

    /* Define the iterators to the dimension lengths and the coarse data */
    auto coarse_predictor_iter = this->data_pyramid_.rbegin();
    auto coarse_dim_lengths_iter = this->dim_length_pyramid_.rbegin();

    /* Set the start data */
    std::vector<T> data = *coarse_predictor_iter;

    /* Perform the iterative compression steps up until the root level */
    for (int32_t lvl_idx{0}; lvl_idx < this->num_compression_lvls_; ++lvl_idx, ++coarse_dim_lengths_iter, ++coarse_predictor_iter)
    {
        /* Define a reference for the ease of notation */
        const std::vector<T>& coarse_data = data;
        const std::array<int32_t, DIM>& coarse_dim_lengths = *coarse_dim_lengths_iter;
        
        const std::vector<T>& next_lvl_init_data = *std::next(coarse_predictor_iter);
        const std::array<int32_t, DIM>& next_lvl_dim_lengths = *std::next(coarse_dim_lengths_iter);

        cmc_debug_msg("A coarsening iteration is initialized (step: ", lvl_idx, ").");
        cmc_debug_msg("Level's data dimensions are: ", coarse_dim_lengths[kLonID]);
        
        /* Allocate the next level's data */
        std::vector<T> next_level_data(next_lvl_init_data.size());

        /* Set up a vector which will hold the encoding for this level */
        this->levelwise_patch_encodings_.emplace_back();
        this->levelwise_patch_encodings_.back().reserve(coarse_data.size());

        /* Iterate over patches */
        for (int32_t lon = 0; lon < coarse_dim_lengths[kLonID]; ++lon)
        {
            /** Now this element will be refined via prediction **/

            /* Gather the face control values */
            const std::array<T, k1DNumControlValues> control_points = GetFaceControlValues<T, DIM>(coarse_data, lon, coarse_dim_lengths);

            /* Reset the access idx for each new patch */
            paccess_idx = 0;

            /* Gather the values for this initial patch */
            for (int32_t lon_idx = 0; lon_idx < kDIMReductionFactor; ++lon_idx)
            {
                if (kDIMReductionFactor * lon + lon_idx >= next_lvl_dim_lengths[kLonID]) [[unlikely]]
                {
                    continue;
                } else
                {
                    /* Retrieve the corresponding value from this patch */
                    init_patch_data[paccess_idx] = GetValue<T>(next_lvl_init_data, kDIMReductionFactor * lon + lon_idx,
                                                               next_lvl_dim_lengths[kLonID]);
                }
                ++paccess_idx;
            }

            /* Perform the prediction onto the next level */
            PatchEncoding<T, DIM> elem_encoding = this->PerformPrediction(control_points, init_patch_data,  this->abs_error_, lon, next_level_data, next_lvl_dim_lengths);

            /* Store the encoding */
            this->levelwise_patch_encodings_.back().push_back(std::move(elem_encoding));
        }

        /* Store the data for the next level */
        data = std::move(next_level_data);

        cmc_debug_msg("The refinement prediction iteration is finished.");
    }

    /* Encode the gathered data */
    this->EncodeData();

    cmc_debug_msg("The lossy compression of is finished.");
}

}

#endif /* !CMC_PATCH_LOSSY_MULTI_RES_RBF_HXX */
