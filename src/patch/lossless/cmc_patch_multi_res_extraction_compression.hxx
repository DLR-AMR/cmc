#ifndef CMC_PATCH_MULTI_RES_EXTRACTION_COMPRESSION_HXX
#define CMC_PATCH_MULTI_RES_EXTRACTION_COMPRESSION_HXX

#include "cmc.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_prefix_entropy_coder.hxx"
#include "utilities/cmc_serialization.hxx"
#include "patch/lossless/cmc_patch_byte_compression_variable.hxx"
#include "utilities/cmc_multi_res_entropy_coder.hxx"
#include "utilities/cmc_multi_res_extraction_util.hxx"
#include "utilities/cmc_lossless_multi_res_extraction_residual_computation.hxx"
#include "utilities/cmc_interpolation_fn.hxx"

#include <utility>
#include <vector>
#include <cstdint>
#include <memory>

namespace cmc::patch::lossless::multi_res
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template<class T, size_t Dim>
class PatchCompressionVariable : public AbstractPatchByteCompressionVariable<T, Dim>
{
public:
    PatchCompressionVariable() = delete;

    PatchCompressionVariable(input::Var& input_variable)
    : AbstractPatchByteCompressionVariable<T, Dim>(input_variable) {}

    CompressionSchema GetCompressionSchema() const override
    {
        return CompressionSchema::PatchMultiResExtraction;
    }

protected:
    ExtractionData<T> PerformExtraction(const std::vector<CompressionValue<T>>& patch) override;
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const override;
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const override;
    
    void InitializeExtractionIteration() override;
    void CompleteExtraction(const std::vector<size_t>& fine_dim_lengths,  const size_t kDimReductionFactor, [[maybe_unused]] const std::vector<CompressionValue<T>>& fine_vals, [[maybe_unused]] const std::vector<CompressionValue<T>>& coarse_vals) override;
private:
void
    CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const;

    std::unique_ptr<entropy_coding::IByteCompressionEntropyCoder> entropy_coder_{std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResEncoder<T>>()}; //!< The entropy coder to use in order to encode information
    std::vector<uint8_t> resdiual_order_indications_;
};

template <typename T>
inline void
SetValue(std::vector<T>& data, const T& value, const int lev, const int lat, const int lon, const int kLonLength, const int kLatLength,  const int kLevLength)
{
    cmc_assert(lev * (kLatLength * kLonLength) + lat * kLonLength + lon < data.size());
    data[lev * (kLatLength * kLonLength) + lat * kLonLength + lon] = value;
}

template<class T, size_t Dim>
void
PatchCompressionVariable<T, Dim>::CompleteExtraction(const std::vector<size_t>& fine_dim_lengths, const size_t kDimReductionFactor, [[maybe_unused]] const std::vector<CompressionValue<T>>& fine_vals, [[maybe_unused]] const std::vector<CompressionValue<T>>& coarse_vals)
{
    cmc_assert(fine_dim_lengths.size() == 3);
    
    std::vector<uint8_t> res_indications(fine_dim_lengths[0] * fine_dim_lengths[1] * fine_dim_lengths[2]);

    size_t res_lin_idx{0};

    /* Iterate over patches */
    for (size_t lev = 0; lev < fine_dim_lengths[0]; lev += kDimReductionFactor)
    {
        for (size_t lat = 0; lat < fine_dim_lengths[1]; lat += kDimReductionFactor)
        {
            for (size_t lon = 0; lon < fine_dim_lengths[2]; lon += kDimReductionFactor)
            {
                /* Gather the values for this patch */
                for (size_t lev_idx = 0; lev_idx < kDimReductionFactor; ++lev_idx)
                {
                    for (size_t lat_idx = 0; lat_idx < kDimReductionFactor; ++lat_idx)
                    {
                        for (size_t lon_idx = 0; lon_idx < kDimReductionFactor; ++lon_idx)
                        {
                            if (lev + lev_idx >= fine_dim_lengths[0] || lat + lat_idx >= fine_dim_lengths[1] || lon + lon_idx >= fine_dim_lengths[2])
                            {
                                continue;
                            } else
                            {
                                SetValue<uint8_t>(res_indications, resdiual_order_indications_[res_lin_idx], lev + lev_idx, lat + lat_idx, lon + lon_idx, fine_dim_lengths[2], fine_dim_lengths[1], fine_dim_lengths[0]);
                                ++res_lin_idx;
                            }
                        }
                    }
                }
            }
        }
    }

    std::swap(res_indications, resdiual_order_indications_);
}

template<class T, size_t Dim>
void
PatchCompressionVariable<T, Dim>::InitializeExtractionIteration()
{
    cmc_debug_msg("initialize patch multi res extraction iteration");
    resdiual_order_indications_ = std::vector<uint8_t>();
}

template<typename T>
T
GetCoarseApproximationMaximizingResidualsLZC(const std::vector<CompressionValue<T>>& values)
{
    cmc_assert(values.size() >= 1);

    T current_best_predictor;
    int current_max_lzc{-1};

    const VectorView<CompressionValue<T>> value_view(values);

    /* Try each value from the view as a predictor */
    for (auto pred_iter = value_view.begin(); pred_iter != value_view.end(); ++pred_iter)
    {
        /* Convert the compression value back to its origianl type */
        const T predictor = pred_iter->template ReinterpretDataAs<T>();

        /* Compute the LZC in the residuals for this approximation */
        const int pred_lzc = cmc::lossless::multi_res::GetCumulativeResidualsLZC<T>(predictor, value_view);

        /* Potentially, update the current predictor */
        if (pred_lzc > current_max_lzc)
        {
            current_best_predictor = predictor;
            current_max_lzc = pred_lzc;
        }
    }

    /* Convert the view to actual values of the underlying data type */
    const std::vector<T> converted_vals = ConvertCompressionValues<T>(value_view);

    /* Try the mid-range as an predictor */
    const T mid_range = InterpolateToMidRange<T>(converted_vals);

    /* Compute the LZC in the residuals for this approximation */
    const int mid_range_lzc = cmc::lossless::multi_res::GetCumulativeResidualsLZC<T>(mid_range, value_view);

    /* Potentially, update the current predictor */
    if (mid_range_lzc > current_max_lzc)
    {
        current_best_predictor = mid_range;
        current_max_lzc = mid_range_lzc;
    }

    /* Try the arithmetic mean as an predictor */
    const T mean = InterpolateToArithmeticMean<T>(converted_vals);

    /* Compute the LZC in the residuals for this approximation */
    const int mean_lzc = cmc::lossless::multi_res::GetCumulativeResidualsLZC<T>(mean, value_view);

    /* Potentially, update the current predictor */
    if (mean_lzc > current_max_lzc)
    {
        current_best_predictor = mean;
        current_max_lzc = mean_lzc;
    }

    return current_best_predictor;
}

template<class T, size_t Dim>
ExtractionData<T>
PatchCompressionVariable<T, Dim>::PerformExtraction(const std::vector<CompressionValue<T>>& patch)
{
    /* Get the coarse approximation for these values which maximizes the cumulative residual LZC */
    const T coarse_approximation = GetCoarseApproximationMaximizingResidualsLZC<T>(patch);

    std::vector<CompressionValue<T>> fine_values;
    fine_values.reserve(patch.size());

    /* Compute the residuals for the given predictor */
    for (auto val_iter = patch.begin(); val_iter != patch.end(); ++val_iter)
    {
        /* Compute the residual between approximation and actual value */
        auto [is_approximation_greater, residual] = cmc::lossless::multi_res::ComputeResidual<T>(coarse_approximation, *val_iter);

        /* Store whether the approximation is greater than the real value */
        resdiual_order_indications_.push_back(is_approximation_greater ? uint8_t{1} : uint8_t{0});

        /* Store the residual */
        fine_values.push_back(residual);
    }

    return ExtractionData<T>(CompressionValue<T>(coarse_approximation), std::move(fine_values));
}

/**
 * @brief  The encodig of the level-wise compression data is handled within the function, we encode the data
 * on the leaf level differently, than all other levels.
 * 
 * @tparam T The original data type of the underlying data (e.g. float)
 * @param level_byte_values The remaining "fine compression" values after an extraction iteration
 * @return std::vector<uint8_t> The encoded data stream
 */
template<class T, size_t Dim>
inline std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
PatchCompressionVariable<T, Dim>::EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    cmc_debug_msg("The encoding of the CompressionValues after the patch-based multi-resolution extraction iteration starts...");
    
    /* The encoded data will be stored in a BitVector */
    cmc::bit_vector::BitVector encoding;
    encoding.Reserve(3 * level_byte_values.size());

    /* Reset the entropy coder and initialize the alphabet */
    entropy_coder_->Reset(std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResCompressionAlphabet<T>>());
    entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Collect the symbols and their frequencies for the entropy coder */
    CollectSymbolFrequenciesForEntropyCoding(level_byte_values);

    /* Setup the interior structure for encoding */
    entropy_coder_->SetupEncoding(MPI_COMM_SELF);

    /* Iterate over all values and encode them */
    size_t res_iter{0};
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter, ++res_iter)
    {
        /* Get the current value */
        CompressionValue<T> val = *val_iter;

        /* Get the encoded LZC */
        uint32_t signum = cmc::lossless::multi_res::util::GetSignumForEncoding(resdiual_order_indications_[res_iter] != 0 ? true : false);
        const uint32_t first_one_bit = val.GetNumberLeadingZeros();
        if (val.GetNumberLeadingZeros() == sizeof(T) * bit_map::kCharBit)
        {
            signum = 0;
        }
        
        /* Encode the LZC and the residual flag together */
        entropy_coder_->EncodeSymbol(signum + first_one_bit);

        /* The first one is implicitly given by the leading zero count; therefore we set the "front bit" in order to discard the leading zeros and the following one */
        val.SetFrontBit(first_one_bit + 1);

        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            encoding.AppendBits(val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
        }
    }

    /* We set an indicaton symbol that the process local end of the values have been reached */
    entropy_coder_->EncodeSymbol(entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte);

    /* Indicate that the encoding has been finished and flush all pending encodings */
    entropy_coder_->FinishEncoding();

    /* Set up the BitVector holding the encoded data for further use */
    encoding.TrimToContent();

    /* Get the local remaining significant bits */
    const uint64_t local_remaining_significant_bits_num_bytes = static_cast<uint64_t>(encoding.size());

    /* Get the local encoded LZC, respectively the encoded first "one-bit" positions */
    cmc::bit_map::BitMap local_encoded_lzc_stream = entropy_coder_->GetEncodedBitStream();
    const uint64_t local_encoded_lzc_stream_num_bytes = static_cast<uint64_t>(local_encoded_lzc_stream.size_bytes());

    /* We need exchange the encoded lengths */
    const std::vector<uint64_t> local_bytes{local_encoded_lzc_stream_num_bytes, local_remaining_significant_bits_num_bytes};

    /* Declare the buffers for the encoded data */
    std::vector<uint8_t> encoded_entropy_codes;
    std::vector<uint8_t> encoded_data;

    /* Get the encoded alphabet */
    cmc::bit_vector::BitVector encoded_alphabet = entropy_coder_->EncodeAlphabet();
    const uint64_t encoded_alphabet_num_bytes = static_cast<uint64_t>(encoded_alphabet.size());

    /* Calculate the overall amount of bytes on the root rank */
    const uint64_t num_locally_encoded_entropy_codes_bytes = 4 * sizeof(uint64_t) + encoded_alphabet_num_bytes + local_encoded_lzc_stream_num_bytes;

    /* Allocate memory for the encoded data */
    encoded_entropy_codes.reserve(num_locally_encoded_entropy_codes_bytes);
    
    /** We store global information about the encoded level **/
    /* Push back the overall byte count for the level */
    const uint64_t num_global_bytes_encoded_level = 4 * sizeof(uint64_t) + encoded_alphabet_num_bytes + local_bytes[0] + local_bytes[1];
    PushBackValueToByteStream<uint64_t>(encoded_entropy_codes, num_global_bytes_encoded_level);

    /* Push back the byte count for the encoded alphabet */
    PushBackValueToByteStream<uint64_t>(encoded_entropy_codes, encoded_alphabet_num_bytes);

    /* Push back the byte count for the encoded first "one-bit" positions */
    PushBackValueToByteStream<uint64_t>(encoded_entropy_codes, local_bytes[0]);

    /* Push back the byte count for the remaining significant bits */
    PushBackValueToByteStream<uint64_t>(encoded_entropy_codes, local_bytes[1]);

    /* Afterwards, we store the encoded alphabet */
    std::copy_n(encoded_alphabet.begin(), encoded_alphabet_num_bytes, std::back_insert_iterator(encoded_entropy_codes));

    /* Finally, copy the entropy codes */
    std::copy_n(local_encoded_lzc_stream.begin_bytes(), local_encoded_lzc_stream_num_bytes, std::back_insert_iterator(encoded_entropy_codes));

    /* Get the encoded remaining significant bits */
    encoding.MoveDataInto(encoded_data);

    cmc_debug_msg("The entropy encoder of the patch-based multi-resolution extraction compression completed the encoding of the CompressionValues of this iteration.");
    
    return std::make_pair(encoded_entropy_codes, encoded_data);
}

template<class T, size_t Dim>
void
PatchCompressionVariable<T, Dim>::CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    size_t res_iter{0};
    /* Collect the symbol frequencies */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter, ++res_iter)
    {
        /* Get the current residual */
        CompressionValue<T> val = *val_iter;

        /* Get the encoded LZC */
        uint32_t signum = cmc::lossless::multi_res::util::GetSignumForEncoding(resdiual_order_indications_[res_iter] != 0 ? true : false);
        if (val.GetNumberLeadingZeros() == sizeof(T) * bit_map::kCharBit)
        {
            signum = 0;
        }
        const uint32_t first_one_bit = val.GetNumberLeadingZeros();

        /* Update this symbol for encoding */
        entropy_coder_->UpdateSymbolFrequency(signum + first_one_bit);
    }
}


template<class T, size_t Dim>
std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
PatchCompressionVariable<T, Dim>::EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const
{
    cmc_debug_msg("The encoding of the root level values of the patch-based multi-resolution compression starts.");

    std::vector<uint8_t> encoded_stream;
    encoded_stream.reserve(root_level_values.size() * sizeof(T));

    for (auto val_iter = root_level_values.begin(); val_iter != root_level_values.end(); ++val_iter)
    {
        const T val = val_iter->template ReinterpretDataAs<T>();
        cmc_debug_msg("Root val to be encoded: ", val);
        PushBackValueToByteStream(encoded_stream, val);
    }

    cmc_debug_msg("The entropy encoder of the patch-based multi-resolution extraction compression stored the root-level CompressionValues within ", encoded_stream.size(), " bytes.");
    return std::make_pair(std::vector<uint8_t>(), std::move(encoded_stream));
}

}



#endif /* !CMC_PATCH_MULTI_RES_EXTRACTION_COMPRESSION_HXX */
