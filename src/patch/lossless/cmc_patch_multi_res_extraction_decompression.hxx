#ifndef CMC_PATCH_MULTI_RES_EXTRACTION_DECOMPRESSION_HXX
#define CMC_PATCH_MULTI_RES_EXTRACTION_DECOMPRESSION_HXX

#include "cmc.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_serialization.hxx"
#include "patch/lossless/cmc_patch_byte_decompression_variable.hxx"
#include "utilities/cmc_multi_res_entropy_coder.hxx"
#include "utilities/cmc_multi_res_extraction_util.hxx"

#include <utility>
#include <vector>
#include <memory>
#include <cstdint>

namespace cmc::patch::decompression::multi_res
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;


template <typename T, size_t Dim>
class DecompressionVariable : public AbstractPatchByteDecompressionVariable<T, Dim>
{
public:
    DecompressionVariable() = delete;

    DecompressionVariable(std::vector<uint8_t>&& encoded_data_byte_stream, const std::string& name, std::vector<std::vector<size_t>>&& dim_lengths_pyramid,
        const GeoDomain& init_domain, const DataLayout init_data_layout, const int num_decompression_lvls)
    : AbstractPatchByteDecompressionVariable<T, Dim>(std::move(encoded_data_byte_stream), name, std::move(dim_lengths_pyramid), init_domain, init_data_layout, num_decompression_lvls)
    {}

protected:
    CompressionValue<T> GetInitDecompressionValue() override;
    CompressionValue<T> RefineNextValue(const CompressionValue<T>& coarse_value) override;
    
    void InitializeDecompressionIteration() override;
    void FinalizeDecompressionIteration() override;

private:
    uint32_t GetNextEncodedResidualLength();
    std::vector<uint8_t> GetNextResidualBitSequence(const size_t num_residual_bits);
    CompressionValue<T> GetNextResidualAppliedValue(const CompressionValue<T>& coarse_value);
    uint32_t ApplyProcessBoundarySymbol();

    size_t level_byte_offset_{0};
    bit_vector::BitVectorView alphabet_;
    bit_map::BitMapView encoded_lzcs_;
    bit_vector::BitVectorView residual_bits_;
    int count_adaptation_step_{0};
    std::unique_ptr<cmc::entropy_coding::arithmetic_coding::Decoder> entropy_decoder_{nullptr};
};

template <typename T, size_t Dim>
void
DecompressionVariable<T, Dim>::FinalizeDecompressionIteration()
{
    cmc_debug_msg("The patch-based multi-resolution decompression iteration (", count_adaptation_step_, ") has been finalized.");
    ++count_adaptation_step_;
}

template <typename T, size_t Dim>
CompressionValue<T>
DecompressionVariable<T, Dim>::GetInitDecompressionValue()
{
    cmc_debug_msg("The setup of the initial value is performed.");

    /* Get the next suffixed value */
    const T val = GetValueFromByteStream<T>(this->GetEncodedDataStream().data());
    cmc_debug_msg("Root level value: ", val);
    
    level_byte_offset_ += sizeof(T);

    return CompressionValue<T>(val);
}

template <typename T, size_t Dim>
void
DecompressionVariable<T, Dim>::InitializeDecompressionIteration()
{
    cmc_debug_msg("A patch-based multi-resolution decompression iteration is initialized.");

    constexpr size_t offset = sizeof(uint64_t);

    size_t processed_bytes = level_byte_offset_;

    /* Get a pointer to the beginning of the encoded data stream */
    const auto data_start_ptr = this->GetEncodedDataStream().data();

    /* Get the amount of relevant bytes for this decompression level */
    const uint64_t current_level_bytes = GetValueFromByteStream<uint64_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;

    cmc_debug_msg("The current refinement level is described by ", current_level_bytes, " bytes.");
    
    /* Get the bytes for the encoded alphabet */
    const uint64_t alphabet_bytes = GetValueFromByteStream<uint64_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;

    /* Get the bytes for the encoded prefix lengths */
    const uint64_t encoded_lzc_bytes = GetValueFromByteStream<uint64_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;

    /* Get the bytes for the remaining bits */
    const uint64_t residual_bytes = GetValueFromByteStream<uint64_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;

    /* Set the view on the alphabet */
    alphabet_ = bit_vector::BitVectorView(data_start_ptr + processed_bytes, alphabet_bytes);
    processed_bytes += alphabet_bytes;

    /* Set the view on the encoded prefix lengths */
    encoded_lzcs_ = bit_map::BitMapView(data_start_ptr + processed_bytes, bit_map::kCharBit * encoded_lzc_bytes);
    processed_bytes += encoded_lzc_bytes;

    /* Set the view on the remaining bits */
    residual_bits_ = bit_vector::BitVectorView(data_start_ptr + processed_bytes, residual_bytes);
    processed_bytes += residual_bytes;

    /* Update the byte count */
    level_byte_offset_ = processed_bytes;

    /* Setup the entropy decoder */
    entropy_decoder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResDecoder<T>>(alphabet_.begin(), encoded_lzcs_);
    entropy_decoder_->SetupDecoding(); 
}

template <typename T, size_t Dim>
inline uint32_t
DecompressionVariable<T, Dim>::ApplyProcessBoundarySymbol()
{
    cmc_debug_msg("\n\nIs this process boundary symnbol happening\n\n");
    bool is_process_boundary = true;
    uint32_t next_symbol;

    /* Iterate until we have cleared (all subsequent) process boundaries */
    while (is_process_boundary)
    {
        /* Reset the decoder */
        entropy_decoder_->ResetAfterProcessBoundary();
        next_symbol = this->GetNextEncodedResidualLength();

        if (next_symbol != entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte)
        {
            is_process_boundary = false;
        }
    }

    /* Move the encoded signifcant bits view to the next byte as well */
    residual_bits_.MoveToNextByte();

    /* Return the newly obtained symbol */
    return next_symbol;
}


template <typename T, size_t Dim>
inline uint32_t
DecompressionVariable<T, Dim>::GetNextEncodedResidualLength()
{
    cmc_assert(entropy_decoder_ != nullptr);
    return entropy_decoder_->DecodeNextSymbol();
}

template <typename T, size_t Dim>
inline std::vector<uint8_t>
DecompressionVariable<T, Dim>::GetNextResidualBitSequence(const size_t num_residual_bits)
{
    return residual_bits_.GetNextBitSequence(num_residual_bits);
}


template <typename T, size_t Dim>
CompressionValue<T>
DecompressionVariable<T, Dim>::GetNextResidualAppliedValue(const CompressionValue<T>& coarse_value)
{
    CompressionValue<T> value = coarse_value;

    /* Get the LZC of the next residual */
    uint32_t encoded_lzc = this->GetNextEncodedResidualLength();

    /* Check if a process-boundary symbol has been encoded */
    if (encoded_lzc == entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte)
    {
        /* Apply the symbol an get the next 'unequal' symbol */
        encoded_lzc = this->ApplyProcessBoundarySymbol();
    }

    /* Get the actual LZC of the residual */
    auto [residual_operation, lzc] = cmc::lossless::multi_res::util::DecodeLZC(encoded_lzc);

    /* Get the maximum length of for the given type */
    const uint32_t max_length_type = sizeof(T) * bit_map::kCharBit;

    /* Check if there is a residual to add/subtract */
    if (lzc < max_length_type)
    {
        /* Determine the length of the encoded residual (excluding the implicit one-bit which has not been stored explicitly) */
        const uint32_t residual_length = max_length_type - lzc - 1;

        /* Get the residual bit sequence */
        const std::vector<uint8_t> residual_bits = this->GetNextResidualBitSequence(residual_length);

        /** Construct a CompressionValue holding the residual **/
        std::array<uint8_t, sizeof(T)> serialized_residual;
        serialized_residual.fill(uint8_t{0});

        CompressionValue<T> residual(serialized_residual);

        /* Set the tail such that the LZC is represented */
        residual.SetTailBit(static_cast<uint8_t>(residual_length + 1));

        /* We add the the implicit one bit */
        residual.ApplySuffix(std::vector<uint8_t>{0x80}, 1);

        cmc_assert((not residual_bits.empty()) || ((residual_length == 0) && residual_bits.empty()));

        if (residual.GetTailBit() > 0)
        {
            /* And finally, we combine it with the actual remaining residual bits */
            residual.ApplySuffix(residual_bits, residual_length);
        }
        if (residual.GetTailBit() != 0)
        {
            cmc_debug_msg("\n\n\nHier ist residual nocht komplett\n\n\n");
        }
        /* Add or subtract the residual */
        if (residual_operation == cmc::lossless::multi_res::util::IntegerAddition)
        {
            value.PerformIntegerAddition(residual);
        } else if (residual_operation == cmc::lossless::multi_res::util::IntegerSubtraction)
        {
            value.PerformIntegerSubtraction(residual);
        }
    }

    return value;
}

template <typename T, size_t Dim>
CompressionValue<T>
DecompressionVariable<T, Dim>::RefineNextValue(const CompressionValue<T>& coarse_value)
{
    const CompressionValue<T> fine_val = this->GetNextResidualAppliedValue(coarse_value);

    return fine_val;
}


}

#endif /* !CMC_PATCH_MULTI_RES_EXTRACTION_DECOMPRESSION_HXX */
