#ifndef CMC_PATCH_PREFIX_EXTRACTION_PLAIN_SUFFIXES_DECOMPRESSION_HXX
#define CMC_PATCH_PREFIX_EXTRACTION_PLAIN_SUFFIXES_DECOMPRESSION_HXX

#include "cmc.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_prefix_entropy_coder.hxx"
#include "patch/lossless/cmc_patch_byte_decompression_variable.hxx"

#include <utility>
#include <vector>
#include <memory>
#include <cstdint>

namespace cmc::patch::decompression::prefix::plain_suffix
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
    uint32_t GetNextPrefixLength();
    std::vector<uint8_t> GetNextPrefixBitSequence(const size_t num_bits);
    CompressionValue<T> GetNextSuffixedValue(const CompressionValue<T>& value);
    uint32_t ApplyProcessBoundarySymbol();
    CompressionValue<T> GetNextInteriorSuffixedValue(const CompressionValue<T>& value);
    CompressionValue<T> GetNextFinalSuffixedValue(const CompressionValue<T>& value);
    bool IsNextIterationLeafLevelDecompression() const;
    void InitializeInteriorDecompressionIteration();
    void InitializeLeafDecompressionIteration();  

    size_t level_byte_offset_{0};
    bit_vector::BitVectorView alphabet_;
    bit_map::BitMapView encoded_prefix_length_;
    bit_vector::BitVectorView prefix_bits_;
    int count_adaptation_step_{0};
    std::unique_ptr<cmc::entropy_coding::arithmetic_coding::Decoder> entropy_decoder_{nullptr};

};

template <typename T, size_t Dim>
void
DecompressionVariable<T, Dim>::FinalizeDecompressionIteration()
{
    ++count_adaptation_step_;
}

template <typename T, size_t Dim>
CompressionValue<T>
DecompressionVariable<T, Dim>::GetInitDecompressionValue()
{
    cmc_debug_msg("The setup of the initial value is performed.");

    this->InitializeDecompressionIteration();

    /* Get the next suffixed value */
    const CompressionValue<T> init_value = GetNextSuffixedValue(CompressionValue<T>());

    this->FinalizeDecompressionIteration();

    return init_value;
}

template <typename T, size_t Dim>
inline uint32_t
DecompressionVariable<T, Dim>::GetNextPrefixLength()
{
    cmc_assert(entropy_decoder_ != nullptr);
    return entropy_decoder_->DecodeNextSymbol();
}

template <typename T, size_t Dim>
inline std::vector<uint8_t>
DecompressionVariable<T, Dim>::GetNextPrefixBitSequence(const size_t num_prefix_bits)
{
    return prefix_bits_.GetNextBitSequence(num_prefix_bits);
}

template <typename T, size_t Dim>
bool
DecompressionVariable<T, Dim>::IsNextIterationLeafLevelDecompression() const 
{
    return (this->GetNumDecompressionIterations() - 1 == count_adaptation_step_);
}

template <typename T, size_t Dim>
void
DecompressionVariable<T, Dim>::InitializeDecompressionIteration()
{
    /* Check if the following decompression iteration is the leaf level decompression */
    cmc_debug_msg("\n\nIs this iteration a leaf level decompression? ", IsNextIterationLeafLevelDecompression(), "\n\n");

    if (IsNextIterationLeafLevelDecompression())
    {
        InitializeLeafDecompressionIteration();
    } else
    {
        InitializeInteriorDecompressionIteration();
    }
}

template <typename T, size_t Dim>
void
DecompressionVariable<T, Dim>::InitializeInteriorDecompressionIteration()
{
    cmc_debug_msg("A prefix decompression iteration is initialized.");
    
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
    const uint64_t encoded_prefix_length_bytes = GetValueFromByteStream<uint64_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;

    /* Get the bytes for the remaining bits */
    const uint64_t prefix_bytes = GetValueFromByteStream<uint64_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;

    /* Set the view on the alphabet */
    alphabet_ = bit_vector::BitVectorView(data_start_ptr + processed_bytes, alphabet_bytes);
    processed_bytes += alphabet_bytes;

    /* Set the view on the encoded prefix lengths */
    encoded_prefix_length_ = bit_map::BitMapView(data_start_ptr + processed_bytes, bit_map::kCharBit * encoded_prefix_length_bytes);
    processed_bytes += encoded_prefix_length_bytes;

    /* Set the view on the remaining bits */
    prefix_bits_ = bit_vector::BitVectorView(data_start_ptr + processed_bytes, prefix_bytes);
    processed_bytes += prefix_bytes;

    /* Update the byte count */
    level_byte_offset_ = processed_bytes;

    /* Setup the entropy decoder */
    entropy_decoder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::PrefixDecoder<T>>(alphabet_.begin(), encoded_prefix_length_);
    entropy_decoder_->SetupDecoding(); 
}

template <typename T, size_t Dim>
void
DecompressionVariable<T, Dim>::InitializeLeafDecompressionIteration()
{
    cmc_debug_msg("A prefix decompression iteration is initialized.");
    
    /* Check if the following decompression iteration is the leaf level decompression */
    cmc_debug_msg("\n\nIs this iteration a leaf level decompression? ", IsNextIterationLeafLevelDecompression(), "\n\n");
    constexpr size_t offset = sizeof(uint64_t);

    size_t processed_bytes = level_byte_offset_;

    /* Get a pointer to the beginning of the encoded data stream */
    const auto data_start_ptr = this->GetEncodedDataStream().data();

    /* Get the amount of relevant bytes for this decompression level */
    const uint64_t current_level_bytes = GetValueFromByteStream<uint64_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;

    cmc_debug_msg("The current refinement level is described by ", current_level_bytes, " bytes.");

    /* Set the view on the remaining bits */
    const uint64_t prefix_bytes = current_level_bytes - sizeof(uint64_t);
    prefix_bits_ = bit_vector::BitVectorView(data_start_ptr + processed_bytes, prefix_bytes);
    processed_bytes += prefix_bytes;

    /* Update the byte count */
    level_byte_offset_ = processed_bytes;
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
        next_symbol = this->GetNextPrefixLength();

        if (next_symbol != entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte)
        {
            is_process_boundary = false;
        }
    }

    /* Move the encoded signifcant bits view to the next byte as well */
    prefix_bits_.MoveToNextByte();

    /* Return the newly obtained symbol */
    return next_symbol;
}


template <typename T, size_t Dim>
CompressionValue<T>
DecompressionVariable<T, Dim>::GetNextSuffixedValue(const CompressionValue<T>& value)
{
    if (this->IsNextIterationLeafLevelDecompression())
    {
        return GetNextFinalSuffixedValue(value);
    } else
    {
        return GetNextInteriorSuffixedValue(value);
    }
}

template <typename T, size_t Dim>
CompressionValue<T>
DecompressionVariable<T, Dim>::GetNextInteriorSuffixedValue(const CompressionValue<T>& value)
{
    CompressionValue<T> suffixed_value = value;

    /* Get the length of the suffix to be appended */
    uint32_t encoded_suffix_length = this->GetNextPrefixLength();

    /* Check if a process-boundary symbol has been encoded */
    if (encoded_suffix_length == entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte)
    {
        /* Apply the symbol an get the next 'unequal' symbol */
        encoded_suffix_length = this->ApplyProcessBoundarySymbol();
    }

    const int suffix_length = static_cast<int>(encoded_suffix_length);

    if (suffix_length > 0)
    {
        if (suffixed_value.GetTailBit() == 0)
        {
            cmc_debug_msg("\n\n\nNextInteriorSuffix: suffix_length: ", suffix_length, ", tail_bit: ", static_cast<int>(suffixed_value.GetTailBit()),"\n\n");

            return suffixed_value;
        }
        /* Get the actual suffix */
        const std::vector<uint8_t> suffix = this->GetNextPrefixBitSequence(suffix_length);

        /* Apply the suffix to the current value */
        suffixed_value.ApplySuffix(suffix, suffix_length);
    }

    return suffixed_value;
}


template <typename T, size_t Dim>
CompressionValue<T>
DecompressionVariable<T, Dim>::GetNextFinalSuffixedValue(const CompressionValue<T>& value)
{
    CompressionValue<T> suffixed_value = value;
    cmc_assert(suffixed_value.GetFrontBit() == 0);

    const int curernt_length = suffixed_value.GetCountOfSignificantBits();
    const int num_missing_bits = suffixed_value.GetSizeBits() - curernt_length;

    if (num_missing_bits > 0)
    {
        /* Get the actual suffix */
        const std::vector<uint8_t> suffix = this->GetNextPrefixBitSequence(num_missing_bits);

        /* Apply the suffix to the current value */
        suffixed_value.ApplySuffix(suffix, num_missing_bits);
    }

    cmc_assert(suffixed_value.GetTailBit() == 0);
    return suffixed_value;
}

template <typename T, size_t Dim>
CompressionValue<T>
DecompressionVariable<T, Dim>::RefineNextValue(const CompressionValue<T>& coarse_value)
{
    const CompressionValue<T> fine_val = this->GetNextSuffixedValue(coarse_value);

    return fine_val;
}


}


#endif /* !CMC_PATCH_PREFIX_EXTRACTION_PLAIN_SUFFIXES_DECOMPRESSION_HXX */
