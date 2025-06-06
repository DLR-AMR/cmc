#ifndef CMC_EMBEDDED_PREFIX_EXTRACTION_DECOMPRESSION_PLAIN_SUFFIXES_HXX
#define CMC_EMBEDDED_PREFIX_EXTRACTION_DECOMPRESSION_PLAIN_SUFFIXES_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_serialization.hxx"
#include "lossless/cmc_embedded_byte_decompression_variable.hxx"
#include "utilities/cmc_prefix_entropy_coder.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "mesh_compression/cmc_embedded_mesh_decoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>
#include <memory>

namespace cmc::lossless::embedded::prefix::plain_suffix
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template<typename T>
class PrefixEmbeddedDecompressionAdaptData : public cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>
{
public:
    PrefixEmbeddedDecompressionAdaptData() = delete;
    PrefixEmbeddedDecompressionAdaptData(cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>* variable)
    : cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>(variable) {};

    bool IsDecompressionProgressing() const override;
    void InitializeDecompressionIteration() override;
    void FinalizeDecompressionIteration() override;
    void CompleteDecompressionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::vector<CompressionValue<T>> DecodeRootLevel(const t8_locidx_t num_local_root_values) override;
protected:
    cmc::decompression::embedded::RefinementData<T> PerformRefinement(const int which_tree, const int lelement_id, const CompressionValue<T> value, const int num_refined_elements) override;
    cmc::decompression::embedded::UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) override;

private:
    uint32_t GetNextPrefixLength();
    std::vector<uint8_t> GetNextPrefixBitSequence(const size_t num_bits);
    CompressionValue<T> GetNextSuffixedValue(const CompressionValue<T>& value);
    uint32_t ApplyProcessBoundarySymbol();
    CompressionValue<T> GetNextInteriorSuffixedValue(const CompressionValue<T>& value);
    CompressionValue<T> GetNextFinalSuffixedValue(const CompressionValue<T>& value);
    bool GetIsLeafLevelDecompression() const;
    bool IsNextIterationLeafLevelDecompression() const;
    void InitializeInteriorDecompressionIteration();
    void InitializeLeafDecompressionIteration();

    size_t level_byte_offset_{0};
    bit_vector::BitVectorView alphabet_;
    bit_map::BitMapView encoded_prefix_length_;
    bit_vector::BitVectorView prefix_bits_;

    std::unique_ptr<cmc::entropy_coding::arithmetic_coding::Decoder> entropy_decoder_{nullptr};

    int count_adaptation_step_{0};
    bool is_leaf_level_decompression_{false};
};

template<typename T>
inline bool
PrefixEmbeddedDecompressionAdaptData<T>::IsDecompressionProgressing() const
{
    return (level_byte_offset_ < cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>::encoded_data_byte_stream_.size());
}

template<typename T>
inline uint32_t
PrefixEmbeddedDecompressionAdaptData<T>::GetNextPrefixLength()
{
    cmc_assert(entropy_decoder_ != nullptr);
    return entropy_decoder_->DecodeNextSymbol();
}

template<typename T>
inline std::vector<uint8_t>
PrefixEmbeddedDecompressionAdaptData<T>::GetNextPrefixBitSequence(const size_t num_prefix_bits)
{
    return prefix_bits_.GetNextBitSequence(num_prefix_bits);
}


template <typename T>
std::vector<CompressionValue<T>>
PrefixEmbeddedDecompressionAdaptData<T>::DecodeRootLevel(const t8_locidx_t num_local_root_values)
{
    cmc_debug_msg("The setup of the root level values is performed.");

    std::vector<CompressionValue<T>> root_values;
    root_values.reserve(num_local_root_values);

    this->InitializeDecompressionIteration();

    for (t8_locidx_t idx = 0; idx < num_local_root_values; ++idx)
    {
        /* Get the next suffixed value */
        const CompressionValue<T> suffixed_value = GetNextSuffixedValue(CompressionValue<T>());
        root_values.push_back(suffixed_value);
    }

    this->FinalizeDecompressionIteration();

    return root_values;
}

template <typename T>
bool
PrefixEmbeddedDecompressionAdaptData<T>::IsNextIterationLeafLevelDecompression() const 
{
    /* We need to check against "count_adaptation_step_ - 1" because the root level decoding is counted as a decomrpession step as well */
    return (this->GetAmrMesh().GetInitialRefinementLevel() == count_adaptation_step_);
}

template <typename T>
void
PrefixEmbeddedDecompressionAdaptData<T>::InitializeDecompressionIteration()
{
    /* Check if the following decompression iteration is the leaf level decompression */
    is_leaf_level_decompression_ = IsNextIterationLeafLevelDecompression();
    cmc_debug_msg("\n\nIs this iteration a leaf level decompression? ", is_leaf_level_decompression_, " and init ref elvel: ", this->GetAmrMesh().GetInitialRefinementLevel(), "\n\n");

    if (is_leaf_level_decompression_ == true)
    {
        InitializeLeafDecompressionIteration();
    } else
    {
        InitializeInteriorDecompressionIteration();
    }
}

template <typename T>
void
PrefixEmbeddedDecompressionAdaptData<T>::InitializeInteriorDecompressionIteration()
{
    cmc_debug_msg("A prefix decompression iteration is initialized.");
    
    constexpr size_t offset = sizeof(uint64_t);

    size_t processed_bytes = level_byte_offset_;

    /* Get a pointer to the beginning of the encoded data stream */
    const auto data_start_ptr = cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>::encoded_data_byte_stream_.data();

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


template <typename T>
void
PrefixEmbeddedDecompressionAdaptData<T>::InitializeLeafDecompressionIteration()
{
    cmc_debug_msg("A prefix decompression iteration is initialized.");
    
    /* Check if the following decompression iteration is the leaf level decompression */
    is_leaf_level_decompression_ = IsNextIterationLeafLevelDecompression();
    cmc_debug_msg("\n\nIs this iteration a leaf level decompression? ", is_leaf_level_decompression_, " and init ref elvel: ", this->GetAmrMesh().GetInitialRefinementLevel(), "\n\n");
    constexpr size_t offset = sizeof(uint64_t);

    size_t processed_bytes = level_byte_offset_;

    /* Get a pointer to the beginning of the encoded data stream */
    const auto data_start_ptr = cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>::encoded_data_byte_stream_.data();

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

template <typename T>
void
PrefixEmbeddedDecompressionAdaptData<T>::FinalizeDecompressionIteration()
{
    ++count_adaptation_step_;
    cmc_debug_msg("The prefix decompression iteration (", count_adaptation_step_, ") has been finalized.");
}

template <typename T>
void
PrefixEmbeddedDecompressionAdaptData<T>::CompleteDecompressionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    //Nothing to be done here!
}

template <typename T>
void
PrefixEmbeddedDecompressionAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

template <typename T>
inline bool
PrefixEmbeddedDecompressionAdaptData<T>::GetIsLeafLevelDecompression() const
{
    return is_leaf_level_decompression_;
}

template <typename T>
inline uint32_t
PrefixEmbeddedDecompressionAdaptData<T>::ApplyProcessBoundarySymbol()
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

template <typename T>
CompressionValue<T>
PrefixEmbeddedDecompressionAdaptData<T>::GetNextSuffixedValue(const CompressionValue<T>& value)
{
    if (this->GetIsLeafLevelDecompression())
    {
        return GetNextFinalSuffixedValue(value);
    } else
    {
        return GetNextInteriorSuffixedValue(value);
    }
}
static int idxx = 0;

template <typename T>
CompressionValue<T>
PrefixEmbeddedDecompressionAdaptData<T>::GetNextInteriorSuffixedValue(const CompressionValue<T>& value)
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

    //if (count_adaptation_step_ < 3)
    //{
    //    cmc_debug_msg("Idx: ", idxx, ", Current val front bit: ", static_cast<int>(suffixed_value.GetFrontBit()), ", tail bit: ", static_cast<int>(suffixed_value.GetTailBit()), ", suff_length: ", suffix_length);
    //    ++idxx;
    //}

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


template <typename T>
CompressionValue<T>
PrefixEmbeddedDecompressionAdaptData<T>::GetNextFinalSuffixedValue(const CompressionValue<T>& value)
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

template <typename T>
cmc::decompression::embedded::RefinementData<T>
PrefixEmbeddedDecompressionAdaptData<T>::PerformRefinement(const int which_tree, const int lelement_id, const CompressionValue<T> value, const int num_refined_elements)
{
    /* Create the refinement data to return */
    cmc::decompression::embedded::RefinementData<T> refinement_data;
    refinement_data.fine_values.reserve(num_refined_elements);

    /* Apply all children suffixes */
    for (int idx = 0; idx < num_refined_elements; ++idx)
    {
        /* Get the next value with the applied residual and store it wihtin the refinement data */
        refinement_data.fine_values.emplace_back(this->GetNextSuffixedValue(value));
    }

    return refinement_data;
}

template <typename T>
cmc::decompression::embedded::UnchangedData<T>
PrefixEmbeddedDecompressionAdaptData<T>::ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value)
{
    /* Get the netx suffixed value */
    const CompressionValue<T> suffixed_value = GetNextSuffixedValue(value);

    return cmc::decompression::embedded::UnchangedData<T>(suffixed_value);
}

template <typename T>
inline cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>*
CreatePrefixDecompressionAdaptationClass(cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>* abstract_var)
{
    return new PrefixEmbeddedDecompressionAdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyPrefixDecompressionAdaptationClass(cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class DecompressionVariable : public cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>
{
public:
    DecompressionVariable() = delete;

    DecompressionVariable(const std::string& name, std::vector<uint8_t>&& encoded_data_byte_stream, std::vector<uint8_t>&& encoded_mesh_byte_stream, VariableAttributes<T>&& attributes, const bool are_refinement_bits_stored, const MPI_Comm comm)
    : cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>(std::move(encoded_data_byte_stream), std::move(encoded_mesh_byte_stream), std::move(attributes), are_refinement_bits_stored, comm)
    {
        this->SetName(name);
        cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>::adaptation_creator_ = CreatePrefixDecompressionAdaptationClass<T>;
        cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>::adaptation_destructor_ = DestroyPrefixDecompressionAdaptationClass<T>;
        cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>::mesh_decoder_ = std::make_unique<mesh_compression::EmbeddedMeshDecoder>(this->GetEncodedMeshStream());
    };
};

}

#endif /* !CMC_EMBEDDED_PREFIX_EXTRACTION_DECOMPRESSION_PLAIN_SUFFIXES_HXX */
