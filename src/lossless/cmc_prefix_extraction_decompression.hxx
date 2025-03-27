#ifndef CMC_PREFIX_EXTRACTION_DECOMPRESSION_HXX
#define CMC_PREFIX_EXTRACTION_DECOMPRESSION_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_serialization.hxx"
#include "decompression/cmc_byte_decompression_variable.hxx"
#include "utilities/cmc_arithmetic_encoding.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "mesh_compression/cmc_mesh_decoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>
#include <memory>

namespace cmc::lossless::prefix
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template<typename T>
class PrefixDecompressionAdaptData : public cmc::decompression::IDecompressionAdaptData<T>
{
public:
    PrefixDecompressionAdaptData() = delete;
    PrefixDecompressionAdaptData(cmc::decompression::AbstractByteDecompressionVariable<T>* variable)
    : cmc::decompression::IDecompressionAdaptData<T>(variable) {};

    bool IsDecompressionProgressing() const override;
    void InitializeDecompressionIteration() override;
    void FinalizeDecompressionIteration() override;
    void CompleteDecompressionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::vector<CompressionValue<T>> DecodeRootLevel(const t8_locidx_t num_local_root_values) override;
protected:
    cmc::decompression::RefinementData<T> PerformRefinement(const int which_tree, const int lelement_id, const CompressionValue<T> value, const int num_refined_elements) override;
    cmc::decompression::UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) override;

private:
    uint32_t GetNextPrefixLength();
    std::vector<uint8_t> GetNextPrefixBitSequence(const size_t num_bits);
    CompressionValue<T> GetNextSuffixedValue(const CompressionValue<T>& value);

    size_t level_byte_offset_{0};
    bit_vector::BitVectorView alphabet_;
    bit_map::BitMapView encoded_prefix_length_;
    bit_vector::BitVectorView prefix_bits_;

    std::unique_ptr<cmc::entropy_coding::arithmetic_coding::Decoder> entropy_decoder_{nullptr};

    int count_adaptation_step_{0};
};

template<typename T>
inline bool
PrefixDecompressionAdaptData<T>::IsDecompressionProgressing() const
{
    return (level_byte_offset_ < cmc::decompression::IDecompressionAdaptData<T>::encoded_data_byte_stream_.size());
}

template<typename T>
inline uint32_t
PrefixDecompressionAdaptData<T>::GetNextPrefixLength()
{
    cmc_assert(entropy_decoder_ != nullptr);
    return entropy_decoder_->DecodeNextSymbol();
}

template<typename T>
inline std::vector<uint8_t>
PrefixDecompressionAdaptData<T>::GetNextPrefixBitSequence(const size_t num_prefix_bits)
{
    return prefix_bits_.GetNextBitSequence(num_prefix_bits);
}


template <typename T>
std::vector<CompressionValue<T>>
PrefixDecompressionAdaptData<T>::DecodeRootLevel(const t8_locidx_t num_local_root_values)
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
void
PrefixDecompressionAdaptData<T>::InitializeDecompressionIteration()
{
    cmc_debug_msg("A prefix decompression iteration is initialized.");

    constexpr size_t offset = sizeof(size_t);

    size_t processed_bytes = level_byte_offset_;

    /* Get a pointer to the beginning of the encoded data stream */
    const auto data_start_ptr = cmc::decompression::IDecompressionAdaptData<T>::encoded_data_byte_stream_.data();

    /* Get the amount of relevant bytes for this decompression level */
    const size_t current_level_bytes = GetValueFromByteStream<size_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;

    cmc_debug_msg("The current refinement level is described by ", current_level_bytes, " bytes.");
    
    /* Get the bytes for the encoded alphabet */
    const size_t alphabet_bytes = GetValueFromByteStream<size_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;
    cmc_debug_msg("alphabet_bytes size: ", alphabet_bytes);

    /* Get the bytes for the encoded prefix lengths */
    const size_t encoded_prefix_length_bytes = GetValueFromByteStream<size_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;
    cmc_debug_msg("encoded_prefix_length_bytes size: ", encoded_prefix_length_bytes);

    /* Get the bytes for the remaining bits */
    const size_t prefix_bytes = GetValueFromByteStream<size_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;
    cmc_debug_msg("prefix_bytes size: ", prefix_bytes);

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

    cmc_debug_msg("Processed bytes: ", level_byte_offset_);
    /* Decode the alphabet */
    [[maybe_unused]] auto [frequency_model, num_alphabet_bytes] = cmc::entropy_coding::arithmetic_coding::DecodeStaticFrequencyAlphabet(alphabet_.begin());

    /* Setup the entropy decoder for the prefix lengths */
    entropy_decoder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::Decoder>(std::make_unique<cmc::entropy_coding::arithmetic_coding::StaticFrequencyModel>(frequency_model), encoded_prefix_length_);
}

template <typename T>
void
PrefixDecompressionAdaptData<T>::FinalizeDecompressionIteration()
{
    ++count_adaptation_step_;
    cmc_debug_msg("The prefix decompression iteration (", count_adaptation_step_, ") has been finalized.");
}

template <typename T>
void
PrefixDecompressionAdaptData<T>::CompleteDecompressionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    //Nothing to be done here!
}

template <typename T>
void
PrefixDecompressionAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

template <typename T>
CompressionValue<T>
PrefixDecompressionAdaptData<T>::GetNextSuffixedValue(const CompressionValue<T>& value)
{
    CompressionValue<T> suffixed_value = value;

    /* Get the length of the suffix to be appended */
    const int suffix_length = static_cast<int>(this->GetNextPrefixLength());

    cmc_debug_msg("Next suffix length: ", suffix_length);
    if (suffix_length > 0)
    {
        /* Get the actual suffix */
        const std::vector<uint8_t> suffix = this->GetNextPrefixBitSequence(suffix_length);

        /* Apply the suffix to the current value */
        suffixed_value.ApplySuffix(suffix, suffix_length);
    }

    return suffixed_value;
}


template <typename T>
cmc::decompression::RefinementData<T>
PrefixDecompressionAdaptData<T>::PerformRefinement(const int which_tree, const int lelement_id, const CompressionValue<T> value, const int num_refined_elements)
{
    /* Create the refinement data to return */
    cmc::decompression::RefinementData<T> refinement_data;
    refinement_data.fine_values.reserve(num_refined_elements);

    /* Apply all children residuals */
    for (int idx = 0; idx < num_refined_elements; ++idx)
    {
        /* Get the next value with the applied residual and store it wihtin the refinement data */
        refinement_data.fine_values.emplace_back(this->GetNextSuffixedValue(value));
    }

    return refinement_data;
}

template <typename T>
cmc::decompression::UnchangedData<T>
PrefixDecompressionAdaptData<T>::ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value)
{
    /* Get the netx suffixed value */
    const CompressionValue<T> suffixed_value = GetNextSuffixedValue(value);

    return cmc::decompression::UnchangedData<T>(suffixed_value);
}

template <typename T>
inline cmc::decompression::IDecompressionAdaptData<T>*
CreatePrefixDecompressionAdaptationClass(cmc::decompression::AbstractByteDecompressionVariable<T>* abstract_var)
{
    return new PrefixDecompressionAdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyPrefixDecompressionAdaptationClass(cmc::decompression::IDecompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class DecompressionVariable : public cmc::decompression::AbstractByteDecompressionVariable<T>
{
public:
    DecompressionVariable() = delete;

    DecompressionVariable(const std::string& name, std::vector<uint8_t>&& encoded_data_byte_stream, std::vector<uint8_t>&& encoded_mesh_byte_stream)
    : cmc::decompression::AbstractByteDecompressionVariable<T>(std::move(encoded_data_byte_stream), std::move(encoded_mesh_byte_stream))
    {
        this->SetName(name);
        cmc::decompression::AbstractByteDecompressionVariable<T>::adaptation_creator_ = CreatePrefixDecompressionAdaptationClass<T>;
        cmc::decompression::AbstractByteDecompressionVariable<T>::adaptation_destructor_ = DestroyPrefixDecompressionAdaptationClass<T>;
        cmc::decompression::AbstractByteDecompressionVariable<T>::mesh_decoder_ = std::make_unique<mesh_compression::MeshDecoder>(this->GetEncodedMeshStream());
    };
};

}

#endif /* !CMC_PREFIX_EXTRACTION_DECOMPRESSION_HXX */
