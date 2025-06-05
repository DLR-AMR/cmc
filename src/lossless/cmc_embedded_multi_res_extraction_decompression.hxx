#ifndef CMC_EMBEDDED_MULTI_RES_EXTRACTION_DECOMPRESSION_HXX
#define CMC_EMBEDDED_MULTI_RES_EXTRACTION_DECOMPRESSION_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_serialization.hxx"
#include "lossless/cmc_embedded_byte_decompression_variable.hxx"
#include "utilities/cmc_multi_res_entropy_coder.hxx"
#include "utilities/cmc_multi_res_extraction_util.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "mesh_compression/cmc_embedded_mesh_decoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>
#include <memory>

namespace cmc::lossless::embedded::multi_res
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template<typename T>
class MultiResEmbeddedDecompressionAdaptData : public cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>
{
public:
    MultiResEmbeddedDecompressionAdaptData() = delete;
    MultiResEmbeddedDecompressionAdaptData(cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>* variable)
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
    uint32_t GetNextEncodedResidualLength();
    std::vector<uint8_t> GetNextResidualBitSequence(const size_t num_bits);
    CompressionValue<T> GetNextResidualAppliedValue(const CompressionValue<T>& value);
    uint32_t ApplyProcessBoundarySymbol();

    size_t level_byte_offset_{0};
    bit_vector::BitVectorView alphabet_;
    bit_map::BitMapView encoded_lzcs_;
    bit_vector::BitVectorView residual_bits_;

    std::unique_ptr<cmc::entropy_coding::arithmetic_coding::Decoder> entropy_decoder_{nullptr};

    int count_adaptation_step_{0};
};

template<typename T>
inline bool
MultiResEmbeddedDecompressionAdaptData<T>::IsDecompressionProgressing() const
{
    return (level_byte_offset_ < cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>::encoded_data_byte_stream_.size());
}

template<typename T>
inline uint32_t
MultiResEmbeddedDecompressionAdaptData<T>::GetNextEncodedResidualLength()
{
    cmc_assert(entropy_decoder_ != nullptr);
    return entropy_decoder_->DecodeNextSymbol();
}

template<typename T>
inline std::vector<uint8_t>
MultiResEmbeddedDecompressionAdaptData<T>::GetNextResidualBitSequence(const size_t num_residual_bits)
{
    return residual_bits_.GetNextBitSequence(num_residual_bits);
}


template <typename T>
std::vector<CompressionValue<T>>
MultiResEmbeddedDecompressionAdaptData<T>::DecodeRootLevel(const t8_locidx_t num_local_root_values)
{
    cmc_debug_msg("The setup of the root level values is performed.");

    std::vector<CompressionValue<T>> root_values;
    root_values.reserve(num_local_root_values);

    const size_t offset = sizeof(T);

    for (t8_locidx_t idx = 0; idx < num_local_root_values; ++idx)
    {
        const T val = GetValueFromByteStream<T>(cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>::encoded_data_byte_stream_.data() + offset * idx);
        cmc_debug_msg("Root level value: ", val, ", fuer idx: ", idx);
        root_values.emplace_back(CompressionValue<T>(val));
    }

    level_byte_offset_ += num_local_root_values * sizeof(T);
    
    return root_values;
}

template <typename T>
void
MultiResEmbeddedDecompressionAdaptData<T>::InitializeDecompressionIteration()
{
    cmc_debug_msg("A multi-resolution decompression iteration is initialized.");

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

template <typename T>
void
MultiResEmbeddedDecompressionAdaptData<T>::FinalizeDecompressionIteration()
{
    ++count_adaptation_step_;
    cmc_debug_msg("The multi-resolution decompression iteration (", count_adaptation_step_, ") has been finalized.");
}

template <typename T>
void
MultiResEmbeddedDecompressionAdaptData<T>::CompleteDecompressionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    //Nothing to be done here!
}

template <typename T>
void
MultiResEmbeddedDecompressionAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}


template <typename T>
inline uint32_t
MultiResEmbeddedDecompressionAdaptData<T>::ApplyProcessBoundarySymbol()
{
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

template <typename T>
CompressionValue<T>
MultiResEmbeddedDecompressionAdaptData<T>::GetNextResidualAppliedValue(const CompressionValue<T>& coarse_value)
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


template <typename T>
cmc::decompression::embedded::RefinementData<T>
MultiResEmbeddedDecompressionAdaptData<T>::PerformRefinement(const int which_tree, const int lelement_id, const CompressionValue<T> value, const int num_refined_elements)
{
    /* Create the refinement data to return */
    cmc::decompression::embedded::RefinementData<T> refinement_data;
    refinement_data.fine_values.reserve(num_refined_elements);

    /* Apply all children suffixes */
    for (int idx = 0; idx < num_refined_elements; ++idx)
    {
        /* Get the next value with the applied residual and store it wihtin the refinement data */
        refinement_data.fine_values.emplace_back(this->GetNextResidualAppliedValue(value));
    }

    return refinement_data;
}

template <typename T>
cmc::decompression::embedded::UnchangedData<T>
MultiResEmbeddedDecompressionAdaptData<T>::ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value)
{
    /* Get the netx suffixed value */
    const CompressionValue<T> residual_applied_value = GetNextResidualAppliedValue(value);

    return cmc::decompression::embedded::UnchangedData<T>(residual_applied_value);
}

template <typename T>
inline cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>*
CreateMultiResDecompressionAdaptationClass(cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>* abstract_var)
{
    return new MultiResEmbeddedDecompressionAdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyMultiResDecompressionAdaptationClass(cmc::decompression::embedded::IEmbeddedDecompressionAdaptData<T>* iadapt_data)
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
        cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>::adaptation_creator_ = CreateMultiResDecompressionAdaptationClass<T>;
        cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>::adaptation_destructor_ = DestroyMultiResDecompressionAdaptationClass<T>;
        cmc::decompression::embedded::AbstractEmbeddedByteDecompressionVariable<T>::mesh_decoder_ = std::make_unique<mesh_compression::EmbeddedMeshDecoder>(this->GetEncodedMeshStream());
    };
};

}


#endif /* !CMC_EMBEDDED_MULTI_RES_EXTRACTION_DECOMPRESSION_HXX */
