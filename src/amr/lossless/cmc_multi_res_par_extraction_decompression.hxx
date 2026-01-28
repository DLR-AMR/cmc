#ifndef CMC_MULTI_RES_PAR_EXTRACTION_DECOMPRESSION_HXX
#define CMC_MULTI_RES_PAR_EXTRACTION_DECOMPRESSION_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_huffman_coder.hxx"
//#include "utilities/cmc_multi_res_extraction_util.hxx"
#include "amr/lossless/cmc_byte_par_decompression_variable.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "mesh_compression/cmc_mesh_par_decoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>
#include <memory>

namespace cmc::lossless::par::multi_res
{
   
/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template<typename T>
class MultiResDecompressionAdaptData : public cmc::decompression::par::IDecompressionAdaptData<T>
{
public:
    MultiResDecompressionAdaptData() = delete;
    MultiResDecompressionAdaptData(cmc::decompression::AbstractByteDecompressionVariable<T>* variable)
    : cmc::decompression::IDecompressionAdaptData<T>(variable) {};
    ~MultiResDecompressionAdaptData() = default;
    
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
MultiResDecompressionAdaptData<T>::IsDecompressionProgressing() const
{
    return (level_byte_offset_ < cmc::decompression::IDecompressionAdaptData<T>::encoded_data_byte_stream_.size());
}


template<typename T>
inline uint32_t
MultiResDecompressionAdaptData<T>::GetNextEncodedResidualLength()
{
    cmc_assert(entropy_decoder_ != nullptr);
    return entropy_decoder_->DecodeNextSymbol();
}

template<typename T>
inline std::vector<uint8_t>
MultiResDecompressionAdaptData<T>::GetNextResidualBitSequence(const size_t num_residual_bits)
{
    return residual_bits_.GetNextBitSequence(num_residual_bits);
}


template <typename T>
std::vector<CompressionValue<T>>
MultiResDecompressionAdaptData<T>::DecodeRootLevel(const t8_locidx_t num_local_root_values)
{
    cmc_debug_msg("The setup of the root level values is performed.");

    std::vector<CompressionValue<T>> root_values;
    root_values.reserve(num_local_root_values);

    const size_t offset = sizeof(T);

    for (t8_locidx_t idx = 0; idx < num_local_root_values; ++idx)
    {
        const T val = GetValueFromByteStream<T>(cmc::decompression::IDecompressionAdaptData<T>::encoded_data_byte_stream_.data() + offset * idx);
        cmc_debug_msg("Root level value: ", val, ", fuer idx: ", idx);
        root_values.emplace_back(CompressionValue<T>(val));
    }

    level_byte_offset_ += num_local_root_values * sizeof(T);
    
    return root_values;
}

template <typename T>
void
MultiResDecompressionAdaptData<T>::InitializeDecompressionIteration()
{
    cmc_debug_msg("A parallel multi-resolution decompression iteration is initialized.");

    constexpr size_t offset = sizeof(uint64_t);

    size_t processed_bytes = level_byte_offset_;

    /* Get a pointer to the beginning of the encoded data stream */
    const auto data_start_ptr = cmc::decompression::IDecompressionAdaptData<T>::encoded_data_byte_stream_.data();

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
    entropy_decoder_ = std::make_unique<typename cmc::entropy_coding::arithmetic_coding::MultiResDecoder<T>>(alphabet_.begin(), encoded_lzcs_);
    entropy_decoder_->SetupDecoding(); 
}

template <typename T>
void
MultiResDecompressionAdaptData<T>::FinalizeDecompressionIteration()
{
    ++count_adaptation_step_;
    cmc_debug_msg("The multi-resolution decompression iteration (", count_adaptation_step_, ") has been finalized.");
}

template <typename T>
void
MultiResDecompressionAdaptData<T>::CompleteDecompressionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    //Nothing to be done here!
}

template <typename T>
void
MultiResDecompressionAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

template <typename T>
inline uint32_t
MultiResDecompressionAdaptData<T>::ApplyProcessBoundarySymbol()
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
MultiResDecompressionAdaptData<T>::GetNextResidualAppliedValue(const CompressionValue<T>& coarse_value)
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
cmc::decompression::RefinementData<T>
MultiResDecompressionAdaptData<T>::PerformRefinement([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, const CompressionValue<T> value, const int num_refined_elements)
{
    /* Create the refinement data to return */
    cmc::decompression::RefinementData<T> refinement_data;
    refinement_data.fine_values.reserve(num_refined_elements);

    /* Apply all children residuals */
    for (int idx = 0; idx < num_refined_elements; ++idx)
    {
        /* Get the next value with the applied residual and store it wihtin the refinement data */
        refinement_data.fine_values.emplace_back(this->GetNextResidualAppliedValue(value));
    }

    return refinement_data;
}

template <typename T>
cmc::decompression::UnchangedData<T>
MultiResDecompressionAdaptData<T>::ElementStaysUnchanged([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, const CompressionValue<T>& value)
{
    if constexpr (kDecodeLessResiduals)
    {
        /* There is no residual to decode and the element remains unchanged */
        return cmc::decompression::UnchangedData<T>(value);
    } else
    {
        /* Get the next suffixed value */
        const CompressionValue<T> residual_applied_value = this->GetNextResidualAppliedValue(value);

        const T residual_applied_value_reinter = residual_applied_value.template ReinterpretDataAs<T>();
        const T init_val = value.template ReinterpretDataAs<T>();
        if (not ApproxCompare(residual_applied_value_reinter, init_val))
        {
            cmc_debug_msg("Hier ist der Wert vor element unchanged ein anderer als nachher: residual applied: ", residual_applied_value_reinter, ", init val: ", init_val);
        }
        return cmc::decompression::UnchangedData<T>(residual_applied_value);
    }
}

template <typename T>
inline cmc::decompression::IDecompressionAdaptData<T>*
CreateMultiResDecompressionAdaptationClass(cmc::decompression::AbstractByteDecompressionVariable<T>* abstract_var)
{
    return new MultiResDecompressionAdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyMultiResDecompressionAdaptationClass(cmc::decompression::IDecompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class DecompressionVariable : public cmc::decompression::par::AbstractByteDecompressionVariable<T>
{
public:
    DecompressionVariable() = delete;

    explciit DecompressionVariable(const std::string& name, std::vector<uint8_t>&& global_level_num_elems, std::vector<uint8_t>&& encoded_mesh_stream, 
        std::vector<uint8_t>&& global_level_data_bytes, const uint64_t file_byte_offset_encoded_data, 
        std::vector<ProcLevelByteStreamOffsets>&& level_offset_hints, const std::string& file_name, const MPI_Comm comm)
    : cmc::decompression::par::AbstractByteParDecompressionVariable(name, std::move(global_level_num_elems), std::move(encoded_mesh_stream), std::move(global_level_data_bytes), file_byte_offset_encoded_data, std::move(level_offset_hints), file_name, comm)
    {
        cmc::decompression::par::AbstractByteDecompressionVariable<T>::adaptation_creator_ = CreateMultiResDecompressionAdaptationClass<T>;
        cmc::decompression::par::AbstractByteDecompressionVariable<T>::adaptation_destructor_ = DestroyMultiResDecompressionAdaptationClass<T>;
        cmc::decompression::par::AbstractByteDecompressionVariable<T>::mesh_decoder_ = std::make_unique<mesh_compression::MeshDecoder>(this->GetEncodedMeshStreamPtr());
    }

    void SetupLevelDecodingStart(const uint8_t* global_var_start, const t8_gloidx_t proc_lvl_elem_offset, const ProcLevelByteStreamOffsets& offset_hints) override;

private:
    std::unique_ptr<cmc::entropy_coding::huffman::HuffmanDecoder<T>> entropy_decoder_{nullptr};

};

template<class T>
void
DecompressionVariable<T>::SetupLevelDecodingStart(const uint8_t* global_var_start, const t8_gloidx_t proc_lvl_entropy_offset, const ProcLevelByteStreamOffsets& hints)
{
    cmc_assert(hints.offset_hints.size() >= 1);

    size_t offset{0};

    /* Get the global number of bytes for this level */
    [[maybe_unused]] const uint64_t global_level_bytes = GetValueFromByteStream<uint64_t>(global_var_start + offset);
    offset += sizeof(uint64_t);

    /* Get the number of bytes for the encoding (without the Huffman symbol frequency table) */
    const uint64_t global_level_encoding_bytes = GetValueFromByteStream<uint64_t>(global_var_start + offset);
    offset += sizeof(uint64_t);

    /* Get the number of symbols in the Huffman symbol requency table */
    [[maybe_unused]] const uint64_t num_huff_symbols = GetValueFromByteStream<uint64_t>(global_var_start + offset);
    offset += sizeof(uint64_t);

    /* Find the best-suited start-point for tht entropy offset based on the hints */
    auto start_iter = std::lower_bound(hints.offset_hints.begin(), hints.offset_hints.end(), proc_lvl_entropy_offset, [](const OffsetHint& hint, const t8_gloidx_t& entropy_offset){
        return hint.entropy_code_id < entropy_offset;
    });

    /* Choose the correct start position for the given entropy offset */
    auto search_start_iter = (start_iter == hints.offset_hints.end() ? std::prev(hints.offset_hints.end()) : start_iter);

    /* Decode the entropy dictionary at the beginning of the encoding of the level */
    entropy_decoder_ = std::make_unique<cmc::entropy_coding::huffman::HuffmanDecoder<uint32_t>>(global_var_start + offset);

    /* Update the start position by offsetting the bytes needed for the reconstruction of the symbol-frequency-table */
    offset += entropy_decoder_->GetNumberOfProcessedBytesForSymbolFrequencyTable();

    cmc_assert(search_start_iter->entropy_code_id <= proc_lvl_entropy_offset); 

    /* Set the vierw to encoded data */
    bit_vector::BitVectorView encoded_level_data(global_var_start + offset, global_level_encoding_bytes);
    entropy_decoder_->StartDecoding(encoded_level_data);

    /* Move to the best-suited start position */
    if (search_start_iter->entropy_code_id != proc_lvl_entropy_offset)
    {
        /* Compute the number of entropy codes that need to be skipped */
        const int num_codes_to_skip = proc_lvl_entropy_offset - search_start_iter->entropy_code_id;

        /* Iterate until the beginning of the process-local start has been reached and leave the variable in this state */
        for (int code_iter{0}; code_iter < num_codes_to_skip; ++code_iter)
        {
            /* Decode the current entropy encoding */
            const int32_t entropy_symbol = entropy_decoder_->DecodeNextSymbol();

            /* Convert the entropy_symbol to an actual symbol (i.e. first one bit position) */
            auto [_, lzc] = ConvertFrequencySymbolToValue(entropy_symbol);

            /* Check if it has been a process boundary */
            if ()
            {
                //How to get this information into here, because the proc offset does not account for process boundary symbols rn
                //....
                //Store actual partitioning in an addtional var for each level?
                Pro Level <ElemOffset; EntropyOffset; ByteOffset> als Tuple 3x 32 Bits 
                Bei EntropyOffsetCalculation: Wenn 
            }

            /* Determine the encoded residual and account for the potential one bit that has been discarded */
            if (lzc < bit_vector::kCharBit * sizeof(T))
            {
                /* Compute the encoded residual */
                const uint32_t residual_length = bit_vector::kCharBit * sizeof(T) - lzc - 1;
                
                /* Skip the corresponding amount of encoded significant bits */
                entropy_decoder_->SkipNextNumberOfBits(residual_length);
            }

        }
    }
} 


}


#endif /* !CMC_MULTI_RES_PAR_EXTRACTION_DECOMPRESSION_HXX */
