#ifndef CMC_EMBEDDED_PREFIX_EXTRACTION_COMPRESSION_HXX
#define CMC_EMBEDDED_PREFIX_EXTRACTION_COMPRESSION_HXX


#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "lossless/cmc_embedded_byte_compression_variable.hxx"
#include "utilities/cmc_prefix_entropy_coder.hxx"
#include "utilities/cmc_serialization.hxx"
#include "mesh_compression/cmc_embedded_mesh_encoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>

namespace cmc::lossless::embedded::prefix
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template<typename T>
class PrefixEmbeddedAdaptData : public IEmbeddedCompressionAdaptData<T>
{
public:
    PrefixEmbeddedAdaptData() = delete;
    PrefixEmbeddedAdaptData(AbstractEmbeddedByteCompressionVariable<T>* variable)
    : IEmbeddedCompressionAdaptData<T>(variable) {
        IEmbeddedCompressionAdaptData<T>::entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::PrefixEncoder<T>>();
    };

    void InitializeExtractionIteration() override;
    void FinalizeExtractionIteration() override;
    void CompleteExtractionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const override;
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const override;

protected:
    ExtractionData<T> PerformExtraction(const int which_tree, const int lelement_id, const int num_elements, const VectorView<CompressionValue<T>> values) override;
    UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) override;

private:
    void CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const;

    int count_adaptation_step_{0};
};


template <typename T>
void
PrefixEmbeddedAdaptData<T>::InitializeExtractionIteration()
{
    //There is currently nothing to be done here!
}

template <typename T>
void
PrefixEmbeddedAdaptData<T>::FinalizeExtractionIteration()
{
    ++count_adaptation_step_;
}

template <typename T>
void
PrefixEmbeddedAdaptData<T>::CompleteExtractionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    //Nothing to be done here!
}

template <typename T>
void
PrefixEmbeddedAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

template<typename T>
std::pair<bool, CompressionValue<T>>
EvaluateEmbeddedCommonPrefix(const VectorView<CompressionValue<T>>& compression_values)
{
    cmc_assert(compression_values.size() >= 2);

    /* Check if all elements are holding an actual prefix */
    for (auto cv_iter = compression_values.begin(); cv_iter != compression_values.end(); ++cv_iter)
    {
        if (cv_iter->IsEmpty())
        {
            /* Since this prefix value is empty, there cannot be a common prefix */
            return std::make_pair(false, CompressionValue<T>());
        }
    }

    /* Determine a common prefix of the first two values */
    CompressionValue<T> prefix = GetCommonPrefix<sizeof(T)>(compression_values[0], compression_values[1]);

    /* Check if there is common prefix between the first two values */
    if (prefix.GetCountOfSignificantBits() == 0)
    {
        /* There is no common prefix */
        return std::make_pair(false, CompressionValue<T>());
    }

    /* Check if there is a common prefix with the other values within the view */
    for (size_t index = 2; index < compression_values.size(); ++index)
    {
        /* Find a common prefix of all variables */
        prefix = GetCommonPrefix(prefix, compression_values[index]);

        if (prefix.GetCountOfSignificantBits() == 0)
        {
            /* There is no common prefix */
            return std::make_pair(false, CompressionValue<T>());
        }
    }

    /* If the function arrives here, we do have found a common prefix which can be extracted from the 'previous prefixes' */
    return std::make_pair(true, prefix);
}

template <typename T>
ExtractionData<T>
PrefixEmbeddedAdaptData<T>::PerformExtraction([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, [[maybe_unused]] const int num_elements, const VectorView<CompressionValue<T>> values)
{
    /* Evaluate whether there is a common prefix to extract */
    auto [is_prefix_found, prefix] = EvaluateEmbeddedCommonPrefix<T>(values);

    if (is_prefix_found)
    {
        std::vector<CompressionValue<T>> fine_values;
        fine_values.reserve(values.size());

        /* We need to trim the previous prefixes by the extracted common prefix */
        for (auto val_iter = values.begin(); val_iter != values.end(); ++val_iter)
        {
            fine_values.push_back(*val_iter);
            fine_values.back().SetFrontBit(sizeof(T) * bit_map::kCharBit - prefix.GetTailBit());
        }

        /* Return the common prefix and the trimmed remaining errors */
        return ExtractionData<T>(prefix, std::move(fine_values));
    } else
    {
        std::vector<CompressionValue<T>> fine_values;
        fine_values.reserve(values.size());

        /* We copy the fine values over and leave them unchanged */
        for (auto val_iter = values.begin(); val_iter != values.end(); ++val_iter)
        {
            fine_values.push_back(*val_iter);
        }

        /* We return an empty prefix and a the unchanged previous prefixes */
        return ExtractionData<T>(CompressionValue<T>(), std::move(fine_values));
    }
}

template <typename T>
UnchangedData<T>
PrefixEmbeddedAdaptData<T>::ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value)
{
    /* We drag this value along to the "coarser level" until it this element is passed with its siblings
     *  as a family of elements into the adaptation callback. Moreover, we leave an empty value 
     * on the "finer level" */
    return UnchangedData<T>(value, CompressionValue<T>());
}

template <typename T>
void
PrefixEmbeddedAdaptData<T>::CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    /* Collect the symbol frequencies */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /* Get the length of this prefix */
        const uint32_t pref_length = static_cast<uint32_t>(val_iter->GetCountOfSignificantBits());

        /* Update the frequency */
        IEmbeddedCompressionAdaptData<T>::entropy_coder_->UpdateSymbolFrequency(pref_length);
    }
}

/**
 * @brief  We use an arithmetic encoder to encode the position of the first "one-bit" in the compression value
 * 
 * @tparam T The original data type of the underlying data (e.g. float)
 * @param level_byte_values The remaining "fine compression" values after an extraction iteration
 * @return std::vector<uint8_t> The encoded data stream
 */
template <typename T>
std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
PrefixEmbeddedAdaptData<T>::EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    cmc_debug_msg("The encoding of the CompressionValues after the prefix extraction iteration starts...");
    
    cmc_assert(IEmbeddedCompressionAdaptData<T>::entropy_coder_ != nullptr);

    /* Get the rank of the mpi process within the communicator */
    int rank{0};
    int ret_val = MPI_Comm_rank(this->GetMPIComm(), &rank);
    MPICheckError(ret_val);

    /* Define the root rank */
    const int root_rank = 0;

    /* The encoded data will be stored in a BitVector */
    cmc::bit_vector::BitVector encoding;
    encoding.Reserve(3 * level_byte_values.size());

    /* Reset the entropy coder and initialize the alphabet */
    IEmbeddedCompressionAdaptData<T>::entropy_coder_->Reset(std::make_unique<cmc::entropy_coding::arithmetic_coding::PrefixCompressionAlphabet<T>>());
    IEmbeddedCompressionAdaptData<T>::entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Collect the symbols and their frequencies for the entropy coder */
    CollectSymbolFrequenciesForEntropyCoding(level_byte_values);

    /* Setup the interior structure for encoding */
    IEmbeddedCompressionAdaptData<T>::entropy_coder_->SetupEncoding(this->GetMPIComm());

    /* Iterate over all values and encode them */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /* Get the current value */
        CompressionValue<T> current_val = *val_iter;

        /* Get the length of the prefix */
        const int pref_length = current_val.GetCountOfSignificantBits();

        /* Encode the prefix length */
        IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(static_cast<uint32_t>(pref_length));

        /* In case there is a prefix, the actual bits of the prefix will be stored */
        if (pref_length > 0)
        {
            encoding.AppendBits(current_val.GetSignificantBitsInBigEndianOrdering(), pref_length);
        }
    }

    /* We set an indication symbol that the process local end of the values have been reached */
    IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte);

    /* Indicate that the encoding has been finished and flush all pending encodings */
    IEmbeddedCompressionAdaptData<T>::entropy_coder_->FinishEncoding();

    /* Set up the BitVector holding the encoded data for further use */
    encoding.TrimToContent();

    /* Get the local remaining significant bits */
    const uint64_t local_remaining_significant_bits_num_bytes = static_cast<uint64_t>(encoding.size());

    /* Get the local encoded prefix lengths */
    cmc::bit_map::BitMap local_encoded_prefix_length_stream = IEmbeddedCompressionAdaptData<T>::entropy_coder_->GetEncodedBitStream();
    const uint64_t local_encoded_prefix_length_stream_num_bytes = static_cast<uint64_t>(local_encoded_prefix_length_stream.size_bytes());

    /* We need exchange the encoded lengths */
    const std::vector<uint64_t> local_bytes{local_encoded_prefix_length_stream_num_bytes, local_remaining_significant_bits_num_bytes};
    std::vector<uint64_t> global_bytes{0, 0};

    ret_val = MPI_Reduce(local_bytes.data(), global_bytes.data(), 2, MPI_UINT64_T, MPI_SUM, root_rank, this->GetMPIComm());
    MPICheckError(ret_val);

    /* Declare the buffers for the encoded data */
    std::vector<uint8_t> encoded_entropy_codes;
    std::vector<uint8_t> encoded_data;

    /* Only the root rank needs to encode the encoded sizes */
    if (rank == root_rank)
    {
        /* Get the encoded alphabet */
        cmc::bit_vector::BitVector encoded_alphabet = IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeAlphabet();
        const uint64_t encoded_alphabet_num_bytes = static_cast<uint64_t>(encoded_alphabet.size());

        /* Calculate the overall amount of bytes on the root rank */
        const uint64_t num_locally_encoded_entropy_codes_bytes = 4 * sizeof(uint64_t) + encoded_alphabet_num_bytes + local_encoded_prefix_length_stream_num_bytes;

        /* Allocate memory for the encoded data */
        encoded_entropy_codes.reserve(num_locally_encoded_entropy_codes_bytes);

        /** We store global information about the encoded level **/
        /* Push back the overall byte count for the level */
        const uint64_t num_global_bytes_encoded_level = 4 * sizeof(uint64_t) + encoded_alphabet_num_bytes + global_bytes[0] + global_bytes[1];
        PushBackValueToByteStream<uint64_t>(encoded_entropy_codes, num_global_bytes_encoded_level);

        /* Push back the byte count for the encoded alphabet */
        PushBackValueToByteStream<uint64_t>(encoded_entropy_codes, encoded_alphabet_num_bytes);

        /* Push back the byte count for the encoded first "one-bit" positions */
        PushBackValueToByteStream<uint64_t>(encoded_entropy_codes, global_bytes[0]);

        /* Push back the byte count for the remaining significant bits */
        PushBackValueToByteStream<uint64_t>(encoded_entropy_codes, global_bytes[1]);

        /* Afterwards, we store the encoded alphabet */
        std::copy_n(encoded_alphabet.begin(), encoded_alphabet_num_bytes, std::back_inserter(encoded_entropy_codes));
    
        /* Finally, copy the entropy codes */
        std::copy_n(local_encoded_prefix_length_stream.begin_bytes(), local_encoded_prefix_length_stream_num_bytes, std::back_insert_iterator(encoded_entropy_codes));
    } else
    {
        /* Otherwise, the rank only hold the entropy codes */
        local_encoded_prefix_length_stream.MoveDataInto(encoded_entropy_codes);
    }

    /* Get the encoded remaining significant bits */
    encoding.MoveDataInto(encoded_data);

    cmc_debug_msg("The entropy encoder of the prefix extraction compression completed the encoding of the CompressionValues of this iteration.");
    
    return std::make_pair(encoded_entropy_codes, encoded_data);
}

template <typename T>
inline std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
PrefixEmbeddedAdaptData<T>::EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const
{
    cmc_debug_msg("The encoding of the root level values of the prefix compression starts.");

    auto encoded_streams = EncodeLevelData(root_level_values);
    
    cmc_debug_msg("The entropy encoder of the prefix exctraction compression completed the encoding of the root-level CompressionValues.");

    return encoded_streams;
}

template <typename T>
inline IEmbeddedCompressionAdaptData<T>*
CreatePrefixEmbeddedExtractionAdaptationClass(AbstractEmbeddedByteCompressionVariable<T>* abstract_var)
{
    return new PrefixEmbeddedAdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyPrefixEmbeddedExtractionAdaptationClass(IEmbeddedCompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class EmbeddedCompressionVariable : public AbstractEmbeddedByteCompressionVariable<T>
{
public:
    EmbeddedCompressionVariable() = delete;

    EmbeddedCompressionVariable(input::Var& input_variable)
    : AbstractEmbeddedByteCompressionVariable<T>(input_variable)
    {
        //cmc_assert(input_variable.IsValid());//TODO: implement 

        this->SetName(input_variable.GetName());
        this->IndicateWhetherMeshRefinementBitsWillBeStored(true);

        AbstractEmbeddedByteCompressionVariable<T>::adaptation_creator_ = CreatePrefixEmbeddedExtractionAdaptationClass<T>;
        AbstractEmbeddedByteCompressionVariable<T>::adaptation_destructor_ = DestroyPrefixEmbeddedExtractionAdaptationClass<T>;
        AbstractEmbeddedByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::EmbeddedMeshEncoder>();
    };

    CompressionSchema GetCompressionSchema() const override
    {
        return CompressionSchema::EmbeddedPrefixExtraction;
    }

private:

};

}

#endif /* !CMC_EMBEDDED_PREFIX_EXTRACTION_COMPRESSION_HXX */
