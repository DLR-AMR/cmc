#ifndef CMC_PREFIX_EXTRACTION_COMPRESSION_HXX
#define CMC_PREFIX_EXTRACTION_COMPRESSION_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "lossless/cmc_byte_compression_variable.hxx"
#include "utilities/cmc_arithmetic_encoding.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>

namespace cmc::lossless::prefix
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = CompressionValue<T>;

template<typename T>
class PrefixAdaptData : public ICompressionAdaptData<T>
{
public:
    PrefixAdaptData() = delete;
    PrefixAdaptData(AbstractByteCompressionVariable<T>* variable)
    : ICompressionAdaptData<T>(variable) {
        ICompressionAdaptData<T>::entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::Encoder>();
    };

    void InitializeExtractionIteration() override;
    void FinalizeExtractionIteration() override;
    void CompleteExtractionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::vector<uint8_t> EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const override;
    
protected:
    ExtractionData<T> PerformExtraction(const int which_tree, const int lelement_id, const int num_elements, const VectorView<CompressionValue<T>> values) override;
    UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) override;

private:
    void IndicateCoarsening();
    void IndicateElementStaysUnchanged();

    bit_map::BitMap refinement_indications_;
    int count_adaptation_step_{0};
};

template<typename T>
inline void
PrefixAdaptData<T>::IndicateCoarsening()
{
    refinement_indications_.AppendSetBit();
}

template<typename T>
inline void
PrefixAdaptData<T>::IndicateElementStaysUnchanged()
{
    refinement_indications_.AppendUnsetBit();
}

template <typename T>
void
PrefixAdaptData<T>::InitializeExtractionIteration()
{
    refinement_indications_ = bit_map::BitMap();
}

template <typename T>
void
PrefixAdaptData<T>::FinalizeExtractionIteration()
{
    ++count_adaptation_step_;
}

template <typename T>
void
PrefixAdaptData<T>::CompleteExtractionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    //Nothing to be done here!
}

template <typename T>
void
PrefixAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

template<typename T>
std::pair<bool, CompressionValue<T>>
EvaluateCommonPrefix(const VectorView<CompressionValue<T>>& compression_values)
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
PrefixAdaptData<T>::PerformExtraction([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, [[maybe_unused]] const int num_elements, const VectorView<CompressionValue<T>> values)
{
    /* Since we perform an extraction, a family of elements is coarsened */
    IndicateCoarsening();

    /* Evaluate whether there is a common prefix to extract */
    auto [is_prefix_found, prefix] = EvaluateCommonPrefix<T>(values);

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
PrefixAdaptData<T>::ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value)
{
    /* In case the element stays unchanged, we indicate that no coarsening is possible */
    IndicateElementStaysUnchanged();

    /* We drag this value along to the "coarser level" until it this element is passed with its siblings
     *  as a family of elements into the adaptation callback. Moreover, we leave an empty value 
     * on the "finer level" */
    return UnchangedData<T>(value, CompressionValue<T>());
}

/**
 * @brief Pushes a value byte-wise in Big-Endian ordering back to a vector
 * 
 * @tparam T The data type of the \a value
 * @param byte_stream The stream to which the value is byte-wise appended
 * @param value The value to be pushed back
 */
template <typename T>
void
PushBackValueToByteStream(std::vector<uint8_t>& byte_stream, const T& value)
{
    /* Serialize the value to a collection of bytes (in big endian order) */
    const std::array<uint8_t, sizeof(T)> serialized_value = SerializeValue(value, Endian::Big);

    /* Append the the bytes to the byte stream */
    std::copy_n(serialized_value.begin(), sizeof(T), std::back_insert_iterator(byte_stream));
}

/**
 * @brief  We use an arithmetic encoder to encode the position of the first "one-bit" in the compression value
 * 
 * @tparam T The original data type of the underlying data (e.g. float)
 * @param level_byte_values The remaining "fine compression" values after an extraction iteration
 * @return std::vector<uint8_t> The encoded data stream
 */
template <typename T>
std::vector<uint8_t>
PrefixAdaptData<T>::EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    //TODO: Rework for parallel usage

    cmc_debug_msg("The encoding of the CompressionValues after the prefix extraction iteration starts...");
    
    cmc_assert(ICompressionAdaptData<T>::entropy_coder_ != nullptr);

    /* The encoded data will be stored in a BitVector */
    cmc::bit_vector::BitVector encoding;
    encoding.Reserve(3 * level_byte_values.size());

    /* If we are on the last level*/
    /* Reset the entropy coder and initialize the alphabet */
    ICompressionAdaptData<T>::entropy_coder_->Reset();
    ICompressionAdaptData<T>::entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Collect the symbol frequencies */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /* Get the leading zero count in the suffix (which is equivalent to the position of the fist "one-bit") */
        const uint32_t lzc = static_cast<uint32_t>(val_iter->GetLeadingZeroCountInSignificantBits());

        /* Update the frequency */
        ICompressionAdaptData<T>::entropy_coder_->UpdateSymbolFrequency(lzc);
    }

    /* Setup the interior structure for encoding */
    ICompressionAdaptData<T>::entropy_coder_->SetupEncoding();

    /* Iterate over all values and encode them */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /* Get the current value */
        CompressionValue<T> current_val = *val_iter;

        /* Get the leading zero count in the suffix (which is equivalent to the position of the fist "one-bit") */
        const int lzc = current_val.GetLeadingZeroCountInSignificantBits();

        /* Encode the leading zero count which equals the first "one-bit" position */
        ICompressionAdaptData<T>::entropy_coder_->EncodeSymbol(static_cast<uint32_t>(lzc));

        /* The remaining bits in the value are encoded/stored */
        if (lzc == 0 && not current_val.IsEmpty())
        {
            /* In case there is no LZC the whol value need to be stored */
            encoding.AppendBits(current_val.GetSignificantBitsInBigEndianOrdering(), current_val.GetCountOfSignificantBits());
        } else
        {
            /* When its greater than zero, we can trim the zeros */
            current_val.SetFrontBit(current_val.GetFrontBit() + lzc);

            if (not current_val.IsEmpty())
            {
                encoding.AppendBits(current_val.GetSignificantBitsInBigEndianOrdering(), current_val.GetCountOfSignificantBits());
            }
        }

    }

    /* Indicate that the encoding has been finished and flush all pending encodings */
    ICompressionAdaptData<T>::entropy_coder_->FinishEncoding();

    /* Get the encoded alphabet */
    cmc::bit_vector::BitVector encoded_alphabet = ICompressionAdaptData<T>::entropy_coder_->EncodeAlphabet();

    /* Get the encoded LZC, respectivel first "one-bit" positions */
    cmc::bit_map::BitMap encoded_lzc_stream = ICompressionAdaptData<T>::entropy_coder_->GetEncodedBitStream();

    /* Collect the overall bytes for the encoding */
    const size_t overall_level_bytes = 5 * sizeof(size_t) + refinement_indications_.size_bytes() + encoded_alphabet.size() + encoded_lzc_stream.size_bytes() + encoding.size();

    /* Now everything is put together to a single stream */
    std::vector<uint8_t> encoded_stream;
    encoded_stream.reserve(overall_level_bytes);
    
    /* Push back the overall byte count for the level */
    PushBackValueToByteStream(encoded_stream, overall_level_bytes);

    /* Push back the refinement indications */
    const size_t refinement_indications_bytes = refinement_indications_.size_bytes();
    PushBackValueToByteStream(encoded_stream, refinement_indications_bytes);

    /* Push back the byte count for the encoded alphabet */
    const size_t encoded_alphabet_bytes = encoded_alphabet.size();
    PushBackValueToByteStream(encoded_stream, encoded_alphabet_bytes);

    /* Push back the byte count for the encoded first "one-bit" positions */
    const size_t encoded_lzc_stream_bytes = encoded_lzc_stream.size_bytes();
    PushBackValueToByteStream(encoded_stream, encoded_lzc_stream_bytes);

    /* Push back the byte count for the remaining bits of the prefixes */
    const size_t reamaining_value_bytes = encoding.size();
    PushBackValueToByteStream(encoded_stream, reamaining_value_bytes);

    /* Push back the refinement_indications, the encoded alphabet, the encoded LZC and the remaining bits in the given order */
    std::copy_n(refinement_indications_.begin_bytes(), refinement_indications_bytes, std::back_insert_iterator(encoded_stream));
    std::copy_n(encoded_alphabet.begin(), encoded_alphabet_bytes, std::back_insert_iterator(encoded_stream));
    std::copy_n(encoded_lzc_stream.begin_bytes(), encoded_lzc_stream_bytes, std::back_insert_iterator(encoded_stream));
    std::copy_n(encoding.begin(), reamaining_value_bytes, std::back_insert_iterator(encoded_stream));

    cmc_debug_msg("The entropy encoder of the prefix ectraction compression stored the CompressionValues of this iteration within ", overall_level_bytes, " bytes.");
    return encoded_stream;
}

template <typename T>
inline ICompressionAdaptData<T>*
CreatePrefixExtractionAdaptationClass(AbstractByteCompressionVariable<T>* abstract_var)
{
    return new PrefixAdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyPrefixExtractionAdaptationClass(ICompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class CompressionVariable : public AbstractByteCompressionVariable<T>
{
public:
    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, const std::vector<T>& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(variable_data);
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreatePrefixExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyPrefixExtractionAdaptationClass<T>;
    };

    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, const std::vector<CompressionValue<T>>& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(variable_data);
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreatePrefixExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyPrefixExtractionAdaptationClass<T>;
    };

    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, std::vector<CompressionValue<T>>&& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(std::move(variable_data));
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreatePrefixExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyPrefixExtractionAdaptationClass<T>;
    };
};



}



#endif /* !CMC_PREFIX_EXTRACTION_COMPRESSION_HXX */
