#ifndef CMC_MULTI_RES_EXTRACTION_COMPRESSION_HXX
#define CMC_MULTI_RES_EXTRACTION_COMPRESSION_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "lossless/cmc_byte_compression_variable.hxx"
#include "utilities/cmc_arithmetic_encoding.hxx"
#include "lossless/cmc_multi_res_extraction_residual_computation.hxx"
#include "utilities/cmc_interpolation_fn.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>

namespace cmc::lossless::multi_res
{

template<typename T>
class MultiResAdaptData : public ICompressionAdaptData<T>
{
public:
    MultiResAdaptData() = delete;
    MultiResAdaptData(AbstractByteCompressionVariable<T>* variable)
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
    bit_map::BitMap resdiual_order_indications_;
    int count_adaptation_step_{0};
};

template<typename T>
inline void
MultiResAdaptData<T>::IndicateCoarsening()
{
    refinement_indications_.AppendSetBit();
}

template<typename T>
inline void
MultiResAdaptData<T>::IndicateElementStaysUnchanged()
{
    refinement_indications_.AppendUnsetBit();
}

template <typename T>
void
MultiResAdaptData<T>::InitializeExtractionIteration()
{
    refinement_indications_ = bit_map::BitMap();
    resdiual_order_indications_ = bit_map::BitMap();
}

template <typename T>
void
MultiResAdaptData<T>::FinalizeExtractionIteration()
{
    ++count_adaptation_step_;
}

template <typename T>
void
MultiResAdaptData<T>::CompleteExtractionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    //Nothing to be done here!
}

template <typename T>
void
MultiResAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

template<typename T>
T
GetCoarseApproximationMaximizingResidualsLZC(const VectorView<CompressionValue<T>> values)
{
    cmc_assert(values.size() >= 1);

    T current_best_predictor;
    int current_max_lzc{-1};

    /* Try each value from the view as a predictor */
    for (auto pred_iter = values.begin(); pred_iter != values.end(); ++pred_iter)
    {
        /* Convert the compression value back to its origianl type */
        const T predictor = pred_iter->template ReinterpretDataAs<T>();

        /* Compute the LZC in the residuals for this approximation */
        const int pred_lzc = GetCumulativeResidualsLZC<T>(predictor, values);

        /* Potentially, update the current predictor */
        if (pred_lzc > current_max_lzc)
        {
            current_best_predictor = predictor;
            current_max_lzc = pred_lzc;
        }
    }

    /* Convert the view to actual values of the underlying data type */
    const std::vector<T> converted_vals = ConvertCompressionValues<T>(values);

    /* Try the mid-range as an predictor */
    const T mid_range = InterpolateToMidRange<T>(converted_vals);

    /* Compute the LZC in the residuals for this approximation */
    const int mid_range_lzc = GetCumulativeResidualsLZC<T>(mid_range, values);

    /* Potentially, update the current predictor */
    if (mid_range_lzc > current_max_lzc)
    {
        current_best_predictor = mid_range;
        current_max_lzc = mid_range_lzc;
    }

    /* Try the arithmetic mean as an predictor */
    const T mean = InterpolateToArithmeticMean<T>(converted_vals);

    /* Compute the LZC in the residuals for this approximation */
    const int mean_lzc = GetCumulativeResidualsLZC<T>(mean, values);

    /* Potentially, update the current predictor */
    if (mean_lzc > current_max_lzc)
    {
        current_best_predictor = mean;
        current_max_lzc = mean_lzc;
    }

    return current_best_predictor;
}

template <typename T>
ExtractionData<T>
MultiResAdaptData<T>::PerformExtraction([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, [[maybe_unused]] const int num_elements, const VectorView<CompressionValue<T>> values)
{
    /* Since we perform an extraction, a family of elements is coarsened */
    IndicateCoarsening();

    /* Get the coarse approximation for these values which maximizes the cumulative residual LZC */
    const T coarse_approximation = GetCoarseApproximationMaximizingResidualsLZC<T>(values);

    std::vector<CompressionValue<T>> fine_values;
    fine_values.reserve(values.size());

    /* Compute the residuals for the given predictor */
    for (auto val_iter = values.begin(); val_iter != values.end(); ++val_iter)
    {
        /* Compute the residual between approximation and actual value */
        auto [is_approximation_greater, residual] = ComputeResidual<T>(coarse_approximation, *val_iter);

        /* Store whether the approximation is greater than the real value */
        resdiual_order_indications_.AppendBit(is_approximation_greater);

        /* Store the residual */
        fine_values.push_back(residual);
    }

    return ExtractionData<T>(CompressionValue<T>(coarse_approximation), std::move(fine_values));
}

template <typename T>
UnchangedData<T>
MultiResAdaptData<T>::ElementStaysUnchanged([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, const CompressionValue<T>& value)
{
    /* In case the element stays unchanged, we indicate that no coarsening is possible */
    IndicateElementStaysUnchanged();

    /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
    resdiual_order_indications_.AppendUnsetBit();

    /* We drag this value along to the "coarser level" until it this element is passed with its siblings
     * as a family of elements into the adaptation callback. Moreover, we leave an empty value 
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

    /* Append the bytes to the byte stream */
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
MultiResAdaptData<T>::EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    //TODO: Rework for parallel usage

    cmc_debug_msg("The encoding of the CompressionValues after the prefix extraction iteration starts...");
    
    cmc_assert(ICompressionAdaptData<T>::entropy_coder_ != nullptr);

    /* The encoded data will be stored in a BitVector */
    cmc::bit_vector::BitVector encoding;
    encoding.Reserve(3 * level_byte_values.size());

    /* Reset the entropy coder and initialize the alphabet */
    ICompressionAdaptData<T>::entropy_coder_->Reset();
    ICompressionAdaptData<T>::entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Get a view on the bitmap storing the residual addition/subtraction flags */
    bit_map::BitMapView residual_flags(resdiual_order_indications_);

    /* Collect the symbol frequencies */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /* Get the current residual */
        CompressionValue<T> val = *val_iter;

        /* Get the encoded LZC */
        const uint32_t signum = (residual_flags.GetNextBit() == true ? 0 : cmc::entropy_coding::arithmetic_coding::kMSBBit);
        const uint32_t first_one_bit = val.GetNumberLeadingZeros();

        /* Update this symbol for encoding */
        ICompressionAdaptData<T>::entropy_coder_->UpdateSymbolFrequency(signum + first_one_bit);
    }

    /* Setup the interior structure for encoding */
    ICompressionAdaptData<T>::entropy_coder_->SetupEncoding();

    /* Reset the residual flag view */
    residual_flags = bit_map::BitMapView(resdiual_order_indications_);
    
    /* Iterate over all values and encode them */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /* Get the current value */
        CompressionValue<T> val = *val_iter;

        /* Get the encoded LZC */
        const uint32_t signum = (residual_flags.GetNextBit() == true ? 0 : cmc::entropy_coding::arithmetic_coding::kMSBBit);
        const uint32_t first_one_bit = val.GetNumberLeadingZeros();

        /* Encode the LZC and the residual fflag together */
        ICompressionAdaptData<T>::entropy_coder_->EncodeSymbol(signum + first_one_bit);

        /* The first one is implicitly given by the leading zero count; therefore we set the "front bit" in order to discard the leading zeros and the following one */
        val.SetFrontBit(first_one_bit + 1);

        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            encoding.AppendBits(val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
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
    return new MultiResAdaptData<T>(abstract_var);
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

#endif /* !CMC_MULTI_RES_EXTRACTION_COMPRESSION_HXX */
