#ifndef CMC_EMBEDDED_PREFIX_QUANTIZATION_COMPRESSION_HXX
#define CMC_EMBEDDED_PREFIX_QUANTIZATION_COMPRESSION_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "lossy/cmc_embedded_byte_compression_variable.hxx"
#include "utilities/cmc_prefix_entropy_coder.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_interpolation_fn.hxx"
#include "mesh_compression/cmc_embedded_mesh_encoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>

namespace cmc::lossy::embedded::prefix
{


/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template<typename T>
class PrefixEmbeddedAdaptData : public IEmbeddedCompressionAdaptData<T>
{
public:
    PrefixEmbeddedAdaptData() = delete;
    PrefixEmbeddedAdaptData(AbstractEmbeddedByteCompressionVariable<T>* variable,  const CompressionSettings& settings)
    : IEmbeddedCompressionAdaptData<T>(variable, settings) {
        IEmbeddedCompressionAdaptData<T>::entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::PrefixEncoder<T>>();
    };

    void InitializeExtractionIteration() override;
    void FinalizeExtractionIteration() override;
    void CompleteExtractionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const override;
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const override;

protected:
    ExtractionData<T> PerformExtraction(t8_forest_t forest, int tree_id, int lelement_id, const t8_scheme_c* ts, const int num_elements,
        const t8_element_t* elements[], const VectorView<CompressionValue<T>> values) override;
    UnchangedData<T> ElementStaysUnchanged(t8_forest_t forest, int tree_id, int lelement_id, const t8_scheme_c* ts, const int num_elements,
        const t8_element_t* elements[], const CompressionValue<T>& value) override;
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
PrefixEmbeddedAdaptData<T>::PerformExtraction(t8_forest_t forest, int tree_id, int lelement_id, const t8_scheme_c* ts, const int num_elements,
    const t8_element_t* elements[], const VectorView<CompressionValue<T>> values)
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
        return ExtractionData<T>(prefix, 0.0, std::move(fine_values));
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
        return ExtractionData<T>(CompressionValue<T>(), 0.0, std::move(fine_values));
    }
}

template <typename T>
UnchangedData<T>
PrefixEmbeddedAdaptData<T>::ElementStaysUnchanged(t8_forest_t forest, int tree_id, int lelement_id, const t8_scheme_c* ts, const int num_elements,
    const t8_element_t* elements[], const CompressionValue<T>& value)
{
    /* We drag this value along to the "coarser level" until it this element is passed with its siblings
     *  as a family of elements into the adaptation callback. Moreover, we leave an empty value 
     * on the "finer level" */
    return UnchangedData<T>(value, 0.0, CompressionValue<T>());
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
CreatePrefixEmbeddedExtractionAdaptationClass(AbstractEmbeddedByteCompressionVariable<T>* abstract_var, const CompressionSettings& settings)
{
    return new PrefixEmbeddedAdaptData<T>(abstract_var, settings);
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

    EmbeddedCompressionVariable(const CompressionSettings& settings, input::Var& input_variable)
    : AbstractEmbeddedByteCompressionVariable<T>(settings, input_variable)
    {
        //cmc_assert(input_variable.IsValid());//TODO: implement 

        this->SetName(input_variable.GetName());
        this->IndicateWhetherMeshRefinementBitsWillBeStored(true);
        this->SetInaccuracyTracking(false);

        AbstractEmbeddedByteCompressionVariable<T>::adaptation_creator_ = CreatePrefixEmbeddedExtractionAdaptationClass<T>;
        AbstractEmbeddedByteCompressionVariable<T>::adaptation_destructor_ = DestroyPrefixEmbeddedExtractionAdaptationClass<T>;
        AbstractEmbeddedByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::EmbeddedMeshEncoder>();
    };

    CompressionSchema GetCompressionSchema() const override
    {
        return CompressionSchema::EmbeddedQuantizedPrefixExtraction;
    }

    void PreCompressionProcessing(std::vector<CompressionValue<T>>& initial_data) override;

private:
    void PerformTailTruncation(CompressionValue<T>& initial_byte_value, const std::vector<PermittedError>& permitted_errors);
};


struct ErrorCompliance
{
    ErrorCompliance() = delete;
    ErrorCompliance(const bool is_error_threshold_fulfilled, const double max_error)
    : is_error_threshold_satisfied{is_error_threshold_fulfilled}, max_introduced_error{max_error}{};

    const bool is_error_threshold_satisfied;
    const double max_introduced_error;
};

template<typename T>
auto
ComputeSingleRelativeDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    /* Calculate the relative deviation */
    return static_cast<double>(std::abs(initial_value - nominal_value)) / static_cast<double>(std::abs(initial_value));
}

template<typename T>
auto
ComputeSingleRelativeDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    /* Calculate the relative deviation */
    return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value)) / static_cast<double>(initial_value);
}

template<typename T>
auto
ComputeSingleAbsoluteDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    /* Calculate the absolute deviation */
    return static_cast<double>(std::abs(initial_value - nominal_value));
}

template<typename T>
auto
ComputeSingleAbsoluteDeviation(const T initial_value, const T nominal_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    /* Calculate the absolute deviation */
    return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value));
}


template<typename T>
auto
ComputeSingleAbsoluteDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the absolute deviation */
        return static_cast<double>(std::abs(initial_value - nominal_value));
    } else
    {
        return 0.0;
    }
}

template<typename T>
auto
ComputeSingleAbsoluteDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the absolute deviation */
        return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value));
    } else
    {
        return 0.0;
    }
}

template<typename T>
auto
ComputeSingleRelativeDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_signed_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the relative deviation */
        return static_cast<double>(std::abs(initial_value - nominal_value)) / static_cast<double>(std::abs(initial_value));
    } else
    {
        return 0.0;
    }
}

template<typename T>
auto
ComputeSingleRelativeDeviationSkipMissingValues(const T initial_value, const T nominal_value, const T& missing_value)
 -> std::enable_if_t<std::is_unsigned_v<T>, double>
{
    if (!ApproxCompare(missing_value, initial_value))
    {
        /* Calculate the relative deviation */
        return static_cast<double>((initial_value > nominal_value ? initial_value - nominal_value : nominal_value - initial_value)) / static_cast<double>(initial_value);
    } else
    {
        return 0.0;
    }
}

template<class T>
ErrorCompliance
IsValueErrorCompliant(const std::vector<PermittedError>& permitted_errors, const T initial_value, const T nominal_value, const T& missing_value)
{
    /* Get the absolute inaccuracy for all values (the absolute error is needed in any case) */
    const double abs_inaccuracy = ComputeSingleAbsoluteDeviationSkipMissingValues(initial_value, nominal_value, missing_value);

    for (auto pe_iter = permitted_errors.begin(); pe_iter != permitted_errors.end(); ++pe_iter)
    {
        switch (pe_iter->criterion)
        {
            case CompressionCriterion::AbsoluteErrorThreshold:
            {
                /* Check if it is compliant with the permitted error */
                if (abs_inaccuracy > pe_iter->error)
                {
                    return ErrorCompliance(false, 0.0);
                }
            }
            break;
            case CompressionCriterion::RelativeErrorThreshold:
            {
                /* Get the relative inaccuracy for all values */
                const double rel_inaccuracy = ComputeSingleRelativeDeviationSkipMissingValues(initial_value, nominal_value, missing_value);
                /* Check if it is compliant with the permitted error */
                if (rel_inaccuracy > pe_iter->error)
                {
                    return ErrorCompliance(false, 0.0);
                }
            }
            break;
                default:
                cmc_err_msg("The error specifications hold an unrecognized criterion.");
                return ErrorCompliance(false, 0.0);
        }
    }

    /* If the funciton reaches this point, the interpolated value complies with the permitted errors */
    return ErrorCompliance(true, abs_inaccuracy);
}

template<typename T>
CompressionValue<T>
GetMaximumTailToggledValue(const CompressionValue<T>& initial_serialized_value, const std::vector<PermittedError>& permitted_errors, const T& missing_value)
{
    bool is_toogling_progressing = true;
    CompressionValue<T> toggled_value = initial_serialized_value;

    const T initial_val = initial_serialized_value.template ReinterpretDataAs<T>();

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_toogling_progressing && iteration_count < max_iteration_count)
    {
        const CompressionValue<T> save_previous_value = toggled_value;

        /* Toggle all ones up until the next unset bit (inclusive) */
        toggled_value.ToggleTailUntilNextUnsetBit();
        const T reinterpreted_value = toggled_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = IsValueErrorCompliant(permitted_errors, initial_val, reinterpreted_value, missing_value);

        if (!error_evaluation.is_error_threshold_satisfied)
        {
            /* Revert the changes to the value */
            toggled_value = save_previous_value;
            is_toogling_progressing = false;
        }

        ++iteration_count;
    }

    return toggled_value;
}


template<typename T>
CompressionValue<T>
GetMaximumTailClearedValue(const CompressionValue<T>& initial_serialized_value, const std::vector<PermittedError>& permitted_errors, const T& missing_value)
{
    bool is_clearing_progressing = true;
    CompressionValue<T> cleared_value = initial_serialized_value;

    const T initial_val = initial_serialized_value.template ReinterpretDataAs<T>();

    int iteration_count = 0;
    const int max_iteration_count = sizeof(T) * CHAR_BIT;

    while (is_clearing_progressing && iteration_count < max_iteration_count)
    {
        const CompressionValue<T> save_previous_value = cleared_value;

        /* Clear the next set bit from the tail */
        cleared_value.ClearNextSetBitFromTail();
        const T reinterpreted_value = cleared_value.template ReinterpretDataAs<T>();

        /* Check if it is error compliant */
        const ErrorCompliance error_evaluation = IsValueErrorCompliant(permitted_errors, initial_val, reinterpreted_value, missing_value);

        if (!error_evaluation.is_error_threshold_satisfied)
        {
            /* Revert the changes to the value */
            cleared_value = save_previous_value;
            is_clearing_progressing = false;
        }

        ++iteration_count;
    }

    return cleared_value;
}

static size_t bit_removal_trunc = 0;
template<typename T>
void
EmbeddedCompressionVariable<T>::PerformTailTruncation(CompressionValue<T>& initial_byte_value, const std::vector<PermittedError>& permitted_errors)
{
    const T missing_value = this->GetVariableAttributes().GetMissingValue();

    const T reinterpreted_val = initial_byte_value.template ReinterpretDataAs<T>();

    if (!ApproxCompare(reinterpreted_val, missing_value))
    {
        /* Get the value which has been transformed by toggling as many ones from the back while setting the succeeding 'zero' bits to one */
        const CompressionValue<T> toggled_value = GetMaximumTailToggledValue(initial_byte_value, permitted_errors, missing_value);

        /* Get the value which has been transformed by clearing as many of the last set bits as possible */
        const CompressionValue<T> cleared_value = GetMaximumTailClearedValue(initial_byte_value, permitted_errors, missing_value);

        /* Check which approach leads to more zero bits at the end */
        const int num_toogled_trailing_zeros = toggled_value.GetNumberTrailingZeros();
        const int num_cleared_trailing_zeros = cleared_value.GetNumberTrailingZeros();

        /* Replace the initial value with the transformed one */
        if (num_toogled_trailing_zeros >= num_cleared_trailing_zeros)
        {
            /* If the toggling approach has been more successfull */
            initial_byte_value = toggled_value;
            bit_removal_trunc += num_toogled_trailing_zeros;
            //initial_byte_value.SetTailBit(num_toogled_trailing_zeros);
        } else
        {
            /* If the clearing approach has been more successfull */
            initial_byte_value = cleared_value;
            bit_removal_trunc += num_cleared_trailing_zeros;
            //initial_byte_value.SetTailBit(num_cleared_trailing_zeros);
        }

        /* Update the trail bit count for the new value */
        initial_byte_value.UpdateTailBitCount();
    } else
    {
        /* In order to not change missing values, we are just able to trim their trailing zeros */
        initial_byte_value.UpdateTailBitCount();
    }

    //Test implicit one before TZC 
    const int tail = initial_byte_value.GetTailBit();
    if (tail < sizeof(T) * CHAR_BIT)
        initial_byte_value.SetTailBit(tail + 1);
}

//TODO: Since we remove the trailing zeros, we can remove the implciit one at the end of the bit-sequence
//!!!!!!
/* Perform the tail truncation as a pre-processing step */
template<class T>
void
EmbeddedCompressionVariable<T>::PreCompressionProcessing(std::vector<CompressionValue<T>>& initial_data) 
{
    t8_forest_t mesh = this->GetAmrMesh().GetMesh();
    const t8_scheme_c* scheme =  t8_forest_get_scheme(mesh);

    const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees(mesh);
    int val_idx = 0;

    /* Iterate over all elements in all trees */
    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx)
    {
        const t8_eclass_t tree_class = t8_forest_get_tree_class (mesh, tree_idx);
        const t8_locidx_t  num_elements_in_tree = t8_forest_get_tree_num_elements (mesh, tree_idx);
        for (t8_locidx_t elem_idx = 0; elem_idx < num_elements_in_tree; ++elem_idx, ++val_idx)
        {
            /* Get the current element */
            const t8_element_t* element = t8_forest_get_element_in_tree (mesh, tree_idx, elem_idx);

            /* Get the permitted error for this element */
            std::vector<PermittedError> permitted_errors = this->GetRestrictingErrors(mesh, tree_idx, elem_idx, scheme, 1, &element);

            /* Perform the tail truncation */
            this->PerformTailTruncation(initial_data[val_idx], permitted_errors);
        }
    }

    cmc_debug_msg("Bit removal tail truncation: ", bit_removal_trunc, " in bytes: ", bit_removal_trunc / 8);
}


}



#endif /* !CMC_EMBEDDED_PREFIX_QUANTIZATION_COMPRESSION_HXX */
