#ifndef CMC_EMBEDDED_PREFIX_QUANTIZATION_COMPRESSION_HXX
#define CMC_EMBEDDED_PREFIX_QUANTIZATION_COMPRESSION_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "embedded/lossy/cmc_embedded_byte_compression_variable.hxx"
#include "utilities/cmc_prefix_entropy_coder.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_interpolation_fn.hxx"
#include "mesh_compression/cmc_embedded_mesh_encoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>

#define kPreferPrefixTruncationOverPrefixElongation 1

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
    ExtractionData<T> PerformExtraction(t8_forest_t forest, int tree_id, const t8_eclass_t tree_class, int lelement_id, const t8_scheme_c* ts, const int num_elements,
        const t8_element_t* elements[], const VectorView<CompressionValue<T>> values) override;
    UnchangedData<T> ElementStaysUnchanged(t8_forest_t forest, int tree_id, const t8_eclass_t tree_class, int lelement_id, const t8_scheme_c* ts, const int num_elements,
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
PrefixEmbeddedAdaptData<T>::PerformExtraction(t8_forest_t forest, int tree_id, [[maybe_unused]] const t8_eclass_t tree_class, int lelement_id, const t8_scheme_c* ts, const int num_elements,
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
PrefixEmbeddedAdaptData<T>::ElementStaysUnchanged([[maybe_unused]] t8_forest_t forest, [[maybe_unused]] int tree_id, [[maybe_unused]] const t8_eclass_t tree_class, [[maybe_unused]] int lelement_id, [[maybe_unused]] const t8_scheme_c* ts, [[maybe_unused]] const int num_elements,
    [[maybe_unused]] const t8_element_t* elements[], const CompressionValue<T>& value)
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
    cmc_debug_msg("Num bytes for prefix length encoding: ", local_encoded_prefix_length_stream_num_bytes);
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

    /* For simplified testing, the Perform TailTruncation function has been made public. Revert when testing stage is done. */
    void PerformTailTruncation(CompressionValue<T>& initial_byte_value, const std::vector<PermittedError>& permitted_errors);
    void PerformTailTruncationOnward(CompressionValue<T>& initial_byte_value, CompressionValue<T>& current_value, const std::vector<PermittedError>& permitted_errors, const int max_num_trunc_iter);
    void PerformTailTruncationOnFamily(std::vector<CompressionValue<T>>& initial_byte_values, const int start_index, const int num_elements, const std::vector<std::vector<PermittedError>>& permitted_errors);

private:
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
GetMaximumTailToggledValue(const CompressionValue<T>& initial_serialized_value, const std::vector<PermittedError>& permitted_errors, const T& missing_value, const int max_iteration_count = sizeof(T) * CHAR_BIT)
{
    bool is_toogling_progressing = true;
    CompressionValue<T> toggled_value = initial_serialized_value;

    const T initial_val = initial_serialized_value.template ReinterpretDataAs<T>();

    int iteration_count = 0;

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
GetMaximumTailClearedValue(const CompressionValue<T>& initial_serialized_value, const std::vector<PermittedError>& permitted_errors, const T& missing_value, const int max_iteration_count = sizeof(T) * CHAR_BIT)
{
    bool is_clearing_progressing = true;
    CompressionValue<T> cleared_value = initial_serialized_value;

    const T initial_val = initial_serialized_value.template ReinterpretDataAs<T>();

    int iteration_count = 0;

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

//TODO: Check for NaN and Infs during tail truncation
static size_t bit_removal_trunc = 0;
template<typename T>
void
EmbeddedCompressionVariable<T>::PerformTailTruncation(CompressionValue<T>& initial_byte_value, const std::vector<PermittedError>& permitted_errors)
{
    #if 1
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
    #endif
}

#if kPreferPrefixTruncationOverPrefixElongation

#if 1

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
        const t8_locidx_t  num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (mesh, tree_idx);
        for (t8_locidx_t elem_idx = 0; elem_idx < num_elements_in_tree; ++elem_idx, ++val_idx)
        {
            /* Get the current element */
            const t8_element_t* element = t8_forest_get_leaf_element_in_tree (mesh, tree_idx, elem_idx);

            /* Get the permitted error for this element */
            std::vector<PermittedError> permitted_errors = this->GetRestrictingErrors(mesh, tree_idx, tree_class, elem_idx, scheme, 1, &element);

            /* Perform the tail truncation */
            this->PerformTailTruncation(initial_data[val_idx], permitted_errors);
        }
    }

    cmc_debug_msg("Bit removal tail truncation: ", bit_removal_trunc, " in bytes: ", bit_removal_trunc / 8);
}
#endif


#else

template<typename T>
CompressionValue<T>
GetMaximumTailToggledValue(const CompressionValue<T>& initial_serialized_value, const CompressionValue<T>& current_value, const std::vector<PermittedError>& permitted_errors, const T& missing_value, const int max_iteration_count = sizeof(T) * CHAR_BIT)
{
    bool is_toogling_progressing = true;
    CompressionValue<T> toggled_value = current_value;

    const T initial_val = initial_serialized_value.template ReinterpretDataAs<T>();

    int iteration_count = 0;

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
GetMaximumTailClearedValue(const CompressionValue<T>& initial_serialized_value, const CompressionValue<T>& current_value, const std::vector<PermittedError>& permitted_errors, const T& missing_value, const int max_iteration_count = sizeof(T) * CHAR_BIT)
{
    bool is_clearing_progressing = true;
    CompressionValue<T> cleared_value = current_value;

    const T initial_val = initial_serialized_value.template ReinterpretDataAs<T>();

    int iteration_count = 0;

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


template<typename T>
void
EmbeddedCompressionVariable<T>::PerformTailTruncationOnward(CompressionValue<T>& initial_byte_value, CompressionValue<T>& current_value, const std::vector<PermittedError>& permitted_errors, const int max_num_trunc_iter)
{
    #if 1
    const T missing_value = this->GetVariableAttributes().GetMissingValue();

    const T reinterpreted_val = initial_byte_value.template ReinterpretDataAs<T>();

    const T reinterpreted_current_val = current_value.template ReinterpretDataAs<T>();
     
    if (!ApproxCompare(reinterpreted_val, missing_value))
    {
        /* Get the value which has been transformed by toggling as many ones from the back while setting the succeeding 'zero' bits to one */
        const CompressionValue<T> toggled_value = GetMaximumTailToggledValue(initial_byte_value, reinterpreted_current_val, permitted_errors, missing_value, max_num_trunc_iter);

        /* Get the value which has been transformed by clearing as many of the last set bits as possible */
        const CompressionValue<T> cleared_value = GetMaximumTailClearedValue(initial_byte_value, reinterpreted_current_val, permitted_errors, missing_value, max_num_trunc_iter);

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
    #endif
}

template <typename T>
struct MaximizeCommonPrefixAdaptData
{
    MaximizeCommonPrefixAdaptData(EmbeddedCompressionVariable<T>* variable, std::vector<CompressionValue<T>>& initial_vals, const T missing_val)
    : var{variable}, compression_vals{initial_vals}, missing_value{missing_val} {};

    EmbeddedCompressionVariable<T>* var;
    std::vector<CompressionValue<T>>& compression_vals;
    const T missing_value;
};

template<typename T>
std::pair<bool, CompressionValue<T>>
EvaluateInitialCommonPrefix(const VectorView<CompressionValue<T>>& compression_values)
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

template<typename T>
std::pair<std::vector<CompressionValue<T>>, int>
AlterToMaximizePrefix(const VectorView<CompressionValue<T>>& values, const CompressionValue<T>& prefix, const std::vector<std::vector<PermittedError>>& permitted_errors, const T& missing_value)
{
    /* Copy the values */
    std::vector<CompressionValue<T>> data;
    data.reserve(values.size());
    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        data.push_back(*iter);
    }

    std::vector<CompressionValue<T>> data_save;
    bool is_altering_progressing = true;

    /* Get the prefix length and the "unequal bit position" */
    int pref_length = prefix.GetCountOfSignificantBits();
    int bit_pos = prefix.GetTypeNumBits() - pref_length - 1;

    while (is_altering_progressing && bit_pos > 0)
    {
        data_save = data;

        bit_pos = prefix.GetTypeNumBits() - pref_length - 1;
        cmc_assert(bit_pos >= 0);

        int count_one_bits{0};

        /* Count the bit difference */
        for (auto iter = data.begin(); iter != data.end(); ++iter)
        {
            if (iter->IsBitSetAtPosition(bit_pos))
            {
                ++count_one_bits;
            }
        }

        if (count_one_bits >= data.size() / 2)
        {
            /* More one bits are present */
            for (auto data_iter = data.begin(); data_iter != data.end(); ++data_iter)
            {
                if (not data_iter->IsBitSetAtPosition(bit_pos))
                {
                    /* In case the bit is not set, we alter the value */
                    data_iter->ToggleBitAtPositionAndResetPreviousBits(bit_pos);
                }
            }
        } else
        {
            /* More zero bits are present */
            for (auto data_iter = data.begin(); data_iter != data.end(); ++data_iter)
            {
                if (data_iter->IsBitSetAtPosition(bit_pos))
                {
                    /* In case the bit is set, we alter the value */
                    data_iter->ToggleBitAtPositionAndResetPreviousBits(bit_pos);
                }
            }
        }

        bool is_error_compliant = true;
        /* Check if the altered values are error-compliant */
        int idx = 0;

        for (auto data_iter = data.begin(); data_iter != data.end(); ++data_iter, ++idx)
        {
            const ErrorCompliance error_eval = IsValueErrorCompliant<T>(permitted_errors[idx], values[idx].template ReinterpretDataAs<T>(), data_iter->template ReinterpretDataAs<T>(), missing_value);

            if (not error_eval.is_error_threshold_satisfied)
            {
                is_error_compliant = false;
                break;
            }
        }

        if (is_error_compliant && pref_length < prefix.GetTypeNumBits())
        {
            ++pref_length;
            data_save = data;

            //cmc_debug_msg("Successfull prefix elongation");
        } else
        {
            is_altering_progressing = false;
            data = data_save;
        }

    }

    return std::make_pair(data, bit_pos);
}


template<typename T>
inline t8_locidx_t
MaximizeCommonPrefix (t8_forest_t forest,
                      t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         [[maybe_unused]] const t8_eclass_t tree_class,
                         t8_locidx_t lelement_id,
                         [[maybe_unused]] const t8_scheme_c * ts,
                         const int is_family,
                         const int num_elements,
                         [[maybe_unused]] t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    MaximizeCommonPrefixAdaptData<T>* adapt_data = static_cast<MaximizeCommonPrefixAdaptData<T>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    const int start_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

    VectorView<CompressionValue<T>> data(&adapt_data->compression_vals[start_index], num_elements);

    /* Collect the permitted errrors  */
    std::vector<std::vector<PermittedError>> permitted_errors;
    permitted_errors.reserve(num_elements);

    for (int elem_id = 0; elem_id < num_elements; ++elem_id)
    {
        const t8_element_t* elem = elements[elem_id];
        permitted_errors.push_back(adapt_data->var->GetRestrictingErrors(forest_from, which_tree, lelement_id + elem_id, ts, 1, &elem));
    }


    // If it is a family we try to maximize the common prefix
    if (is_family)
    {
        auto [is_prefix_present, prefix] = EvaluateInitialCommonPrefix<T>(data);

        /* If there is a common prefix, we try to maximize it, if there is no common prefix, we just perform the tail truncation */
        if (is_prefix_present)
        {
            if (prefix.GetCountOfSignificantBits() == prefix.GetTypeNumBits())
            {
                /* If the prefix is full, nothing has to be done */
                return -1;
            }

            const int pref_length = prefix.GetCountOfSignificantBits();
            const int bit_pos = prefix.GetTypeNumBits() - pref_length - 1;

            cmc_assert(bit_pos >= 0);

            /* Maximize the prefixes */
            auto [altered_values, new_bit_pos] = AlterToMaximizePrefix<T>(data, prefix, permitted_errors, adapt_data->missing_value);

            /* Perform tail truncation */
            int alt_idx = 0;
            for (auto alt_data_iter = altered_values.begin(); alt_data_iter != altered_values.end(); ++alt_data_iter, ++alt_idx)
            {
                adapt_data->var->PerformTailTruncationOnward(adapt_data->compression_vals[start_index + alt_idx], *alt_data_iter, permitted_errors[alt_idx], new_bit_pos);
            }

        } else
        {
            //Perform tailt runcaiton on all family members here 
            int d_idx = 0;
            
            for (auto data_iter = data.begin(); data_iter != data.end(); ++data_iter, ++d_idx)
            {
                adapt_data->var->PerformTailTruncation(adapt_data->compression_vals[start_index + d_idx], permitted_errors[d_idx]);
            }

        }
        return -1;
    } else
    {
        // If it is no family, we just perform the tail truncation 
        //In the adaptive case, the element values needs to be checked for prefix maximizing the first time they would be coarsened and if it is not applicable tail truncation can be used

        adapt_data->var->PerformTailTruncation(adapt_data->compression_vals[start_index], permitted_errors[0]);
    }

    return 0;
}

/* Try to maximize common prefixes */
template<class T>
void
EmbeddedCompressionVariable<T>::PreCompressionProcessing(std::vector<CompressionValue<T>>& initial_data) 
{
    t8_forest_t mesh = this->GetAmrMesh().GetMesh();
    const T missing_value = this->GetVariableAttributes().GetMissingValue();

    t8_forest_ref(mesh);
    MaximizeCommonPrefixAdaptData<T> adapt_data(this, initial_data, missing_value);

    t8_forest_t mesh2 = t8_forest_new_adapt(mesh, MaximizeCommonPrefix<T>, 0, 0, &adapt_data);
    t8_forest_unref(&mesh2);
}



#endif



#if 0

template <typename T>
struct BlockTruncPrefixAdaptData
{
    BlockTruncPrefixAdaptData(EmbeddedCompressionVariable<T>* variable, std::vector<CompressionValue<T>>& initial_vals, const T missing_val)
    : var{variable}, compression_vals{initial_vals}, missing_value{missing_val} {};

    EmbeddedCompressionVariable<T>* var;
    std::vector<CompressionValue<T>>& compression_vals;
    const T missing_value;
};


template<typename T>
void
EmbeddedCompressionVariable<T>::PerformTailTruncationOnFamily(std::vector<CompressionValue<T>>& initial_byte_values, const int start_index, const int num_elements, const std::vector<std::vector<PermittedError>>& permitted_errors)
{
    #if 1
    const T missing_value = this->GetVariableAttributes().GetMissingValue();

    int min_num_tzc = 1000;
    std::vector<CompressionValue<T>> truncated_vals;
    truncated_vals.reserve(num_elements);

    for (int idx = 0; idx < num_elements; ++idx)
    {
        CompressionValue<T> initial_byte_value = initial_byte_values[start_index + idx];

        /* Get the value which has been transformed by toggling as many ones from the back while setting the succeeding 'zero' bits to one */
        const CompressionValue<T> toggled_value = GetMaximumTailToggledValue(initial_byte_value, permitted_errors[idx], missing_value);
        /* Get the value which has been transformed by clearing as many of the last set bits as possible */
        const CompressionValue<T> cleared_value = GetMaximumTailClearedValue(initial_byte_value, permitted_errors[idx], missing_value);

        /* Check which approach leads to more zero bits at the end */
        const int num_toogled_trailing_zeros = toggled_value.GetNumberTrailingZeros();
        const int num_cleared_trailing_zeros = cleared_value.GetNumberTrailingZeros();

        /* Replace the initial value with the transformed one */
        if (num_toogled_trailing_zeros >= num_cleared_trailing_zeros)
        {
            /* If the toggling approach has been more successfull */
            truncated_vals.push_back(toggled_value);
            if (min_num_tzc > num_toogled_trailing_zeros)
            {
                min_num_tzc = num_toogled_trailing_zeros;
            }
        } else
        {
            /* If the clearing approach has been more successfull */
            truncated_vals.push_back(cleared_value);
            if (min_num_tzc > num_cleared_trailing_zeros)
            {
                min_num_tzc = num_cleared_trailing_zeros;
            }
        }
    }

    /* Set the tail bit for this family */
    for (auto iter = truncated_vals.begin(); iter != truncated_vals.end(); ++iter)
    {
        iter->SetTailBit(min_num_tzc);
    }

    /* Copy truncated values over */
    for (int idx = 0; idx < num_elements; ++idx)
    {
        initial_byte_values[start_index + idx] = truncated_vals[idx];
    }

    #endif
}


template<typename T>
inline t8_locidx_t
BlockTruncatePrefixes (t8_forest_t forest,
                      t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         [[maybe_unused]] const t8_eclass_t tree_class,
                         t8_locidx_t lelement_id,
                         [[maybe_unused]] const t8_scheme_c * ts,
                         const int is_family,
                         const int num_elements,
                         [[maybe_unused]] t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    BlockTruncPrefixAdaptData<T>* adapt_data = static_cast<BlockTruncPrefixAdaptData<T>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    const int start_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

    VectorView<CompressionValue<T>> data(&adapt_data->compression_vals[start_index], num_elements);

    /* Collect the permitted errrors  */
    std::vector<std::vector<PermittedError>> permitted_errors;
    permitted_errors.reserve(num_elements);

    for (int elem_id = 0; elem_id < num_elements; ++elem_id)
    {
        const t8_element_t* elem = elements[elem_id];
        permitted_errors.push_back(adapt_data->var->GetRestrictingErrors(forest_from, which_tree, lelement_id + elem_id, ts, 1, &elem));
    }


    // If it is a family we try to maximize the common prefix
    if (is_family)
    {
        adapt_data->var->PerformTailTruncationOnFamily(adapt_data->compression_vals, start_index, num_elements, permitted_errors);

        return -1;
    } else
    {
        // If it is no family, we just perform the tail truncation 
        //In the adaptive case, the element values needs to be checked for prefix maximizing the first time they would be coarsened and if it is not applicable tail truncation can be used

        adapt_data->var->PerformTailTruncation(adapt_data->compression_vals[start_index], permitted_errors[0]);
    }

    return 0;
}


/* Try to maximize common prefixes */
template<class T>
void
EmbeddedCompressionVariable<T>::PreCompressionProcessing(std::vector<CompressionValue<T>>& initial_data) 
{
    t8_forest_t mesh = this->GetAmrMesh().GetMesh();
    const T missing_value = this->GetVariableAttributes().GetMissingValue();

    t8_forest_ref(mesh);
    BlockTruncPrefixAdaptData<T> adapt_data(this, initial_data, missing_value);

    t8_forest_t mesh2 = t8_forest_new_adapt(mesh, BlockTruncatePrefixes<T>, 0, 0, &adapt_data);
    t8_forest_unref(&mesh2);
}

#endif


#if 0

template <typename T>
struct TreeQuantAdaptData
{
    TreeQuantAdaptData(EmbeddedCompressionVariable<T>* variable, std::vector<CompressionValue<T>>& initial_vals, const T missing_val)
    : var{variable}, compression_vals{initial_vals}, missing_value{missing_val} {};

    EmbeddedCompressionVariable<T>* var;
    std::vector<CompressionValue<T>>& compression_vals;
    const T missing_value;
};



static int el_counter = 0;

static size_t bit_counter = 0;

template <typename T>
T CalculateMidRange(const VectorView<CompressionValue<T>>& values)
{
    #if 1
    //Mid Range
    cmc_assert(!values.empty());
    cmc_assert(values.size() >= 2);

    T min = std::numeric_limits<T>::max();
    T max = std::numeric_limits<T>::min();

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        const T val = iter->template ReinterpretDataAs<T>();
        if (min > val)
        {
            min = val;
        }
        if (max < val)
        {
            max = val;
        }
    }

    return ((max / 2) + (min / 2));

    #else
    //Arithmetic Mean

    double sum = 0;

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        const T val = iter->template ReinterpretDataAs<T>();
        sum += static_cast<double>(val);
    }

    return static_cast<T>(sum / static_cast<double>(values.size()));

    #endif
}

template <typename T>
std::pair<int, int>
FindSFCEncoding(const T mid_val, const T init_val, const T init_res, const std::vector<PermittedError>& permitted_error)
{
    int lvl_count = 0;
    bool encoding_continues = true;
    T res_range = init_res;
    T quant_val = mid_val;

    int encoding{0};

    while (encoding_continues)
    {
        /* Check if error compliance is fullfilled */
        ErrorCompliance err_compl = IsValueErrorCompliant(permitted_error, init_val, quant_val, T(10000000000));

        if (err_compl.is_error_threshold_satisfied)
        {
            encoding_continues = false;
        } else
        {
            res_range = res_range / 2.0;

            if (quant_val >= init_val)
            {
                quant_val -= res_range;
                encoding <<= 1;
            } else
            {
                quant_val += res_range;
                encoding <<= 1;
                encoding |= int{0x00000001};
            }

            ++lvl_count;
        }
    }

    return std::make_pair(encoding, lvl_count);
}

template<typename T>
inline t8_locidx_t
TreeQuant (t8_forest_t forest,
                      t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         [[maybe_unused]] const t8_eclass_t tree_class,
                         t8_locidx_t lelement_id,
                         [[maybe_unused]] const t8_scheme_c * ts,
                         const int is_family,
                         const int num_elements,
                         [[maybe_unused]] t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    TreeQuantAdaptData<T>* adapt_data = static_cast<TreeQuantAdaptData<T>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    const int start_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

    /* Collect the permitted errrors  */
    std::vector<std::vector<PermittedError>> permitted_errors;
    permitted_errors.reserve(num_elements);

    for (int elem_id = 0; elem_id < num_elements; ++elem_id)
    {
        const t8_element_t* elem = elements[elem_id];
        permitted_errors.push_back(adapt_data->var->GetRestrictingErrors(forest_from, which_tree, lelement_id + elem_id, ts, 1, &elem));
    }

    if (is_family)
    {
        VectorView<CompressionValue<T>> data(&adapt_data->compression_vals[start_index], num_elements);

        /* Find the mid range */
        const T mid_range = CalculateMidRange<T>(data);

        if (el_counter >= 5000000 && el_counter < 5001000)
        {
            T max_res = T(0.0);

            for (auto iter = data.begin(); iter != data.end(); ++iter)
            {
                if (std::abs(mid_range - iter->template ReinterpretDataAs<T>()) > max_res)
                {
                    max_res = std::abs(mid_range - iter->template ReinterpretDataAs<T>());
                }
            }

            int pe_idx = 0;
            for (auto iter = data.begin(); iter != data.end(); ++iter, ++pe_idx)
            {
                cmc_debug_msg("Idx: ", el_counter, ", mid range: ", mid_range, ", init_val: ", iter->template ReinterpretDataAs<T>(), ", res: ", mid_range - iter->template ReinterpretDataAs<T>());
                auto [encoding, lvl_count] = FindSFCEncoding(mid_range, iter->template ReinterpretDataAs<T>(), max_res, permitted_errors[pe_idx]);
                cmc_debug_msg("Lvls for encoding: ", lvl_count, ", Encodiung: ", encoding, ", max_res was: ", max_res);
                cmc_debug_msg("\n");

            }
        }



        el_counter += num_elements;
        return -1;
    } else
    {
        ++el_counter;
        return 0;
    }
}

template<class T>
void
EmbeddedCompressionVariable<T>::PreCompressionProcessing(std::vector<CompressionValue<T>>& initial_data) 
{
    t8_forest_t mesh = this->GetAmrMesh().GetMesh();
    const T missing_value = this->GetVariableAttributes().GetMissingValue();

    t8_forest_ref(mesh);
    TreeQuantAdaptData<T> adapt_data(this, initial_data, missing_value);

    t8_forest_t mesh2 = t8_forest_new_adapt(mesh, TreeQuant<T>, 0, 0, &adapt_data);
    t8_forest_unref(&mesh2);

    cmc_err_msg("We stop here for testing purposes.");
}



#endif


}



#endif /* !CMC_EMBEDDED_PREFIX_QUANTIZATION_COMPRESSION_HXX */
