#ifndef LOSSY_CMC_MULTI_RES_TRIMMED_RESIDUALS_HXX
#define LOSSY_CMC_MULTI_RES_TRIMMED_RESIDUALS_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_interpolation_fn.hxx"
#include "utilities/cmc_serialization.hxx"
//#include "utilities/cmc_trimmed_multi_res_entropy_coder.hxx"
#include "utilities/cmc_multi_res_entropy_coder.hxx"
#include "utilities/cmc_multi_res_extraction_util.hxx"
#include "lossy/cmc_byte_compression_variable.hxx"
#include "utilities/cmc_lossy_multi_res_extraction_residual_computation.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>

namespace cmc::lossy::multi_res
{

template<typename T>
class MultiResAdaptData : public IByteCompressionAdaptData<T>
{
public:
    MultiResAdaptData() = delete;
    MultiResAdaptData(AbstractByteCompressionVariable<T>* variable, const CompressionSettings& settings)
    : IByteCompressionAdaptData<T>(variable, settings) {
        IByteCompressionAdaptData<T>::entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResEncoder<T>>();
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

    bit_map::BitMap resdiual_order_indications_;
    int count_adaptation_step_{0};
};

template <typename T>
void
MultiResAdaptData<T>::InitializeExtractionIteration()
{
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

    T current_best_predictor{T()};
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
MultiResAdaptData<T>::PerformExtraction(t8_forest_t forest, int tree_id, int lelement_id, const t8_scheme_c* ts, const int num_elements,
                                        const t8_element_t* elements[], const VectorView<CompressionValue<T>> values)
{
    cmc_assert(num_elements == values.size());

    /* Get the coarse approximation for these values which maximizes the cumulative residual LZC */
    const T coarse_approximation = GetCoarseApproximationMaximizingResidualsLZC<T>(values);

    std::vector<CompressionValue<T>> fine_values;
    fine_values.reserve(num_elements);

    double max_abs_inaccuracy{0.0};

    for (int idx = 0; idx < num_elements; ++idx)
    {
        /* Get the value in the correct type */
        T current_value = values[idx].template ReinterpretDataAs<T>();

        /* Get the restricting error for this element */
        std::vector<PermittedError> permitted_errors = this->GetRestrictingErrors(forest, tree_id, lelement_id, ts, 1, &elements[idx]);

        /* Get the previous absolute deviation */
        const double previous_absolute_error = this->GetPreviousAbsoluteError(tree_id, lelement_id);

        /* Compute the residual between approximation and actual value */
        auto [is_approximation_greater, residual, absolute_inaccuracy] = ComputeMinimalResidual<T>(coarse_approximation, values[idx], permitted_errors, previous_absolute_error);
        //cmc_debug_msg("Residual: LZC: ", residual.GetNumberLeadingZeros(), ", Length: ", residual.GetCountOfSignificantBits(), ", TZC: ", residual.GetNumberTrailingZeros());
        /* Store whether the approximation is greater than the real value */
        resdiual_order_indications_.AppendBit(is_approximation_greater);

        /* Store the residual */
        fine_values.push_back(residual);

        /* Update the maximal inaccuracy for this family */
        if (max_abs_inaccuracy < absolute_inaccuracy)
        {
            max_abs_inaccuracy = absolute_inaccuracy;
        }
    }

    return ExtractionData<T>(CompressionValue<T>(coarse_approximation), max_abs_inaccuracy, std::move(fine_values));
}

template <typename T>
UnchangedData<T>
MultiResAdaptData<T>::ElementStaysUnchanged([[maybe_unused]] t8_forest_t forest, int tree_id, int lelement_id, [[maybe_unused]] const t8_scheme_c* ts, [[maybe_unused]] const int num_elements,
                                            [[maybe_unused]] const t8_element_t* elements[], const CompressionValue<T>& value)
{
    /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
    resdiual_order_indications_.AppendUnsetBit();

    /* Get the previous absolute deviation */
    const double previous_absolute_error = this->GetPreviousAbsoluteError(tree_id, lelement_id);

    /* We drag this value along to the "coarser level" until it this element is passed with its siblings
     * as a family of elements into the adaptation callback. Moreover, we leave an empty value 
     * on the "finer level" */
    return UnchangedData<T>(value, previous_absolute_error, CompressionValue<T>());
}

template <typename T>
void
MultiResAdaptData<T>::CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    /* Get a view on the bitmap storing the residual addition/subtraction flags */
    bit_map::BitMapView residual_flags(resdiual_order_indications_);

    /* Collect the symbol frequencies */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /** First, we collect the signum with the LZC as a symbol **/

        /* Get the current residual */
        CompressionValue<T> val = *val_iter;

        /* Get the encoded LZC */
        uint32_t signum = cmc::lossless::multi_res::util::GetSignumForEncoding(residual_flags.GetNextBit());
        const uint32_t first_one_bit = val.GetNumberLeadingZeros();
        if (val.GetNumberLeadingZeros() == sizeof(T) * bit_map::kCharBit)
        {
            signum = 0;
        }

        /* Update this symbol for encoding */
        IByteCompressionAdaptData<T>::entropy_coder_->UpdateSymbolFrequency(signum + first_one_bit);


        /** Afterwards, we collect the trailing zero count (TZC) **/
        /* The collection of the TZC is only done, when the LZC does not fully describe the residual */
        if (first_one_bit >= sizeof(T) * bit_map::kCharBit - 1)
        {
            continue;
        }

        /* Otherwise, compute the TZC */
        const uint32_t last_one_bit = val.GetNumberTrailingZeros();

        /* Update this symbol for encoding */
        IByteCompressionAdaptData<T>::entropy_coder_->UpdateSymbolFrequency(last_one_bit);       
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
MultiResAdaptData<T>::EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    cmc_debug_msg("The encoding of the CompressionValues after the multi-resolution extraction iteration starts...");
    
    cmc_assert(IByteCompressionAdaptData<T>::entropy_coder_ != nullptr);

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
    IByteCompressionAdaptData<T>::entropy_coder_->Reset(std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResCompressionAlphabet<T>>());
    IByteCompressionAdaptData<T>::entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Collect the symbols and their frequencies for the entropy coder */
    CollectSymbolFrequenciesForEntropyCoding(level_byte_values);

    /* Setup the interior structure for encoding */
    IByteCompressionAdaptData<T>::entropy_coder_->SetupEncoding(this->GetMPIComm());

    /* Get a view on the residual flags */
    bit_map::BitMapView residual_flags(resdiual_order_indications_);
    
    /* Iterate over all values and encode them */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /* Get the current value */
        CompressionValue<T> val = *val_iter;

        /* Get the encoded LZC */
        uint32_t signum = cmc::lossless::multi_res::util::GetSignumForEncoding(residual_flags.GetNextBit());
        const uint32_t first_one_bit = val.GetNumberLeadingZeros();
        if (val.GetNumberLeadingZeros() == sizeof(T) * bit_map::kCharBit)
        {
            signum = 0;
        }
        
        /* Encode the LZC and the residual flag together */
        IByteCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(signum + first_one_bit);

        /* The first one is implicitly given by the leading zero count; therefore we set the "front bit" in order to discard the leading zeros and the following one */
        val.SetFrontBit(first_one_bit + 1);

        /* Pontentially, encode the trailing zero count */
        if (first_one_bit >= sizeof(T) * bit_map::kCharBit - 1)
        {
            continue;
        }

        /* Otherwise, compute the TZC */
        const uint32_t last_one_bit = val.GetNumberTrailingZeros();

        /* Check, if we can make use of the last implicit one */
        if (first_one_bit + 1 + last_one_bit < sizeof(T) * bit_map::kCharBit)
        {
            /* In this case, the residual has at least two significant one bits, and therefore, we could make use of the 
             * implicit one bit at the end of the residual implied by the TZC */
            val.SetTailBit(last_one_bit + 1);
        } else
        {
            /* Otherwise, we could not make use of the implicit one bit by the TZC */
            val.SetTailBit(last_one_bit);
        }

        /* Encode the TZC */
        IByteCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(last_one_bit);

        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            encoding.AppendBits(val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
        }
    }

    /* We set an indicaton symbol that the process local end of the values have been reached */
    IByteCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte);

    /* Indicate that the encoding has been finished and flush all pending encodings */
    IByteCompressionAdaptData<T>::entropy_coder_->FinishEncoding();

    /* Set up the BitVector holding the encoded data for further use */
    encoding.TrimToContent();

    /* Get the local remaining significant bits */
    const uint64_t local_remaining_significant_bits_num_bytes = static_cast<uint64_t>(encoding.size());

    /* Get the local encoded LZC, respectively the encoded first "one-bit" positions */
    cmc::bit_map::BitMap local_encoded_lzc_stream = IByteCompressionAdaptData<T>::entropy_coder_->GetEncodedBitStream();
    const uint64_t local_encoded_lzc_stream_num_bytes = static_cast<uint64_t>(local_encoded_lzc_stream.size_bytes());

    /* We need exchange the encoded lengths */
    const std::vector<uint64_t> local_bytes{local_encoded_lzc_stream_num_bytes, local_remaining_significant_bits_num_bytes};
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
        cmc::bit_vector::BitVector encoded_alphabet = IByteCompressionAdaptData<T>::entropy_coder_->EncodeAlphabet();
        const uint64_t encoded_alphabet_num_bytes = static_cast<uint64_t>(encoded_alphabet.size());

        /* Calculate the overall amount of bytes on the root rank */
        const uint64_t num_locally_encoded_entropy_codes_bytes = 4 * sizeof(uint64_t) + encoded_alphabet_num_bytes + local_encoded_lzc_stream_num_bytes;

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
        std::copy_n(encoded_alphabet.begin(), encoded_alphabet_num_bytes, std::back_insert_iterator(encoded_entropy_codes));

        /* Finally, copy the entropy codes */
        std::copy_n(local_encoded_lzc_stream.begin_bytes(), local_encoded_lzc_stream_num_bytes, std::back_insert_iterator(encoded_entropy_codes));
    } else
    {
        /* Otherwise, the rank only hold the entropy codes */
        local_encoded_lzc_stream.MoveDataInto(encoded_entropy_codes);
    }

    /* Get the encoded remaining significant bits */
    encoding.MoveDataInto(encoded_data);

    cmc_debug_msg("The entropy encoder of the multi-resolution extraction compression completed the encoding of the CompressionValues of this iteration.");
    
    return std::make_pair(encoded_entropy_codes, encoded_data);
}

template <typename T>
std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
MultiResAdaptData<T>::EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const
{
    cmc_debug_msg("The encoding of the root level values of the multi-resolution compression starts.");

    std::vector<uint8_t> encoded_stream;
    encoded_stream.reserve(root_level_values.size() * sizeof(T));

    for (auto val_iter = root_level_values.begin(); val_iter != root_level_values.end(); ++val_iter)
    {
        const T val = val_iter->template ReinterpretDataAs<T>();
        cmc_debug_msg("Root val to be encoded: ", val);
        PushBackValueToByteStream(encoded_stream, val);
    }

    cmc_debug_msg("The entropy encoder of the multi-resolution extraction compression stored the root-level CompressionValues within ", encoded_stream.size(), " bytes.");
    return std::make_pair(std::vector<uint8_t>(), std::move(encoded_stream));
}


template <typename T>
inline IByteCompressionAdaptData<T>*
CreateMultiResExtractionAdaptationClass(AbstractByteCompressionVariable<T>* abstract_var, const CompressionSettings& settings)
{
    return new MultiResAdaptData<T>(abstract_var, settings);
}

template <typename T>
inline void
DestroyMultiResExtractionAdaptationClass(IByteCompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class CompressionVariable : public AbstractByteCompressionVariable<T>
{
public:
    CompressionVariable() = delete;

    CompressionVariable(const CompressionSettings& settings, const std::string& name, t8_forest_t initial_mesh, const std::vector<T>& variable_data)
    : AbstractByteCompressionVariable<T>(settings)
    {
        if (static_cast<size_t>(t8_forest_get_local_num_leaf_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(variable_data);
        StoreMeshMPIComm();
        this->SetupInaccuracyTracking();
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreateMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::MeshEncoder>();  
    };

    CompressionVariable(const CompressionSettings& settings, const std::string& name, t8_forest_t initial_mesh, const std::vector<CompressionValue<T>>& variable_data)
    : AbstractByteCompressionVariable<T>(settings)
    {
        if (static_cast<size_t>(t8_forest_get_local_num_leaf_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(variable_data);
        StoreMeshMPIComm();
        this->SetupInaccuracyTracking();
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreateMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::MeshEncoder>();  
    };

    CompressionVariable(const CompressionSettings& settings, const std::string& name, t8_forest_t initial_mesh, std::vector<CompressionValue<T>>&& variable_data)
    : AbstractByteCompressionVariable<T>(settings)
    {
        if (static_cast<size_t>(t8_forest_get_local_num_leaf_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(std::move(variable_data));
        StoreMeshMPIComm();
        this->SetupInaccuracyTracking();
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreateMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::MeshEncoder>();  
    };

    CompressionSchema GetCompressionSchema() const override
    {
        return CompressionSchema::TrimmedMultiResResiduals;
    }

private:
    void StoreMeshMPIComm(){this->SetMPIComm(t8_forest_get_mpicomm(this->GetAmrMesh().GetMesh()));};

};

}

#endif /* !LOSSY_CMC_MULTI_RES_TRIMMED_RESIDUALS_HXX */
