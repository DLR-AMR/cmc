#ifndef CMC_TEST_PCP4_EMBEDDED_COMPRESSION_HXX
#define CMC_TEST_PCP4_EMBEDDED_COMPRESSION_HXX


#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "embedded/lossless/cmc_embedded_byte_compression_variable.hxx"
#include "utilities/cmc_multi_res_entropy_coder.hxx"
#include "utilities/cmc_multi_res_extraction_util.hxx"
#include "utilities/cmc_lossless_multi_res_extraction_residual_computation.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_interpolation_fn.hxx"
#include "mesh_compression/cmc_embedded_mesh_encoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>

namespace cmc::lossless::embedded::test_pcp4
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template<typename T>
class TestPcp4EmbeddedAdaptData : public IEmbeddedCompressionAdaptData<T>
{
public:
    TestPcp4EmbeddedAdaptData() = delete;
    TestPcp4EmbeddedAdaptData(AbstractEmbeddedByteCompressionVariable<T>* variable)
    : IEmbeddedCompressionAdaptData<T>(variable) {
        IEmbeddedCompressionAdaptData<T>::entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResEncoder<T>>();
    };

    void InitializeExtractionIteration() override;
    void FinalizeExtractionIteration() override;
    void CompleteExtractionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const override;
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const override;
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLeafLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const;
protected:
    ExtractionData<T> PerformExtraction(const int which_tree, const int lelement_id, const int num_elements, const VectorView<CompressionValue<T>> values) override;
    UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) override;

private:
    void CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const;
    void CollectLeafLevelSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const;
    
    bit_vector::BitVector fixed_four_byte_family_codes_;

    int count_adaptation_step_{0};
};


template <typename T>
void
TestPcp4EmbeddedAdaptData<T>::InitializeExtractionIteration()
{
    fixed_four_byte_family_codes_ = bit_vector::BitVector();
    //fixed_four_byte_family_codes_.Reserve();
}

template <typename T>
void
TestPcp4EmbeddedAdaptData<T>::FinalizeExtractionIteration()
{
    ++count_adaptation_step_;
}

template <typename T>
void
TestPcp4EmbeddedAdaptData<T>::CompleteExtractionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    fixed_four_byte_family_codes_.TrimToContent();
}

template <typename T>
void
TestPcp4EmbeddedAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}


template <typename T>
std::vector<T>
GetCompressionValuesAs(const VectorView<CompressionValue<T>> values)
{
    std::vector<T> vals;
    vals.reserve(values.size());

    for (auto iter = values.begin(); iter != values.end(); ++iter)
    {
        vals.push_back(iter->template ReinterpretDataAs<T>());
    }

    return vals;
}


template<typename T>
T
GetInterpoaltionMaximizingLZCInXOR(const VectorView<CompressionValue<T>> values)
{
    if constexpr (std::is_same_v<float, T>)
    {
    const std::vector<T> vals = GetCompressionValuesAs<T>(values);

    /* Get a view on the values */
    const VectorView<T> value_view(vals);

    int max_lzc = -1;
    int index = 0;

    int pred_idx = 0;
    for (auto pred_iter = value_view.begin(); pred_iter != value_view.end(); ++pred_iter, ++pred_idx)
    {
        T current_predictor = *pred_iter;
        uint32_t cp;
        std::memcpy(&cp, &current_predictor, 4);

        uint32_t xor_result{0};

        for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
        {
            T val = *val_iter;
            uint32_t uiv;
            std::memcpy(&uiv, &val, 4);

            const uint32_t diff = cp ^ uiv;

            xor_result |= diff;
        }        

        CompressionValue<T> crv(xor_result);
        const int lzc_count = crv.GetNumberLeadingZeros();

        if (lzc_count > max_lzc)
        {
            max_lzc = lzc_count;
            index = pred_idx;
        }
    }

    //Try midrange and arithmetic mean as well
    const T mid_range = InterpolateToMidRange<T>(value_view);
    uint32_t ui_mid_range;
    std::memcpy(&ui_mid_range, &mid_range, 4);

    uint32_t mr_xor_result{0};

    for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
    {
        T val = *val_iter;
        uint32_t uiv;
        std::memcpy(&uiv, &val, 4);

        const uint32_t diff = ui_mid_range ^ uiv;

        mr_xor_result |= diff;
    }    

    CompressionValue<T> mr_crv(mr_xor_result);
    const int mr_lzc_count = mr_crv.GetNumberLeadingZeros();

    /* Check the arithmetic mean as well */
    const T mean = InterpolateToArithmeticMean<T>(value_view);
    uint32_t ui_mean;
    std::memcpy(&ui_mean, &mean, 4);

    uint32_t mean_xor_result{0};

    for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
    {
        T val = *val_iter;
        uint32_t uiv;
        std::memcpy(&uiv, &val, 4);

        const uint32_t diff = ui_mean ^ uiv;

        mean_xor_result |= diff;
    }    

    CompressionValue<T> mean_crv(mean_xor_result);
    const int mean_lzc_count = mean_crv.GetNumberLeadingZeros();

    /* Check which value to return */
    if (max_lzc >= mean_lzc_count && max_lzc >= mr_lzc_count)
    {
        return value_view[index];
    }  else if (mean_lzc_count >= max_lzc && mean_lzc_count >= mr_lzc_count)
    {
        return mean;
    } else if (mr_lzc_count >= max_lzc && mr_lzc_count >= mean_lzc_count)
    {
        return mid_range;
    } else
    {
        cmc_err_msg("One of the above cases should have been selected");
        return T(0.0);
    }

    } else
    {
        cmc_err_msg("Only floating point data compression is enabled for this comparison test.");
        return T();
    }
}

template <typename T>
ExtractionData<T>
TestPcp4EmbeddedAdaptData<T>::PerformExtraction([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, [[maybe_unused]] const int num_elements, const VectorView<CompressionValue<T>> values)
{
    /* Get the coarse approximation for these values which maximizes the cumulative residual LZC */
    const T coarse_approximation = GetInterpoaltionMaximizingLZCInXOR<T>(values);

    uint32_t interpolation_result;
    std::memcpy(&interpolation_result, &coarse_approximation, sizeof(T));

    std::vector<CompressionValue<T>> fine_values;
    fine_values.reserve(values.size());

    const std::vector<T> vals = GetCompressionValuesAs<T>(values);

    uint32_t xor_result{0};

    /* Compute the residuals for the given predictor */
    for (size_t idx = 0; idx < vals.size(); ++idx)
    {
        uint32_t uiv;
        std::memcpy(&uiv, &vals[idx], sizeof(T));

        const uint32_t xor_residual = interpolation_result ^ uiv;

        fine_values.push_back(CompressionValue<T>(xor_residual));
     
        xor_result |= xor_residual;
    }

    CompressionValue<T> crv(xor_result);
    int lzc = crv.GetNumberLeadingZeros();

    /* Check whether the LZC exceeds 15 */
    if (lzc > 15)
    {
        lzc = 15;
    }

    const uint8_t final_lzc = static_cast<uint8_t>(lzc);

    /* Store the LZC */
    fixed_four_byte_family_codes_.AppendFourBits(final_lzc);

    for (auto fiter = fine_values.begin(); fiter != fine_values.end(); ++fiter)
    {
        fiter->SetFrontBit(lzc);
    }

    return ExtractionData<T>(CompressionValue<T>(coarse_approximation), std::move(fine_values));
}

template <typename T>
UnchangedData<T>
TestPcp4EmbeddedAdaptData<T>::ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value)
{
    /* We drag this value along to the "coarser level" until it this element is passed with its siblings
     * as a family of elements into the adaptation callback. Moreover, we leave an empty value 
     * on the "finer level" */
    CompressionValue<T> coarse_val = value;
    fixed_four_byte_family_codes_.AppendFourBits(uint8_t(15));
    CompressionValue<T> fine_val;
    fine_val.SetFrontBit(15);
    fine_val.SetTailBit(0);

    return UnchangedData<T>(coarse_val, fine_val);
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
TestPcp4EmbeddedAdaptData<T>::EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    cmc_debug_msg("The encoding of the CompressionValues after the multi-resolution pcp4 extraction iteration starts...");
    
    /* Get the rank of the mpi process within the communicator */
    int rank{0};
    int ret_val = MPI_Comm_rank(this->GetMPIComm(), &rank);
    MPICheckError(ret_val);

    /* Define the root rank */
    const int root_rank = 0;

    /* The encoded data will be stored in a BitVector */
    cmc::bit_vector::BitVector encoding;
    encoding.Reserve(3 * level_byte_values.size());
    
    /* Iterate over all values and encode them */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /* Get the current value */
        CompressionValue<T> val = *val_iter;

        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            encoding.AppendBits(val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
        }
    }

    /* Set up the BitVector holding the encoded data for further use */
    encoding.TrimToContent();

    /* Get the local remaining significant bits */
    const uint64_t local_remaining_significant_bits_num_bytes = static_cast<uint64_t>(encoding.size());

    /* Get the local encoded LZC */
    const uint64_t local_encoded_lzc_stream_num_bytes = static_cast<uint64_t>(fixed_four_byte_family_codes_.size());
    cmc_debug_msg("Size of LZC for this level: ", local_encoded_lzc_stream_num_bytes);
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
        /* Calculate the overall amount of bytes on the root rank */
        const uint64_t num_locally_encoded_entropy_codes_bytes = 3 * sizeof(uint64_t) + local_encoded_lzc_stream_num_bytes;

        /* Allocate memory for the encoded data */
        encoded_entropy_codes.reserve(num_locally_encoded_entropy_codes_bytes);
        
        /** We store global information about the encoded level **/
        /* Push back the overall byte count for the level */
        const uint64_t num_global_bytes_encoded_level = 3 * sizeof(uint64_t) + global_bytes[0] + global_bytes[1];
        PushBackValueToByteStream<uint64_t>(encoded_entropy_codes, num_global_bytes_encoded_level);

        /* Push back the byte count for the encoded first "one-bit" positions */
        PushBackValueToByteStream<uint64_t>(encoded_entropy_codes, global_bytes[0]);

        /* Push back the byte count for the remaining significant bits */
        PushBackValueToByteStream<uint64_t>(encoded_entropy_codes, global_bytes[1]);

        /* Finally, copy the entropy codes */
        std::copy_n(fixed_four_byte_family_codes_.begin(), local_encoded_lzc_stream_num_bytes, std::back_inserter(encoded_entropy_codes));
    } else
    {
        /* Otherwise, the rank only hold the entropy codes */
        std::copy_n(fixed_four_byte_family_codes_.begin(), local_encoded_lzc_stream_num_bytes, std::back_inserter(encoded_entropy_codes));
    }

    /* Get the encoded remaining significant bits */
    encoding.MoveDataInto(encoded_data);

    cmc_debug_msg("The entropy encoder of the multi-resolution pcp4 extraction compression completed the encoding of the CompressionValues of this iteration.");
    
    return std::make_pair(encoded_entropy_codes, encoded_data);
}

template <typename T>
inline std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
TestPcp4EmbeddedAdaptData<T>::EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const
{
    cmc_debug_msg("The encoding of the root level values of the multi-resolution pcp4 compression starts.");

    std::vector<uint8_t> encoded_stream;
    encoded_stream.reserve(root_level_values.size() * sizeof(T));

    for (auto val_iter = root_level_values.begin(); val_iter != root_level_values.end(); ++val_iter)
    {
        const T val = val_iter->template ReinterpretDataAs<T>();
        cmc_debug_msg("Root val to be encoded: ", val);
        PushBackValueToByteStream(encoded_stream, val);
    }

    cmc_debug_msg("The entropy encoder of the multi-resolution pcp4 extraction compression stored the root-level CompressionValues within ", encoded_stream.size(), " bytes.");
    return std::make_pair(std::vector<uint8_t>(), std::move(encoded_stream));
}

template <typename T>
inline IEmbeddedCompressionAdaptData<T>*
CreateTestPcp4EmbeddedExtractionAdaptationClass(AbstractEmbeddedByteCompressionVariable<T>* abstract_var)
{
    return new TestPcp4EmbeddedAdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyTestPcp4EmbeddedExtractionAdaptationClass(IEmbeddedCompressionAdaptData<T>* iadapt_data)
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
        //this->SetMPIComm(input_variable.GetMPIComm());
        //cmc_debug_msg("Here is comm: ", input_variable.GetMPIComm());
        this->IndicateWhetherMeshRefinementBitsWillBeStored(true);

        AbstractEmbeddedByteCompressionVariable<T>::adaptation_creator_ = CreateTestPcp4EmbeddedExtractionAdaptationClass<T>;
        AbstractEmbeddedByteCompressionVariable<T>::adaptation_destructor_ = DestroyTestPcp4EmbeddedExtractionAdaptationClass<T>;
        AbstractEmbeddedByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::EmbeddedMeshEncoder>();
    };

    CompressionSchema GetCompressionSchema() const override
    {
        return CompressionSchema::_TestEmbeddedPCP4Extraction;
    }

private:

};



}


#endif /* !CMC_TEST_PCP4_EMBEDDED_COMPRESSION_HXX */
