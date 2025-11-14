#ifndef CMC_TEST_PCP4_COMPRESSION_HXX
#define CMC_TEST_PCP4_COMPRESSION_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_interpolation_fn.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_multi_res_entropy_coder.hxx"
#include "lossless/cmc_byte_compression_variable.hxx"


#include <utility>
#include <vector>
#include <array>
#include <algorithm>

namespace cmc::lossless::test_pcp4
{

template<typename T>
class TestPcp4AdaptData : public ICompressionAdaptData<T>
{
public:
    TestPcp4AdaptData() = delete;
    TestPcp4AdaptData(AbstractByteCompressionVariable<T>* variable)
    : ICompressionAdaptData<T>(variable) {
        ICompressionAdaptData<T>::entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResEncoder<T>>();
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

    bit_vector::BitVector fixed_four_byte_family_codes_;

    int count_adaptation_step_{0};
};

template <typename T>
void
TestPcp4AdaptData<T>::InitializeExtractionIteration()
{
    fixed_four_byte_family_codes_ = bit_vector::BitVector();
}

template <typename T>
void
TestPcp4AdaptData<T>::FinalizeExtractionIteration()
{
    ++count_adaptation_step_;
}

template <typename T>
void
TestPcp4AdaptData<T>::CompleteExtractionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    fixed_four_byte_family_codes_.TrimToContent();
}

template <typename T>
void
TestPcp4AdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
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
std::pair<T, int>
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
        return std::make_pair(value_view[index], max_lzc);
    }  else if (mean_lzc_count >= max_lzc && mean_lzc_count >= mr_lzc_count)
    {
        return std::make_pair(mean, mean_lzc_count);
    } else if (mr_lzc_count >= max_lzc && mr_lzc_count >= mean_lzc_count)
    {
        return std::make_pair(mid_range, mr_lzc_count);
    } else
    {
        cmc_err_msg("One of the above cases should have been selected");
        return std::make_pair(T(0.0), 0);
    }

    } else if constexpr (std::is_same_v<double, T>)
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
            uint64_t cp;
            std::memcpy(&cp, &current_predictor, 8);
    
            uint64_t xor_result{0};
    
            for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
            {
                T val = *val_iter;
                uint64_t uiv;
                std::memcpy(&uiv, &val, 8);
    
                const uint64_t diff = cp ^ uiv;
    
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
        uint64_t ui_mid_range;
        std::memcpy(&ui_mid_range, &mid_range, 8);
    
        uint64_t mr_xor_result{0};
    
        for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
        {
            T val = *val_iter;
            uint64_t uiv;
            std::memcpy(&uiv, &val, 8);
    
            const uint64_t diff = ui_mid_range ^ uiv;
    
            mr_xor_result |= diff;
        }    
    
        CompressionValue<T> mr_crv(mr_xor_result);
        const int mr_lzc_count = mr_crv.GetNumberLeadingZeros();
    
        /* Check the arithmetic mean as well */
        const T mean = InterpolateToArithmeticMean<T>(value_view);
        uint64_t ui_mean;
        std::memcpy(&ui_mean, &mean, 8);
    
        uint64_t mean_xor_result{0};
    
        for (auto val_iter = value_view.begin(); val_iter != value_view.end(); ++val_iter)
        {
            T val = *val_iter;
            uint64_t uiv;
            std::memcpy(&uiv, &val, 8);
    
            const uint64_t diff = ui_mean ^ uiv;
    
            mean_xor_result |= diff;
        }
    
        CompressionValue<T> mean_crv(mean_xor_result);
        const int mean_lzc_count = mean_crv.GetNumberLeadingZeros();
    
        /* Check which value to return */
        if (max_lzc >= mean_lzc_count && max_lzc >= mr_lzc_count)
        {
            return std::make_pair(value_view[index], max_lzc);
        }  else if (mean_lzc_count >= max_lzc && mean_lzc_count >= mr_lzc_count)
        {
            return std::make_pair(mean, mean_lzc_count);
        } else if (mr_lzc_count >= max_lzc && mr_lzc_count >= mean_lzc_count)
        {
            return std::make_pair(mid_range, mr_lzc_count);
        } else
        {
            cmc_err_msg("One of the above cases should have been selected");
            return std::make_pair(T(0.0), 0);
        }
    } else
    {
        cmc_err_msg("Only floating point data compression is enabled for this comparison test.");
        return std::make_pair(T(),0);
    }
}

template <typename T>
ExtractionData<T>
TestPcp4AdaptData<T>::PerformExtraction([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, [[maybe_unused]] const int num_elements, const VectorView<CompressionValue<T>> values)
{
    /* Get the coarse approximation for these values which maximizes the cumulative residual LZC */
    auto [coarse_approximation, lzc] = GetInterpoaltionMaximizingLZCInXOR<T>(values);

    std::vector<CompressionValue<T>> fine_values;
    fine_values.reserve(values.size());

    const std::vector<T> vals = GetCompressionValuesAs<T>(values);

    if constexpr (std::is_same_v<float, T>)
    {
        uint32_t xor_result{0};
        uint32_t interpolation_result;
        std::memcpy(&interpolation_result, &coarse_approximation, sizeof(T));

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
        cmc_assert(crv.GetNumberLeadingZeros() == lzc);
    } else if constexpr (std::is_same_v<double, T>)
    {
        uint64_t xor_result{0};
        uint64_t interpolation_result;
        std::memcpy(&interpolation_result, &coarse_approximation, sizeof(T));

        /* Compute the residuals for the given predictor */
        for (size_t idx = 0; idx < vals.size(); ++idx)
        {
            uint64_t uiv;
            std::memcpy(&uiv, &vals[idx], sizeof(T));

            const uint64_t xor_residual = interpolation_result ^ uiv;

            fine_values.push_back(CompressionValue<T>(xor_residual));
        
            xor_result |= xor_residual;
        }
        CompressionValue<T> crv(xor_result);
        cmc_assert(crv.GetNumberLeadingZeros() == lzc);
    } else
    {
        cmc_err_msg("Only double and float compression is enabled for the testing purposes of the PCP4 scheme.");
    }

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
TestPcp4AdaptData<T>::ElementStaysUnchanged([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, const CompressionValue<T>& value)
{
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
TestPcp4AdaptData<T>::EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    cmc_debug_msg("The encoding of the CompressionValues after the multi-resolution pcp4 extraction iteration starts...");
    
    cmc_assert(ICompressionAdaptData<T>::entropy_coder_ != nullptr);

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

    /* Get the local encoded LZC, respectively the encoded first "one-bit" positions */
    const uint64_t local_encoded_lzc_stream_num_bytes = static_cast<uint64_t>(fixed_four_byte_family_codes_.size());

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
std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
TestPcp4AdaptData<T>::EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const
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
inline ICompressionAdaptData<T>*
CreateTestPcp4ExtractionAdaptationClass(AbstractByteCompressionVariable<T>* abstract_var)
{
    return new TestPcp4AdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyTestPcp4ExtractionAdaptationClass(ICompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class CompressionVariable : public AbstractByteCompressionVariable<T>
{
public:
    CompressionVariable() = delete;

    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, const std::vector<T>& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_leaf_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(variable_data);
        StoreMeshMPIComm();
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreateTestPcp4ExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyTestPcp4ExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::MeshEncoder>();  
    };

    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, const std::vector<CompressionValue<T>>& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_leaf_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(variable_data);
        StoreMeshMPIComm();
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreateTestPcp4ExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyTestPcp4ExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::MeshEncoder>();  
    };

    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, std::vector<CompressionValue<T>>&& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_leaf_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(std::move(variable_data));
        StoreMeshMPIComm();
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreateTestPcp4ExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyTestPcp4ExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::MeshEncoder>();  
    };

    CompressionSchema GetCompressionSchema() const override
    {
        return CompressionSchema::_TestPCP4Extraction;
    }

private:
    void StoreMeshMPIComm(){this->SetMPIComm(t8_forest_get_mpicomm(this->GetAmrMesh().GetMesh()));};

};

}

#endif /* !CMC_TEST_PCP4_COMPRESSION_HXX */
