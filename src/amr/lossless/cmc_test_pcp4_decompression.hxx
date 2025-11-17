#ifndef CMC_TEST_PCP4_DECOMPRESSION_HXX
#define CMC_TEST_PCP4_DECOMPRESSION_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_multi_res_entropy_coder.hxx"
#include "utilities/cmc_multi_res_extraction_util.hxx"
#include "amr/lossless/cmc_byte_decompression_variable.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "mesh_compression/cmc_mesh_decoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>
#include <memory>

namespace cmc::lossless::test_pcp4
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template<typename T>
class TestPcp4DecompressionAdaptData : public cmc::decompression::IDecompressionAdaptData<T>
{
public:
    TestPcp4DecompressionAdaptData() = delete;
    TestPcp4DecompressionAdaptData(cmc::decompression::AbstractByteDecompressionVariable<T>* variable)
    : cmc::decompression::IDecompressionAdaptData<T>(variable) {};

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
    CompressionValue<T> GetNextResidualAppliedValue(const int lzc, const CompressionValue<T>& value);

    size_t level_byte_offset_{0};
    bit_vector::BitVectorView encoded_lzcs_;
    bit_vector::BitVectorView residual_bits_;

    std::unique_ptr<cmc::entropy_coding::arithmetic_coding::Decoder> entropy_decoder_{nullptr};

    int count_adaptation_step_{0};
};

template<typename T>
inline bool
TestPcp4DecompressionAdaptData<T>::IsDecompressionProgressing() const
{
    return (level_byte_offset_ < cmc::decompression::IDecompressionAdaptData<T>::encoded_data_byte_stream_.size());
}


template<typename T>
inline uint32_t
TestPcp4DecompressionAdaptData<T>::GetNextEncodedResidualLength()
{
    std::vector<uint8_t> encoded_lzc = encoded_lzcs_.GetNextBitSequence(4);
    cmc_assert(encoded_lzc.size() == 1);

    return static_cast<uint32_t>(encoded_lzc.front() >> 4);
}

template<typename T>
inline std::vector<uint8_t>
TestPcp4DecompressionAdaptData<T>::GetNextResidualBitSequence(const size_t num_residual_bits)
{
    return residual_bits_.GetNextBitSequence(num_residual_bits);
}


template <typename T>
std::vector<CompressionValue<T>>
TestPcp4DecompressionAdaptData<T>::DecodeRootLevel(const t8_locidx_t num_local_root_values)
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
TestPcp4DecompressionAdaptData<T>::InitializeDecompressionIteration()
{
    cmc_debug_msg("A multi-resolution pcp4 decompression iteration is initialized.");

    constexpr size_t offset = sizeof(uint64_t);

    size_t processed_bytes = level_byte_offset_;

    /* Get a pointer to the beginning of the encoded data stream */
    const auto data_start_ptr = cmc::decompression::IDecompressionAdaptData<T>::encoded_data_byte_stream_.data();

    /* Get the amount of relevant bytes for this decompression level */
    const uint64_t current_level_bytes = GetValueFromByteStream<uint64_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;

    cmc_debug_msg("The current refinement level is described by ", current_level_bytes, " bytes.");

    /* Get the bytes for the encoded prefix lengths */
    const uint64_t encoded_lzc_bytes = GetValueFromByteStream<uint64_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;

    /* Get the bytes for the remaining bits */
    const uint64_t residual_bytes = GetValueFromByteStream<uint64_t>(data_start_ptr + processed_bytes);
    processed_bytes += offset;


    /* Set the view on the encoded prefix lengths */
    encoded_lzcs_ = bit_vector::BitVectorView(data_start_ptr + processed_bytes, encoded_lzc_bytes);
    processed_bytes += encoded_lzc_bytes;

    /* Set the view on the remaining bits */
    residual_bits_ = bit_vector::BitVectorView(data_start_ptr + processed_bytes, residual_bytes);
    processed_bytes += residual_bytes;

    /* Update the byte count */
    level_byte_offset_ = processed_bytes;
}

template <typename T>
void
TestPcp4DecompressionAdaptData<T>::FinalizeDecompressionIteration()
{
    ++count_adaptation_step_;
    cmc_debug_msg("The multi-resolution decompression iteration (", count_adaptation_step_, ") has been finalized.");

    //if constexpr (std::is_same_v<T, double>)
    //{
    //    std::vector<double> double_data;
    //    double_data.reserve(this->GetDecompressedData().size());
//
    //    for (auto val_iter = this->GetDecompressedData().begin(); val_iter != this->GetDecompressedData().end(); ++val_iter)
    //    {
    //        double_data.push_back(val_iter->template ReinterpretDataAs<double>());
    //    }
//
    //    t8_vtk_data_field_t vtk_data[1];
    //    snprintf (vtk_data[0].description, BUFSIZ, "ArtificialData");
    //    vtk_data[0].type = T8_VTK_SCALAR;
    //    vtk_data[0].data = double_data.data();
//
    //    std::string filename = "HybridCircleSquareReconstruction_Step_" + std::to_string(count_adaptation_step_ - 1);
    //    t8_forest_write_vtk_ext (this->GetAmrMesh().GetMesh(), filename.c_str(), 0, 0, 0, 0, 0, 0, 0, 1, vtk_data);
    //}
}

template <typename T>
void
TestPcp4DecompressionAdaptData<T>::CompleteDecompressionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    //Nothing to be done here!
}

template <typename T>
void
TestPcp4DecompressionAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

template <typename T>
CompressionValue<T>
TestPcp4DecompressionAdaptData<T>::GetNextResidualAppliedValue(const int lzc, const CompressionValue<T>& coarse_value)
{
    CompressionValue<T> value = coarse_value;

    /* Get the maximum length of for the given type */
    const uint32_t max_length_type = sizeof(T) * bit_map::kCharBit;

    /* Check if there is a residual to add/subtract */
    if (lzc < max_length_type)
    {
        /* Determine the length of the encoded residual */
        const uint32_t residual_length = max_length_type - lzc;

        /* Get the residual bit sequence */
        const std::vector<uint8_t> residual_bits = this->GetNextResidualBitSequence(residual_length);

        /** Construct a CompressionValue holding the residual **/
        std::array<uint8_t, sizeof(T)> serialized_residual;
        serialized_residual.fill(uint8_t{0});

        CompressionValue<T> residual(serialized_residual);

        /* Set the tail such that the LZC is represented */
        residual.SetTailBit(static_cast<uint8_t>(residual_length));

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
        value._TestPCP4ReverseXORResidual(residual);
    }

    return value;
}

template <typename T>
cmc::decompression::RefinementData<T>
TestPcp4DecompressionAdaptData<T>::PerformRefinement([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, const CompressionValue<T> value, const int num_refined_elements)
{
   /* Get the LZC of the next residual */
   uint32_t lzc = this->GetNextEncodedResidualLength();

   /* Create the refinement data to return */
   cmc::decompression::RefinementData<T> refinement_data;
   refinement_data.fine_values.reserve(num_refined_elements);

   /* Apply all children suffixes */
   for (int idx = 0; idx < num_refined_elements; ++idx)
   {
       /* Get the next value with the applied residual and store it wihtin the refinement data */
       refinement_data.fine_values.emplace_back(this->GetNextResidualAppliedValue(lzc, value));
   }

   return refinement_data;
}

template <typename T>
cmc::decompression::UnchangedData<T>
TestPcp4DecompressionAdaptData<T>::ElementStaysUnchanged([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, const CompressionValue<T>& value)
{
   /* Get the LZC of the next residual */
   uint32_t lzc = this->GetNextEncodedResidualLength();

   /* Get the netx suffixed value */
   const CompressionValue<T> residual_applied_value = GetNextResidualAppliedValue(lzc, value);

   return cmc::decompression::UnchangedData<T>(residual_applied_value);
}

template <typename T>
inline cmc::decompression::IDecompressionAdaptData<T>*
CreateTestPcp4DecompressionAdaptationClass(cmc::decompression::AbstractByteDecompressionVariable<T>* abstract_var)
{
    return new TestPcp4DecompressionAdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyTestPcp4DecompressionAdaptationClass(cmc::decompression::IDecompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class DecompressionVariable : public cmc::decompression::AbstractByteDecompressionVariable<T>
{
public:
    DecompressionVariable() = delete;

    DecompressionVariable(const std::string& name, std::vector<uint8_t>&& encoded_data_byte_stream, std::vector<uint8_t>&& encoded_mesh_byte_stream, const int max_num_decompression_iterations)
    : cmc::decompression::AbstractByteDecompressionVariable<T>(std::move(encoded_data_byte_stream), std::move(encoded_mesh_byte_stream), max_num_decompression_iterations)
    {
        this->SetName(name);
        cmc::decompression::AbstractByteDecompressionVariable<T>::adaptation_creator_ = CreateTestPcp4DecompressionAdaptationClass<T>;
        cmc::decompression::AbstractByteDecompressionVariable<T>::adaptation_destructor_ = DestroyTestPcp4DecompressionAdaptationClass<T>;
        cmc::decompression::AbstractByteDecompressionVariable<T>::mesh_decoder_ = std::make_unique<mesh_compression::MeshDecoder>(this->GetEncodedMeshStream());
    };
};


}



#endif /* !CMC_TEST_PCP4_DECOMPRESSION_HXX */
