#ifndef CMC_EMBEDDED_MULTI_RES_EXTRACTION_LIN_RECONSTRUCTION_COMPRESSION_HXX
#define CMC_EMBEDDED_MULTI_RES_EXTRACTION_LIN_RECONSTRUCTION_COMPRESSION_HXX


#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "lossless/cmc_embedded_byte_compression_ext_variable.hxx"
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

namespace cmc::lossless::embedded::multi_res::ext_reconstruction
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

/* A container for the collection of element face values */
template <typename T>
struct FaceDataValues
{
    std::vector<CompressionValue<T>> values;
    int num_values_within_family{0};
};

template<typename T>
class MultiResEmbeddedAdaptData : public cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>
{
public:
    MultiResEmbeddedAdaptData() = delete;
    MultiResEmbeddedAdaptData(cmc::lossless::embedded::ext::AbstractEmbeddedByteCompressionVariable<T>* variable)
    : cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>(variable) {
        cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>::entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResEncoder<T>>();
    };

    void InitializeExtractionIteration() override;
    void FinalizingExtractionIteration() override;
    void CompleteExtractionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const override;
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const override;

    protected:
    cmc::lossless::embedded::ext::ExtractionData<T> PerformExtraction(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                        const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[]) override;
    cmc::lossless::embedded::ext::UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) override;

private:
    void CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const;
    std::vector<FaceDataValues<T>> GatherElementFaceValues(t8_forest_t forest, t8_locidx_t ltreeid, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                      const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[]) const;

    bit_map::BitMap resdiual_order_indications_;
    int count_adaptation_step_{0};
};


template <typename T>
void
MultiResEmbeddedAdaptData<T>::InitializeExtractionIteration()
{
    resdiual_order_indications_ = bit_map::BitMap();
}

template <typename T>
void
MultiResEmbeddedAdaptData<T>::FinalizingExtractionIteration()
{
    ++count_adaptation_step_;
}

template <typename T>
void
MultiResEmbeddedAdaptData<T>::CompleteExtractionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    //Nothing to be done here!
}

template <typename T>
void
MultiResEmbeddedAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

//TODO: parallelize
template <typename T>
std::vector<FaceDataValues<T>>
MultiResEmbeddedAdaptData<T>::GatherElementFaceValues(t8_forest_t forest, t8_locidx_t ltreeid, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                      const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[]) const
{
    /* Only a family of elements is passed to this function, therefore all elements share the same parent element  */
  
    /* Allocate two elements in order to compare whether values belong to same family */
    std::array<t8_element_t*, 2> parent_elems{nullptr, nullptr};
    scheme->element_new (tree_class, 2, parent_elems.data());
    
    /* Create the parent element for this family */
    scheme->element_get_parent (tree_class, elements[0], parent_elems[0]); 

    /* Get the global tree id */
    const t8_gloidx_t global_tree_id = t8_forest_global_tree_id (forest, ltreeid);

    std::vector<FaceDataValues<T>> face_values;
    face_values.reserve(num_elements);

    for (int idx = 0; idx < num_elements; ++idx)
    {
        /* Get the number of faces for this element */
        const int num_faces = scheme->element_get_num_faces(tree_class, elements[idx]);

        /* Push back a face value collector for this element */
        face_values.emplace_back();
        face_values.back().values.reserve(num_faces);

        /* Iterate over all faces and gather the corresponding values */
        for (int face_idx = 0; face_idx < num_faces; ++face_idx)
        {
            t8_element_t** neighbor_leaves;
            int* dual_faces;
            int num_neighbors{0};
            t8_locidx_t* neighbor_element_indices;
            t8_eclass_t neighbor_tree_class;
            t8_gloidx_t global_neighbor_tree_id{0};
            int orientation{0};

            /* Get the face neighbors */
            t8_forest_leaf_face_neighbors_ext (forest, ltreeid, elements[idx], &neighbor_leaves, face_idx, &dual_faces, &num_neighbors,
                                               &neighbor_element_indices, &neighbor_tree_class, 1, &global_neighbor_tree_id, &orientation);
            
            //t8_forest_leaf_face_neighbors (forest, ltreeid, elements[idx], &neighbor_leaves, face_idx, &dual_faces, &num_neighbors,
            //                              &neighbor_element_indices, &neighbor_tree_class, 1);
            /* Check if there is a neighboring element at the face */
            if (num_neighbors > 0)
            {
                /* Since the forest is balanced, there should be exactly one face neighbor element */
                cmc_assert(num_neighbors == 1);

                /* Check if the element is from the same tree */
                if (global_neighbor_tree_id == global_tree_id)
                //if (neighbor_tree_class == tree_class)
                {
                    /* Only in this case, we need to check whether the neighboring element is from the same family */
                    /* Construct the parent element */
                    scheme->element_get_parent (tree_class, neighbor_leaves[0], parent_elems[1]);

                    /* Check if the parent elements are equal */
                    const bool have_same_parent_element = scheme->element_is_equal(tree_class, parent_elems[0], parent_elems[1]);
                    if (have_same_parent_element)
                    {
                        /* In case the elements are equal, they belong to the same family */
                        ++(face_values.back().num_values_within_family);
                    } else
                    {
                        /* In case they do not belong to the same family, we gather the value */
                        /* We obtain the local index of all local elements if the element is not a ghost.
                         * Therfore, we can directly access the corresponding value  */
                        cmc_assert(neighbor_element_indices[0] < t8_forest_get_local_num_leaf_elements(forest));
                        face_values.back().values.push_back(this->GetDataValueAtIndex(neighbor_element_indices[0]));
                    }
                } else
                {
                    /* In case the neighboring element is from a different tree, they cannot be part of the same family */
                    /* In case they do not belong to the same family, we gather the value */
                    /* We obtain the local index of all local elements if the element is not a ghost.
                     * Therfore, we can directly access the corresponding value  */
                    cmc_assert(neighbor_element_indices[0] < t8_forest_get_local_num_leaf_elements(forest));
                    face_values.back().values.push_back(this->GetDataValueAtIndex(neighbor_element_indices[0]));
                }

                /* Deallocate the memory for the face neighbor construction */
                scheme->element_destroy (neighbor_tree_class, num_neighbors, neighbor_leaves);
                T8_FREE (neighbor_leaves);
                T8_FREE (neighbor_element_indices);
                T8_FREE (dual_faces);
           
            }
        }
    }

    /* Deallocate the parent elements */
    scheme->element_destroy(tree_class, 2, parent_elems.data());

    /* Return the collected information per face */
    return face_values;
}

template <typename T>
inline T
ComputeReconstructionViaFaceValues(const T prediction_value, const FaceDataValues<T>& face_values)
{   
    const int num_face_vals = face_values.values.size();
    
    if (num_face_vals > 0)
    {
        /* We construct the mean with the neighboring face values */
        T prediction(0);
        for (auto face_iter = face_values.values.begin(); face_iter != face_values.values.end(); ++face_iter)
        {
            prediction += face_iter->template ReinterpretDataAs<T>();
        }
        prediction += (face_values.num_values_within_family * prediction_value);
        prediction /= (face_values.num_values_within_family + num_face_vals);
        
        return prediction;
    } else 
    {
        /* In case there are no face values, we return the prediction */
        return prediction_value;
    }
}

template <typename T>
auto
ComputeCollectiveLZCForPredictor(const CompressionValue<T>& predictor, const VectorView<CompressionValue<T>>& initial_values, const std::vector<FaceDataValues<T>>& face_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint8_t), std::tuple<int, std::vector<CompressionValue<T>>, bit_map::BitMap>>
{
    cmc_assert(initial_values.size() == face_values.size());

    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Collect the residuals */
    std::vector<CompressionValue<T>> residuals;
    residuals.reserve(initial_values.size());

    /* Collect the residual indications */
    bit_map::BitMap residual_indications;
    residual_indications.Reserve(initial_values.size());

    /* Convert the compression value back to its origianl type */
    const T predictor_value = predictor.template ReinterpretDataAs<T>();

    int elem_idx = 0;
    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = initial_values.begin(); val_iter != initial_values.end(); ++val_iter, ++elem_idx)
    {
        /* Get an unsigned representation of the approximation value */
        const T prediction = ComputeReconstructionViaFaceValues<T>(predictor_value, face_values[elem_idx]);
        
        /* Get an unsigned representation of the prediction */
        uint8_t unsigned_prediction;
        std::memcpy(&unsigned_prediction, &prediction, N);

        /* Get an unsigned representation of the inital value */
        uint8_t unsigned_initial_value;
        std::memcpy(&unsigned_initial_value, val_iter->GetMemoryForReading().data(), N);

        /* Check whether the approximation is greater than the initial value */
        const bool is_approximation_greater = (unsigned_prediction >= unsigned_initial_value);
        residual_indications.AppendBit(is_approximation_greater);

        /* Compute the unsigned residual and determine the LZC */
        const uint8_t diff = (unsigned_initial_value >= unsigned_prediction ? unsigned_initial_value - unsigned_prediction : unsigned_prediction - unsigned_initial_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        /* Add the number of the LZC */
        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();

        /* Store the residual */
        residuals.push_back(diff_compression_val);
    }

    return std::make_tuple(cumulative_lzc, residuals, residual_indications);
}


template <typename T>
auto
ComputeCollectiveLZCForPredictor(const CompressionValue<T>& predictor, const VectorView<CompressionValue<T>>& initial_values, const std::vector<FaceDataValues<T>>& face_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint16_t), std::tuple<int, std::vector<CompressionValue<T>>, bit_map::BitMap>>
{
    cmc_assert(initial_values.size() == face_values.size());

    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Collect the residuals */
    std::vector<CompressionValue<T>> residuals;
    residuals.reserve(initial_values.size());

    /* Collect the residual indications */
    bit_map::BitMap residual_indications;
    residual_indications.Reserve(initial_values.size());

    /* Convert the compression value back to its origianl type */
    const T predictor_value = predictor.template ReinterpretDataAs<T>();

    int elem_idx = 0;
    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = initial_values.begin(); val_iter != initial_values.end(); ++val_iter, ++elem_idx)
    {
        /* Get an unsigned representation of the approximation value */
        const T prediction = ComputeReconstructionViaFaceValues<T>(predictor_value, face_values[elem_idx]);
        
        /* Get an unsigned representation of the prediction */
        uint16_t unsigned_prediction;
        std::memcpy(&unsigned_prediction, &prediction, N);

        /* Get an unsigned representation of the inital value */
        uint16_t unsigned_initial_value;
        std::memcpy(&unsigned_initial_value, val_iter->GetMemoryForReading().data(), N);

        /* Check whether the approximation is greater than the initial value */
        const bool is_approximation_greater = (unsigned_prediction >= unsigned_initial_value);
        residual_indications.AppendBit(is_approximation_greater);

        /* Compute the unsigned residual and determine the LZC */
        const uint16_t diff = (unsigned_initial_value >= unsigned_prediction ? unsigned_initial_value - unsigned_prediction : unsigned_prediction - unsigned_initial_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        /* Add the number of the LZC */
        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();

        /* Store the residual */
        residuals.push_back(diff_compression_val);
    }

    return std::make_tuple(cumulative_lzc, residuals, residual_indications);
}


template <typename T>
auto
ComputeCollectiveLZCForPredictor(const CompressionValue<T>& predictor, const VectorView<CompressionValue<T>>& initial_values, const std::vector<FaceDataValues<T>>& face_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint32_t), std::tuple<int, std::vector<CompressionValue<T>>, bit_map::BitMap>>
{
    cmc_assert(initial_values.size() == face_values.size());

    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Collect the residuals */
    std::vector<CompressionValue<T>> residuals;
    residuals.reserve(initial_values.size());

    /* Collect the residual indications */
    bit_map::BitMap residual_indications;
    residual_indications.Reserve(initial_values.size());

    /* Convert the compression value back to its origianl type */
    const T predictor_value = predictor.template ReinterpretDataAs<T>();

    int elem_idx = 0;
    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = initial_values.begin(); val_iter != initial_values.end(); ++val_iter, ++elem_idx)
    {
        /* Get an unsigned representation of the approximation value */
        const T prediction = ComputeReconstructionViaFaceValues<T>(predictor_value, face_values[elem_idx]);
        
        /* Get an unsigned representation of the prediction */
        uint32_t unsigned_prediction;
        std::memcpy(&unsigned_prediction, &prediction, N);

        /* Get an unsigned representation of the inital value */
        uint32_t unsigned_initial_value;
        std::memcpy(&unsigned_initial_value, val_iter->GetMemoryForReading().data(), N);

        /* Check whether the approximation is greater than the initial value */
        const bool is_approximation_greater = (unsigned_prediction >= unsigned_initial_value);
        residual_indications.AppendBit(is_approximation_greater);

        /* Compute the unsigned residual and determine the LZC */
        const uint32_t diff = (unsigned_initial_value >= unsigned_prediction ? unsigned_initial_value - unsigned_prediction : unsigned_prediction - unsigned_initial_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        /* Add the number of the LZC */
        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();

        /* Store the residual */
        residuals.push_back(diff_compression_val);
    }

    return std::make_tuple(cumulative_lzc, residuals, residual_indications);
}


template <typename T>
auto
ComputeCollectiveLZCForPredictor(const CompressionValue<T>& predictor, const VectorView<CompressionValue<T>>& initial_values, const std::vector<FaceDataValues<T>>& face_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint64_t), std::tuple<int, std::vector<CompressionValue<T>>, bit_map::BitMap>>
{
    cmc_assert(initial_values.size() == face_values.size());

    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Collect the residuals */
    std::vector<CompressionValue<T>> residuals;
    residuals.reserve(initial_values.size());

    /* Collect the residual indications */
    bit_map::BitMap residual_indications;
    residual_indications.Reserve(initial_values.size());

    /* Convert the compression value back to its origianl type */
    const T predictor_value = predictor.template ReinterpretDataAs<T>();

    int elem_idx = 0;
    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = initial_values.begin(); val_iter != initial_values.end(); ++val_iter, ++elem_idx)
    {
        /* Get an unsigned representation of the approximation value */
        const T prediction = ComputeReconstructionViaFaceValues<T>(predictor_value, face_values[elem_idx]);
        
        /* Get an unsigned representation of the prediction */
        uint64_t unsigned_prediction;
        std::memcpy(&unsigned_prediction, &prediction, N);

        /* Get an unsigned representation of the inital value */
        uint64_t unsigned_initial_value;
        std::memcpy(&unsigned_initial_value, val_iter->GetMemoryForReading().data(), N);

        /* Check whether the approximation is greater than the initial value */
        const bool is_approximation_greater = (unsigned_prediction >= unsigned_initial_value);
        residual_indications.AppendBit(is_approximation_greater);

        /* Compute the unsigned residual and determine the LZC */
        const uint64_t diff = (unsigned_initial_value >= unsigned_prediction ? unsigned_initial_value - unsigned_prediction : unsigned_prediction - unsigned_initial_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        /* Add the number of the LZC */
        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();

        /* Store the residual */
        residuals.push_back(diff_compression_val);
    }

    return std::make_tuple(cumulative_lzc, residuals, residual_indications);
}

template <typename T>
std::tuple<T, std::vector<CompressionValue<T>>, bit_map::BitMap>
GetCoarseApproximationMaximizingResidualsLZC(const VectorView<CompressionValue<T>>& initial_values, const std::vector<FaceDataValues<T>>& face_values)
{
    cmc_assert(initial_values.size() >= 1);

    std::vector<CompressionValue<T>> prediction_residuals;
    bit_map::BitMap residual_idications;
    T current_best_predictor{T()};
    int current_max_lzc{-1};

    /* Try each value from the view as a predictor */
    for (auto pred_iter = initial_values.begin(); pred_iter != initial_values.end(); ++pred_iter)
    {
        /* Compute the LZC for this predictor */
        auto [prediction_lzc, residuals, residual_op_indications] = ComputeCollectiveLZCForPredictor<T>(*pred_iter, initial_values, face_values);
    
        /* Potentially, update the current best predictor */
        if (prediction_lzc > current_max_lzc)
        {
            /* Convert the compression value back to its origianl type */
            current_best_predictor = pred_iter->template ReinterpretDataAs<T>();
            current_max_lzc = prediction_lzc;
            prediction_residuals = std::move(residuals);
            residual_idications = std::move(residual_op_indications);
        }
    }

    /* Convert the view to actual values of the underlying data type */
    const std::vector<T> converted_values = ConvertCompressionValues<T>(initial_values);

    /* Try the mid-range as an predictor */
    const T mid_range = InterpolateToMidRange<T>(converted_values);

    /* Compute the LZC in the residuals for this approximation */
    auto [mid_range_lzc, mid_range_residuals, mid_range_residual_op_indications]  = ComputeCollectiveLZCForPredictor<T>(CompressionValue<T>(mid_range), initial_values, face_values);

    /* Potentially, update the current best predictor */
    if (mid_range_lzc > current_max_lzc)
    {
        current_best_predictor = mid_range;
        current_max_lzc = mid_range_lzc;
        prediction_residuals = std::move(mid_range_residuals);
        residual_idications = std::move(mid_range_residual_op_indications);
    }

    /* Try the arithmetic mean as an predictor */
    const T mean = InterpolateToArithmeticMean<T>(converted_values);

    /* Compute the LZC in the residuals for this approximation */
    auto [mean_lzc, mean_residuals, mean_residual_op_indications] = ComputeCollectiveLZCForPredictor<T>(CompressionValue<T>(mean), initial_values, face_values);

    /* Potentially, update the current best predictor */
    if (mean_lzc > current_max_lzc)
    {
        current_best_predictor = mean;
        current_max_lzc = mean_lzc;
        prediction_residuals = std::move(mean_residuals);
        residual_idications = std::move(mean_residual_op_indications);
    }

    return std::make_tuple(current_best_predictor, std::move(prediction_residuals), std::move(residual_idications));
}

template <typename T>
cmc::lossless::embedded::ext::ExtractionData<T>
MultiResEmbeddedAdaptData<T>::PerformExtraction(t8_forest_t forest, t8_locidx_t ltreeid, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[])
{
    /* Get the view on the current values of this family of elements */
    const VectorView<CompressionValue<T>> values = this->GetView(ltreeid, lelement_id, num_elements);

    /* Gather the information about the face neighbors */
    const std::vector<FaceDataValues<T>> face_values = GatherElementFaceValues(forest, ltreeid, tree_class, lelement_id, scheme, num_elements, elements);

    /* Get the coarse approximation for these values which maximizes the cumulative residual LZC */
    auto [coarse_approximation, residuals, residual_indications] = GetCoarseApproximationMaximizingResidualsLZC<T>(values, face_values);

    /* We store the residual order indication */
    resdiual_order_indications_.AppendBits(residual_indications);

    /* And return the extarcted data */
    return cmc::lossless::embedded::ext::ExtractionData<T>(CompressionValue<T>(coarse_approximation), std::move(residuals));
}

template <typename T>
cmc::lossless::embedded::ext::UnchangedData<T>
MultiResEmbeddedAdaptData<T>::ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value)
{
    /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
    resdiual_order_indications_.AppendUnsetBit();

    /* We drag this value along to the "coarser level" until it this element is passed with its siblings
     * as a family of elements into the adaptation callback. Moreover, we leave an empty value 
     * on the "finer level" */
    return cmc::lossless::embedded::ext::UnchangedData<T>(value, CompressionValue<T>());
}

template <typename T>
void
MultiResEmbeddedAdaptData<T>::CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    /* Get a view on the bitmap storing the residual addition/subtraction flags */
    bit_map::BitMapView residual_flags(resdiual_order_indications_);

    /* Collect the symbol frequencies */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /* Get the current residual */
        CompressionValue<T> val = *val_iter;

        /* Get the encoded LZC */
        const uint32_t signum = cmc::lossless::multi_res::util::GetSignumForEncoding(residual_flags.GetNextBit());
        const uint32_t first_one_bit = val.GetNumberLeadingZeros();

        /* Update this symbol for encoding */
        cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>::entropy_coder_->UpdateSymbolFrequency(signum + first_one_bit);
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
MultiResEmbeddedAdaptData<T>::EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    cmc_debug_msg("The encoding of the CompressionValues after the multi-resolution extraction iteration starts...");
    
    cmc_assert(cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>::entropy_coder_ != nullptr);

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
    cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>::entropy_coder_->Reset(std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResCompressionAlphabet<T>>());
    cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>::entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Collect the symbols and their frequencies for the entropy coder */
    CollectSymbolFrequenciesForEntropyCoding(level_byte_values);

    /* Setup the interior structure for encoding */
    cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>::entropy_coder_->SetupEncoding(this->GetMPIComm());

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
        cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(signum + first_one_bit);

        /* The first one is implicitly given by the leading zero count; therefore we set the "front bit" in order to discard the leading zeros and the following one */
        val.SetFrontBit(first_one_bit + 1);

        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            encoding.AppendBits(val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
        }
    }

    /* We set an indicaton symbol that the process local end of the values have been reached */
    cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte);

    /* Indicate that the encoding has been finished and flush all pending encodings */
    cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>::entropy_coder_->FinishEncoding();

    /* Set up the BitVector holding the encoded data for further use */
    encoding.TrimToContent();

    /* Get the local remaining significant bits */
    const uint64_t local_remaining_significant_bits_num_bytes = static_cast<uint64_t>(encoding.size());

    /* Get the local encoded LZC, respectively the encoded first "one-bit" positions */
    cmc::bit_map::BitMap local_encoded_lzc_stream = cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>::entropy_coder_->GetEncodedBitStream();
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
        cmc::bit_vector::BitVector encoded_alphabet = cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeAlphabet();
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
inline std::pair<std::vector<uint8_t>, std::vector<uint8_t>>
MultiResEmbeddedAdaptData<T>::EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const
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
inline cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>*
CreateMultiResEmbeddedExtractionAdaptationClass(cmc::lossless::embedded::ext::AbstractEmbeddedByteCompressionVariable<T>* abstract_var)
{
    return new MultiResEmbeddedAdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyMultiResEmbeddedExtractionAdaptationClass(cmc::lossless::embedded::ext::IEmbeddedCompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class EmbeddedCompressionVariable : public cmc::lossless::embedded::ext::AbstractEmbeddedByteCompressionVariable<T>
{
public:
    EmbeddedCompressionVariable() = delete;

    EmbeddedCompressionVariable(input::Var& input_variable)
    : cmc::lossless::embedded::ext::AbstractEmbeddedByteCompressionVariable<T>(input_variable)
    {
        //cmc_assert(input_variable.IsValid());//TODO: implement 

        this->SetName(input_variable.GetName());
        this->IndicateWhetherMeshRefinementBitsWillBeStored(true);

        cmc::lossless::embedded::ext::AbstractEmbeddedByteCompressionVariable<T>::adaptation_creator_ = CreateMultiResEmbeddedExtractionAdaptationClass<T>;
        cmc::lossless::embedded::ext::AbstractEmbeddedByteCompressionVariable<T>::adaptation_destructor_ = DestroyMultiResEmbeddedExtractionAdaptationClass<T>;
        cmc::lossless::embedded::ext::AbstractEmbeddedByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::EmbeddedMeshEncoder>();
    };

    CompressionSchema GetCompressionSchema() const override
    {
        return CompressionSchema::EmbeddedMultiResExtractionLinReconstruction;
    }

private:

};



}


#endif /* !CMC_EMBEDDED_MULTI_RES_EXTRACTION_LIN_RECONSTRUCTION_COMPRESSION_HXX */
