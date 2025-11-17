#ifndef LOSSY_CMC_EMBEDDED_MULTI_RES_EXTRACTION_NEAREST_NEIGHBOR_RECONSTRUCTION_COMPRESSION_HXX
#define LOSSY_CMC_EMBEDDED_MULTI_RES_EXTRACTION_NEAREST_NEIGHBOR_RECONSTRUCTION_COMPRESSION_HXX


#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "lossy/cmc_embedded_byte_compression_scanline_variable.hxx"
#include "utilities/cmc_multi_res_entropy_coder.hxx"
#include "utilities/cmc_multi_res_extraction_util.hxx"
//#include "lossless/cmc_multi_res_extraction_residual_computation.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_interpolation_fn.hxx"
#include "mesh_compression/cmc_embedded_mesh_encoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>
#include <cmath>

namespace cmc::lossy::embedded::multi_res::scanline
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

using err_type_t = cmc::lossy::embedded::scanline::err_type_t;

/* A container for the collection of element face values */
template <typename T>
struct FaceDataValues
{
    std::vector<CompressionValue<T>> values;
    int num_values_within_family{0};
};

template <typename T>
struct FaceValues
{
    FaceValues(const t8_scheme_c *scheme, const t8_eclass_t tree_class, const t8_element_t* elem)
    {
        values.insert(values.end(), scheme->element_get_num_faces(tree_class, elem), CompressionValue<T>());
    }

    std::vector<CompressionValue<T>> values;
};

template<typename T>
class MultiResEmbeddedAdaptData : public cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>
{
public:
    MultiResEmbeddedAdaptData() = delete;
    MultiResEmbeddedAdaptData(cmc::lossy::embedded::scanline::AbstractEmbeddedByteCompressionVariable<T>* variable, const CompressionSettings& settings)
    : cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>(variable, settings) {
        cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>::entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResEncoder<T>>();
    };

    void InitializeExtractionIteration() override;
    void FinalizingExtractionIteration() override;
    std::pair<cmc::lossy::embedded::scanline::ResidualData<T>, err_type_t> ComputeResiduals(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
        const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
        const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
        const t8_locidx_t first_incoming) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const override;
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const override;

    protected:
    cmc::lossy::embedded::scanline::ExtractionData<T> PerformExtraction(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                        const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[]) override;
    cmc::lossy::embedded::scanline::UnchangedData<T> ElementStaysUnchanged(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                                           const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[]) override;

private:
    void CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const;
    std::vector<FaceDataValues<T>> GatherElementFaceValues(t8_forest_t forest, t8_locidx_t ltreeid, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                      const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[]) const;
    std::pair<bool, CompressionValue<T>> GetFaceValue(t8_forest_t forest, const int tree_id, const t8_scheme_c *scheme, const t8_element_t* elem, const int face_idx) const;
    FaceValues<T> ConstructHexFaceData(t8_forest_t forest_new, const int which_tree, const t8_eclass_t tree_class, const t8_scheme_c *scheme, const t8_element_t* elem) const;
    VectorView<cmc::lossy::embedded::scanline::ElementData<T>> GetScanline(t8_forest_t forest_new, t8_locidx_t which_tree,
                                                                           const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int first_incoming, const int scanline_length);

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
std::pair<bool, CompressionValue<T>>
MultiResEmbeddedAdaptData<T>::GetFaceValue(t8_forest_t forest, const int tree_id, const t8_scheme_c *scheme, const t8_element_t* elem, const int face_idx) const 
{
    t8_element_t** neighbor_leaves;
    int* dual_faces;
    int num_neighbors{0};
    t8_locidx_t* neighbor_element_indices;
    t8_eclass_t neighbor_tree_class;

    /* Gathe the face neighbor via this face */
    t8_forest_leaf_face_neighbors (forest, tree_id, elem, &neighbor_leaves, face_idx, &dual_faces, &num_neighbors,
                                   &neighbor_element_indices, &neighbor_tree_class, 1);

    /* Check if there is a neighboring element at the face */
    if (num_neighbors > 0)
    {
        /* Since the forest is balanced, there should be exactly one face neighbor element */
        cmc_assert(num_neighbors == 1);
        cmc_assert(neighbor_element_indices[0] < t8_forest_get_local_num_leaf_elements(forest));

        CompressionValue<T> coarse_face_value = this->GetAdaptedDataValueAtIndex(neighbor_element_indices[0]);

        /* Deallocate the memory for the face neighbor construction */
        scheme->element_destroy (neighbor_tree_class, num_neighbors, neighbor_leaves);
        T8_FREE (neighbor_leaves);
        T8_FREE (neighbor_element_indices);
        T8_FREE (dual_faces);

        /* Return that there is a face neighbor and the value the face neighboring element holds */
        return std::make_pair(true, coarse_face_value);
    }
    
    /* In case there is no face neighbor at the requested face */
    return std::make_pair(false, CompressionValue<T>());
}

template<typename T>
FaceValues<T> 
MultiResEmbeddedAdaptData<T>::ConstructHexFaceData(t8_forest_t forest_new, const int which_tree, const t8_eclass_t tree_class, const t8_scheme_c *scheme, const t8_element_t* elem) const
{
    /* Allocate data for the face values */
    FaceValues<T> face_values(scheme, tree_class, elem);

    /* Check if there is face value at the face with index 0 */
    int face_idx = 0;
    auto [is_face_value_present0, face_value0] = this->GetFaceValue(forest_new, which_tree, scheme, elem, face_idx);
    if (is_face_value_present0)
    {
        face_values.values[face_idx] = face_value0;
    }

    /* Check if there is face value at the face with index 2 */
    face_idx = 2;
    auto [is_face_value_present2, face_value2] = this->GetFaceValue(forest_new, which_tree, scheme, elem, face_idx);
    if (is_face_value_present2)
    {
        face_values.values[face_idx] = face_value2;
    }

    /* Check if there is face value at the face with index 4 */
    face_idx = 4;
    auto [is_face_value_present4, face_value4] = this->GetFaceValue(forest_new, which_tree, scheme, elem, face_idx);
    if (is_face_value_present4)
    {
        face_values.values[face_idx] = face_value4;
    }

    /* Now, we check the opposite faces and pontentially mirror the values if no face value exists */

    /* Check if there is face value at the face with index 1 */
    face_idx = 1;
    int opposite_face_idx = 0;
    auto [is_face_value_present1, face_value1] = this->GetFaceValue(forest_new, which_tree, scheme, elem, face_idx);
    if (is_face_value_present1)
    {
        face_values.values[face_idx] = face_value1;

        /* Check if face zero has a neighbor, if not, we mirror this face */
        if (face_values.values[opposite_face_idx].IsEmpty())
        {
            face_values.values[opposite_face_idx] = face_value1;
        }
    } else
    {
        /* If there is no face, we check whether the opposite face has a value and mirror it */
        if (not face_values.values[opposite_face_idx].IsEmpty())
        {
            face_values.values[face_idx] = face_values.values[opposite_face_idx];
        } 
    }

    /* Check if there is face value at the face with index 3 */
    face_idx = 3;
    opposite_face_idx = 2;
    auto [is_face_value_present3, face_value3] = this->GetFaceValue(forest_new, which_tree, scheme, elem, face_idx);
    if (is_face_value_present3)
    {
        face_values.values[face_idx] = face_value3;

        /* Check if face zero has a neighbor, if not, we mirror this face */
        if (face_values.values[opposite_face_idx].IsEmpty())
        {
            face_values.values[opposite_face_idx] = face_value3;
        }
    } else
    {
        /* If there is no face, we check whether the opposite face has a value and mirror it */
        if (not face_values.values[opposite_face_idx].IsEmpty())
        {
            face_values.values[face_idx] = face_values.values[opposite_face_idx];
        } 
    }

    /* Check if there is face value at the face with index 3 */
    face_idx = 5;
    opposite_face_idx = 4;
    auto [is_face_value_present5, face_value5] = this->GetFaceValue(forest_new, which_tree, scheme, elem, face_idx);
    if (is_face_value_present5)
    {
        face_values.values[face_idx] = face_value5;

        /* Check if face zero has a neighbor, if not, we mirror this face */
        if (face_values.values[opposite_face_idx].IsEmpty())
        {
            face_values.values[opposite_face_idx] = face_value5;
        }
    } else
    {
        /* If there is no face, we check whether the opposite face has a value and mirror it */
        if (not face_values.values[opposite_face_idx].IsEmpty())
        {
            face_values.values[face_idx] = face_values.values[opposite_face_idx];
        } 
    }

    return face_values;
}

#if 0
//constexpr double weight_coarse_approx = 0.25824569976124335585176;
//constexpr double weight_face_neighbor = 0.24725143341291888138275;

//constexpr double weight_coarse_approx = 0.25;
//constexpr double weight_face_neighbor = 0.25;

//constexpr double weight_coarse_approx = 1.0;
//constexpr double weight_face_neighbor = 0.0;

template<typename T>
T
ComputeHexPrediction(const CompressionValue<T>& coarse_approximation, const FaceValues<T>& face_values, const int child_id)
{
    /* Get the actual value of the coarse approximation */
    const T approximation = coarse_approximation.template ReinterpretDataAs<T>();

    switch (child_id)
    {
        case 0:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();

            const T prediction = (weight_coarse_approx * approximation + weight_face_neighbor * f0 + weight_face_neighbor * f2 + weight_face_neighbor * f4);
            return prediction;
        }
        break;
        case 1:
        {
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();

            const T prediction = (weight_coarse_approx * approximation + weight_face_neighbor * f1 + weight_face_neighbor * f2 + weight_face_neighbor * f4);
            return prediction;
        }
        break;
        case 2:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();

            const T prediction = (weight_coarse_approx * approximation + weight_face_neighbor * f0 + weight_face_neighbor * f3 + weight_face_neighbor * f4);
            return prediction;
        }
        break;
        case 3:
        {
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();

            const T prediction = (weight_coarse_approx * approximation + weight_face_neighbor * f1 + weight_face_neighbor * f3 + weight_face_neighbor * f4);
            return prediction;
        }
        break;
        case 4:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (weight_coarse_approx * approximation + weight_face_neighbor * f0 + weight_face_neighbor * f2 + weight_face_neighbor * f5);
            return prediction;
        }
        break;
        case 5:
        {
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (weight_coarse_approx * approximation + weight_face_neighbor * f1 + weight_face_neighbor * f2 + weight_face_neighbor * f5);
            return prediction;
        }
        break;
        case 6:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (weight_coarse_approx * approximation + weight_face_neighbor * f0 + weight_face_neighbor * f3 + weight_face_neighbor * f5);
            return prediction;
        }
        break;
        case 7:
        {
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (weight_coarse_approx * approximation + weight_face_neighbor * f1 + weight_face_neighbor * f3 + weight_face_neighbor * f5);
            return prediction;
        }
        break;
        default:
            cmc_err_msg("The child id (", child_id, ") does not coincide with a refinement of an hex element.");
            return T();
        break;
    }
}

#endif


constexpr double weight_coarse_approx = 1.0 / 5.0;
constexpr double weight_face_neighbor = 1.0 / 6.0;
constexpr double weight_far_face_neighbor = 1.0 / 10.0;

constexpr double wca = weight_coarse_approx;
constexpr double wfn = weight_face_neighbor;
constexpr double wffn = weight_far_face_neighbor;

template<typename T>
T
ComputeHexPrediction(const CompressionValue<T>& coarse_approximation, const FaceValues<T>& face_values, const int child_id)
{
    /* Get the actual value of the coarse approximation */
    const T approximation = coarse_approximation.template ReinterpretDataAs<T>();

    switch (child_id)
    {
        case 0:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (wca * approximation + wfn * f0 + wfn * f2 + wfn * f4 + wffn * f1 + wffn * f3 + wffn * f5);
            return prediction;
        }
        break;
        case 1:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (wca * approximation + wfn * f1 + wfn * f2 + wfn * f4 + wffn * f0 + wffn * f3 + wffn * f5);
            return prediction;
        }
        break;
        case 2:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (wca * approximation + wfn * f0 + wfn * f3 + wfn * f4 + wffn * f1 + wffn * f2 + wffn * f5);
            return prediction;
        }
        break;
        case 3:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (wca * approximation + wfn * f1 + wfn * f3 + wfn * f4 + wffn * f0 + wffn * f2 + wffn * f5);
            return prediction;
        }
        break;
        case 4:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (wca * approximation + wfn * f0 + wfn * f2 + wfn * f5 + wffn * f1 + wffn * f3 + wffn * f4);
            return prediction;
        }
        break;
        case 5:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (wca * approximation + wfn * f1 + wfn * f2 + wfn * f5 + wffn * f0 + wffn * f3 + wffn * f4);
            return prediction;
        }
        break;
        case 6:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (wca * approximation + wfn * f0 + wfn * f3 + wfn * f5 + wffn * f1 + wffn * f2 + wffn * f4);
            return prediction;
        }
        break;
        case 7:
        {
            const T f0 = face_values.values[0].template ReinterpretDataAs<T>();
            const T f1 = face_values.values[1].template ReinterpretDataAs<T>();
            const T f2 = face_values.values[2].template ReinterpretDataAs<T>();
            const T f3 = face_values.values[3].template ReinterpretDataAs<T>();
            const T f4 = face_values.values[4].template ReinterpretDataAs<T>();
            const T f5 = face_values.values[5].template ReinterpretDataAs<T>();

            const T prediction = (wca * approximation + wfn * f1 + wfn * f3 + wfn * f5 + wffn * f0 + wffn * f2 + wffn * f4);
            return prediction;
        }
        break;
        default:
            cmc_err_msg("The child id (", child_id, ") does not coincide with a refinement of an hex element.");
            return T();
        break;
    }
}



template<typename T>
auto
PredictAndComputeResidual(const CompressionValue<T>& coarse_approximation, const FaceValues<T>& face_values, const CompressionValue<T>& initial_value, const err_type_t remaining_permitted_error, const int child_id)
 -> std::enable_if_t<sizeof(T) == sizeof(uint8_t), std::tuple<bool, CompressionValue<T>, err_type_t>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Predict the value as it would have been done by the reconstruction */
    const T prediction = ComputeHexPrediction<T>(coarse_approximation, face_values, child_id);

    /* Reinterpret the initial value as the underlying data type */
    const T real_init_value = initial_value.template ReinterpretDataAs<T>();

    /* We will Check whether the residual is smaller than the remaining error */
    err_type_t residual_abs_err{std::numeric_limits<err_type_t>::max()};

    /* Compute the residual in its underlying type */
    if constexpr (std::is_signed_v<T>)
    {
        residual_abs_err = std::abs(static_cast<err_type_t>(prediction - real_init_value));
    } else
    {
        residual_abs_err = static_cast<err_type_t>(prediction >= real_init_value ? prediction - real_init_value : real_init_value - prediction);
    }

    /* Check whether the residual is smaller than the remaining error */
    if (residual_abs_err <= remaining_permitted_error)
    {
        /* If so, we do not need to store the residual since its prediction is within the error bound */
        const err_type_t new_remaining_permitted_error = remaining_permitted_error - residual_abs_err;
        return std::make_tuple(true, CompressionValue<T>(),  new_remaining_permitted_error);
    }

    /* In case the absolute residual error is greater than the remaining error, we need to store the residual */
    //TODO: Trim the residual compliant to the error bound 

    /* Get an unsigned representation of the prediction */
    uint8_t unsigned_prediction;
    std::memcpy(&unsigned_prediction, &prediction, N);

    /* Get an unsigned representation of the inital value */
    uint8_t unsigned_initial_value;
    std::memcpy(&unsigned_initial_value, initial_value.GetMemoryForReading().data(), N);

    /* Check whether the approximation is greater than the initial value */
    const bool is_approximation_greater = (unsigned_prediction >= unsigned_initial_value);

    /* Compute the unsigned residual and determine the LZC */
    const uint8_t diff = (unsigned_initial_value >= unsigned_prediction ? unsigned_initial_value - unsigned_prediction : unsigned_prediction - unsigned_initial_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);

    //The remaining error stays the same, since we do not trim the residual currently
    return std::make_tuple(is_approximation_greater, diff_compression_val, remaining_permitted_error);
}


template<typename T>
auto
PredictAndComputeResidual(const CompressionValue<T>& coarse_approximation, const FaceValues<T>& face_values, const CompressionValue<T>& initial_value, const err_type_t remaining_permitted_error, const int child_id)
 -> std::enable_if_t<sizeof(T) == sizeof(uint16_t), std::tuple<bool, CompressionValue<T>, err_type_t>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Predict the value as it would have been done by the reconstruction */
    const T prediction = ComputeHexPrediction<T>(coarse_approximation, face_values, child_id);

    /* Reinterpret the initial value as the underlying data type */
    const T real_init_value = initial_value.template ReinterpretDataAs<T>();

    /* We will Check whether the residual is smaller than the remaining error */
    err_type_t residual_abs_err{std::numeric_limits<err_type_t>::max()};

    /* Compute the residual in its underlying type */
    if constexpr (std::is_signed_v<T>)
    {
        residual_abs_err = std::abs(static_cast<err_type_t>(prediction - real_init_value));
    } else
    {
        residual_abs_err = static_cast<err_type_t>(prediction >= real_init_value ? prediction - real_init_value : real_init_value - prediction);
    }

    /* Check whether the residual is smaller than the remaining error */
    if (residual_abs_err <= remaining_permitted_error)
    {
        /* If so, we do not need to store the residual since its prediction is within the error bound */
        const err_type_t new_remaining_permitted_error = remaining_permitted_error - residual_abs_err;
        return std::make_tuple(true, CompressionValue<T>(),  new_remaining_permitted_error);
    }

    /* In case the absolute residual error is greater than the remaining error, we need to store the residual */
    //TODO: Trim the residual compliant to the error bound 

    /* Get an unsigned representation of the prediction */
    uint16_t unsigned_prediction;
    std::memcpy(&unsigned_prediction, &prediction, N);

    /* Get an unsigned representation of the inital value */
    uint16_t unsigned_initial_value;
    std::memcpy(&unsigned_initial_value, initial_value.GetMemoryForReading().data(), N);

    /* Check whether the approximation is greater than the initial value */
    const bool is_approximation_greater = (unsigned_prediction >= unsigned_initial_value);

    /* Compute the unsigned residual and determine the LZC */
    const uint16_t diff = (unsigned_initial_value >= unsigned_prediction ? unsigned_initial_value - unsigned_prediction : unsigned_prediction - unsigned_initial_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);

    //The remaining error stays the same, since we do not trim the residual currently
    return std::make_tuple(is_approximation_greater, diff_compression_val, remaining_permitted_error);
}

static size_t iiiiiii = 0;
static size_t zero_residuals = 0;

template<typename T>
auto
PredictAndComputeResidual(const CompressionValue<T>& coarse_approximation, const FaceValues<T>& face_values, const CompressionValue<T>& initial_value, const err_type_t remaining_permitted_error, const int child_id)
 -> std::enable_if_t<sizeof(T) == sizeof(uint32_t), std::tuple<bool, CompressionValue<T>, err_type_t>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Predict the value as it would have been done by the reconstruction */
    const T prediction = ComputeHexPrediction<T>(coarse_approximation, face_values, child_id);

    /* Reinterpret the initial value as the underlying data type */
    const T real_init_value = initial_value.template ReinterpretDataAs<T>();

    /* We will Check whether the residual is smaller than the remaining error */
    err_type_t residual_abs_err{std::numeric_limits<err_type_t>::max()};

    /* Compute the residual in its underlying type */
    if constexpr (std::is_signed_v<T>)
    {
        residual_abs_err = std::abs(static_cast<err_type_t>(prediction - real_init_value));
    } else
    {
        residual_abs_err = static_cast<err_type_t>(prediction >= real_init_value ? prediction - real_init_value : real_init_value - prediction);
    }

    /* Check whether the residual is smaller than the remaining error */
    if (residual_abs_err <= remaining_permitted_error)
    {
        ++zero_residuals;
        /* If so, we do not need to store the residual since its prediction is within the error bound */
        const err_type_t new_remaining_permitted_error = remaining_permitted_error - residual_abs_err;
        return std::make_tuple(true, CompressionValue<T>(),  new_remaining_permitted_error);
    }

    /* In case the absolute residual error is greater than the remaining error, we need to store the residual */
    //TODO: Trim the residual compliant to the error bound 

    /* Get an unsigned representation of the prediction */
    uint32_t unsigned_prediction;
    std::memcpy(&unsigned_prediction, &prediction, N);

    /* Get an unsigned representation of the inital value */
    uint32_t unsigned_initial_value;
    std::memcpy(&unsigned_initial_value, initial_value.GetMemoryForReading().data(), N);

    /* Check whether the approximation is greater than the initial value */
    const bool is_approximation_greater = (unsigned_prediction >= unsigned_initial_value);

    /* Compute the unsigned residual and determine the LZC */
    const uint32_t diff = (unsigned_initial_value >= unsigned_prediction ? unsigned_initial_value - unsigned_prediction : unsigned_prediction - unsigned_initial_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);

    //++iiiiiii;
    //if (iiiiiii >= 8000000 && iiiiiii < 8001000)
    //{
    //    const T initial_val = initial_value.template ReinterpretDataAs<T>();
    //    cmc_debug_msg("Index: ", iiiiiii, ", Residual: ", prediction - initial_val);
    //}

    //The remaining error stays the same, since we do not trim the residual currently
    return std::make_tuple(is_approximation_greater, diff_compression_val, remaining_permitted_error);
}


template<typename T>
auto
PredictAndComputeResidual(const CompressionValue<T>& coarse_approximation, const FaceValues<T>& face_values, const CompressionValue<T>& initial_value, const err_type_t remaining_permitted_error, const int child_id)
 -> std::enable_if_t<sizeof(T) == sizeof(uint64_t), std::tuple<bool, CompressionValue<T>, err_type_t>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Predict the value as it would have been done by the reconstruction */
    const T prediction = ComputeHexPrediction<T>(coarse_approximation, face_values, child_id);

    /* Reinterpret the initial value as the underlying data type */
    const T real_init_value = initial_value.template ReinterpretDataAs<T>();

    /* We will Check whether the residual is smaller than the remaining error */
    err_type_t residual_abs_err{std::numeric_limits<err_type_t>::max()};

    /* Compute the residual in its underlying type */
    if constexpr (std::is_signed_v<T>)
    {
        residual_abs_err = std::abs(static_cast<err_type_t>(prediction - real_init_value));
    } else
    {
        residual_abs_err = static_cast<err_type_t>(prediction >= real_init_value ? prediction - real_init_value : real_init_value - prediction);
    }

    /* Check whether the residual is smaller than the remaining error */
    if (residual_abs_err <= remaining_permitted_error)
    {
        /* If so, we do not need to store the residual since its prediction is within the error bound */
        const err_type_t new_remaining_permitted_error = remaining_permitted_error - residual_abs_err;
        return std::make_tuple(true, CompressionValue<T>(),  new_remaining_permitted_error);
    }

    /* In case the absolute residual error is greater than the remaining error, we need to store the residual */
    //TODO: Trim the residual compliant to the error bound 

    /* Get an unsigned representation of the prediction */
    uint64_t unsigned_prediction;
    std::memcpy(&unsigned_prediction, &prediction, N);

    /* Get an unsigned representation of the inital value */
    uint64_t unsigned_initial_value;
    std::memcpy(&unsigned_initial_value, initial_value.GetMemoryForReading().data(), N);

    /* Check whether the approximation is greater than the initial value */
    const bool is_approximation_greater = (unsigned_prediction >= unsigned_initial_value);

    /* Compute the unsigned residual and determine the LZC */
    const uint64_t diff = (unsigned_initial_value >= unsigned_prediction ? unsigned_initial_value - unsigned_prediction : unsigned_prediction - unsigned_initial_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);

    //The remaining error stays the same, since we do not trim the residual currently
    return std::make_tuple(is_approximation_greater, diff_compression_val, remaining_permitted_error);
}

#if 0
template <typename T>
std::pair<cmc::lossy::embedded::scanline::ResidualData<T>, err_type_t>
MultiResEmbeddedAdaptData<T>::ComputeResiduals(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
    const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
    const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
    const t8_locidx_t first_incoming)
{
    /* The computation of residuals takes place here. Therefore, we do need to coarse approximation of the data 
     * and compute the reconstruction of the fine values */
    /* We only need to compute the residuals, if a coarsening has actually been taken place, otherwise the values are just dragged along the compression 
     * until they are a part of a family that is passed to the adaptation callback */
    if (refine != cmc::lossy::embedded::scanline::kIndicationOfCoarseningDuringAdaptation)
    {
        /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
        resdiual_order_indications_.AppendUnsetBit();

        /* Therfore there is no residual associated with this element and we return an empty residual */
        cmc::lossy::embedded::scanline::ResidualData<T> residual_data;
        residual_data.values.push_back(CompressionValue<T>());

        /* Get the remaining error that is associated with this element */
        const err_type_t previous_remaining_permitted_error = this->GetRemainingPermittedError(which_tree, first_outgoing);

        return std::make_pair(residual_data, previous_remaining_permitted_error);
    }

    // This hold only for hex elements currently 
    cmc_assert(tree_class == T8_ECLASS_HEX);

    /* Get the adapted coarse element */
    const t8_element_t* elem = t8_forest_get_leaf_element_in_tree(forest_new, which_tree, first_incoming);

    /* Compute the id of the corresponding data value in the local contiguous array */
    const int local_coarse_elem_index = t8_forest_get_tree_element_offset (forest_new, which_tree) + first_incoming;

    /* Get the coarse approximation of this element */
    const CompressionValue<T> coarse_approximation = this->GetAdaptedDataValueAtIndex(local_coarse_elem_index);

    /* Gather the data of the face values */
    const FaceValues<T> face_values = this->ConstructHexFaceData(forest_new, which_tree, tree_class, scheme, elem);

    /* Get a view on the initial values */
    const VectorView<CompressionValue<T>> initial_values = this->GetView(which_tree, first_outgoing, num_outgoing);

    /* Get a view on the remaining permitted errors */
    const VectorView<err_type_t> remaining_permitted_errors = this->GetViewOnRemainingPermittedErrors(which_tree, first_outgoing, num_outgoing);

    /* Allocate the residuals */
    std::vector<CompressionValue<T>> residuals;
    residuals.reserve(num_outgoing);

    err_type_t min_remaining_permitted_error{std::numeric_limits<err_type_t>::max()};

    /* We iterate over all outgoing fine values predict their value and store the remaining residual */
    for (int elem_child_idx = 0; elem_child_idx < num_outgoing; ++elem_child_idx)
    {
        auto [is_approximation_greater, residual, remaining_error] = PredictAndComputeResidual<T>(coarse_approximation, face_values, initial_values[elem_child_idx], remaining_permitted_errors[elem_child_idx], elem_child_idx);
        //if (elem_child_idx < num_outgoing - 1)
        //{
        resdiual_order_indications_.AppendBit(is_approximation_greater);
        //}
        residuals.push_back(residual);

        /* Find the minimum remaining permitted error for the coarse element */
        if (min_remaining_permitted_error > remaining_error)
        {
            min_remaining_permitted_error = remaining_error;
        }
    }

    cmc_assert(not cmc::ApproxCompare(min_remaining_permitted_error, std::numeric_limits<err_type_t>::max()));

    return std::make_pair(cmc::lossy::embedded::scanline::ResidualData<T>(std::move(residuals)), min_remaining_permitted_error);
}

#endif

template <typename T>
VectorView<cmc::lossy::embedded::scanline::ElementData<T>> 
MultiResEmbeddedAdaptData<T>::GetScanline(t8_forest_t forest_new, t8_locidx_t which_tree,
    const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int first_incoming, const int scanline_length)
{
    const t8_locidx_t offset = t8_forest_get_tree_element_offset (forest_new, which_tree);
    const int start_index = offset + first_incoming;
    cmc_assert(start_index < t8_forest_get_local_num_leaf_elements(forest_new));

    const int half_scan = (scanline_length - 1) / 2;
    cmc_assert(half_scan >= 0);


    const int start_scan = std::max(0, start_index - half_scan);

    const int end_scan = std::min(t8_forest_get_global_num_leaf_elements(forest_new), static_cast<t8_gloidx_t>(start_scan + scanline_length));

    const int count  = end_scan - start_scan;

    cmc_assert(count > 0);
    cmc_assert(count <= scanline_length);

    return this->GetViewOnAdaptedData(start_scan, count);
}


#if 0
static int zero_res_per_def = 0;

template <typename T>
T
PerformScanlinePrediction(const VectorView<cmc::lossy::embedded::scanline::ElementData<T>>& scanline, const std::vector<double>& elem_midpoint)
{
    std::vector<double> weights;
    weights.reserve(scanline.size() * 2);

    for (auto iter = scanline.begin(); iter != scanline.end(); ++iter)
    {
        /* Compute the distance */
        double max_weight = std::sqrt((iter->max_x - elem_midpoint[0]) * (iter->max_x - elem_midpoint[0]) + (iter->max_y - elem_midpoint[1]) * (iter->max_y - elem_midpoint[1]) + (iter->max_z - elem_midpoint[2]) * (iter->max_z - elem_midpoint[2]));
        
        /*  Check whether the disatnce equals zero, if so, we return the actual value at this position as a prediction */
        if (ApproxCompare(max_weight, 0.0))
        {
            ++zero_res_per_def;
            return iter->max;
        }

        /* Store the inverse weight */
        weights.push_back(1.0 / (max_weight * max_weight * max_weight * max_weight * max_weight * max_weight * max_weight * max_weight));

        /* Check if a minimum is also present */
        if(not ApproxCompare(iter->min, std::numeric_limits<T>::max()))
        {
            /* Compute the distance */
            double min_weight = std::sqrt((iter->min_x - elem_midpoint[0]) * (iter->min_x - elem_midpoint[0]) + (iter->min_y - elem_midpoint[1]) * (iter->min_y - elem_midpoint[1]) + (iter->min_z - elem_midpoint[2]) * (iter->min_z - elem_midpoint[2]));
            
            /*  Check whether the disatnce equals zero, if so, we return the actual value at this position as a prediction */
            if (ApproxCompare(min_weight, 0.0))
            {
                ++zero_res_per_def;
                return iter->min;
            }

            /* Store the inverse weight */
            weights.push_back(1.0 / (min_weight * min_weight * min_weight));
        }
    }

    /* Compute the denominator */
    const double denominator = std::accumulate(weights.begin(), weights.end(), 0.0);

    /* Compte the numerator */
    double numerator = 0.0;

    int count = 0;
    for (auto iter = scanline.begin(); iter != scanline.end(); ++iter)
    {
        numerator += weights[count] * iter->max;
        ++count;

        if(not ApproxCompare(iter->min, std::numeric_limits<T>::max()))
        {
            numerator += weights[count] * iter->min;
        ++count;
        }
    }

    /* Compute and return the prediction */
    return numerator / denominator;
}

#else


static int zero_res_per_def = 0;

template <typename T>
T
PerformScanlinePrediction(const VectorView<cmc::lossy::embedded::scanline::ElementData<T>>& scanline, const std::vector<double>& elem_midpoint)
{
    //const int kNumNearestWeights = 27;
    const int kNumNearestWeights = 6;
    cmc_assert(scanline.size() > static_cast<size_t>(kNumNearestWeights));

    std::vector<double> dist;
    dist.reserve(scanline.size() * 2);

    for (auto iter = scanline.begin(); iter != scanline.end(); ++iter)
    {
        /* Compute the distance */
        double dist_to_max = std::sqrt((iter->max_x - elem_midpoint[0]) * (iter->max_x - elem_midpoint[0]) + (iter->max_y - elem_midpoint[1]) * (iter->max_y - elem_midpoint[1]) + (iter->max_z - elem_midpoint[2]) * (iter->max_z - elem_midpoint[2]));
        
        /*  Check whether the disatnce equals zero, if so, we return the actual value at this position as a prediction */
        if (ApproxCompare(dist_to_max, 0.0))
        {
            ++zero_res_per_def;
            return iter->max;
        }
        dist.push_back(dist_to_max);

        /* Check if a minimum is also present */
        if(not ApproxCompare(iter->min, std::numeric_limits<T>::max()))
        {
            /* Compute the distance */
            double dist_to_min = std::sqrt((iter->min_x - elem_midpoint[0]) * (iter->min_x - elem_midpoint[0]) + (iter->min_y - elem_midpoint[1]) * (iter->min_y - elem_midpoint[1]) + (iter->min_z - elem_midpoint[2]) * (iter->min_z - elem_midpoint[2]));
            
            /*  Check whether the disatnce equals zero, if so, we return the actual value at this position as a prediction */
            if (ApproxCompare(dist_to_min, 0.0))
            {
                ++zero_res_per_def;
                return iter->min;
            }

            /* Store the inverse weight */
            dist.push_back(dist_to_min);
        }
    }

    std::vector<double> weight_selector = dist;
    std::sort(weight_selector.begin(), weight_selector.end());

    const double max_radii = weight_selector[kNumNearestWeights];

    std::vector<double> weights = dist;

    for (auto weight_iter = weights.begin(); weight_iter != weights.end(); ++weight_iter)
    {
        if (*weight_iter >= max_radii)
        {
            *weight_iter = 0.0;
        } else
        {
            #if 0
            *weight_iter = (max_radii - *weight_iter) / (max_radii * (*weight_iter));
            *weight_iter = (*weight_iter) * (*weight_iter);
            #else
            *weight_iter = 1.0 / ((*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter));
            //*weight_iter = 1.0 / ((*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter) * (*weight_iter));
            #endif
        }
    }

    /* Compute the denominator */
    const double denominator = std::accumulate(weights.begin(), weights.end(), 0.0);

    /* Compte the numerator */
    double numerator = 0.0;

    int count = 0;
    for (auto iter = scanline.begin(); iter != scanline.end(); ++iter)
    {
        numerator += weights[count] * iter->max;
        ++count;

        if(not ApproxCompare(iter->min, std::numeric_limits<T>::max()))
        {
            numerator += weights[count] * iter->min;
        ++count;
        }
    }

    /* Compute and return the prediction */
    return numerator / denominator;
}

#endif

template <typename T>
T
ChooseMedianValuePrediction(const VectorView<CompressionValue<T>>& values)
{
    /* Convert the view to actual values of the underlying data type */
    std::vector<T> converted_values = ConvertCompressionValues<T>(values);

    std::sort(converted_values.begin(), converted_values.end());

    const int median_accessor = converted_values.size() / 2;

    return converted_values[median_accessor];
}



static int count_res_in_bins = 0;
static int count_bin_1 = 0;
static int count_bin_2 = 0;
static int count_bin_3 = 0;
static int count_bin_4 = 0;
static int count_bin_5 = 0;
static int count_bin_6 = 0;
static int count_bin_7 = 0;
static int count_bin_8 = 0;
static int count_bin_9 = 0;
static int count_bin_10 = 0;
static int outlier = 0;

static int smaller_than_128 = 0;

inline
void
CountResidualBins(const double residual, const double remaining_err)
{
    if (std::abs(residual) <= 250 * remaining_err)
    {
        ++smaller_than_128;
    }

    if (std::abs(residual) <= remaining_err)
    {
        ++count_bin_1;
    } else if (std::abs(residual) <= 2 * remaining_err)
    {
        ++count_bin_2;
    } else if (std::abs(residual) <= 3 * remaining_err)
    {
        ++count_bin_3;   
    } else if (std::abs(residual) <= 4 * remaining_err)
    {
        ++count_bin_4;
    } else if (std::abs(residual) <= 5 * remaining_err)
    {
        ++count_bin_5;
    } else if (std::abs(residual) <= 6 * remaining_err)
    {
        ++count_bin_6;
    } else if (std::abs(residual) <= 7 * remaining_err)
    {
        ++count_bin_7;
    } else if (std::abs(residual) <= 8 * remaining_err)
    {
        ++count_bin_8;
    } else if (std::abs(residual) <= 9 * remaining_err)
    {
        ++count_bin_9;
    } else if (std::abs(residual) <= 10 * remaining_err)
    {
        ++count_bin_10;
    } else
    {
        ++outlier;
    }

    ++count_res_in_bins;
}

static int count_zero_residuals = 0;
static int iiiiiiiiii = 0;
static double mse = 0.0;
static int mes_sample_size = 0;
static bool print_rmse = false;
static int ov_count = 0;
static int whole_family_zero = 0;
static int all_fam_elements_in_bins = 0;

template <typename T>
std::pair<cmc::lossy::embedded::scanline::ResidualData<T>, err_type_t>
MultiResEmbeddedAdaptData<T>::ComputeResiduals(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
    const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
    const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
    const t8_locidx_t first_incoming)
{
    ov_count += num_outgoing;
    /* The computation of residuals takes place here. Therefore, we do need to coarse approximation of the data 
     * and compute the reconstruction of the fine values */
    /* We only need to compute the residuals, if a coarsening has actually been taken place, otherwise the values are just dragged along the compression 
     * until they are a part of a family that is passed to the adaptation callback */
    if (refine != cmc::lossy::embedded::scanline::kIndicationOfCoarseningDuringAdaptation)
    {
        /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
        resdiual_order_indications_.AppendUnsetBit();

        /* Therfore there is no residual associated with this element and we return an empty residual */
        cmc::lossy::embedded::scanline::ResidualData<T> residual_data;
        residual_data.values.push_back(CompressionValue<T>());

        /* Get the remaining error that is associated with this element */
        const err_type_t previous_remaining_permitted_error = this->GetRemainingPermittedError(which_tree, first_outgoing);

        return std::make_pair(residual_data, previous_remaining_permitted_error);
    }

    // This hold only for hex elements currently 
    cmc_assert(tree_class == T8_ECLASS_HEX);
    cmc::lossy::embedded::scanline::ResidualData<T> residual_data;
    residual_data.values.reserve(num_outgoing);

    const int kScanlineLength = 53;

    /* Get a view on the remaining permitted errors */
    const VectorView<err_type_t> remaining_permitted_errors = this->GetViewOnRemainingPermittedErrors(which_tree, first_outgoing, num_outgoing);

    auto min_iter = std::min_element(remaining_permitted_errors.begin(), remaining_permitted_errors.end());
    const double min_family_remaining_err = *min_iter;


    /* Get all elements from the scanline */
    VectorView<cmc::lossy::embedded::scanline::ElementData<T>> scanline = this->GetScanline(forest_new, which_tree, tree_class, scheme, first_incoming, kScanlineLength);

    err_type_t min_remaining_permitted_error{std::numeric_limits<err_type_t>::max()};

    bool are_all_elements_zero_res = true;

    int are_alle_elements_in_encoded_bins = 0;

    for (int elem_child_idx = 0; elem_child_idx < num_outgoing; ++elem_child_idx)
    {
        /* get the child element from the previous forest */
        const t8_element_t* element = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing + elem_child_idx);

        /* Get the midpoint */
        std::vector<double> elem_midpoint(3);
        t8_forest_element_centroid (forest_old, which_tree, element, elem_midpoint.data());

        /* Get the current value corresponding to the child element */
        const CompressionValue<T> init_val = this->GetDataValueAtIndex(which_tree, first_outgoing + elem_child_idx);
        const T reinterpreted_init_val = init_val.template ReinterpretDataAs<T>();

        const T prediction = PerformScanlinePrediction<T>(scanline, elem_midpoint);
        //const T prediction = (this->GetDataValueAtIndex(which_tree, first_outgoing)).template ReinterpretDataAs<T>(); //Test one family value without costly prediction
        //const VectorView<CompressionValue<T>> previous_vals = this->GetView(which_tree, first_outgoing, num_outgoing);
        //const T prediction = ChooseMedianValuePrediction<T>(previous_vals);

        const double residual = std::abs(static_cast<double>(reinterpreted_init_val) - static_cast<double>(prediction));

        //CountResidualBins(residual, remaining_permitted_errors[elem_child_idx]);
        CountResidualBins(residual, min_family_remaining_err);

        if (residual <= remaining_permitted_errors[elem_child_idx])
        {
            ++count_zero_residuals;
            residual_data.values.push_back(CompressionValue<T>());

            if (min_remaining_permitted_error > remaining_permitted_errors[elem_child_idx] - residual)
            {
                min_remaining_permitted_error = remaining_permitted_errors[elem_child_idx] - residual;
            }
        } else
        {
            /* Compute usnigned residual and */
            residual_data.values.push_back(CompressionValue<T>());
            are_all_elements_zero_res = false;
        }
        
        if (residual > 250 * min_family_remaining_err)
        {
            are_alle_elements_in_encoded_bins += 1;
        }
        
        #if 0
        ++iiiiiiiiii;
        if (iiiiiiiiii >= 8000000 && iiiiiiiiii < 8001000)
        {
            cmc_debug_msg("Index: ", iiiiiiiiii, ", first outgoing: ", first_outgoing, ", Residual: ", reinterpreted_init_val - prediction, ", init val: ", reinterpreted_init_val, ", prediciton: ", prediction);
        }
        #endif
        
        if (first_outgoing >= 8048968 && first_outgoing <= 8049968)
        {
            cmc_debug_msg("first outgoing: ", first_outgoing, ", Residual: ", reinterpreted_init_val - prediction, ", init val: ", reinterpreted_init_val, ", prediciton: ", prediction);
        }
        mse += residual * residual;
        ++mes_sample_size;
    }
    
    if (are_all_elements_zero_res)
    {
        ++whole_family_zero;
    }

    if (are_alle_elements_in_encoded_bins <= 0)
    {
        ++all_fam_elements_in_bins;
    }
    if (not print_rmse && ov_count > 25422260)
    {
        print_rmse = true;
        const double rmse = std::sqrt(mse / static_cast<double>(mes_sample_size));
        cmc_debug_msg("RMSE is: ", rmse, " mit mse = ", mse, " und sample size: ", mes_sample_size);
        cmc_debug_msg("Num zero residuals: ", count_zero_residuals);
        cmc_debug_msg("Num zero residuals per def: ", zero_res_per_def, "; Correcte zero residuals: ", count_zero_residuals - zero_res_per_def);
        cmc_debug_msg("the number of whole familes that have zero residuals is: ", whole_family_zero);
        cmc_debug_msg("Count of families whose residuals are all in bins: ", all_fam_elements_in_bins);
        cmc_debug_msg("Residuals in bins: ", count_res_in_bins, ", Bin1: ", count_bin_1, ", Bin2: ", count_bin_2, ", Bin3: ", count_bin_3, ", Bin4: ", count_bin_4, ", Bin5: ", count_bin_5, ", Bin6: ", count_bin_6, ", Bin7: ", count_bin_7, ", Bin8: ", count_bin_8, ", Bin9: ", count_bin_9, ", Bin10: ", count_bin_10, ", Outlier: ", outlier, ", In 128-er bins: ", smaller_than_128);
    }
    return std::make_pair(std::move(residual_data), min_remaining_permitted_error);




    #if 0
    /* Get the adapted coarse element */
    const t8_element_t* elem = t8_forest_get_leaf_element_in_tree(forest_new, which_tree, first_incoming);

    /* Compute the id of the corresponding data value in the local contiguous array */
    const int local_coarse_elem_index = t8_forest_get_tree_element_offset (forest_new, which_tree) + first_incoming;

    /* Get the coarse approximation of this element */
    const CompressionValue<T> coarse_approximation = this->GetAdaptedDataValueAtIndex(local_coarse_elem_index);

    /* Gather the data of the face values */
    const FaceValues<T> face_values = this->ConstructHexFaceData(forest_new, which_tree, tree_class, scheme, elem);

    /* Get a view on the initial values */
    const VectorView<CompressionValue<T>> initial_values = this->GetView(which_tree, first_outgoing, num_outgoing);

    /* Get a view on the remaining permitted errors */
    const VectorView<err_type_t> remaining_permitted_errors = this->GetViewOnRemainingPermittedErrors(which_tree, first_outgoing, num_outgoing);

    /* Allocate the residuals */
    std::vector<CompressionValue<T>> residuals;
    residuals.reserve(num_outgoing);

    err_type_t min_remaining_permitted_error{std::numeric_limits<err_type_t>::max()};

    /* We iterate over all outgoing fine values predict their value and store the remaining residual */
    for (int elem_child_idx = 0; elem_child_idx < num_outgoing; ++elem_child_idx)
    {
        auto [is_approximation_greater, residual, remaining_error] = PredictAndComputeResidual<T>(coarse_approximation, face_values, initial_values[elem_child_idx], remaining_permitted_errors[elem_child_idx], elem_child_idx);
        //if (elem_child_idx < num_outgoing - 1)
        //{
        resdiual_order_indications_.AppendBit(is_approximation_greater);
        //}
        residuals.push_back(residual);

        /* Find the minimum remaining permitted error for the coarse element */
        if (min_remaining_permitted_error > remaining_error)
        {
            min_remaining_permitted_error = remaining_error;
        }
    }

    cmc_assert(not cmc::ApproxCompare(min_remaining_permitted_error, std::numeric_limits<err_type_t>::max()));

    return std::make_pair(cmc::lossy::embedded::scanline::ResidualData<T>(std::move(residuals)), min_remaining_permitted_error);
    #endif
}


template <typename T>
void
MultiResEmbeddedAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

#if 0

template <typename T>
cmc::lossy::embedded::scanline::ExtractionData<T>
MultiResEmbeddedAdaptData<T>::PerformExtraction(t8_forest_t forest, t8_locidx_t ltreeid, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[])
{
    /* Only calculate the arithmetic mean of the values from the family of elements */
    
    /* Get the view on the current values of this family of elements */
    const VectorView<CompressionValue<T>> values = this->GetView(ltreeid, lelement_id, num_elements);

    /* Convert the view to actual values of the underlying data type */
    const std::vector<T> converted_values = ConvertCompressionValues<T>(values);

    cmc::lossy::embedded::scanline::ElementData<T> extraction;

    int min_idx{-1};
    int max_idx{-1};

    /* Find the minimum and maximum values and ensure that min and max are not the same element */
    for (int idx = 0; idx < num_elements; ++idx)
    {
        if (extraction.min > converted_values[idx])
        {
            extraction.min = converted_values[idx];
            min_idx = idx;
        }

        if (extraction.max <= converted_values[idx])
        {
            extraction.max = converted_values[idx];
            max_idx = idx;
        }
    }

    cmc_assert(min_idx != max_idx);
    cmc_assert(min_idx >= 0 && min_idx < num_elements);
    cmc_assert(max_idx >= 0 && max_idx < num_elements);

    std::vector<double> midpoint(3);

    /* Get the midpoint of the min element */
    t8_forest_element_centroid (forest, ltreeid, elements[min_idx], midpoint.data());
    extraction.min_x = static_cast<float>(midpoint[0]);
    extraction.min_y = static_cast<float>(midpoint[1]);
    extraction.min_z = static_cast<float>(midpoint[2]);
    
    /* Get the midpoint of the max element */
    t8_forest_element_centroid (forest, ltreeid, elements[max_idx], midpoint.data());
    extraction.max_x = static_cast<float>(midpoint[0]);
    extraction.max_y = static_cast<float>(midpoint[1]);
    extraction.max_z = static_cast<float>(midpoint[2]);

    /* Store the encoding of the child id */
    extraction.min_encoding = (extraction.min_encoding << 3) || static_cast<uint16_t>(min_idx);
    extraction.min_encoding = (extraction.min_encoding << 3) || static_cast<uint16_t>(min_idx);

    return cmc::lossy::embedded::scanline::ExtractionData<T>(std::move(extraction));
}


#else

static int counter_evenodd = 0;

template <typename T>
cmc::lossy::embedded::scanline::ExtractionData<T>
MultiResEmbeddedAdaptData<T>::PerformExtraction(t8_forest_t forest, t8_locidx_t ltreeid, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[])
{
    #if 1
    /* Only calculate the arithmetic mean of the values from the family of elements */
    #if 0
    /* Get the view on the current values of this family of elements */
    const VectorView<CompressionValue<T>> values = this->GetView(ltreeid, lelement_id, num_elements);

    /* Convert the view to actual values of the underlying data type */
    const std::vector<T> converted_values = ConvertCompressionValues<T>(values);

    const T mean = InterpolateToArithmeticMean<T>(converted_values);
    //const T mean = InterpolateToMidRange<T>(converted_values);
    
    cmc::lossy::embedded::scanline::ElementData<T> extraction;

    extraction.max = mean;//converted_values[std::rand() % 8];//mean;

    std::array<t8_element_t*, 1> parent_elem{nullptr};
    scheme->element_new (tree_class, 1, parent_elem.data());
    
    /* Create the parent element for this family */
    scheme->element_get_parent (tree_class, elements[0], parent_elem[0]); 

    std::vector<double> midpoint(3);

    /* Get the midpoint of the min element */
    t8_forest_element_centroid (forest, ltreeid, parent_elem[0], midpoint.data());

    /* Get the midpoint of the max element */
    extraction.max_x = static_cast<float>(midpoint[0]);
    extraction.max_y = static_cast<float>(midpoint[1]);
    extraction.max_z = static_cast<float>(midpoint[2]);

    scheme->element_destroy(tree_class, 1, parent_elem.data());

    return cmc::lossy::embedded::scanline::ExtractionData<T>(std::move(extraction));

    #else

     /* Get the view on the current values of this family of elements */
     const VectorView<CompressionValue<T>> values = this->GetView(ltreeid, lelement_id, num_elements);

     /* Convert the view to actual values of the underlying data type */
     const std::vector<T> converted_values = ConvertCompressionValues<T>(values);

     cmc::lossy::embedded::scanline::ElementData<T> extraction;

    extraction.max = converted_values[0];


     std::vector<double> midpoint(3);

    /* Get the midpoint of the min element */
    t8_forest_element_centroid (forest, ltreeid, elements[0], midpoint.data());

    /* Get the midpoint of the max element */
    extraction.max_x = static_cast<float>(midpoint[0]);
    extraction.max_y = static_cast<float>(midpoint[1]);
    extraction.max_z = static_cast<float>(midpoint[2]);


    return cmc::lossy::embedded::scanline::ExtractionData<T>(std::move(extraction));

    #endif

    #else

    /* Get the view on the current values of this family of elements */
    const VectorView<CompressionValue<T>> values = this->GetView(ltreeid, lelement_id, num_elements);

    /* Convert the view to actual values of the underlying data type */
    const std::vector<T> converted_values = ConvertCompressionValues<T>(values);

    cmc::lossy::embedded::scanline::ElementData<T> extraction;

    int min_idx{-1};
    int max_idx{-1};

    float min_val = extraction.min ;
    float max_val = extraction.max;

    /* Find the minimum and maximum values and ensure that min and max are not the same element */
    for (int idx = 0; idx < num_elements; ++idx)
    {
        if (min_val > converted_values[idx])
        {
            min_val = converted_values[idx];
            min_idx = idx;
        }

        if (max_val <= converted_values[idx])
        {
            max_val = converted_values[idx];
            max_idx = idx;
        }
    }

    cmc_assert(min_idx != max_idx);
    cmc_assert(min_idx >= 0 && min_idx < num_elements);
    cmc_assert(max_idx >= 0 && max_idx < num_elements);

    std::vector<double> midpoint(3);

    if (counter_evenodd % 2 == 0)
    {
        /* Get the midpoint of the max element */
        t8_forest_element_centroid (forest, ltreeid, elements[max_idx], midpoint.data());
        extraction.max_x = static_cast<float>(midpoint[0]);
        extraction.max_y = static_cast<float>(midpoint[1]);
        extraction.max_z = static_cast<float>(midpoint[2]);

        extraction.max = max_val;

    } else
    {
        /* Get the midpoint of the min element */
        t8_forest_element_centroid (forest, ltreeid, elements[min_idx], midpoint.data());
        extraction.max_x = static_cast<float>(midpoint[0]);
        extraction.max_y = static_cast<float>(midpoint[1]);
        extraction.max_z = static_cast<float>(midpoint[2]);

        extraction.max = min_val;
    }

    ++counter_evenodd;

    return cmc::lossy::embedded::scanline::ExtractionData<T>(std::move(extraction));
    #endif
}

#endif


template <typename T>
cmc::lossy::embedded::scanline::UnchangedData<T>
MultiResEmbeddedAdaptData<T>::ElementStaysUnchanged(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                    const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[])
{
    cmc::lossy::embedded::scanline::ElementData<T> extraction;

    const CompressionValue<T> value = this->GetDataValueAtIndex(which_tree, lelement_id);
    extraction.max = value.template ReinterpretDataAs<T>();

    std::vector<double> midpoint(3);
    /* Get the midpoint of the current element */
    t8_forest_element_centroid (forest, which_tree, elements[0], midpoint.data());
    extraction.max_x = static_cast<float>(midpoint[0]);
    extraction.max_y = static_cast<float>(midpoint[1]);
    extraction.max_z = static_cast<float>(midpoint[2]);

    //What to do with encoding of child ids???

    return cmc::lossy::embedded::scanline::UnchangedData<T>(std::move(extraction));
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
        uint32_t signum = cmc::lossless::multi_res::util::GetSignumForEncoding(residual_flags.GetNextBit());
        if (val.GetNumberLeadingZeros() == sizeof(T) * bit_map::kCharBit)
        {
            signum = 0;
        }
        const uint32_t first_one_bit = val.GetNumberLeadingZeros();

        /* Update this symbol for encoding */
        cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>::entropy_coder_->UpdateSymbolFrequency(signum + first_one_bit);
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
    
    cmc_debug_msg("This iteration had ", zero_residuals, " zero residuals.");
    zero_residuals = 0;
    
    cmc_assert(cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>::entropy_coder_ != nullptr);

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
    cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>::entropy_coder_->Reset(std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResCompressionAlphabet<T>>());
    cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>::entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Collect the symbols and their frequencies for the entropy coder */
    CollectSymbolFrequenciesForEntropyCoding(level_byte_values);

    /* Setup the interior structure for encoding */
    cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>::entropy_coder_->SetupEncoding(this->GetMPIComm());

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
        cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(signum + first_one_bit);

        /* The first one is implicitly given by the leading zero count; therefore we set the "front bit" in order to discard the leading zeros and the following one */
        val.SetFrontBit(first_one_bit + 1);

        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            encoding.AppendBits(val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
        }
    }

    /* We set an indicaton symbol that the process local end of the values have been reached */
    cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte);

    /* Indicate that the encoding has been finished and flush all pending encodings */
    cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>::entropy_coder_->FinishEncoding();

    /* Set up the BitVector holding the encoded data for further use */
    encoding.TrimToContent();

    /* Get the local remaining significant bits */
    const uint64_t local_remaining_significant_bits_num_bytes = static_cast<uint64_t>(encoding.size());

    /* Get the local encoded LZC, respectively the encoded first "one-bit" positions */
    cmc::bit_map::BitMap local_encoded_lzc_stream = cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>::entropy_coder_->GetEncodedBitStream();
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
        cmc::bit_vector::BitVector encoded_alphabet = cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeAlphabet();
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
    cmc_debug_msg("The encoding of the entropy codes amounts to ", encoded_entropy_codes.size(), " bytes and the significant bits amount to ", encoded_data.size(), " bytes.");
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
inline cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>*
CreateMultiResEmbeddedExtractionAdaptationClass(cmc::lossy::embedded::scanline::AbstractEmbeddedByteCompressionVariable<T>* abstract_var, const CompressionSettings& settings)
{
    return new MultiResEmbeddedAdaptData<T>(abstract_var, settings);
}

template <typename T>
inline void
DestroyMultiResEmbeddedExtractionAdaptationClass(cmc::lossy::embedded::scanline::IEmbeddedCompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class EmbeddedCompressionVariable : public cmc::lossy::embedded::scanline::AbstractEmbeddedByteCompressionVariable<T>
{
public:
    EmbeddedCompressionVariable() = delete;

    EmbeddedCompressionVariable(const CompressionSettings& settings, input::Var& input_variable)
    : cmc::lossy::embedded::scanline::AbstractEmbeddedByteCompressionVariable<T>(settings, input_variable)
    {
        //cmc_assert(input_variable.IsValid());//TODO: implement 

        this->SetName(input_variable.GetName());
        this->IndicateWhetherMeshRefinementBitsWillBeStored(true);

        cmc::lossy::embedded::scanline::AbstractEmbeddedByteCompressionVariable<T>::adaptation_creator_ = CreateMultiResEmbeddedExtractionAdaptationClass<T>;
        cmc::lossy::embedded::scanline::AbstractEmbeddedByteCompressionVariable<T>::adaptation_destructor_ = DestroyMultiResEmbeddedExtractionAdaptationClass<T>;
        cmc::lossy::embedded::scanline::AbstractEmbeddedByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::EmbeddedMeshEncoder>();
    };

    CompressionSchema GetCompressionSchema() const override
    {
        return CompressionSchema::LossyEmbeddedMultiResExtractionNearestNeighborReconstruction;
    }

private:

};



}


#endif /* !LOSSY_CMC_EMBEDDED_MULTI_RES_EXTRACTION_NEAREST_NEIGHBOR_RECONSTRUCTION_COMPRESSION_HXX */
