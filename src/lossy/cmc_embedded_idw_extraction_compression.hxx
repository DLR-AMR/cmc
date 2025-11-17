#ifndef LOSSY_CMC_EMBEDDED_IDW_EXTRACTION_COMPRESSION_HXX
#define LOSSY_CMC_EMBEDDED_IDW_EXTRACTION_COMPRESSION_HXX


#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "lossy/cmc_embedded_idw_compression_variable.hxx"
#include "utilities/cmc_residual_quantization_encoder.hxx"
#include "utilities/cmc_multi_res_extraction_util.hxx"
#include "utilities/cmc_lossy_multi_res_extraction_residual_computation.hxx"
#include "utilities/cmc_multi_res_entropy_coder.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_interpolation_fn.hxx"
#include "mesh_compression/cmc_embedded_mesh_encoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>
#include <functional>
#include <type_traits>
#include <cmath>

namespace cmc::lossy::embedded::idw
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

using err_type_t = cmc::lossy::embedded::idw::err_type_t;

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
class MultiResEmbeddedAdaptData : public cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>
{
public:
    MultiResEmbeddedAdaptData() = delete;
    MultiResEmbeddedAdaptData(cmc::lossy::embedded::idw::AbstractEmbeddedByteCompressionVariable<T>* variable, const CompressionSettings& settings)
    : cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>(variable, settings) {
        cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>::entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::ResidualQuantizationEncoder<T>>();
    };

    void InitializeExtractionIteration() override;
    void FinalizingExtractionIteration() override;
    std::pair<cmc::lossy::embedded::idw::ResidualData<T>, err_type_t> ComputeResiduals(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
        const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
        const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
        const t8_locidx_t first_incoming) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const override;
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const override;

    void TestQuantizationEntropyCoding() override;

    protected:
    cmc::lossy::embedded::idw::ExtractionData<T> PerformExtraction(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                        const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[]) override;
    cmc::lossy::embedded::idw::UnchangedData<T> ElementStaysUnchanged(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
        const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[]) override;

private:
    void CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const;
    std::vector<FaceDataValues<T>> GatherElementFaceValues(t8_forest_t forest, t8_locidx_t ltreeid, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                      const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[]) const;
    std::pair<bool, CompressionValue<T>> GetFaceValue(t8_forest_t forest, const int tree_id, const t8_scheme_c *scheme, const t8_element_t* elem, const int face_idx) const;
    FaceValues<T> ConstructHexFaceData(t8_forest_t forest_new, const int which_tree, const t8_eclass_t tree_class, const t8_scheme_c *scheme, const t8_element_t* elem) const;

    std::vector<ElementData<T>> GetFaceValues(t8_forest_t forest, const int tree_id, const t8_eclass_t tree_class, const t8_scheme_c *scheme, const t8_element_t* elem) const;
    T PerformTailTruncation(const T init_value, const int which_tree, const int lelement_idx) const;

    bit_map::BitMap resdiual_order_indications_;
    int count_adaptation_step_{0};

    std::vector<uint16_t> residual_codes_;

    bit_map::BitMap resdiual_order_indications_test;
    std::vector<CompressionValue<T>> compr_residuals;

    std::vector<uint8_t> test_entropy_codes;
    std::vector<uint8_t> extracted_child_ids;
    int extract_child_id_iter{0};

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
std::vector<ElementData<T>>
MultiResEmbeddedAdaptData<T>::GetFaceValues(t8_forest_t forest, const int tree_id, const t8_eclass_t tree_class, const t8_scheme_c *scheme, const t8_element_t* elem) const
{
    const int num_faces = scheme->element_get_num_faces(tree_class, elem);

    std::vector<ElementData<T>> face_values;
    face_values.reserve(num_faces);

    t8_element_t** neighbor_leaves;
    int* dual_faces;
    int num_neighbors{0};
    t8_locidx_t* neighbor_element_indices;
    t8_eclass_t neighbor_tree_class;

    for (int face_idx = 0; face_idx < num_faces; ++face_idx)
    {
        /* Gather the face neighbor via this face */
        t8_forest_leaf_face_neighbors (forest, tree_id, elem, &neighbor_leaves, face_idx, &dual_faces, &num_neighbors,
            &neighbor_element_indices, &neighbor_tree_class, 1);

        /* Check if there is a neighboring element at the face */
        if (num_neighbors > 0)
        {
            /* Since the forest is balanced, there should be exactly one face neighbor element */
            cmc_assert(num_neighbors == 1);
            cmc_assert(neighbor_element_indices[0] < t8_forest_get_local_num_leaf_elements(forest));

            face_values.push_back(this->GetAdaptedDataValueAtIndex(neighbor_element_indices[0]));

            /* Deallocate the memory for the face neighbor construction */
            scheme->element_destroy (neighbor_tree_class, num_neighbors, neighbor_leaves);
            T8_FREE (neighbor_leaves);
            T8_FREE (neighbor_element_indices);
            T8_FREE (dual_faces);
        }
    }

    return face_values;
}

template <typename T>
T
ComputeIDWPredictionViaFaceValues(const std::vector<ElementData<T>>& face_values, const ElementData<T>& coarse_approximation, t8_forest_t forest_old, const int which_tree, const t8_element_t* child_element)
{
    /* Allocate the weights */
    std::vector<double> weights;
    weights.reserve(face_values.size());

    /* Get the midpoint of the element to predict */
    std::vector<double> elem_midpoint(3);
    t8_forest_element_centroid (forest_old, which_tree, child_element, elem_midpoint.data());

    /* The first weight corresponds to the coarse approximation of the element */
    const double dist_coarse_approx = Dist(coarse_approximation.coordinates, elem_midpoint);
    #if 1
    // exponent 16
    const double weight_coarse_approx = 1.0 / (dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx
                                               * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx);
    #else
    //exponent 10
    const double weight_coarse_approx = 1.0 / (dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx * dist_coarse_approx
        * dist_coarse_approx * dist_coarse_approx);

    #endif
    /* Gather the weights of the face values */
    for (auto face_iter = face_values.begin(); face_iter != face_values.end(); ++face_iter)
    {
        /* Get the distance to the face value's point */
        const double dist = Dist(face_iter->coordinates, elem_midpoint);

        #if 1
        // exponent 16
        /* Compute the weight */
        weights.push_back(1.0 / (dist * dist * dist * dist * dist * dist * dist * dist
                                 * dist * dist * dist * dist * dist * dist * dist * dist));
        #else
        //exponent 10
        weights.push_back(1.0 / (dist * dist * dist * dist * dist * dist * dist * dist
            * dist * dist));
        #endif
    }

    /* Compute the denominator */
    const double denominator = std::accumulate(weights.begin(), weights.end(), weight_coarse_approx);

    /* Compute the numerator */
    double numerator = weight_coarse_approx * coarse_approximation.value;

    int face_value_idx = 0;
    for (auto weight_iter = weights.begin(); weight_iter != weights.end(); ++weight_iter, ++face_value_idx)
    {
        numerator += *weight_iter * static_cast<double>(face_values[face_value_idx].value);
    }

    /* Compute and return the prediction */
    return static_cast<T>(numerator / denominator);
}









#if 0
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

#if 1
constexpr double weight_coarse_approx = 1.0 / 5.0;
constexpr double weight_face_neighbor = 1.0 / 6.0;
constexpr double weight_far_face_neighbor = 1.0 / 10.0;
#else
constexpr double weight_coarse_approx = 1.0;
constexpr double weight_face_neighbor = 0.0;
constexpr double weight_far_face_neighbor = 0.0;
#endif

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
    if (std::abs(residual) <= 512 * remaining_err)
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


#if 0
static int zero_res_per_def = 0;

template <typename T>
T
PerformScanlinePrediction(const VectorView<cmc::lossy::embedded::scanline::ElementData<T>>& scanline, const std::vector<double>& elem_midpoint)
{
    const int kNumNearestWeights = 27;
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

static size_t iiiiiii = 0;
static size_t zero_residuals = 0;

static int count_zero_residuals = 0;
static int iiiiiiiiii = 0;
static double mse = 0.0;
static int mes_sample_size = 0;
static bool print_rmse = false;
static int ov_count = 0;

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

    mse += residual_abs_err * residual_abs_err;
    ++mes_sample_size;

    CountResidualBins(residual_abs_err, remaining_permitted_error);
    
    /* Check whether the residual is smaller than the remaining error */
    if (residual_abs_err <= remaining_permitted_error)
    {
        ++zero_residuals;
        ++count_zero_residuals;
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

    ++iiiiiii;
    if (iiiiiii >= 8000000 && iiiiiii < 8001000)
    {
        const T initial_val = initial_value.template ReinterpretDataAs<T>();
        cmc_debug_msg("Index: ", iiiiiii, ", Residual: ", prediction - initial_val);
    }

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

#endif


















// Continue below! 

const double test_rel_err = 0.001;
const double test_abs_err = 0.00005;
const int num_sym_bins = 16382;
//const int num_sym_bins = 254;

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
    if (std::abs(residual) < num_sym_bins * remaining_err)
    {
        ++smaller_than_128;
    }

    //if (std::abs(residual) <= remaining_err)
    //{
    //    ++count_bin_1;
    //} else 
    if (std::abs(residual) <= 2 * remaining_err)
    {
        ++count_bin_2;
    } else if (std::abs(residual) <= 4 * remaining_err)
    {
        ++count_bin_3;   
    } else if (std::abs(residual) <= 6 * remaining_err)
    {
        ++count_bin_4;
    } else if (std::abs(residual) <= 8 * remaining_err)
    {
        ++count_bin_5;
    } else if (std::abs(residual) <= 10 * remaining_err)
    {
        ++count_bin_6;
    } else if (std::abs(residual) <= 12 * remaining_err)
    {
        ++count_bin_7;
    } else if (std::abs(residual) <= 14 * remaining_err)
    {
        ++count_bin_8;
    } else if (std::abs(residual) <= 16 * remaining_err)
    {
        ++count_bin_9;
    } else if (std::abs(residual) <= 18 * remaining_err)
    {
        ++count_bin_10;
    } else
    {
        ++outlier;
    }

    //++count_res_in_bins;
}

static int count_zero_residuals = 0;
static int iiiiiiiiii = 0;
static double mse = 0.0;
static int mes_sample_size = 0;
static bool print_rmse = false;
static int ov_count = 0;
static int whole_family_zero = 0;
static int all_fam_elements_in_bins = 0;
static int zero_res_per_def = 0;

static int count_within_bin = 0;
static int outlier_count = 0;
static size_t zero_residuals = 0;

static int count_only_zero_fams = 0;

static int smaller_than_64_bound = 0;
#if 1

template <typename T>
std::pair<cmc::lossy::embedded::idw::ResidualData<T>, err_type_t>
MultiResEmbeddedAdaptData<T>::ComputeResiduals(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
    const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
    const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
    const t8_locidx_t first_incoming)
{
    ov_count += num_outgoing;

    if (not print_rmse && ov_count > 25422260)
    {
        print_rmse = true;
        const double rmse = std::sqrt(mse / static_cast<double>(mes_sample_size));
        cmc_debug_msg("RMSE is: ", rmse, " mit mse = ", mse, " und sample size: ", mes_sample_size);
        cmc_debug_msg("Num zero residuals: ", count_zero_residuals);
        cmc_debug_msg("Num zero residuals per def: ", zero_res_per_def, "; Correcte zero residuals: ", count_zero_residuals - zero_res_per_def);
        //cmc_debug_msg("the number of whole familes that have zero residuals is: ", whole_family_zero);
        //cmc_debug_msg("Count of families whose residuals are all in bins: ", all_fam_elements_in_bins);
        //cmc_debug_msg("Residuals in bins: ", count_res_in_bins, ", In ", num_sym_bins,"-er bins: ", smaller_than_128, ", Outlier: ", count_res_in_bins - smaller_than_128);
        cmc_debug_msg("Residuals within bins are: ", count_within_bin, ", In ", num_sym_bins,"-er bins: ", ", Outlier: ", outlier_count);
        //cmc_err_msg("End here currently for tests ");
        cmc_debug_msg("Residuals in bins: ", count_res_in_bins, ", Bin1: ", count_bin_1, ", Bin2: ", count_bin_2, ", Bin3: ", count_bin_3, ", Bin4: ", count_bin_4, ", Bin5: ", count_bin_5, ", Bin6: ", count_bin_6, ", Bin7: ", count_bin_7, ", Bin8: ", count_bin_8, ", Bin9: ", count_bin_9, ", Bin10: ", count_bin_10, ", Outlier: ", outlier, ", In 128-er bins: ", smaller_than_128);
        cmc_debug_msg("Count only zero fams: ", count_only_zero_fams);
        cmc_debug_msg("Smaller than 64 bound: ", smaller_than_64_bound);
    }

    /* The computation of residuals takes place here. Therefore, we do need to coarse approximation of the data 
     * and compute the reconstruction of the fine values */
    /* We only need to compute the residuals, if a coarsening has actually been taken place, otherwise the values are just dragged along the compression 
     * until they are a part of a family that is passed to the adaptation callback */
    if (refine != cmc::lossy::embedded::idw::kIndicationOfCoarseningDuringAdaptation)
    {
        /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
        resdiual_order_indications_.AppendUnsetBit();

        /* Therefore there is no residual associated with this element and we return an empty residual */
        cmc::lossy::embedded::idw::ResidualData<T> residual_data;
        residual_data.values.push_back(CompressionValue<T>());

        /* Get the remaining error that is associated with this element */
        const err_type_t previous_remaining_permitted_error = this->GetRemainingPermittedError(which_tree, first_outgoing);

        return std::make_pair(residual_data, previous_remaining_permitted_error);
    }

    cmc::lossy::embedded::idw::ResidualData<T> residual_data;

    /* Get the adapted coarse element */
    const t8_element_t* coarse_elem = t8_forest_get_leaf_element_in_tree(forest_new, which_tree, first_incoming);

    /* Compute the id of the corresponding data value in the local contiguous array */
    const int local_coarse_elem_index = t8_forest_get_tree_element_offset (forest_new, which_tree) + first_incoming;

    /* Get the coarse approximation of this element */
    const ElementData<T>& coarse_approximation = this->GetAdaptedDataValueAtIndex(local_coarse_elem_index);

    /* Get the face neighbors */
    std::vector<ElementData<T>> face_values = this->GetFaceValues(forest_new, which_tree, tree_class, scheme, coarse_elem);

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (forest_old, which_tree) + first_outgoing;

    /* Get a view on the initial values */
    const VectorView<T> initial_values = this->GetViewOnInitialData(start_index, num_outgoing);

    /* Get a view on the remaining permitted errors */
    const VectorView<err_type_t> remaining_permitted_errors = this->GetViewOnRemainingPermittedErrors(which_tree, first_outgoing, num_outgoing);

    auto min_iter = std::min_element(std::next(remaining_permitted_errors.begin()), remaining_permitted_errors.end());
    const double min_family_remaining_err = *min_iter;
    bool are_all_elements_zero_res = true; 
    int are_alle_elements_in_encoded_bins = 0;

    /* Allocate the residuals */
    std::vector<CompressionValue<T>> residuals;
    residuals.reserve(num_outgoing);

    //err_type_t min_remaining_permitted_error{std::numeric_limits<err_type_t>::max()};
    err_type_t min_remaining_permitted_error{0.0};

    #if 1

    const int extracted_child_id = extracted_child_ids[extract_child_id_iter];
    ++extract_child_id_iter;

    /* Compute Predictions */
    std::vector<T> predictions;
    //predictions.push_back(initial_values[0]);

    std::vector<double> comp_residuals;
    //comp_residuals.push_back(0.0);
    //++zero_res_per_def;

    /* Compute the predictions */
    for (int elem_child_idx = 0; elem_child_idx < num_outgoing; ++elem_child_idx)
    {
        if (extracted_child_id == elem_child_idx)
        {
            predictions.push_back(initial_values[extracted_child_id]);
            comp_residuals.push_back(0.0);
            //residual_codes_.push_back(8191);
            ++zero_res_per_def;
            continue;
        }

        /* Get the child element from the previous forest */
        const t8_element_t* element = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing + elem_child_idx);

        /* Predict the child's value */
        predictions.push_back(ComputeIDWPredictionViaFaceValues<T>(face_values, coarse_approximation, forest_old, which_tree, element));

        /* Compute Residual */
        comp_residuals.push_back(static_cast<double>(initial_values[elem_child_idx]) - static_cast<double>(predictions.back()));

        ++iiiiiiiiii;
        if (iiiiiiiiii >= 8000000 && iiiiiiiiii < 8001000)
        {
            cmc_debug_msg("Index: ", iiiiiiiiii, ", first outgoing: ", first_outgoing, ", Residual: ", comp_residuals.back(), ", init val: ", initial_values[elem_child_idx], ", prediciton: ", predictions.back(), ", extracted cid: ", extracted_child_id);
        }

        /* Test how the prediction changes, when we add the reconstructed value */
        #if 0
        //NOt sure if correct and/or worthy

        /* Get the midpoint */
        std::vector<double> elem_midpoint(3);
        t8_forest_element_centroid (forest_old, which_tree, element, elem_midpoint.data());

        ElementData<T> new_value;
        new_value.value = predictions.back();
        new_value.coordinates = PointData(static_cast<float>(elem_midpoint[0]),static_cast<float>(elem_midpoint[1]),static_cast<float>(elem_midpoint[2]));
        
        T correction;
        const double quantization_step = min_family_remaining_err;
        if (std::abs(comp_residuals[elem_child_idx]) <= remaining_permitted_errors[elem_child_idx])
        {
            face_values.push_back(new_value);
        } else if (std::abs(comp_residuals[elem_child_idx]) <= num_sym_bins * quantization_step)
        {
            const bool is_negative = comp_residuals[elem_child_idx] < 0.0;
            const uint16_t bin = (static_cast<uint16_t>(std::floor(std::abs(comp_residuals[elem_child_idx]) / (2 * quantization_step))) + 1); 
            correction = (is_negative ? -1 : 1) * ((bin * 2 - 1) * quantization_step);
            new_value.value += correction;
        face_values.push_back(new_value);
        } else
        {
            //const bool is_negative = comp_residuals[elem_child_idx] < 0.0;
            //correction = std::abs(predictions.back() - initial_values[elem_child_idx]) - test_rel_err * initial_values[elem_child_idx];
        }

        #endif

    }

    /* Convert the prediciton to absolute values */
    for (auto pred_iter = predictions.begin(); pred_iter != predictions.end(); ++pred_iter)
    {
        *pred_iter = std::abs(*pred_iter);
    }

    /* Check the min elem quantization step */
    //auto min_elem = std::min_element(predictions.begin(), predictions.end());
    //const double quantization_step = (*min_elem) * test_rel_err;
    ////const double quantization_step = min_family_remaining_err;
    //if (iiiiiiiiii >= 8000000 && iiiiiiiiii < 8001000)
    //cmc_debug_msg("Quantization step is: ", quantization_step);

    //const double quantization_step = initial_values[extracted_child_id] * test_rel_err;
    //ABS error below
    const double quantization_step = test_abs_err;
    bool only_zero_residuals = true;


    /* Sort the residuals into the quantization bins */
    for (int elem_child_idx = 0; elem_child_idx < num_outgoing; ++elem_child_idx)
    {
        //const double quantization_step = test_rel_err * std::abs(initial_values[elem_child_idx]); //Test when we know the exact permitted error

        #if 0
        //Test SZ quantization computaion
        uint16_t bin = std::floor((std::log10(std::abs(initial_values[elem_child_idx])) - std::log10(std::abs(predictions[elem_child_idx]))) / (2*std::log10(1+test_rel_err)) + 0.5);
        if (bin < 8192)
        {
            ++count_within_bin;
        } else
        {
            ++outlier_count;
            bin = 8192;
        }
        residual_codes_.push_back(bin);

        #endif

        #if 0
        if (std::abs(comp_residuals[elem_child_idx]) <= remaining_permitted_errors[elem_child_idx])
        {
            ++count_zero_residuals;
            ++count_within_bin;

            ++count_bin_1;
            ++count_bin_2;

            if (elem_child_idx != extracted_child_id)
            {
                residual_codes_.push_back(0);
            }
            if (iiiiiiiiii >= 8000000 && iiiiiiiiii < 8001000)
            cmc_debug_msg("Residual is in bin zero");
        } else if (std::abs(comp_residuals[elem_child_idx]) >= num_sym_bins * quantization_step)
        {
            /* If it is an outlier */
            residual_codes_.push_back(8192);
            ++outlier_count;
            if (iiiiiiiiii >= 8000000 && iiiiiiiiii < 8001000)
            cmc_debug_msg("Outlier Residual is in bin: ", 8192);
        } else
        {
            ++count_within_bin;
            const bool is_negative = comp_residuals[elem_child_idx] < 0.0;
            const uint16_t bin = (is_negative ? 8192 : 0) + (static_cast<uint16_t>(std::floor(std::abs(comp_residuals[elem_child_idx]) / (2 * quantization_step))) + 1);
            //const uint16_t bin = (static_cast<uint16_t>(std::floor(std::abs(comp_residuals[elem_child_idx]) / (2 * quantization_step))) + 1); //Only positive bins
            residual_codes_.push_back(bin);
            if (iiiiiiiiii >= 8000000 && iiiiiiiiii < 8001000)
            cmc_debug_msg("Residual is in bin: ", bin);
        }

        if (elem_child_idx != extracted_child_id)
        CountResidualBins(comp_residuals[elem_child_idx], quantization_step);

        /* just for testing, push back anything */
        residual_data.values.push_back(CompressionValue<T>());
        ++count_res_in_bins;

        mse += comp_residuals[elem_child_idx] * comp_residuals[elem_child_idx];
        ++mes_sample_size;

        #endif


        #if 1

        if (extracted_child_id == elem_child_idx)
        {
            ++zero_res_per_def;
            ++count_zero_residuals;
            continue;
        }

        if (std::abs(comp_residuals[elem_child_idx]) <= 64 * test_abs_err)
        {
            ++smaller_than_64_bound;
        }
        if (std::abs(comp_residuals[elem_child_idx]) <= remaining_permitted_errors[elem_child_idx])
        {
            ++count_zero_residuals;
        } else
        {
            only_zero_residuals = false;
            const bool is_negative = comp_residuals[elem_child_idx] < 0.0;
            //const uint16_t bin = static_cast<uint16_t>(std::floor(std::abs(comp_residuals[elem_child_idx]) / (2 * quantization_step) + 0.5));
            //ABS error below
            const uint16_t bin = static_cast<uint16_t>(std::floor(std::abs(comp_residuals[elem_child_idx]) / (2 * test_abs_err) + 0.5));
            if (bin < 2045)
            {
                residual_codes_.push_back((is_negative ? 2045 : 0)  + bin);
            } else
            {
                ++outlier_count;
                residual_codes_.push_back(2045);
            }
        }

        /* just for testing, push back anything */
        residual_data.values.push_back(CompressionValue<T>());
        ++count_res_in_bins;

        mse += comp_residuals[elem_child_idx] * comp_residuals[elem_child_idx];
        ++mes_sample_size;

        #endif
    }

    if (only_zero_residuals)
    {
        ++count_only_zero_fams;
    }
    #endif


    #if 0
    /* We iterate over all outgoing fine values predict their value and store the remaining residual */
    for (int elem_child_idx = 0; elem_child_idx < num_outgoing; ++elem_child_idx)
    {
        T prediction;
        double residual;

        /* get the child element from the previous forest */
        const t8_element_t* element = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing + elem_child_idx);

        /* Get the current value corresponding to the child element */
        const T init_val = initial_values[elem_child_idx];

        if (elem_child_idx > 0)
        {
            prediction = ComputeIDWPredictionViaFaceValues<T>(face_values, coarse_approximation, forest_old, which_tree, element);

            residual = std::abs(static_cast<double>(init_val) - static_cast<double>(prediction));
        } else
        {
            /* The first child is always taken as it is, therefore it is exact */
            prediction = initial_values[0];
            residual = 0.0;
            ++zero_res_per_def;
        }

        //const double family_err = std::abs(test_rel_err * coarse_approximation.value);
        //const double family_err = std::abs(test_rel_err * static_cast<double>(prediction));
        const double family_err = min_family_remaining_err; 

        /* Store residual bin */
        if (elem_child_idx > 0)
        {
            if (std::abs(residual) <= remaining_permitted_errors[elem_child_idx])
            {
                /* If the prediction is good enough */
                residual_codes_.push_back(0);
            } else if (std::abs(residual) >= num_sym_bins * family_err)
            {
                /* If it is an outlier */
                residual_codes_.push_back(8192);
            } else
            {
                const bool is_negative = (static_cast<double>(init_val) - static_cast<double>(prediction) < 0.0);
                const uint16_t bin = (is_negative ? 8192 : 0) + (static_cast<uint16_t>(std::floor(std::abs(residual) / (2 * family_err))) + 1); 
                residual_codes_.push_back(bin);
            }
        }

        if (std::abs(residual) <= remaining_permitted_errors[elem_child_idx])
        {
            ++count_bin_1;
        } else
        {
            //CountResidualBins(residual, std::abs(test_rel_err * coarse_approximation.value));
            //CountResidualBins(residual, min_family_remaining_err);
            CountResidualBins(residual, min_remaining_permitted_error);
        }

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
        
        if (residual >= num_sym_bins * min_family_remaining_err)
        {
            are_alle_elements_in_encoded_bins += 1;
        }
        
        ++iiiiiiiiii;
        if (iiiiiiiiii >= 8000000 && iiiiiiiiii < 8001000)
        {
            cmc_debug_msg("Index: ", iiiiiiiiii, ", first outgoing: ", first_outgoing, ", Residual: ", init_val - prediction, ", init val: ", init_val, ", prediciton: ", prediction);
        }
    
        mse += residual * residual;
        ++mes_sample_size;
    }

    //cmc_assert(not cmc::ApproxCompare(min_remaining_permitted_error, std::numeric_limits<err_type_t>::max()));

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

    #endif


    return std::make_pair(std::move(residuals), min_remaining_permitted_error);
}


#else



template <typename T>
std::pair<cmc::lossy::embedded::idw::ResidualData<T>, err_type_t>
MultiResEmbeddedAdaptData<T>::ComputeResiduals(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
    const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
    const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
    const t8_locidx_t first_incoming)
{
    ov_count += num_outgoing;

    if (not print_rmse && ov_count > 25422260)
    {
        print_rmse = true;
        const double rmse = std::sqrt(mse / static_cast<double>(mes_sample_size));
        cmc_debug_msg("RMSE is: ", rmse, " mit mse = ", mse, " und sample size: ", mes_sample_size);
        cmc_debug_msg("Num zero residuals: ", count_zero_residuals);
        cmc_debug_msg("Num zero residuals per def: ", zero_res_per_def, "; Correcte zero residuals: ", count_zero_residuals - zero_res_per_def);
        //cmc_debug_msg("the number of whole familes that have zero residuals is: ", whole_family_zero);
        //cmc_debug_msg("Count of families whose residuals are all in bins: ", all_fam_elements_in_bins);
        //cmc_debug_msg("Residuals in bins: ", count_res_in_bins, ", In ", num_sym_bins,"-er bins: ", smaller_than_128, ", Outlier: ", count_res_in_bins - smaller_than_128);
        cmc_debug_msg("Residuals within bins are: ", count_within_bin, ", In ", num_sym_bins,"-er bins: ", ", Outlier: ", outlier_count);
        //cmc_err_msg("End here currently for tests ");
        cmc_debug_msg("Residuals in bins: ", count_res_in_bins, ", Bin1: ", count_bin_1, ", Bin2: ", count_bin_2, ", Bin3: ", count_bin_3, ", Bin4: ", count_bin_4, ", Bin5: ", count_bin_5, ", Bin6: ", count_bin_6, ", Bin7: ", count_bin_7, ", Bin8: ", count_bin_8, ", Bin9: ", count_bin_9, ", Bin10: ", count_bin_10, ", Outlier: ", outlier, ", In 128-er bins: ", smaller_than_128);
    }

    /* The computation of residuals takes place here. Therefore, we do need to coarse approximation of the data 
     * and compute the reconstruction of the fine values */
    /* We only need to compute the residuals, if a coarsening has actually been taken place, otherwise the values are just dragged along the compression 
     * until they are a part of a family that is passed to the adaptation callback */
    if (refine != cmc::lossy::embedded::idw::kIndicationOfCoarseningDuringAdaptation)
    {
        /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
        resdiual_order_indications_.AppendUnsetBit();

        /* Therefore there is no residual associated with this element and we return an empty residual */
        cmc::lossy::embedded::idw::ResidualData<T> residual_data;
        residual_data.values.push_back(CompressionValue<T>());

        /* Get the remaining error that is associated with this element */
        const err_type_t previous_remaining_permitted_error = this->GetRemainingPermittedError(which_tree, first_outgoing);

        return std::make_pair(residual_data, previous_remaining_permitted_error);
    }

    cmc::lossy::embedded::idw::ResidualData<T> residual_data;

    /* Get the adapted coarse element */
    const t8_element_t* coarse_elem = t8_forest_get_leaf_element_in_tree(forest_new, which_tree, first_incoming);

    /* Compute the id of the corresponding data value in the local contiguous array */
    const int local_coarse_elem_index = t8_forest_get_tree_element_offset (forest_new, which_tree) + first_incoming;

    /* Get the coarse approximation of this element */
    const ElementData<T>& coarse_approximation = this->GetAdaptedDataValueAtIndex(local_coarse_elem_index);

    /* Get the face neighbors */
    std::vector<ElementData<T>> face_values = this->GetFaceValues(forest_new, which_tree, tree_class, scheme, coarse_elem);

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (forest_old, which_tree) + first_outgoing;

    /* Get a view on the initial values */
    const VectorView<T> initial_values = this->GetViewOnInitialData(start_index, num_outgoing);

    /* Get a view on the remaining permitted errors */
    const VectorView<err_type_t> remaining_permitted_errors = this->GetViewOnRemainingPermittedErrors(which_tree, first_outgoing, num_outgoing);

    auto min_iter = std::min_element(std::next(remaining_permitted_errors.begin()), remaining_permitted_errors.end());
    const double min_family_remaining_err = *min_iter;
    bool are_all_elements_zero_res = true;
    int are_alle_elements_in_encoded_bins = 0;

    /* Allocate the residuals */
    std::vector<CompressionValue<T>> residuals;
    residuals.reserve(num_outgoing);

    //err_type_t min_remaining_permitted_error{std::numeric_limits<err_type_t>::max()};
    err_type_t min_remaining_permitted_error{0.0};

    #if 1
    /* Compute Predictions */
    std::vector<T> predictions;
    predictions.push_back(initial_values[0]);

    std::vector<double> comp_residuals;
    comp_residuals.push_back(0.0);
    ++zero_res_per_def;

    /* Compute the predictions */
    for (int elem_child_idx = 1; elem_child_idx < num_outgoing; ++elem_child_idx)
    {
        /* Get the child element from the previous forest */
        const t8_element_t* element = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing + elem_child_idx);

        /* Predict the child's value */
        predictions.push_back(ComputeIDWPredictionViaFaceValues<T>(face_values, coarse_approximation, forest_old, which_tree, element));

        /* Compute Residual */
        comp_residuals.push_back(static_cast<double>(initial_values[elem_child_idx]) - static_cast<double>(predictions.back()));


        /* Test previous multi res encoding */
        std::vector<PermittedError> permitted_errors{PermittedError{cmc::CompressionCriterion::AbsoluteErrorThreshold, remaining_permitted_errors[elem_child_idx]}};
        const double previous_absolute_error = 0.0;
        auto [is_approximation_greater, residual, absolute_inaccuracy] = cmc::lossy::multi_res::ComputeMinimalResidual<T>(predictions.back(), CompressionValue<T>(initial_values[elem_child_idx]), permitted_errors, previous_absolute_error);
        resdiual_order_indications_test.AppendBit(is_approximation_greater);
        compr_residuals.push_back(residual);
        
        
        
        residuals.push_back(residual);

        ++iiiiiiiiii;
        if (iiiiiiiiii >= 8000000 && iiiiiiiiii < 8001000)
        {
            cmc_debug_msg("Index: ", iiiiiiiiii, ", first outgoing: ", first_outgoing, ", Residual: ", comp_residuals.back(), ", init val: ", initial_values[elem_child_idx], ", prediciton: ", predictions.back());
        }

    }

    #endif


    return std::make_pair(std::move(residuals), min_remaining_permitted_error);
}



#endif



template <typename T>
void
CollectSymbolFrequenciesForEntropyCodingTest(std::unique_ptr<cmc::entropy_coding::IByteCompressionEntropyCoder>& entropy_coder_, bit_map::BitMap& resdiual_order_indications_, std::vector<CompressionValue<T>>& level_byte_values) 
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
        entropy_coder_->UpdateSymbolFrequency(signum + first_one_bit);


        /** Afterwards, we collect the trailing zero count (TZC) **/
        /* The collection of the TZC is only done, when the LZC does not fully describe the residual */
        if (first_one_bit >= sizeof(T) * bit_map::kCharBit - 1)
        {
            continue;
        }

        /* Otherwise, compute the TZC */
        const uint32_t last_one_bit = val.GetNumberTrailingZeros();

        /* Update this symbol for encoding */
        entropy_coder_->UpdateSymbolFrequency(last_one_bit);       
    }
}

template <typename T>
void
MultiResEmbeddedAdaptData<T>::TestQuantizationEntropyCoding()
{

    #if 1
    /* Set up an entropy coder */
    std::unique_ptr<cmc::entropy_coding::IByteCompressionEntropyCoder> entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::ResidualQuantizationEncoder<T>>();
    entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Setup the dictionary */
    for (auto res_iter = residual_codes_.begin(); res_iter != residual_codes_.end(); ++res_iter)
    {
        entropy_coder_->UpdateSymbolFrequency(*res_iter);
    }

    /* Indicate to start the encoding process */
    entropy_coder_->SetupEncoding(MPI_COMM_SELF);

    for (auto res_iter = residual_codes_.begin(); res_iter != residual_codes_.end(); ++res_iter)
    {
        entropy_coder_->EncodeSymbol(*res_iter);
    }

    entropy_coder_->FinishEncoding();

    cmc::bit_map::BitMap local_encoded_lzc_stream = entropy_coder_->GetEncodedBitStream();

    //FILE* file_out = fopen("example_quantized_res_no_entropy_coding.bin", "wb");
    //fwrite(residual_codes_.data(), sizeof(uint16_t), residual_codes_.size(), file_out);
    //fclose(file_out);

    cmc_debug_msg("Local encoded residuals bins size: ", local_encoded_lzc_stream.size_bytes());

    cmc_err_msg("end here for testing");



    #else

    cmc::bit_vector::BitVector encoding;
    encoding.Reserve(3 * compr_residuals.size());

    /* Reset the entropy coder and initialize the alphabet */
    std::unique_ptr<cmc::entropy_coding::IByteCompressionEntropyCoder> entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResEncoder<T>>();
    entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Collect the symbols and their frequencies for the entropy coder */
    CollectSymbolFrequenciesForEntropyCodingTest<T>(entropy_coder_, resdiual_order_indications_test, compr_residuals);

    entropy_coder_->SetupEncoding(MPI_COMM_SELF);

    /* Get a view on the residual flags */
    bit_map::BitMapView residual_flags(resdiual_order_indications_test);
    
    std::vector<uint8_t> entropy_codes_front;
    std::vector<uint8_t> entropy_codes_end;
    
    for (auto val_iter = compr_residuals.begin(); val_iter != compr_residuals.end(); ++val_iter)
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
        entropy_coder_->EncodeSymbol(signum + first_one_bit);
        entropy_codes_front.push_back(first_one_bit);

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
        entropy_coder_->EncodeSymbol(last_one_bit);
        entropy_codes_end.push_back(last_one_bit);

        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            encoding.AppendBits(val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
        }
    }

    /* We set an indicaton symbol that the process local end of the values have been reached */
    entropy_coder_->EncodeSymbol(entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte);

    /* Indicate that the encoding has been finished and flush all pending encodings */
    entropy_coder_->FinishEncoding();

    /* Set up the BitVector holding the encoded data for further use */
    encoding.TrimToContent();

     /* Get the local remaining significant bits */
     const uint64_t local_remaining_significant_bits_num_bytes = static_cast<uint64_t>(encoding.size());

     /* Get the local encoded LZC, respectively the encoded first "one-bit" positions */
     cmc::bit_map::BitMap local_encoded_lzc_stream = entropy_coder_->GetEncodedBitStream();
     const uint64_t local_encoded_lzc_stream_num_bytes = static_cast<uint64_t>(local_encoded_lzc_stream.size_bytes());

    FILE* file_out1 = fopen("lzc_encoding_front.bin", "wb");
    fwrite(entropy_codes_front.data(), sizeof(uint8_t), entropy_codes_front.size(), file_out1);
    fclose(file_out1);

    FILE* file_out2 = fopen("lzc_encoding_end.bin", "wb");
    fwrite(entropy_codes_end.data(), sizeof(uint8_t), entropy_codes_end.size(), file_out2);
    fclose(file_out2);

    cmc_debug_msg("Size Entropy Codes: ", local_encoded_lzc_stream_num_bytes, ", size of remaining significant bits: ", local_remaining_significant_bits_num_bytes);

    cmc_err_msg("End here for testing");

    #endif
}












#if 0

#if 1

template <typename T>
std::pair<cmc::lossy::embedded::idw::ResidualData<T>, err_type_t>
MultiResEmbeddedAdaptData<T>::ComputeResiduals(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
    const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
    const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
    const t8_locidx_t first_incoming)
{
    ov_count += num_outgoing;

    if (not print_rmse && ov_count > 25422260)
    {
        print_rmse = true;
        const double rmse = std::sqrt(mse / static_cast<double>(mes_sample_size));
        cmc_debug_msg("RMSE is: ", rmse, " mit mse = ", mse, " und sample size: ", mes_sample_size);
        cmc_debug_msg("Num zero residuals: ", count_zero_residuals);
        cmc_debug_msg("Residuals in bins: ", count_res_in_bins, ", Bin1: ", count_bin_1, ", Bin2: ", count_bin_2, ", Bin3: ", count_bin_3, ", Bin4: ", count_bin_4, ", Bin5: ", count_bin_5, ", Bin6: ", count_bin_6, ", Bin7: ", count_bin_7, ", Bin8: ", count_bin_8, ", Bin9: ", count_bin_9, ", Bin10: ", count_bin_10, ", Outlier: ", outlier, ", In 128-er bins: ", smaller_than_128);

        cmc_err_msg("End here currently for tests ");
    }

    /* The computation of residuals takes place here. Therefore, we do need to coarse approximation of the data 
     * and compute the reconstruction of the fine values */
    /* We only need to compute the residuals, if a coarsening has actually been taken place, otherwise the values are just dragged along the compression 
     * until they are a part of a family that is passed to the adaptation callback */
    if (refine != cmc::lossy::embedded::idw::kIndicationOfCoarseningDuringAdaptation)
    {
        /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
        resdiual_order_indications_.AppendUnsetBit();

        /* Therfore there is no residual associated with this element and we return an empty residual */
        cmc::lossy::embedded::idw::ResidualData<T> residual_data;
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

    return std::make_pair(cmc::lossy::embedded::idw::ResidualData<T>(std::move(residuals)), min_remaining_permitted_error);
}

#else


template <typename T>
std::pair<cmc::lossy::embedded::idw::ResidualData<T>, err_type_t>
MultiResEmbeddedAdaptData<T>::ComputeResiduals(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
    const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
    const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
    const t8_locidx_t first_incoming)
{
    ov_count += num_outgoing;

    if (not print_rmse && ov_count > 25422260)
    {
        print_rmse = true;
        const double rmse = std::sqrt(mse / static_cast<double>(mes_sample_size));
        cmc_debug_msg("RMSE is: ", rmse, " mit mse = ", mse, " und sample size: ", mes_sample_size);
        cmc_debug_msg("Num zero residuals: ", count_zero_residuals);
        cmc_debug_msg("Residuals in bins: ", count_res_in_bins, ", Bin1: ", count_bin_1, ", Bin2: ", count_bin_2, ", Bin3: ", count_bin_3, ", Bin4: ", count_bin_4, ", Bin5: ", count_bin_5, ", Bin6: ", count_bin_6, ", Bin7: ", count_bin_7, ", Bin8: ", count_bin_8, ", Bin9: ", count_bin_9, ", Bin10: ", count_bin_10, ", Outlier: ", outlier, ", In 128-er bins: ", smaller_than_128);

        cmc_err_msg("End here currently for tests ");
    }

    /* The computation of residuals takes place here. Therefore, we do need to coarse approximation of the data 
     * and compute the reconstruction of the fine values */
    /* We only need to compute the residuals, if a coarsening has actually been taken place, otherwise the values are just dragged along the compression 
     * until they are a part of a family that is passed to the adaptation callback */
    if (refine != cmc::lossy::embedded::idw::kIndicationOfCoarseningDuringAdaptation)
    {
        /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
        resdiual_order_indications_.AppendUnsetBit();

        /* Therfore there is no residual associated with this element and we return an empty residual */
        cmc::lossy::embedded::idw::ResidualData<T> residual_data;
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
         /* get the child element from the previous forest */
         const t8_element_t* element = t8_forest_get_leaf_element_in_tree (forest_old, which_tree, first_outgoing + elem_child_idx);

         /* Get the midpoint */
         std::vector<double> elem_midpoint(3);
         t8_forest_element_centroid (forest_old, which_tree, element, elem_midpoint.data());

        const T prediction = PerformIDWPrediction(face_values, elem_midpoint);


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
        
        
        ++iiiiiiiiii;
        if (iiiiiiiiii >= 8000000 && iiiiiiiiii < 8001000)
        {
            cmc_debug_msg("Index: ", iiiiiiiiii, ", Residual: ", reinterpreted_init_val - prediction, ", init val: ", reinterpreted_init_val, ", prediciton: ", prediction);
        }
    
        mse += residual * residual;
        ++mes_sample_size;
    }


    return std::make_pair(cmc::lossy::embedded::idw::ResidualData<T>(std::move(residuals)), min_remaining_permitted_error);
}


#endif

#endif

template <typename T>
void
MultiResEmbeddedAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

template <typename T>
T
ChooseMedianValuePrediction(const VectorView<T>& values)
{
    std::vector<T> converted_values(8);
    int idx = 0;
    for (auto iter = values.begin(); iter != values.end(); ++iter, ++idx)
    {
        converted_values[idx] = *iter;
    }

    std::sort(converted_values.begin(), converted_values.end());

    const int median_accessor = converted_values.size() / 2;

    return converted_values[median_accessor];
}

static int counter_evenodd = 0;
static bool switch_min_max = false;

#if 0

template <typename T>
cmc::lossy::embedded::idw::ExtractionData<T>
MultiResEmbeddedAdaptData<T>::PerformExtraction(t8_forest_t forest, t8_locidx_t ltreeid, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[])
{
    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (forest, ltreeid) + lelement_id;

    /* Get the view on the current values of this family of elements */
    const VectorView<T> values = this->GetViewOnInitialData(start_index, num_elements);

    cmc::lossy::embedded::idw::ElementData<T> extraction;

    #if 0
    
    //extraction.value = ChooseMedianValuePrediction<T>(values);
    //std::vector<T> converted_values(8);
    //int idx = 0;
    //for (auto iter = values.begin(); iter != values.end(); ++iter, ++idx)
    //{
    //    converted_values[idx] = *iter;
    //}
    //
    //const T mean = InterpolateToArithmeticMean<T>(converted_values);
    //extraction.value = mean;

    if (counter_evenodd % 3 == 0)
    {
        switch_min_max = not switch_min_max;
    }

    if (switch_min_max)
    {
        auto max_elem = std::max_element(values.begin(), values.end());
        extraction.value = *max_elem;
    } else
    {
        auto min_elem = std::min_element(values.begin(), values.end());
        extraction.value = *min_elem;
    }

    ++counter_evenodd;

    #else
    extraction.value = values[0];

    #endif


    #if 0
    /* Test value at parent element */
    std::array<t8_element_t*, 1> parent_elem{nullptr};
    scheme->element_new (tree_class, 1, parent_elem.data());
    scheme->element_get_parent (tree_class, elements[0], parent_elem[0]); 

    std::vector<double> midpoint(3);
   /* Get the midpoint of the min element */
   t8_forest_element_centroid (forest, ltreeid, parent_elem[0], midpoint.data());

   scheme->element_destroy(tree_class, 1, parent_elem.data());
   #else

   std::vector<double> midpoint(3);
   /* Get the midpoint of the min element */
   t8_forest_element_centroid (forest, ltreeid, elements[0], midpoint.data());
   #endif

   /* Get the midpoint of the max element */
   extraction.coordinates.x = static_cast<float>(midpoint[0]);
   extraction.coordinates.y = static_cast<float>(midpoint[1]);
   extraction.coordinates.z = static_cast<float>(midpoint[2]);


   return cmc::lossy::embedded::idw::ExtractionData<T>(std::move(extraction));
}

#else


template<typename T>
CompressionValue<T>
GetMaximumTailToggledValue(const std::vector<PermittedError>& permitted_errors, const CompressionValue<T>& initial_serialized_value, const T& missing_value)
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
        const cmc::lossy::multi_res::ErrorCompliance error_evaluation = cmc::lossy::multi_res::IsValueErrorCompliant(permitted_errors, initial_val, reinterpreted_value);

        if (!error_evaluation.is_error_threshold_satisfied)
        {
            /* Revert the changes to the value */
            toggled_value = save_previous_value;
            is_toogling_progressing = false;
        }

        if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>)
        {
            if (not std::isfinite(reinterpreted_value))
            {
                toggled_value = save_previous_value;
                is_toogling_progressing = false;
            }
        }

        ++iteration_count;
    }

    return toggled_value;
}


template<typename T>
CompressionValue<T>
GetMaximumTailClearedValue(const std::vector<PermittedError>& permitted_errors, const CompressionValue<T>& initial_serialized_value, const T& missing_value)
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
        const cmc::lossy::multi_res::ErrorCompliance error_evaluation = cmc::lossy::multi_res::IsValueErrorCompliant(permitted_errors, initial_val, reinterpreted_value);

        if (!error_evaluation.is_error_threshold_satisfied)
        {
            /* Revert the changes to the value */
            cleared_value = save_previous_value;
            is_clearing_progressing = false;
        }

        if constexpr (std::is_same_v<T, float> || std::is_same_v<T, double>)
        {
            if (not std::isfinite(reinterpreted_value))
            {
                cleared_value = save_previous_value;
                is_clearing_progressing = false;
            }
        }

        ++iteration_count;
    }

    return cleared_value;
}

static int jjjjj = 0;

template<typename T>
T
MultiResEmbeddedAdaptData<T>::PerformTailTruncation(const T init_value, const int which_tree, const int lelement_idx) const
{
    const T missing_value = this->GetVariableAttributes().GetMissingValue();

    CompressionValue<T> initial_value(init_value);

    /* Iterate through the serialized values and try to emplace as many zeros at the tail as possible (compliant to the error threshold) */
    if (!ApproxCompare(init_value, missing_value))
    {
        /* Get the permitted error for the current values */
        const std::vector<PermittedError> permitted_errors{PermittedError(CompressionCriterion::AbsoluteErrorThreshold, this->GetRemainingPermittedError(which_tree, lelement_idx))};

        /* Get the value which has been transformed by toggling as many ones from the back while setting the succeeding 'zero' bits to one */
        const CompressionValue<T> toggled_value = GetMaximumTailToggledValue<T>(permitted_errors, initial_value, missing_value);

        /* Get the value which has been transformed by clearing as many of the last set bits as possible */
        const CompressionValue<T> cleared_value = GetMaximumTailClearedValue<T>(permitted_errors, initial_value, missing_value);

        /* Check which approach leads to more zero bits at the end */
        const int num_toogled_trailing_zeros = toggled_value.GetNumberTrailingZeros();
        const int num_cleared_trailing_zeros = cleared_value.GetNumberTrailingZeros();

        /* Replace the initial value with the transformed one */
        if (num_toogled_trailing_zeros >= num_cleared_trailing_zeros)
        {
            /* If the toggling approach has been more successfull */
            initial_value = toggled_value;
        } else
        {
            /* If the clearing approach has been more successfull */
            initial_value = cleared_value;
        }

        /* Update the trail bit count for the new value */
        initial_value.UpdateTailBitCount();
        //cmc_debug_msg("Trailing Zeros: Toggled: ", num_toogled_trailing_zeros, ", Cleared: ", num_cleared_trailing_zeros);
    } else
    {
        /* In order to not change missing values, we are just able to trim their trailing zeros */
        initial_value.UpdateTailBitCount();
    }

    if (jjjjj >= 8000000 && jjjjj < 8001000)
    cmc_debug_msg("Control value has been trimmed by num bits: ", static_cast<int>(initial_value.GetTailBit()));

    jjjjj += 8;

    return initial_value.template ReinterpretDataAs<T>();
}



template <typename T>
cmc::lossy::embedded::idw::ExtractionData<T>
MultiResEmbeddedAdaptData<T>::PerformExtraction(t8_forest_t forest, t8_locidx_t ltreeid, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[])
{
    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (forest, ltreeid) + lelement_id;

    /* Get the view on the current values of this family of elements */
    const VectorView<T> values = this->GetViewOnInitialData(start_index, num_elements);

    cmc::lossy::embedded::idw::ElementData<T> extraction;

    #if 1

    #if 0
    std::vector<T> converted_values(8);
    int idx = 0;
    for (auto iter = values.begin(); iter != values.end(); ++iter, ++idx)
    {
        converted_values[idx] = std::abs(*iter);
    }

    auto min_elem = std::min_element(values.begin(), values.end());
    extraction.value = *min_elem;

    const int min_idx = std::distance(values.begin(), min_elem);
    
    #else
    const int min_idx = 0;
    //extraction.value = values[0];

    extraction.value = this->PerformTailTruncation(values[min_idx], ltreeid, lelement_id);

    #endif

    std::vector<double> midpoint(3);
    /* Get the midpoint of the min element */
    t8_forest_element_centroid (forest, ltreeid, elements[min_idx], midpoint.data());

    /* Get the midpoint of the max element */
    extraction.coordinates.x = static_cast<float>(midpoint[0]);
    extraction.coordinates.y = static_cast<float>(midpoint[1]);
    extraction.coordinates.z = static_cast<float>(midpoint[2]);

    extracted_child_ids.push_back(min_idx);

    return cmc::lossy::embedded::idw::ExtractionData<T>(std::move(extraction));


    #else


    /* Test value at parent element */
    std::array<t8_element_t*, 1> parent_elem{nullptr};
    scheme->element_new (tree_class, 1, parent_elem.data());
    scheme->element_get_parent (tree_class, elements[0], parent_elem[0]); 

    std::vector<double> midpoint(3);
    /* Get the midpoint of the min element */
    t8_forest_element_centroid (forest, ltreeid, parent_elem[0], midpoint.data());

    scheme->element_destroy(tree_class, 1, parent_elem.data());

    std::vector<T> converted_values(8);
    int idx = 0;
    for (auto iter = values.begin(); iter != values.end(); ++iter, ++idx)
    {
        converted_values[idx] = std::abs(*iter);
    }

    extraction.value = InterpolateToMidRange<T>(converted_values);

    int chosen_child_idx{-1};
    //for (auto iter = values.begin(); iter != values.end(); ++iter)
    //{
    //    if (cmc::ApproxCompare(*iter, extraction.value))
    //    {
    //        chosen_child_idx = std::distance(values.begin(), iter);
    //        break;
    //    }
    //}

    extracted_child_ids.push_back(chosen_child_idx);

    /* Get the midpoint of the max element */
    extraction.coordinates.x = static_cast<float>(midpoint[0]);
    extraction.coordinates.y = static_cast<float>(midpoint[1]);
    extraction.coordinates.z = static_cast<float>(midpoint[2]);


    return cmc::lossy::embedded::idw::ExtractionData<T>(std::move(extraction));

    #endif

}

#endif

template <typename T>
cmc::lossy::embedded::idw::UnchangedData<T>
MultiResEmbeddedAdaptData<T>::ElementStaysUnchanged(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
    const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[])
{
    /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
    //resdiual_order_indications_.AppendUnsetBit();

    //Test something with IDW
    //t8_element_t* element = t8_forest_get_leaf_element_in_tree (this->GetAmrMesh().GetMesh(), which_tree, lelement_id);
    //std::vector<double> midpoint(3);
    /* Get the midpoint of the min element */
    //t8_forest_element_centroid (forest, ltreeid, element, midpoint.data());
    //point_data_.push_back(PointData(static_cast<float>(midpoint[0]), static_cast<float>(midpoint[1]), static_cast<float>(midpoint[2])));

    #if 1
    cmc::lossy::embedded::idw::ElementData<T> extraction;

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (forest, which_tree) + lelement_id;

    const T value = this->GetInitialDataValueAtIndex(start_index);
    extraction.value = value;

    std::vector<double> midpoint(3);
    /* Get the midpoint of the current element */
    t8_forest_element_centroid (forest, which_tree, elements[0], midpoint.data());
    extraction.coordinates.x = static_cast<float>(midpoint[0]);
    extraction.coordinates.y = static_cast<float>(midpoint[1]);
    extraction.coordinates.z = static_cast<float>(midpoint[2]);

    /* We drag this value along to the "coarser level" until it this element is passed with its siblings
     * as a family of elements into the adaptation callback. Moreover, we leave an empty value 
     * on the "finer level" */
    return cmc::lossy::embedded::idw::UnchangedData<T>(extraction);

    #else

    cmc::lossy::embedded::idw::ElementData<T> extraction;

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (forest, which_tree) + lelement_id;

    const T value = this->GetInitialDataValueAtIndex(start_index);
    extraction.value = value;

    /* Construct first child */
    std::array<t8_element_t*, 1> child_elem{nullptr};
    scheme->element_new (tree_class, 1, child_elem.data());
    scheme->element_get_child (tree_class, elements[0], 0, child_elem[0]);

    std::vector<double> midpoint(3);
    /* Get the midpoint of the element */
    t8_forest_element_centroid (forest, which_tree, child_elem[0], midpoint.data());

    scheme->element_destroy(tree_class, 1, child_elem.data());

    t8_forest_element_centroid (forest, which_tree, child_elem[0], midpoint.data());
    extraction.coordinates.x = static_cast<float>(midpoint[0]);
    extraction.coordinates.y = static_cast<float>(midpoint[1]);
    extraction.coordinates.z = static_cast<float>(midpoint[2]);

    return cmc::lossy::embedded::idw::UnchangedData<T>(extraction);

    #endif
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
        cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>::entropy_coder_->UpdateSymbolFrequency(signum + first_one_bit);
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
    
    cmc_assert(cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>::entropy_coder_ != nullptr);

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
    cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>::entropy_coder_->Reset(std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResCompressionAlphabet<T>>());
    cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>::entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Collect the symbols and their frequencies for the entropy coder */
    CollectSymbolFrequenciesForEntropyCoding(level_byte_values);

    /* Setup the interior structure for encoding */
    cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>::entropy_coder_->SetupEncoding(this->GetMPIComm());

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
        cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(signum + first_one_bit);

        /* The first one is implicitly given by the leading zero count; therefore we set the "front bit" in order to discard the leading zeros and the following one */
        val.SetFrontBit(first_one_bit + 1);

        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            encoding.AppendBits(val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
        }
    }

    /* We set an indicaton symbol that the process local end of the values have been reached */
    cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte);

    /* Indicate that the encoding has been finished and flush all pending encodings */
    cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>::entropy_coder_->FinishEncoding();

    /* Set up the BitVector holding the encoded data for further use */
    encoding.TrimToContent();

    /* Get the local remaining significant bits */
    const uint64_t local_remaining_significant_bits_num_bytes = static_cast<uint64_t>(encoding.size());

    /* Get the local encoded LZC, respectively the encoded first "one-bit" positions */
    cmc::bit_map::BitMap local_encoded_lzc_stream = cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>::entropy_coder_->GetEncodedBitStream();
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
        cmc::bit_vector::BitVector encoded_alphabet = cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeAlphabet();
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
inline cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>*
CreateMultiResEmbeddedExtractionAdaptationClass(cmc::lossy::embedded::idw::AbstractEmbeddedByteCompressionVariable<T>* abstract_var, const CompressionSettings& settings)
{
    return new MultiResEmbeddedAdaptData<T>(abstract_var, settings);
}

template <typename T>
inline void
DestroyMultiResEmbeddedExtractionAdaptationClass(cmc::lossy::embedded::idw::IEmbeddedCompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class EmbeddedCompressionVariable : public cmc::lossy::embedded::idw::AbstractEmbeddedByteCompressionVariable<T>
{
public:
    EmbeddedCompressionVariable() = delete;

    EmbeddedCompressionVariable(const CompressionSettings& settings, input::Var& input_variable)
    : cmc::lossy::embedded::idw::AbstractEmbeddedByteCompressionVariable<T>(settings, input_variable)
    {
        //cmc_assert(input_variable.IsValid());//TODO: implement 

        this->SetName(input_variable.GetName());
        this->IndicateWhetherMeshRefinementBitsWillBeStored(true);

        cmc::lossy::embedded::idw::AbstractEmbeddedByteCompressionVariable<T>::adaptation_creator_ = CreateMultiResEmbeddedExtractionAdaptationClass<T>;
        cmc::lossy::embedded::idw::AbstractEmbeddedByteCompressionVariable<T>::adaptation_destructor_ = DestroyMultiResEmbeddedExtractionAdaptationClass<T>;
        cmc::lossy::embedded::idw::AbstractEmbeddedByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::EmbeddedMeshEncoder>();
    };

    CompressionSchema GetCompressionSchema() const override
    {
        return CompressionSchema::LossyEmbeddedMultiResExtractionNearestNeighborReconstruction;
    }

private:

};



}


#endif /* !LOSSY_CMC_EMBEDDED_IDW_EXTRACTION_COMPRESSION_HXX */
