#ifndef LOSSY_CMC_EMBEDDED_MULTI_RES_EXTRACTION_NEAREST_NEIGHBOR_RECONSTRUCTION_TRIMMED_RESIDUALS_HXX
#define LOSSY_CMC_EMBEDDED_MULTI_RES_EXTRACTION_NEAREST_NEIGHBOR_RECONSTRUCTION_TRIMMED_RESIDUALS_HXX


#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "lossy/cmc_embedded_byte_compression_extended_variable.hxx"
#include "utilities/cmc_multi_res_entropy_coder.hxx"
#include "utilities/cmc_multi_res_extraction_util.hxx"
#include "lossy/cmc_multi_res_extraction_residual_computation.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_interpolation_fn.hxx"
#include "mesh_compression/cmc_embedded_mesh_encoder.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>

namespace cmc::lossy::embedded::multi_res::nearest_neighbor_reconstruction
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

using err_type_t = cmc::lossy::embedded::extended::err_type_t;

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
class MultiResEmbeddedAdaptData : public cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>
{
public:
    MultiResEmbeddedAdaptData() = delete;
    MultiResEmbeddedAdaptData(cmc::lossy::embedded::extended::AbstractEmbeddedByteCompressionVariable<T>* variable, const CompressionSettings& settings)
    : cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>(variable, settings) {
        cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_ = std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResEncoder<T>>();
    };

    void InitializeExtractionIteration() override;
    void FinalizingExtractionIteration() override;
    std::pair<cmc::lossy::embedded::extended::ResidualData<T>, err_type_t> ComputeResiduals(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
        const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
        const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
        const t8_locidx_t first_incoming) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const override;
    std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const override;

    protected:
    cmc::lossy::embedded::extended::ExtractionData<T> PerformExtraction(t8_forest_t forest, t8_locidx_t which_tree, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                        const t8_scheme_c * ts, const int num_elements, t8_element_t * elements[]) override;
    cmc::lossy::embedded::extended::UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) override;

private:
    void CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const;
    std::vector<FaceDataValues<T>> GatherElementFaceValues(t8_forest_t forest, t8_locidx_t ltreeid, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                      const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[]) const;
    std::pair<bool, CompressionValue<T>> GetFaceValue(t8_forest_t forest, const int tree_id, const t8_scheme_c *scheme, const t8_element_t* elem, const int face_idx) const;
    FaceValues<T> ConstructHexFaceData(t8_forest_t forest_new, const int which_tree, const t8_eclass_t tree_class, const t8_scheme_c *scheme, const t8_element_t* elem) const;

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
        cmc_assert(neighbor_element_indices[0] < t8_forest_get_local_num_elements(forest));

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

constexpr double weight_coarse_approx = 0.25;
constexpr double weight_face_neighbor = 0.25;

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

#else


#if 1


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

#else


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

            const T prediction = (approximation + f0 + f2 + f4 - f1 - f3 - f5);
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

            const T prediction = (approximation + f1 + f2 + f4 - f0 - f3 - f5);
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

            const T prediction = (approximation + f0 + f3 + f4 - f1 - f2 - f5);
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

            const T prediction = (approximation + f1 + f3 + f4 - f0 - f2 - f5);
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

            const T prediction = (approximation + f0 + f2 + f5 - f1 - f3 - f4);
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

            const T prediction = (approximation + f1 + f2 + f5 - f0 - f3 - f4);
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

            const T prediction = (approximation + f0 + f3 + f5 - f1 - f2 - f4);
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

            const T prediction = (approximation + f1 + f3 + f5 - f0 - f2 - f4);
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

#endif

static size_t iiiiiii = 0;
static size_t zero_residuals = 0;

template<typename T>
std::tuple<bool, CompressionValue<T>, err_type_t>
PredictAndComputeResidual(const CompressionValue<T>& coarse_approximation, const FaceValues<T>& face_values, const CompressionValue<T>& initial_value, const err_type_t remaining_permitted_error, const int child_id)
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
    //TODO: Trim the residual compliant to the error bound: First try implemented below

    /* Setup the error criterion */
    std::vector<PermittedError> error_criterion{PermittedError(CompressionCriterion::AbsoluteErrorThreshold, static_cast<double>(remaining_permitted_error))};

    auto [is_approximation_greater, residual] = cmc::lossy::multi_res::ComputeResidual<T>(prediction, initial_value);

    /* We trim the residual by the means of clearing and toggeling the least significant bits from the residual */
    auto [tail_cleared_residual, tail_cleared_inaccurcy] = cmc::lossy::multi_res::GetMaximumTailClearedResidual(error_criterion, real_init_value, prediction, residual, is_approximation_greater);
    auto [tail_toggled_residual, tail_toggled_inaccuracy] = cmc::lossy::multi_res::GetMaximumTailToggledResidual(error_criterion, real_init_value, prediction, residual, is_approximation_greater);

    /* Check which approach lead to less significant bits */
    if (tail_cleared_residual.GetTailBit() >= tail_toggled_residual.GetTailBit())
    {
        //if (tail_cleared_inaccurcy > static_cast<double>(remaining_permitted_error))
        //{
        //    cmc_debug_msg("Hier ist tail cut inaccuracy: ", tail_cleared_inaccurcy, " und remaining permitetd error: ", remaining_permitted_error);
        //}
        if (std::isnan(tail_cleared_inaccurcy))
        {
            tail_cleared_inaccurcy = remaining_permitted_error;
        }

        cmc_assert(tail_cleared_inaccurcy <= static_cast<double>(remaining_permitted_error));

        const err_type_t new_remaining_permitted_error = remaining_permitted_error - tail_cleared_inaccurcy;

        /* If the tail cleared residual has less significant bits */
        return std::make_tuple(is_approximation_greater, tail_cleared_residual, new_remaining_permitted_error);
    } else
    {
        //if (tail_toggled_inaccuracy > static_cast<double>(remaining_permitted_error))
        //{
        //    cmc_err_msg("Hier ist tail cut inaccuracy: ", tail_toggled_inaccuracy, " und remaining permitetd error: ", remaining_permitted_error);
        //}
        //cmc_debug_msg("Hier ist tail cut inaccuracy: ", tail_toggled_inaccuracy, " und remaining permitetd error: ", remaining_permitted_error);
        //cmc_assert(tail_toggled_inaccuracy <= static_cast<double>(remaining_permitted_error));
        if (std::isnan(tail_toggled_inaccuracy))
        {
            tail_toggled_inaccuracy = remaining_permitted_error;
        }
        const err_type_t new_remaining_permitted_error = remaining_permitted_error - tail_toggled_inaccuracy;

        /* If the tail toggled residual has less significant bits */
        return std::make_tuple(is_approximation_greater, tail_toggled_residual, new_remaining_permitted_error);
    }
}

template <typename T>
std::pair<cmc::lossy::embedded::extended::ResidualData<T>, err_type_t>
MultiResEmbeddedAdaptData<T>::ComputeResiduals(t8_forest_t forest_old, t8_forest_t forest_new, t8_locidx_t which_tree,
    const t8_eclass_t tree_class, const t8_scheme_c *scheme, const int refine,
    const int num_outgoing, const t8_locidx_t first_outgoing, const int num_incoming,
    const t8_locidx_t first_incoming)
{
    /* The computation of residuals takes place here. Therefore, we do need to coarse approximation of the data 
     * and compute the reconstruction of the fine values */
    /* We only need to compute the residuals, if a coarsening has actually been taken place, otherwise the values are just dragged along the compression 
     * until they are a part of a family that is passed to the adaptation callback */
    if (refine != cmc::lossy::embedded::extended::kIndicationOfCoarseningDuringAdaptation)
    {
        /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
        resdiual_order_indications_.AppendUnsetBit();

        /* Therfore there is no residual associated with this element and we return an empty residual */
        cmc::lossy::embedded::extended::ResidualData<T> residual_data;
        residual_data.values.push_back(CompressionValue<T>());

        /* Get the remaining error that is associated with this element */
        const err_type_t previous_remaining_permitted_error = this->GetRemainingPermittedError(which_tree, first_outgoing);

        return std::make_pair(residual_data, previous_remaining_permitted_error);
    }

    // This hold only for hex elements currently 
    cmc_assert(tree_class == T8_ECLASS_HEX);

    /* Get the adapted coarse element */
    const t8_element_t* elem = t8_forest_get_element_in_tree(forest_new, which_tree, first_incoming);

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
        resdiual_order_indications_.AppendBit(is_approximation_greater);
        residuals.push_back(residual);

        /* Find the minimum remaining permitted error for the coarse element */
        if (min_remaining_permitted_error > remaining_error)
        {
            min_remaining_permitted_error = remaining_error;
        }
    }

    cmc_assert(not cmc::ApproxCompare(min_remaining_permitted_error, std::numeric_limits<err_type_t>::max()));

    return std::make_pair(cmc::lossy::embedded::extended::ResidualData<T>(std::move(residuals)), min_remaining_permitted_error);
}

template <typename T>
void
MultiResEmbeddedAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

template <typename T>
cmc::lossy::embedded::extended::ExtractionData<T>
MultiResEmbeddedAdaptData<T>::PerformExtraction(t8_forest_t forest, t8_locidx_t ltreeid, const t8_eclass_t tree_class, t8_locidx_t lelement_id,
                                                const t8_scheme_c * scheme, const int num_elements, t8_element_t * elements[])
{
    /* Only calculate the arithmetic mean of the values from the family of elements */
    
    /* Get the view on the current values of this family of elements */
    const VectorView<CompressionValue<T>> values = this->GetView(ltreeid, lelement_id, num_elements);

    /* Convert the view to actual values of the underlying data type */
    const std::vector<T> converted_values = ConvertCompressionValues<T>(values);

    /* Try the arithmetic mean as an predictor */
    const T mean = InterpolateToArithmeticMean<T>(converted_values);
    //Try the mid-range for comparison
    //const T mean = InterpolateToMidRange<T>(VectorView<T>(converted_values));
    return cmc::lossy::embedded::extended::ExtractionData<T>(CompressionValue<T>(mean));
}

template <typename T>
cmc::lossy::embedded::extended::UnchangedData<T>
MultiResEmbeddedAdaptData<T>::ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value)
{
    /* Since the residual is zero, we indicate that with an unset bit (although it is not of relevance, since the residual will be empty) */
    //resdiual_order_indications_.AppendUnsetBit();

    /* We drag this value along to the "coarser level" until it this element is passed with its siblings
     * as a family of elements into the adaptation callback. Moreover, we leave an empty value 
     * on the "finer level" */
    return cmc::lossy::embedded::extended::UnchangedData<T>(value);
}

#if 0
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
        cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->UpdateSymbolFrequency(signum + first_one_bit);
    }
}
#endif

//Aus trim residuals
template <typename T>
void
MultiResEmbeddedAdaptData<T>::CollectSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values) const
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
        cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->UpdateSymbolFrequency(signum + first_one_bit);


        /** Afterwards, we collect the trailing zero count (TZC) **/
        /* The collection of the TZC is only done, when the LZC does not fully describe the residual */
        if (first_one_bit >= sizeof(T) * bit_map::kCharBit - 1)
        {
            continue;
        }

        /* Otherwise, compute the TZC */
        const uint32_t last_one_bit = val.GetNumberTrailingZeros();

        /* Update this symbol for encoding */
        cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->UpdateSymbolFrequency(last_one_bit);       
    }
}

//Aus trim residuals
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
    
    cmc_assert(cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_ != nullptr);

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
    cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->Reset(std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResCompressionAlphabet<T>>());
    cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Collect the symbols and their frequencies for the entropy coder */
    CollectSymbolFrequenciesForEntropyCoding(level_byte_values);

    /* Setup the interior structure for encoding */
    cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->SetupEncoding(this->GetMPIComm());

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
        cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(signum + first_one_bit);

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
        cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(last_one_bit);

        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            encoding.AppendBits(val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
        }
    }

    /* We set an indicaton symbol that the process local end of the values have been reached */
    cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte);

    /* Indicate that the encoding has been finished and flush all pending encodings */
    cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->FinishEncoding();

    /* Set up the BitVector holding the encoded data for further use */
    encoding.TrimToContent();

    /* Get the local remaining significant bits */
    const uint64_t local_remaining_significant_bits_num_bytes = static_cast<uint64_t>(encoding.size());

    /* Get the local encoded LZC, respectively the encoded first "one-bit" positions */
    cmc::bit_map::BitMap local_encoded_lzc_stream = cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->GetEncodedBitStream();
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
        cmc::bit_vector::BitVector encoded_alphabet = cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeAlphabet();
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


#if 0
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
    
    cmc_assert(cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_ != nullptr);

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
    cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->Reset(std::make_unique<cmc::entropy_coding::arithmetic_coding::MultiResCompressionAlphabet<T>>());
    cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->InitializeAlphabet(sizeof(T));

    /* Collect the symbols and their frequencies for the entropy coder */
    CollectSymbolFrequenciesForEntropyCoding(level_byte_values);

    /* Setup the interior structure for encoding */
    cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->SetupEncoding(this->GetMPIComm());

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
        cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(signum + first_one_bit);

        /* The first one is implicitly given by the leading zero count; therefore we set the "front bit" in order to discard the leading zeros and the following one */
        val.SetFrontBit(first_one_bit + 1);

        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            encoding.AppendBits(val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
        }
    }

    /* We set an indicaton symbol that the process local end of the values have been reached */
    cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeSymbol(entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte);

    /* Indicate that the encoding has been finished and flush all pending encodings */
    cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->FinishEncoding();

    /* Set up the BitVector holding the encoded data for further use */
    encoding.TrimToContent();

    /* Get the local remaining significant bits */
    const uint64_t local_remaining_significant_bits_num_bytes = static_cast<uint64_t>(encoding.size());

    /* Get the local encoded LZC, respectively the encoded first "one-bit" positions */
    cmc::bit_map::BitMap local_encoded_lzc_stream = cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->GetEncodedBitStream();
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
        cmc::bit_vector::BitVector encoded_alphabet = cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>::entropy_coder_->EncodeAlphabet();
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

#endif

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
inline cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>*
CreateMultiResEmbeddedExtractionAdaptationClass(cmc::lossy::embedded::extended::AbstractEmbeddedByteCompressionVariable<T>* abstract_var, const CompressionSettings& settings)
{
    return new MultiResEmbeddedAdaptData<T>(abstract_var, settings);
}

template <typename T>
inline void
DestroyMultiResEmbeddedExtractionAdaptationClass(cmc::lossy::embedded::extended::IEmbeddedCompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class EmbeddedCompressionVariable : public cmc::lossy::embedded::extended::AbstractEmbeddedByteCompressionVariable<T>
{
public:
    EmbeddedCompressionVariable() = delete;

    EmbeddedCompressionVariable(const CompressionSettings& settings, input::Var& input_variable)
    : cmc::lossy::embedded::extended::AbstractEmbeddedByteCompressionVariable<T>(settings, input_variable)
    {
        //cmc_assert(input_variable.IsValid());//TODO: implement 

        this->SetName(input_variable.GetName());
        this->IndicateWhetherMeshRefinementBitsWillBeStored(true);

        cmc::lossy::embedded::extended::AbstractEmbeddedByteCompressionVariable<T>::adaptation_creator_ = CreateMultiResEmbeddedExtractionAdaptationClass<T>;
        cmc::lossy::embedded::extended::AbstractEmbeddedByteCompressionVariable<T>::adaptation_destructor_ = DestroyMultiResEmbeddedExtractionAdaptationClass<T>;
        cmc::lossy::embedded::extended::AbstractEmbeddedByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::EmbeddedMeshEncoder>();
    };

    CompressionSchema GetCompressionSchema() const override
    {
        return CompressionSchema::LossyEmbeddedMultiResExtractionNearestNeighborReconstruction;
    }

private:

};



}


#endif /* !LOSSY_CMC_EMBEDDED_MULTI_RES_EXTRACTION_NEAREST_NEIGHBOR_RECONSTRUCTION_TRIMMED_RESIDUALS_HXX */
