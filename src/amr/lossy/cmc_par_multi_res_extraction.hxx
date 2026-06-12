#ifndef CMC_AMR_LOSSY_PAR_MULTI_RES_EXTRACTION_HXX
#define CMC_AMR_LOSSY_PAR_MULTI_RES_EXTRACTION_HXX

#include "cmc.hxx"
#include "amr/lossy/cmc_par_multi_res_extraction_util.hxx"
#include "amr/lossy/cmc_par_multi_res_error_mesh.hxx"
#include "amr/lossy/cmc_par_multi_res_interpolation_util.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_compression_schema.hxx"

#include <t8_forest/t8_forest_partition.h>

#include <string>
#include <span>
#include <filesystem>

namespace cmc::par::lossy::multi_res
{

/* Forward declaration of the general compression variable */
template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
class CompressionVariableMultiData;

/* Typedef for the general case with one data point per variable */
template<ArithmeticType T, int32_t DIM>
using CompressionVariable = CompressionVariableMultiData<T, DIM, int32_t{1}>;

/* Forward declaration of the DataOuptut struct holding the compressed data to be writetn to disk */
template<ArithmeticType T>
struct DataOutput;

/* Actual class definition of the compression variable */
template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
class CompressionVariableMultiData
{
public:
    CompressionVariableMultiData() = delete;
    CompressionVariableMultiData(std::string name, t8_forest_t forest, const std::span<T> data, const std::vector<ErrorDomain>& error_domains)
    : name_(name), mesh_(forest), init_data_(data), error_domains_(error_domains) {
        /* Get the MPI communicator from the mesh */
        this->comm_ = t8_forest_get_mpicomm(forest);
        /* Get the rank and the size of the communicator */
        const int rv_size = MPI_Comm_size(this->comm_, &(this->comm_size_));
        MPICheckError(rv_size);
        const int rv_rank = MPI_Comm_rank(this->comm_, &(this->comm_rank_));
        MPICheckError(rv_rank);
        /* Increment the reference count, since we do not have ownership of the forest mesh */
        t8_forest_ref(forest);
    }

    void Compress();
    void WriteCompressedData(const std::string& file_name);
    
    void SetMaximumInitialElementLevel(const int max_init_elem_level);
    void SetIsMeshAlreadyPartitionedForCoarsening(const bool is_partitioned_for_coarsening);

    //TODO:
    std::vector<DataOutput<T>> __GetLocalOutputDataInMemory() const;
private:
    void DetermineMaxInitElementLevel();
    void PerformIntraElementCompression();
    void EncodeData();
    std::vector<uint64_t> EncodeRootLevelData() const;
    bool IsMeshCompressionProgressing() const;
    std::pair<t8_forest_t, std::vector<T>> Repartition(t8_forest_t adapted_mesh, std::vector<T>& adapted_data);
    bool HasIntraElementCompression() const;
    bool IsAlreadyPartitionedForCoarsening() const;
    void GenerateOuputStreams();
    void AppendVariableHeaderToStream(std::vector<uint64_t>& stream, const SizeType global_byte_count, const SizeType start_root_level_encoding_offset,
                                      const int mesh_compression_levels, const std::vector<SizeType>& global_bytes_per_level);
    std::vector<std::vector<uint64_t>> AdjustLevelMeshEncodings(const std::vector<PartitionInfo>& partition_info, const int num_global_mesh_encoding_steps,
                                                               const int num_global_encoding_steps);
    void CollectMeshEncodingOnTheRootRank(const std::vector<PartitionInfo>& global_partition_info, const SizeType num_global_encoding_steps, const SizeType num_global_mesh_encoding_steps, const std::vector<std::vector<uint64_t>>& offseted_level_mesh_encodings);
    std::pair<t8_forest_t, std::vector<T>> RepartitionForCompressionIteration(t8_forest_t mesh, std::vector<T>& data, const SizeType previous_bound);
    bool IsCoarsePredictorExtractionProgressing() const;
    void CollectCoarseLevelPredictionPyramid();

    /* The name of the variable */
    std::string name_;

    /* The input data for the compression */
    AmrMesh mesh_;
    const std::span<T> init_data_;
    int32_t max_init_elem_level_{kMaxPresentElementLevelUnknown};

    /* The MPI communicator to use (which is extracted from the mesh) */
    MPI_Comm comm_{MPI_COMM_NULL};
    int comm_size_{1}, comm_rank_{0};

    /* The current data during the extraction */
    std::vector<T> data_;


    const std::vector<ErrorDomain> error_domains_;

    std::vector<std::vector<T>> data_pyramid_;
    int mesh_coarsening_steps_{0};
    int prediction_step_{0};
    std::vector<SizeType> level_partition_offset_;
    //ErrorMesh error_indicator_;
    std::unique_ptr<ErrorMesh> error_mesh_;

    int compression_step_{0};


    /* Members filled by the mesh coarsening compression */
    std::vector<cmc::bits::vector> coarsening_indications_;
    std::vector<std::vector<ElemEncodingData<T, DIM>>> level_elements_encodings_;
    std::vector<PartitionInfo> level_partitioning_info_;

    /* Members filled by the intra element compression */
    std::vector<uint64_t> serialized_intra_element_entropy_dictionary_;
    std::vector<uint64_t> intra_element_encoded_data_;

    /* Encoded data after the compression */
    std::vector<uint64_t> serialized_entropy_dictionary_;
    std::vector<std::vector<uint64_t>> levelwise_encoded_data_;
    bool is_initially_partitioned_for_coarsening_{false};
    
    /* Mesh Encoding after the compression */
    std::vector<SizeType> num_elements_per_level_;
    std::vector<uint64_t> global_mesh_encoding_; //Only filled by the root rank

    /* Error Mesh Encoding */
    uint64_t error_mesh_encoding_size_bytes_{0};
    std::vector<uint64_t> error_mesh_encoding_; //Only filled by the root rank

    /* Data to be written out generated after compression */
    uint64_t global_compressed_byte_count_{0};
    std::vector<DataOutput<T>> output_streams_;
};

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
inline void
CompressionVariableMultiData<T, DIM, N>::SetMaximumInitialElementLevel(const int max_init_elem_level)
{
    if (max_init_elem_level <= 0 || max_init_elem_level > kMaxPossibleInitialRefinementLevel)
    {
        cmc_err_msg("The supplied refinemenet level ", max_init_elem_level, " is not in line with the supported features!");
    }
    max_init_elem_level_ = max_init_elem_level;
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
inline void
CompressionVariableMultiData<T, DIM, N>::SetIsMeshAlreadyPartitionedForCoarsening(const bool is_partitioned_for_coarsening)
{
    is_initially_partitioned_for_coarsening_ = is_partitioned_for_coarsening;
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
inline bool
CompressionVariableMultiData<T, DIM, N>::HasIntraElementCompression() const
{
    return (N > 1);
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
inline bool
CompressionVariableMultiData<T, DIM, N>::IsAlreadyPartitionedForCoarsening() const
{
    return is_initially_partitioned_for_coarsening_;
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
inline bool
CompressionVariableMultiData<T, DIM, N>::IsCoarsePredictorExtractionProgressing() const
{
    return (mesh_.GetNumberGlobalElements() > mesh_.GetNumberGlobalTrees());
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
struct CoarsePredictorIterationData
{
    CoarsePredictorIterationData(const std::span<T> current_data, const int32_t init_max_level, const int32_t coarsening_step)
    : data(current_data), current_coarsening_level{init_max_level - coarsening_step}
    {
        cmc_assert(init_max_level >= coarsening_step);

        coarse_level_data.reserve(current_data.size() / 4 + 1);
        coarsening_indications.Reserve(current_data.size() / 4 + 1);
    }

    void LeaveElementUnchanged(const int local_idx)
    {
        /* Indicate that corsening has been performed */
        this->coarsening_indications.AppendUnsetBit();
        /* The element's value is dragged along until we are able to caorsen it */
        this->coarse_level_data.push_back(this->data[local_idx]);
    }

    void PerformExtraction(const int local_idx)
    {
        /* Indicate that corsening has been performed */
        this->coarsening_indications.AppendSetBit();
        /* We extract the first family element's value as a coarse predictor */
        this->coarse_level_data.push_back(this->data[local_idx]);
    }

    const std::span<T> data;
    const int64_t current_coarsening_level;
    std::vector<T> coarse_level_data;
    cmc::bits::vector coarsening_indications;

};

constexpr bool
CheckIfElementIsEligibleForCoarsening(const int current_coarsening_level, const int element_level)
{
    return (current_coarsening_level == element_level);
}

template<typename T, int32_t DIM>
requires Dimension<DIM>
inline t8_locidx_t
CollectCoarsePredictors (t8_forest_t forest,
                         t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         const t8_eclass_t tree_class,
                         t8_locidx_t lelement_id,
                         const t8_scheme_c * ts,
                         const int is_family,
                         [[maybe_unused]] const int num_elements,
                         t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    CoarsePredictorIterationData<T, DIM>* adapt_data = static_cast<CoarsePredictorIterationData<T, DIM>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Compute the start offset in the local contiguous array of the data*/
    const int local_start_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

    /* Check if a family is supplied to the adaptation function */
    if (is_family == 0 || (not CheckIfElementIsEligibleForCoarsening(adapt_data->current_coarsening_level, ts->element_get_level(tree_class, elements[0]))))
    {
        /* If there is no family, the element stays unchanged */
        adapt_data->LeaveElementUnchanged(local_start_index);
        return cmc::t8::kLeaveElementUnchanged;
    } else
    {
        /* Extract a value of the family and coarsen it */
        adapt_data->PerformExtraction(local_start_index);
        return cmc::t8::kCoarsenElements;
    }
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
inline std::pair<t8_forest_t, std::vector<T>>
CompressionVariableMultiData<T, DIM, N>::Repartition(t8_forest_t adapted_mesh, std::vector<T>& adapted_data)
{
    /** Partition the mesh **/
    /* Keep the not-partitioned forest */
    t8_forest_ref(adapted_mesh);

    /* Allocate a forest */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    constexpr int partition_for_coarsening = 1;
    t8_forest_set_partition(partitioned_forest, adapted_mesh, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    /** Partition the adapted data **/
    /* Create an sc_array_t wrapper of the variable's data */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(adapted_data.data()), sizeof(T), adapted_data.size());

    /* Allocate an output vector for the partitioned data */
    std::vector<T> partitioned_data(t8_forest_get_local_num_leaf_elements(partitioned_forest));

    /* Create a wrapper for the freshly allocated partitioned data */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(partitioned_data.data()), sizeof(T), partitioned_data.size());

    /* Partition the variables data */
    t8_forest_partition_data(adapted_mesh, partitioned_forest, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Free the former forest and store the adapted/repartitioned mesh */
    t8_forest_unref(&adapted_mesh);

    return std::make_pair(partitioned_forest, partitioned_data);
}

static int step3{0};

inline void
WriteDataToVTKTest3(t8_forest_t mesh, const std::vector<float>& data)
{
    std::vector<double> double_data1;

    for (int idx{0}; idx < t8_forest_get_local_num_leaf_elements(mesh); ++idx)
    {
        double_data1.push_back(data[idx]);
    }

    t8_vtk_data_field_t vtk_data[1];
    snprintf (vtk_data[0].description, BUFSIZ, "GeneralData");
    vtk_data[0].type = T8_VTK_SCALAR;
    vtk_data[0].data = double_data1.data();

    const std::string file_name = std::string("cmc_lossy_mr_control_vals_") + std::to_string(step3);
    ++step3;
    t8_forest_write_vtk_ext (mesh, file_name.c_str(), 1, 1, 1, 1, 0, 0, 0, 1, vtk_data);
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
CompressionVariableMultiData<T, DIM, N>::CollectCoarseLevelPredictionPyramid()
{

    /* Allocate for the expected levelwise data-streams */
    this->coarsening_indications_.reserve(this->max_init_elem_level_ + 1);
    this->level_partition_offset_.reserve(this->max_init_elem_level_ + 1);

    // Coarsen until the root level predictors have been reached
    while(this->IsCoarsePredictorExtractionProgressing())
    {
        cmc_debug_msg("A coarsening iteration is initialized in order to extract the predictors.");

        /* Allocate a coarse predictor extraction struct */
        CoarsePredictorIterationData<T, DIM> adapt_data(std::span(this->data_pyramid_.back()), this->max_init_elem_level_, this->mesh_coarsening_steps_);

        /* Coarsen the mesh and the data */
        t8_forest_t coarse_mesh = t8_forest_new_adapt(this->mesh_.GetMesh(), CollectCoarsePredictors<T, DIM>, 0, 0, static_cast<void*>(&adapt_data));

        /* Store the new offsets of the coarser mesh (this is needed in order to perform the prediction later on correctly) */
        const SizeType mesh_offset = static_cast<SizeType>(t8_forest_get_first_local_leaf_element_id(coarse_mesh));
        this->level_partition_offset_.push_back(mesh_offset);

        /* Store the coarsening/refinement indications */
        this->coarsening_indications_.push_back(std::move(adapt_data.coarsening_indications));

        /* Partition the mesh and the data for the next coarsening step */
        //TODO: Repartition would not be needed at the last extraction step, since we need to re-partition it to the manual bound directly
        //when the compression starts (but then, the repartition in the first compression iterations needs to be put out as well)
        auto [partitioned_mesh, partitioned_data] = this->Repartition(coarse_mesh, adapt_data.coarse_level_data);

        /* Set the partitioned mesh */
        this->mesh_.SetMesh(partitioned_mesh);

        /* Store the extracted coarse level predictors */
        this->data_pyramid_.push_back(std::move(partitioned_data));

        /* Increment the step counter */
        ++(this->mesh_coarsening_steps_);

        WriteDataToVTKTest3(this->mesh_.GetMesh(), this->data_pyramid_.back());
    }
}


template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
inline bool
CompressionVariableMultiData<T, DIM, N>::IsMeshCompressionProgressing() const
{
    return (this->compression_step_ < this->mesh_coarsening_steps_);
}


template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
struct CompressionIterationData
{
    CompressionIterationData() = delete;
    CompressionIterationData(const std::span<T> current_data, const std::span<T> data_to_predict, const cmc::bits::vector& refinement_indications, const std::unique_ptr<ErrorMesh>& error_mesh_, const int32_t init_max_level, const int32_t coarsening_step)
    : data(current_data), data_predict(data_to_predict), refinement_indications(refinement_indications), error_mesh{error_mesh_}, current_coarsening_level{init_max_level - coarsening_step}
    {
        cmc_assert(init_max_level >= coarsening_step);

        elem_encodings.reserve(current_data.size());
    }

    bool WillNextElementBeRefined()
    {
        return refinement_indications.GetNextBit();
    }

    T GetData(const int idx) const {return data[idx];}
    
    inline float GetPermittedAbsError(const t8_scheme_c* scheme, const int global_tree_id, const t8_eclass_t tree_class, const t8_element_t* element) const
    {
        return error_mesh->GetPermittedAbsError(scheme, global_tree_id, tree_class, element);
    }

    void LeaveElementUnchanged(const int local_idx);

    void PerformRefinementPrediction(const t8_eclass_t tree_class, const std::span<T> initial_values, const std::vector<T>& control_values, const std::vector<float>& permitted_abs_error);
    void PerformQuadCompression(const std::span<T> initial_values, const std::vector<T>& control_values, const std::vector<float>& permitted_abs_error);
    void PerformHexCompression(const std::span<T> initial_values, const std::vector<T>& control_values, const std::vector<float>& permitted_abs_error);


    const std::span<T> data;
    const std::span<T> data_predict;
    cmc::bits::vector_view_in_memory refinement_indications;
    const std::unique_ptr<ErrorMesh>& error_mesh;
    const int64_t current_coarsening_level;
    std::vector<T> next_level_data_;
    std::vector<ElemEncodingData<T, DIM>> elem_encodings;
    int local_data_accessor_{0};
};

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void
CompressionIterationData<T, DIM>::LeaveElementUnchanged(const int local_index)
{
    /* Just copy the data over to the next level */
    this->next_level_data_.push_back(this->data[local_index]);
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline void 
CompressionIterationData<T, DIM>::PerformRefinementPrediction(const t8_eclass_t tree_class, const std::span<T> initial_values, const std::vector<T>& control_values, const std::vector<float>& permitted_abs_error)
{
    if constexpr (DIM == 1)
    {
        /* 1D Elements */
        cmc_err_msg("The element type t8_eclass_t = ", static_cast<int>(tree_class), " is currently not supported!");
    } else if constexpr (DIM == 2)
    {
        /* 2D ELements */
        switch (tree_class)
        {
            case t8_eclass::T8_ECLASS_QUAD:
                this->PerformQuadCompression(initial_values, control_values, permitted_abs_error);
            break;
            case t8_eclass::T8_ECLASS_TRIANGLE:
                cmc_err_msg("The compression for triangle elements is not yet implemented!");
            break;
            default:
                cmc_err_msg("The element type t8_eclass_t = ", static_cast<int>(tree_class), " is currently not supported!");
        }
    } else if constexpr (DIM == 2)
    {
        /* 3D Elements */
        switch (tree_class)
        {
            case t8_eclass::T8_ECLASS_HEX:
                cmc_err_msg("The compression for hex elements is not yet implemented!");
                this->PerformHexCompression();
            break;
            case t8_eclass::T8_ECLASS_TET:
                cmc_err_msg("The compression for tetrahedral elements is not yet implemented!");
            break;
            case t8_eclass::T8_ECLASS_PRISM:
                cmc_err_msg("The compression for prismatic elements is not yet implemented!");
            break;
            case t8_eclass::T8_ECLASS_PYRAMID:
                cmc_err_msg("The compression for pyramidal elements is not yet implemented!");
            break;
            default:
                cmc_err_msg("The element type t8_eclass_t = ", static_cast<int>(tree_class), " is currently not supported!");
        }
    } else
    {
        cmc_err_msg("An unsupported element type has been supplied!");
    }
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline void 
CompressionIterationData<T, DIM>::PerformQuadCompression(const std::span<T> initial_values, const std::vector<T>& control_values, const std::vector<float>& permitted_abs_error)
{
    cmc_assert(control_values.size() == cmc::par::lossy::rbf::util::kNumQuadControlPoints);
    
    /* Copy the control values */
    const std::array<T, cmc::par::lossy::rbf::util::kNumQuadControlPoints> control_vals{control_values[0], control_values[1], control_values[2], control_values[3], control_values[4]};

    /* Perform the quad prediction */
    const std::array<T, cmc::par::lossy::rbf::util::kNumQuadPredictionPoints> predictions = cmc::par::lossy::rbf::util::PerformQuadPrediction<T>(control_vals);

    /* Create the predictions vector for this family of elements */
    /** The first child element is always equal to the first control point **/
    std::vector<T> fam_predictions;
    fam_predictions.reserve(1 + cmc::par::lossy::rbf::util::kNumQuadPredictionPoints);
    fam_predictions.push_back(control_vals[0]);
    std::copy_n(predictions.begin(), cmc::par::lossy::rbf::util::kNumQuadPredictionPoints, std::back_inserter(fam_predictions));

    this->elem_encodings.emplace_back();

    /* Check the deviation between the predicted and intial values */
    const int num_elems = 1 + cmc::par::lossy::rbf::util::kNumQuadPredictionPoints;
    cmc_assert(static_cast<size_t>(num_elems) == initial_values.size());

    for (int elem_idx{0}; elem_idx < num_elems; ++elem_idx)
    {
        /* Compute the residual between initial value and predicted value */
        const T residual = GetAbsResidual<T>(fam_predictions[elem_idx], initial_values[elem_idx]);

        /* Determine whether the prediction is in line with the permitted error */
        if (residual <= permitted_abs_error[elem_idx])
        {
            /* If the residual is within the permitted error bound */
            this->elem_encodings.back().quantization_bins[elem_idx] = kPredictionWithinBound;

            /* We store the prediction */
            this->next_level_data_.push_back(fam_predictions[elem_idx]);

        } else if (residual <= kResidualMaxDeviationFactor * permitted_abs_error[elem_idx])
        {
            /* If the residual is within the permitted error interval such that the error can be met by quantization */
            const SymbolType bin = static_cast<SymbolType>(std::floor(residual / (2 * permitted_abs_error[elem_idx]) + 0.5));

            /* Check in which direction the quantization goes */
            const bool is_prediction_greater = (fam_predictions[elem_idx] >= initial_values[elem_idx]);

            /* Create the entropy symbol from the information above */
            const SymbolType entropy_symbol = CreateEntropySymbolFromQuantizationBin(bin, is_prediction_greater); 

            /* In this case, we do not need to store anything apart the quantization bin */
            this->elem_encodings.back().quantization_bins[elem_idx] = entropy_symbol;

            /* De-Quantize the value */
            const T decompressed_value = fam_predictions[elem_idx] + (is_prediction_greater ? -2.0 : +2.0) * permitted_abs_error[elem_idx] * bin;

            /* We store the de-quantized value */
            this->next_level_data_.push_back(decompressed_value);

        } else
        {
            /* In this case, the value is unpredictable and we store it, as it is */
            /* We flag the value as unpredictable */
            this->elem_encodings.back().quantization_bins[elem_idx] = kFlagUnpredictable;
            
            /* And we store the actual value */
            this->elem_encodings.back().unpredictable_values[elem_idx] = initial_values[elem_idx];
            this->next_level_data_.push_back(initial_values[elem_idx]);
        }
    }

    /* Store the number of elements */
    this->elem_encodings.back().num_elements = num_elems;
}


template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
void inline 
CompressionIterationData<T, DIM>::PerformHexCompression(const std::span<T> initial_values, const std::vector<T>& control_values, const std::vector<float>& permitted_abs_error)
{
    cmc_assert(control_values.size() == cmc::par::lossy::rbf::util::kNumHexControlPoints);
    
    /* Copy the control values */
    const std::array<T, cmc::par::lossy::rbf::util::kNumHexControlPoints> control_vals{control_values[0], control_values[1], control_values[2], control_values[3], control_values[4], control_values[5], control_values[6]};

    /* Perform the hex predcition */
    const std::array<T, cmc::par::lossy::rbf::util::kNumHexPredictionPoints> predictions = cmc::par::lossy::rbf::util::PerformHexPrediction<T>(control_vals);

    /* Create the predictions vector for this family of elements */
    /** The first child element is always equal to the first control point **/
    std::vector<T> fam_predictions;
    fam_predictions.reserve(1 + cmc::par::lossy::rbf::util::kNumHexPredictionPoints);
    fam_predictions.push_back(control_vals[0]);
    std::copy_n(predictions.begin(), cmc::par::lossy::rbf::util::kNumHexPredictionPoints, std::back_inserter(fam_predictions));

    this->elem_encodings.emplace_back();

    /* Check the deviation between the predicted and intial values */
    const int num_elems = 1 + cmc::par::lossy::rbf::util::kNumHexPredictionPoints;
    cmc_assert(static_cast<size_t>(num_elems) == initial_values.size());
    cmc_assert(static_cast<size_t>(num_elems) == fam_predictions.size());

    for (int elem_idx{0}; elem_idx < num_elems; ++elem_idx)
    {
        /* Compute the residual between initial value and predicted value */
        const T residual = GetAbsResidual<T>(fam_predictions[elem_idx], initial_values[elem_idx]);

        /* Determine whether the prediction is in line with the permitted error */
        if (residual <= permitted_abs_error[elem_idx])
        {
            /* If the residual is within the permitted error bound */
            this->elem_encodings.back().quantization_bins[elem_idx] = kPredictionWithinBound;

            /* We store the prediction */
            this->next_level_data_.push_back(fam_predictions[elem_idx]);

        } else if (residual <= kResidualMaxDeviationFactor * permitted_abs_error[elem_idx])
        {
            /* If the residual is within the permitted error interval such that the error can be met by quantization */
            const SymbolType bin = static_cast<SymbolType>(std::floor(residual / (2 * permitted_abs_error[elem_idx]) + 0.5));

            /* Check in which direction the quantization goes */
            const bool is_prediction_greater = (fam_predictions[elem_idx] >= initial_values[elem_idx]);

            /* Create the entropy symbol from the information above */
            const SymbolType entropy_symbol = CreateEntropySymbolFromQuantizationBin(bin, is_prediction_greater); 

            /* In this case, we do not need to store anything apart the quantization bin */
            this->elem_encodings.back().quantization_bins[elem_idx] = entropy_symbol;

            /* De-Quantize the value */
            const T decompressed_value = fam_predictions[elem_idx] + (is_prediction_greater ? -2.0 : +2.0) * permitted_abs_error[elem_idx] * bin;

            /* We store the de-quantized value */
            this->next_level_data_.push_back(decompressed_value);

        } else
        {
            /* In this case, the value is unpredictable and we store it, as it is */
            /* We flag the value as unpredictable */
            this->elem_encodings.back().quantization_bins[elem_idx] = kFlagUnpredictable;
            
            /* And we store the actual value */
            this->elem_encodings.back().unpredictable_values[elem_idx] = initial_values[elem_idx];
            this->next_level_data_.push_back(initial_values[elem_idx]);
        }
    }

    /* Store the number of elements */
    this->elem_encodings.back().num_elements = num_elems;
}

template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
inline t8_locidx_t
LossyMultiResCompression (t8_forest_t forest,
                          t8_forest_t forest_from,
                          t8_locidx_t which_tree,
                          const t8_eclass_t tree_class,
                          t8_locidx_t lelement_id,
                          const t8_scheme_c * ts,
                          const int is_family,
                          const int num_elements,
                          t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    CompressionIterationData<T, DIM>* adapt_data = static_cast<CompressionIterationData<T, DIM>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);
    
    /* Compute the start offset in the local contiguous array of the data*/
    const int local_start_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

    const int global_elem_idx = t8_forest_get_first_local_leaf_element_id(forest_from) + local_start_index;

    /* Check if the current element will be refined/performs a prediction */
    if (adapt_data->WillNextElementBeRefined())
    {
        /** The element will be refined and the children elements need to be predicted **/

        /* Convert the local to a global tree id */
        const int global_tree_idx = t8_forest_global_tree_id (forest_from, which_tree);

        /* Get the number of children elements that need to be predicted */
        const int num_children = t8_element_get_num_children(ts, tree_class, elements[0]);
        cmc_assert(num_children <= kNumMaxChildrenElements<DIM>);

        /* Get the number of faces of this element */
        const int num_faces = t8_element_get_num_faces(ts, tree_class, elements[0]);

        /* Get this families control value */
        const T fam_control_value = adapt_data->GetData(local_start_index);

        /* Create a vector for all face_values (and fill it with the default value (this families control value)) */
        std::vector<T> control_points(num_faces + 1, fam_control_value);

        /* Allocate some variables to be filled by the face neighbor calls */
        const t8_element_t** neighbor_leaves;
        int* dual_faces;
        int num_neighbors{0};
        t8_locidx_t* neighbor_element_indices;
        t8_eclass_t neighbor_tree_class;

        /* Iterate over all faces and fill the control points */
        for (int face_idx{0}; face_idx < num_faces; ++face_idx)
        {
            /* Gather the face neighbor via this face */
            t8_forest_leaf_face_neighbors (forest_from, which_tree, elements[0], &neighbor_leaves, face_idx, &dual_faces, &num_neighbors,
                &neighbor_element_indices, &neighbor_tree_class);

            /* Check if there is a neighboring element at the face */
            if (num_neighbors == 1) [[likely]]
            {
                /* Get the control point */
                const T control_point = adapt_data->GetData(neighbor_element_indices[0]);
                /* Store the control point */
                control_points[1 + face_idx] = control_point;

                /* Deallocate the memory for the face neighbor construction */
                T8_FREE (neighbor_leaves);
                T8_FREE (neighbor_element_indices);
                T8_FREE (dual_faces);
            }
            else if (num_neighbors > 0) [[unlikely]]
            {
                cmc_err_msg("only one neighbor is supported");
                /* In theory, this case should not happen, since we start on the finest refinement level currently */
                /** If there are some, we use their values as a control point (otherwise we use the family's own control value) **/
                T control_point{};

                for (int neigh_idx{0}; neigh_idx < num_neighbors; ++neigh_idx)
                {
                    /* Sum up over all face neighbors */
                    control_point += adapt_data->GetData(neighbor_element_indices[neigh_idx]);
                }

                /* Build the mean value */
                control_point = control_point / static_cast<T>(num_neighbors);

                /* Store the control point */
                control_points[1 + face_idx] = control_point;

                /* Deallocate the memory for the face neighbor construction */
                T8_FREE (neighbor_leaves);
                T8_FREE (neighbor_element_indices);
                T8_FREE (dual_faces);
            }
        }

        /* We need to get the permitted errors for the children elements that will be constructed */
        std::vector<float> permitted_abs_errors;
        permitted_abs_errors.reserve(num_children);
        
        /* Allocate storage for the children elements */
        std::vector<t8_element_t*> child_elements(num_children, nullptr);
        t8_element_new (ts, tree_class, num_children, child_elements.data());

        /* Create the children and ask for their permitted absolute errors */
        t8_element_get_children(ts, tree_class, elements[0], num_children, child_elements.data());

        /* Get the permitted errors for the child elements */
        for (int child_idx{0}; child_idx < num_children; ++child_idx)
        {
            permitted_abs_errors.emplace_back(adapt_data->GetPermittedAbsError(ts, global_tree_idx, tree_class, child_elements[child_idx]));
        }

        /* Destroy the constructed elements */
        t8_element_destroy(ts, tree_class, num_children, child_elements.data());

        /* Create a span on the initial data */
        const std::span<T> init_data(&(adapt_data->data_predict[adapt_data->local_data_accessor_]), num_children);

        /* Perform the prediction and get the data for encoding */
        adapt_data->PerformRefinementPrediction(tree_class, init_data, control_points, permitted_abs_errors);

        /* Update the local accessor index */
        adapt_data->local_data_accessor_ += num_children;

        return cmc::t8::kRefineElement;
    } else
    {
        /**  The element remains unchanged, therefore no prediction needs to be made **/
        adapt_data->LeaveElementUnchanged(local_start_index);
        
        #ifdef CMC_DEBUG
        /* Compare whether the permitted error still holds */

        #endif
        
        /* Update the accessor for the prediction data */
        ++(adapt_data->local_data_accessor_);

        return cmc::t8::kLeaveElementUnchanged;
    }
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
std::pair<t8_forest_t, std::vector<T>>
CompressionVariableMultiData<T, DIM, N>::RepartitionForCompressionIteration(t8_forest_t mesh, std::vector<T>& data, const SizeType previous_bound)
{
    const t8_locidx_t current_local_elems = t8_forest_get_local_num_leaf_elements(mesh);
    [[maybe_unused]] const t8_locidx_t current_ghost_elems = t8_forest_get_num_ghosts (mesh);

    cmc_assert(current_local_elems + current_ghost_elems == static_cast<t8_locidx_t>(data.size()));

    /* Keep the not-partitioned forest */
    t8_forest_ref(mesh);

    /** Partition the forest correctly and build a halo layer **/
    t8_forest_t partitioned_ghost_mesh;
    t8_forest_init (&partitioned_ghost_mesh);

    /* Set the forest for partitioning */
    t8_forest_set_partition (partitioned_ghost_mesh, mesh, 0);

    /* Set the partition bound explicitly */
    t8_forest_set_partition_offset (partitioned_ghost_mesh, static_cast<t8_gloidx_t>(previous_bound));

    /* Set the forest for creating a face ghost layer */
    t8_forest_set_ghost (partitioned_ghost_mesh, 1, T8_GHOST_FACES);

    /* Commit the forest, this step will perform the partitioning and ghost layer creation. */
    t8_forest_commit (partitioned_ghost_mesh);

    /** Exchange the data corectly to set up the halo layer **/
    /* Get the number of local elements */
    const t8_locidx_t num_local_elements = t8_forest_get_local_num_leaf_elements(partitioned_ghost_mesh);
    
    /* Get the number of ghost elements of forest. */
    const t8_locidx_t num_ghost_elements = t8_forest_get_num_ghosts (partitioned_ghost_mesh);

    /* Create an sc_array_t wrapper of the variable's local element data */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(data.data()), sizeof(T), current_local_elems);

    /* Allocate an output vector for the partitioned data */
    std::vector<T> partitioned_ghost_data(num_local_elements + num_ghost_elements);

    /* Create a wrapper for the freshly allocated partitioned data */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(partitioned_ghost_data.data()), sizeof(T), num_local_elements);

    /* Partition the variables local element data */
    t8_forest_partition_data(mesh, partitioned_ghost_mesh, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Free the former forest */
    t8_forest_unref(&mesh);

    /* Create a wrapper for the ghost exchange */
    sc_array_t* ghost_exchange_data = sc_array_new_data (partitioned_ghost_data.data(), sizeof(T), num_local_elements + num_ghost_elements);

    /* Exchange the ghost data for the newly partitioned data */
    t8_forest_ghost_exchange_data (partitioned_ghost_mesh, ghost_exchange_data);

    /* Destroy the array wrapper */
    sc_array_destroy(ghost_exchange_data);

    return std::make_pair(partitioned_ghost_mesh, std::move(partitioned_ghost_data));
}

static int step2{0};

inline void
WriteDataToVTKTest2(t8_forest_t mesh, const std::vector<float>& data)
{
    std::vector<double> double_data1;

    for (int idx{0}; idx < t8_forest_get_local_num_leaf_elements(mesh); ++idx)
    {
        double_data1.push_back(data[idx]);
    }

    t8_vtk_data_field_t vtk_data[1];
    snprintf (vtk_data[0].description, BUFSIZ, "GeneralData");
    vtk_data[0].type = T8_VTK_SCALAR;
    vtk_data[0].data = double_data1.data();

    const std::string file_name = std::string("cmc_test_lossy_mr_") + std::to_string(step2);
    ++step2;
    t8_forest_write_vtk_ext (mesh, file_name.c_str(), 1, 1, 1, 1, 0, 0, 0, 1, vtk_data);
}



template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
CompressionVariableMultiData<T, DIM, N>::Compress()
{
    cmc_debug_msg("The lossy multi-resolution compression of variable ", this->name_, " is performed.");
    /* Potentially, gather the maximum present element level */
    this->DetermineMaxInitElementLevel();

    /* Store the initial data in the data pyramid */
    std::vector<T> initial_data;
    initial_data.reserve(this->init_data_.size());
    std::copy_n(this->init_data_.begin(), this->init_data_.size(), std::back_inserter(initial_data));
    this->data_pyramid_.push_back(std::move(initial_data));

    /* Potentially, perform partition for coarsening */
    if (not this->IsAlreadyPartitionedForCoarsening())
    {
        auto [partitioned_mesh, partitioned_data] = this->Repartition(this->mesh_.GetMesh(), this->data_pyramid_.back());

        /* Set the partitioned mesh and data */
        this->mesh_.SetMesh(partitioned_mesh);
        this->data_pyramid_.back() = std::move(partitioned_data);
    }

    /* Create the error mesh */
    this->error_mesh_ = std::make_unique<ErrorMesh>(this->mesh_.GetMesh(), this->error_domains_, this->data_pyramid_.back());

    //Perfom Intra Element Coarse-Predictor Extraction?

    /* Collect the coarse level predictors from the mesh coarsening steps */
    this->CollectCoarseLevelPredictionPyramid();
 
    this->compression_step_ = 0;

    /* Start on the coarsest level of the data pyramid */
    auto coarse_data_iter = this->data_pyramid_.rbegin();
    auto level_partition_iter = this->level_partition_offset_.rbegin();
    auto refinement_indicator_iter = this->coarsening_indications_.rbegin();

    /* We start with the prediction on the root level */
    this->data_ = *coarse_data_iter;
    ++coarse_data_iter;


    // Compress until the leaf level
    while(this->IsMeshCompressionProgressing())
    {
        cmc_debug_msg("A coarsening iteration is initialized.");

        /* Repartition the data to the manual bound in order to perform the process-local compression correctly */
        auto [partitioned_ghost_mesh, partitioned_ghost_data] = this->RepartitionForCompressionIteration(this->mesh_.GetMesh(), this->data_, *level_partition_iter);

        // We need to store the number of local elements before the adaptation/refinement/prediction for the partition info 
        this->level_partitioning_info_.emplace_back(static_cast<uint64_t>(t8_forest_get_local_num_leaf_elements(partitioned_ghost_mesh)));

        /* Set the mesh and the data */
        this->mesh_.SetMesh(partitioned_ghost_mesh);
        this->data_ = std::move(partitioned_ghost_data);

        // 1) Allocate a compression iteration 
        CompressionIterationData<T, DIM> adapt_data(std::span(this->data_), *coarse_data_iter, *refinement_indicator_iter, this->error_mesh_, this->max_init_elem_level_, this->compression_step_);

        // 2)  Adapt
        /* Perform a refinement/prediction iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(this->mesh_.GetMesh(), LossyMultiResCompression<T, DIM>, 0, 0, static_cast<void*>(&adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements");

        // 3) Store data to be encoded later
        this->level_elements_encodings_.push_back(std::move(adapt_data.elem_encodings));

        // 4) Set the data/control values for the next iteration 
        this->data_ = std::move(adapt_data.next_level_data_);

        // Set the adapted forest for the next iteration 
        this->mesh_.SetMesh(adapted_forest);

        /* Update the counter and iterators */
        ++(this->compression_step_);
        ++(coarse_data_iter);
        ++(level_partition_iter);
        ++(refinement_indicator_iter);

        //if (this->compression_step_ > 0) {return;}
        WriteDataToVTKTest2(mesh_.GetMesh(), this->data_);
    }

    //Perfom Intra Element Compression

    /* Encode the data that has been collected */
    this->EncodeData();

    /* Clean-up the error mesh MPI-allocations after the encoding of the error mesh */
    this->error_mesh_->DestructErrorMeshCollectively();

    cmc_debug_msg("The lossly multi-resolution compression of variable ", this->name_, " has been completed.");
}


template<ArithmeticType T, int32_t DIM>
requires Dimension<DIM>
std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>>
ExchangeEntropySymbols(const std::vector<std::vector<ElemEncodingData<T, DIM>>>& levelwise_encodings, const MPI_Comm comm)
{
    /* Get the number of all possible entropy symbols */
    constexpr int num_entropy_symbols = GetNumEntropySymbols();

    /* Set the array and zero intialiaze the frequencies */
    std::array<uint64_t, num_entropy_symbols> entropy_symbol_frequencies{};

    /* Iterate through all entropy codes and accumulate their frequencies */
    for (size_t lvl_idx{0}; lvl_idx < levelwise_encodings.size(); ++lvl_idx)
    {
        /* Iterate through all coarsening data on this level */
        for (size_t coarsening_idx{0}; coarsening_idx < levelwise_encodings[lvl_idx].size(); ++coarsening_idx)
        {
            /* Iterate over all entropy codes from this coarsening data */
            for (unsigned entropy_sym_idx{0}; entropy_sym_idx < levelwise_encodings[lvl_idx][coarsening_idx].num_elements; ++entropy_sym_idx)
            {
                /* Convert the symbol to the corresponding array index */
                const int array_idx = MapEntropySymbolToArrayIndex<T>(levelwise_encodings[lvl_idx][coarsening_idx].quantization_bins[entropy_sym_idx]);
                /* Update the frequency */
                ++entropy_symbol_frequencies[array_idx];
            }
        }
    }

    /* Add the process end symbol */
    AddProcessEndSymbol(entropy_symbol_frequencies, levelwise_encodings.size());

    /* After all entropy symbol frequencies have been collected, we exchange them */
    std::array<uint64_t, num_entropy_symbols> exchanged_entropy_symbol_frequencies{};

    /* Exchange the frequencies */
    const int rv_allreduce = MPI_Allreduce(entropy_symbol_frequencies.data(), exchanged_entropy_symbol_frequencies.data(), num_entropy_symbols, MPI_UINT64_T, MPI_SUM, comm);
    MPICheckError(rv_allreduce);

    /* Replicate the global entropy frequencies locally */
    std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>> global_symbol_frequencies;
    global_symbol_frequencies.reserve(num_entropy_symbols);

    for (int idx{0}; idx < num_entropy_symbols; ++idx)
    {
        /* Convert the index back to the entropy symbol */
        const SymbolType entropy_symbol = MapArrayIndexToEntropySymbol<T>(idx);

        /* Store the symbol with the global frequency */
        global_symbol_frequencies.emplace_back(entropy_symbol, exchanged_entropy_symbol_frequencies[idx]);
    }

    return global_symbol_frequencies;
}

template<OneByteArithmeticType T>
inline std::vector<uint64_t>
PerformRootLevelEncoding(const std::vector<T>& data_)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(data_.size() * sizeof(T));
    
    for (size_t idx{0}; idx < data_.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<OneByteResidualType>(std::bit_cast<OneByteResidualType>(data_[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}

template<TwoByteArithmeticType T>
inline std::vector<uint64_t>
PerformRootLevelEncoding(const std::vector<T>& data_)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(data_.size() * sizeof(T));
    
    for (size_t idx{0}; idx < data_.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<TwoByteResidualType>(std::bit_cast<TwoByteResidualType>(data_[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}

template<FourByteArithmeticType T>
inline std::vector<uint64_t>
PerformRootLevelEncoding(const std::vector<T>& data_)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(data_.size() * sizeof(T));
    
    for (size_t idx{0}; idx < data_.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<FourByteResidualType>(std::bit_cast<FourByteResidualType>(data_[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}

template<EightByteArithmeticType T>
inline std::vector<uint64_t>
PerformRootLevelEncoding(const std::vector<T>& data_)
{
    cmc::bits::vector root_lvl_encoding;
    root_lvl_encoding.Reserve(data_.size() * sizeof(T));
    
    for (size_t idx{0}; idx < data_.size(); ++idx)
    {
        /* Store the value */
        root_lvl_encoding.AppendBits<EightByteResidualType>(std::bit_cast<EightByteResidualType>(data_[idx]), 0, 0);
    }
    return root_lvl_encoding.GetSerializedByteStreamBE();
}
//TODO: IF uneven num root elements per process, we have gaps in the encoding if the type is not 64-bit
template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
std::vector<uint64_t>
CompressionVariableMultiData<T, DIM, N>::EncodeRootLevelData() const
{   
    /* The lastly added vector to the data pyramid resembles the root level */
    return PerformRootLevelEncoding<T>(this->data_pyramid_.back());
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
CompressionVariableMultiData<T, DIM, N>::EncodeData()
{
    /* Store the encoding of the error mesh */
    this->error_mesh_encoding_ = this->error_mesh_->GetSerializedErrorMesh();
    this->error_mesh_encoding_size_bytes_ = this->error_mesh_->GetSerializedErrorMeshEncodingSizeBytes();

    /* Collect and exchange all entropy symbols */
    const std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>> entropy_symbols = ExchangeEntropySymbols<T>(this->level_elements_encodings_, this->comm_);
    
    /* Create a Huffman encoder */
    cmc::entropy_coding::huffman::HuffmanCoder<SymbolType> entropy_coder(entropy_symbols);

    /* Store the serialized Huffman coder */
    this->serialized_entropy_dictionary_ = entropy_coder.SerializeHuffmanCodesBEPadded();

    /* Number of overall encoding steps */
    const int num_encoding_steps = this->coarsening_indications_.size() + 1;

    /* Allocate an output vector for the levelwise encoding */
    this->levelwise_encoded_data_.reserve(num_encoding_steps);

    /* We need to encode the data resididng on the root level */
    this->levelwise_encoded_data_.push_back(this->EncodeRootLevelData());

    /* We insert a partition table on the root level at the beginning as well for completeness */
    this->level_partitioning_info_.insert(this->level_partitioning_info_.begin(), PartitionInfo(this->data_pyramid_.back().size(), this->levelwise_encoded_data_.back().size() * sizeof(uint64_t)));

    /* We encode the data from the root level to the leaf level (therefore, we access the indications in reverse, since they are collected during the coarsening) */
    auto lvl_iter = this->coarsening_indications_.rbegin();
    auto enc_iter = this->level_elements_encodings_.begin();

    /* Iterate over all compression levels (Without the root level and the intra element level) */
    for (size_t step_idx{1}; step_idx <= this->coarsening_indications_.size(); ++step_idx, ++lvl_iter, ++enc_iter)
    {
        /* Allocate a bits::vector to store this level's encoded data */
        cmc::bits::vector lvl_data;
        lvl_data.Reserve((enc_iter->size() * sizeof(ElemEncodingData<T, DIM>) * cmc::bits::kCharBit) / 2);

        /* We define a view onto refinement indications */
        cmc::bits::vector_view_in_memory lvl_view(*lvl_iter);

        /* Define a reference on the coarsening data for the ease of notation */
        const std::vector<ElemEncodingData<T, DIM>>& lvl_encoding_data = *enc_iter;

        /* Number of refinement indication bits on this level */
        const size_t num_elems = lvl_iter->size();
        int data_enc_idx{0};

        /* Iterate over this level's refinement indications */
        for (size_t elem_idx{0}; elem_idx < num_elems; ++elem_idx)
        {
            if (lvl_view.GetNextBit() == true)
            {
                /* Define a reference for the ease of notation */
                const ElemEncodingData<T, DIM>& enc_data = lvl_encoding_data[data_enc_idx];

                /* Encode the quantization bins and potentially interleave the unpredicted values */
                for (int child_elem_idx{0}; child_elem_idx < static_cast<int>(enc_data.num_elements); ++child_elem_idx)
                {
                    /* Encode the entropy symbol */
                    const cmc::entropy_coding::huffman::HuffmanCode code = entropy_coder.EncodeSymbol(enc_data.quantization_bins[child_elem_idx]);

                    /* Serialize the encoded entropy symbol */
                    lvl_data.AppendBits(code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - code.code_length), 0);
                
                    /* If the data is not predicatble, we need to store the actual value */
                    if (enc_data.quantization_bins[child_elem_idx] == kFlagUnpredictable) [[unlikely]]
                    {
                        /* Append the not-predictable value fully */
                        lvl_data.AppendBits(TransformToUInteger<T>(enc_data.unpredictable_values[child_elem_idx]), 0, 0);
                    }
                }

                /* Update the data encoding accessing index */
                ++data_enc_idx;
            }
        }

        /* At the end of the local encoding of the level, we append the process-end symbol */
        const cmc::entropy_coding::huffman::HuffmanCode process_lvl_end_code = entropy_coder.EncodeSymbol(kProcessEndSymbol);

        /* Serialize the encoded process end symbol */
        lvl_data.AppendBits(process_lvl_end_code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - process_lvl_end_code.code_length), 0);

        /* We store this level's encoding in the variable's buffer */
        levelwise_encoded_data_.push_back(lvl_data.GetSerializedByteStreamBE());

        /* We store the encoding length in the partition info */
        this->level_partitioning_info_[step_idx].num_bytes_encoding = levelwise_encoded_data_.back().size() * sizeof(uint64_t);
    }
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
CompressionVariableMultiData<T, DIM, N>::DetermineMaxInitElementLevel()
{
    /* If the maximum initial element level is not known, we need to gather it */
    if (this->max_init_elem_level_ == kMaxPresentElementLevelUnknown)
    {
        t8_forest_t mesh = mesh_.GetMesh();
        const t8_scheme_c* scheme =  t8_forest_get_scheme(mesh);

        const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees(mesh);

        int32_t max_elem_lvl{0};

        /* Iterate over all elements in all trees */
        for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx)
        {
            const t8_eclass_t tree_class = t8_forest_get_tree_class (mesh, tree_idx);
            const t8_locidx_t  num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (mesh, tree_idx);
            for (t8_locidx_t elem_idx = 0; elem_idx < num_elements_in_tree; ++elem_idx)
            {
                /* Get the current element */
                const t8_element_t* element = t8_forest_get_leaf_element_in_tree (mesh, tree_idx, elem_idx);

                /* Get the level of the element */
                const int32_t elem_level = scheme->element_get_level(tree_class, element);

                /* Check if the level is larger than the maximum previous level */
                if (max_elem_lvl < elem_level) [[unlikely]]
                {
                    max_elem_lvl = elem_level;
                }
            }
        }

        /* Exchange the maximum present refinement level */
        const int rv_allredc = MPI_Allreduce(&max_elem_lvl, &(this->max_init_elem_level_), 1, MPI_INT32_T, MPI_MAX, this->comm_);
        MPICheckError(rv_allredc);

        cmc_debug_msg("The maximum present element refinement level is ", this->max_init_elem_level_);
    }
}

template<ArithmeticType T>
struct DataOutput
{
    explicit DataOutput(const MPI_Offset byte_offset_, std::vector<uint64_t>&& data_stream_)
    : offset{byte_offset_}, data_stream(std::move(data_stream_)) {}

    MPI_Offset offset;
    std::vector<uint64_t> data_stream;
};

inline SizeType
GetVarHeaderSize(const int mesh_compression_levels, const bool has_intra_element_compression)
{
    static_assert(sizeof(uint8_t) == sizeof(char));
    cmc_assert(mesh_compression_levels > 0);
    return (sizeof(SizeType) //Global Number of Bytes (this value included)
            + sizeof(SizeType) //Num bytes until root level encoding starts (this value included)
            + kNumCharsVariableName * sizeof(char) //Variable Name
            + sizeof(SizeType) //CmcType = Data Type (e.g. float)
            + sizeof(SizeType) //Dimensionality
            + sizeof(SizeType) //Num Data per element
            + sizeof(SizeType) //Size of MPI communicator
            + sizeof(SizeType) //Compression Scheme
            + sizeof(SizeType) //Num Mesh Compression Levels
            + sizeof(SizeType) //Num Intra Element Compression Levels
            + sizeof(SizeType) //Intra Element Compression PackSize
            + (mesh_compression_levels + (has_intra_element_compression ? 1 : 0)) * sizeof(SizeType) //Global Bytes per Level
            + ((mesh_compression_levels - 1) * sizeof(SizeType)) //Global Element Indications per level
            + sizeof(SizeType) //Num Bytes Partition Table
            + sizeof(SizeType) //Num Bytes Serialized Huffman Codes Mesh Encoding Steps
            + sizeof(SizeType) //Num Bytes Serialized Huffman Codes Intra Element Encoding
            + sizeof(SizeType) //Num Bytes global mesh encoding
            + sizeof(SizeType) //Num Bytes Encoding of ErrorMesh
           );
}

inline SizeType
GetPartitionTableSize(const int mesh_compression_levels, const bool has_intra_element_compression, const int comm_size)
{
    cmc_assert(mesh_compression_levels > 0);
    /* We have an offset for each process on each level storing the mesh_elem_offset and the intra level coding byte offset */
    return sizeof(SizeType) * 2 * (mesh_compression_levels + (has_intra_element_compression ? 1 : 0)) * comm_size;
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
CompressionVariableMultiData<T, DIM, N>::AppendVariableHeaderToStream(std::vector<uint64_t>& stream, const SizeType global_byte_count, const SizeType start_root_level_encoding_offset,
                                                                      const int mesh_compression_levels, const std::vector<SizeType>& global_bytes_per_level)
{
    /* Number of mesh encoding steps */
    const int num_global_mesh_encoding_steps = this->levelwise_encoded_data_.size();

    /* Compute the number of bytes needed for the variable header */
    [[maybe_unused]] const SizeType num_bytes_var_header = GetVarHeaderSize(num_global_mesh_encoding_steps, this->HasIntraElementCompression());

    /** Fill the header **/
    /* Global Bytes compressed variable */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(global_byte_count));

    /* Store the offset from the stream start to the start of the root level encoding */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(start_root_level_encoding_offset));

    /* Variable Name */
    static_assert(8 * sizeof(char) == sizeof(uint64_t));
    cmc_assert(kNumCharsVariableName == sizeof(SizeType) * 32);
    if (this->name_.size() > kNumCharsVariableName)
    {
        /* Copy on the the number of bytes that are within the permitted range */
        cmc_warn_msg("The variable name '", this->name_, "' is too long; it gets trimmed to ", kNumCharsVariableName, " characters.");
    } 
    std::array<uint64_t, 32> var_name{}; //Filled with zeros, i.e. chars of '\0'
    int name_char_count{7}, val_idx{0};
    for (size_t idx{0}; idx < this->name_.size() && idx < kNumCharsVariableName; ++idx)
    {
        var_name[val_idx] |= (static_cast<SizeType>(this->name_[idx]) << (cmc::bits::kCharBit * name_char_count));
        --name_char_count;
        if (name_char_count < 0)
        {
            ++val_idx;
            name_char_count = 7;
        }
    }
    /* Copy the name to the stream */
    for (int idx{0}; idx < 32; ++idx)
    {
        stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(var_name[idx]));
    }

    /* Store the data type */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(ConvertToCmcType<T>())));

    /* Store the dimensionality */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(DIM)));

    /* Store the number of data per element */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(N)));

    /* Store the size of the MPI Communicator */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(this->comm_size_)));

    /* Store the compression scheme */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(CompressionSchema::ParallelMultiResExtraction)));

    /* Store the number of mesh compresion levels */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(mesh_compression_levels)));

    /* Store the number of intra element compression levels */
    //stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(ComputeNumIntraCompressionLevels<DIM, N>())));
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(0))); //TODO: Just a placeholder atm

    /* Store the utilized PackSize for intra element compression */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(kPackSize<DIM>)));

    /* Append the global level bytes */
    cmc_assert(static_cast<int>(global_bytes_per_level.size()) == mesh_compression_levels + (this->HasIntraElementCompression() ? 1 : 0));
    for (auto lvl_iter = global_bytes_per_level.begin(); lvl_iter != global_bytes_per_level.end(); ++lvl_iter)
    {
        stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(*lvl_iter)));
    }

    /* Apped the global num elements per level */
    cmc_assert(static_cast<int>(this->num_elements_per_level_.size()) == num_global_mesh_encoding_steps - 1);
    for (auto lvl_elem_iter = this->num_elements_per_level_.begin(); lvl_elem_iter != this->num_elements_per_level_.end(); ++lvl_elem_iter)
    {
        stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(*lvl_elem_iter)));
    }

    /* Store the number of bytes for the partition table */
    const SizeType partition_table_bytes = GetPartitionTableSize(mesh_compression_levels, this->HasIntraElementCompression(), this->comm_size_);
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(partition_table_bytes));

    /* Store the Huffman Code Lengths for the mesh encoding steps */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(this->serialized_entropy_dictionary_.size() * sizeof(uint64_t)));

    /* Store the Huffman Code Lengths for the intra element coding */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(this->serialized_intra_element_entropy_dictionary_.size() * sizeof(uint64_t)));

    /* Store the length of the global mesh encoding */
    cmc_assert(not this->global_mesh_encoding_.empty());
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(this->global_mesh_encoding_.size() * sizeof(uint64_t))));

    /* Store the length if the global error mesh encoding */
    stream.push_back(cmc::bits::ConvertToBigEndian<SizeType>(static_cast<SizeType>(this->error_mesh_encoding_size_bytes_)));
}

inline void
AppendPartitionTableToStream(std::vector<uint64_t>& root_rank_start_stream_no_frac, const std::vector<PartitionInfo>& global_partition_info, const int num_global_encoding_steps, const int comm_size)
{
    cmc_assert(num_global_encoding_steps * comm_size == static_cast<int>(global_partition_info.size()));

    /* We iterate over all levels and store the computed offset originating from the partitioning */
    for (int lvl_idx{0}; lvl_idx < num_global_encoding_steps; ++lvl_idx)
    {
        SizeType intra_level_elem_count_offset{0};
        SizeType intra_level_coding_byte_offset{0};

        /* We iterate through all ranks */
        for (int rank_id{0}; rank_id < comm_size; ++rank_id)
        {
            /* We store the current offset for this rank */
            root_rank_start_stream_no_frac.push_back(cmc::bits::ConvertToBigEndian<SizeType>(intra_level_elem_count_offset));
            root_rank_start_stream_no_frac.push_back(cmc::bits::ConvertToBigEndian<SizeType>(intra_level_coding_byte_offset));

            /* Compute the access index */
            const int access_idx = rank_id * num_global_encoding_steps + lvl_idx;

            /* And add the lengths from this rank to the offsets */
            intra_level_elem_count_offset += global_partition_info[access_idx].num_elems;
            intra_level_coding_byte_offset += global_partition_info[access_idx].num_bytes_encoding;
        }
    }
}

inline std::vector<SizeType>
ComputeGlobalBytesPerLevel(const std::vector<PartitionInfo>& global_partition_info, const int num_global_encoding_steps, const int comm_size)
{
    std::vector<SizeType> bytes_per_lvl;
    bytes_per_lvl.reserve(global_partition_info.size());

    for (int lvl_idx{0}; lvl_idx < num_global_encoding_steps; ++lvl_idx)
    {
        SizeType num_bytes{0};
        for (int rank_id{0}; rank_id < comm_size; ++rank_id)
        {
            const int access_idx = rank_id * num_global_encoding_steps + lvl_idx;
            num_bytes += global_partition_info[access_idx].num_bytes_encoding;
        }

        bytes_per_lvl.push_back(num_bytes);
    }

    return bytes_per_lvl;
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
std::vector<std::vector<uint64_t>>
CompressionVariableMultiData<T, DIM, N>::AdjustLevelMeshEncodings(const std::vector<PartitionInfo>& partition_info, const int num_global_mesh_encoding_steps , const int num_global_encoding_steps)
{
    cmc_assert(num_global_mesh_encoding_steps >= 1);
    const int num_coarsening_encodings = num_global_mesh_encoding_steps - 1;
    cmc_assert(num_coarsening_encodings == static_cast<int>(this->coarsening_indications_.size()));

    std::vector<std::vector<uint64_t>> offseted_level_mesh_encodings;
    offseted_level_mesh_encodings.reserve(num_coarsening_encodings);

    for (int lvl_idx{0}; lvl_idx < num_coarsening_encodings; ++lvl_idx)
    {
        SizeType lvl_elem_offset{0};

        /* Count this levels offset */
        for (int rank_idx{0}; rank_idx < this->comm_rank_; ++rank_idx)
        {
            const int access_idx = rank_idx * num_global_encoding_steps + lvl_idx + 1;
            lvl_elem_offset += partition_info[access_idx].num_elems;
        }

        /* We only store the mesh level encoding if there are actual signficant bits, otherwise we store an empty vector */
        if (this->coarsening_indications_[num_coarsening_encodings - 1 - lvl_idx].size() == 0) [[unlikely]]
        {
           offseted_level_mesh_encodings.emplace_back(std::vector<uint64_t>());
           continue;
        }

        /* Compute the local offset */
        const int shift = lvl_elem_offset % (sizeof(uint64_t) * cmc::bits::kCharBit);

        /* Get the offseted byte stream */
        offseted_level_mesh_encodings.emplace_back(this->coarsening_indications_[num_coarsening_encodings - 1 - lvl_idx].GetSerializedOffsetByteStreamBE(shift));
    }

    return offseted_level_mesh_encodings;
}

inline SizeType
ComputeGlobalMeshEncodingLength(const std::vector<SizeType> global_elems_per_level)
{
    SizeType num_bytes{0};
    for (auto lvl_iter = global_elems_per_level.begin(); lvl_iter != global_elems_per_level.end(); ++lvl_iter)
    {
        num_bytes += (*lvl_iter) / (sizeof(uint64_t) * cmc::bits::kCharBit) + ((*lvl_iter) % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0);
    }

    cmc::cmc_debug_msg("Number of bytes for global mesh encoding: ", num_bytes);
    
    return num_bytes;
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
CompressionVariableMultiData<T, DIM, N>::CollectMeshEncodingOnTheRootRank(const std::vector<PartitionInfo>& global_partition_info, const SizeType num_global_encoding_steps, const SizeType num_global_mesh_encoding_steps,  const std::vector<std::vector<uint64_t>>& offseted_level_mesh_encodings)
{
    /* Compute the number of mesh encoding levels */
    const int num_levels = this->coarsening_indications_.size();
    cmc_assert(num_levels == static_cast<int>(num_global_mesh_encoding_steps) - 1);

    /* Allocate an output vector */
    std::vector<SizeType> elems_per_level;
    elems_per_level.reserve(this->coarsening_indications_.size());

    cmc_assert(num_levels >= 1);

    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        SizeType global_elems_level{0};
        /* Count this levels offset */
        for (int rank_idx{0}; rank_idx < this->comm_size_; ++rank_idx)
        {
            const int access_idx = rank_idx * num_global_encoding_steps + lvl_idx + 1;
            global_elems_level += global_partition_info[access_idx].num_elems;
        }

        elems_per_level.push_back(global_elems_level);
    }

    if (this->comm_rank_ != kRootRank)
    {
        /* Count overall length of accumulated message */
        SizeType num_bytes_msg{0};
        for (const auto& lvl_mesh_encoding : offseted_level_mesh_encodings)
        {
            num_bytes_msg += lvl_mesh_encoding.size();
        }

        std::vector<uint64_t> message;
        message.reserve(num_bytes_msg);

        /* Setup the message accordingly */
        for (const auto& lvl_mesh_encoding : offseted_level_mesh_encodings)
        {
            std::copy_n(lvl_mesh_encoding.begin(), lvl_mesh_encoding.size(), std::back_inserter(message));
        }

        /* Send the message to the root rank */
        const int rv_send = MPI_Send(message.data(), message.size(), MPI_UINT64_T, kRootRank, kTagMeshEncoding, this->comm_);
        MPICheckError(rv_send);
    } else
    {
        /* Allocate a vector collecting the messages from the other ranks */
        const int num_expected_msgs = this->comm_size_ - 1;

        std::vector<std::vector<uint64_t>> recv_messages(num_expected_msgs);
        recv_messages.reserve(num_expected_msgs);


        /* Receive all messages from the other ranks collecting their process-local mesh encoding for all levels */
        for (int msg_idx{0}; msg_idx < num_expected_msgs; ++msg_idx)
        {
            /* Wait for a message to be received  */
            MPI_Status probe_status;
            const int rv_probe = MPI_Probe(MPI_ANY_SOURCE, kTagMeshEncoding, this->comm_, &probe_status);
            MPICheckError(rv_probe);

            /* Get the sending rank */
            const int source_rank = probe_status.MPI_SOURCE;

            /* Get the count */
            int num_elements{0};
            const int rv_count = MPI_Get_count(&probe_status, MPI_UINT64_T, &num_elements);
            MPICheckError(rv_count);

            /* Allocate the vector */
            recv_messages[source_rank - 1] = std::vector<uint64_t>(num_elements);

            /* Receive the message */
            const int rv_recv = MPI_Recv(recv_messages[source_rank - 1].data(), num_elements, MPI_UINT64_T, source_rank, kTagMeshEncoding, this->comm_, MPI_STATUS_IGNORE);
            MPICheckError(rv_recv);
        }
        
        /* After all expeceted messages have been received, we process them and create the global level-wise mesh encoding */

        /* Compute the number of bytes for the global level-wise mesh encoding */
        SizeType mesh_encoding_vals{0};
        for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
        {
            mesh_encoding_vals += (elems_per_level[lvl_idx] / (sizeof(uint64_t) * cmc::bits::kCharBit) + (elems_per_level[lvl_idx] % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0));
        }

        /* Allocate mesh encoding */
        std::vector<uint64_t> global_mesh_encoding;
        global_mesh_encoding.reserve(mesh_encoding_vals);

        std::vector<SizeType> level_rank_offsets(this->comm_size_ - 1, 0);

        for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
        {
            SizeType current_elem_offset{0};

            /* Copy the root-local part into the global mesh encoding */
            std::copy_n(offseted_level_mesh_encodings[lvl_idx].begin(), offseted_level_mesh_encodings[lvl_idx].size(), std::back_inserter(global_mesh_encoding));

            cmc_assert(not global_mesh_encoding.empty());
            if (global_mesh_encoding.empty()) [[unlikely]] {global_mesh_encoding.push_back(uint64_t{0});}
            
            /* Update the offset by the root local elements */
            current_elem_offset += global_partition_info[lvl_idx + 1].num_elems;

            /* Iterate over the root level from all other ranks and append their encodings */
            for (int rank_idx{1}; rank_idx < this->comm_size_; ++rank_idx)
            {
                const int access_idx = rank_idx * num_global_encoding_steps + lvl_idx + 1;

                const SizeType lvl_rank_num_elems = global_partition_info[access_idx].num_elems;

                if (lvl_rank_num_elems > 0) [[likely]]
                {
                    const SizeType lvl_rank_bits = lvl_rank_num_elems + (current_elem_offset % (sizeof(uint64_t) * cmc::bits::kCharBit));
                    //const SizeType lvl_rank_bits = lvl_rank_num_elems + ((sizeof(uint64_t) * cmc::bits::kCharBit) - (current_elem_offset % (sizeof(uint64_t) * cmc::bits::kCharBit))); //Test1

                    SizeType start_idx{0};
                    SizeType num_vals_to_copy = lvl_rank_bits / (sizeof(uint64_t) * cmc::bits::kCharBit) + (lvl_rank_bits % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0);

                    /* Check whether we need to interleave the first value */
                    if (current_elem_offset % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0) [[likely]]
                    {
                        /* Interleave the first byte */
                        global_mesh_encoding.back() |= recv_messages[rank_idx - 1].operator[](level_rank_offsets[rank_idx - 1]);
                        start_idx = 1;
                        --num_vals_to_copy;
                    }

                    if (num_vals_to_copy > 0) [[likely]]
                    {
                        /* Copy the remaining bytes over */
                        std::copy_n(recv_messages[rank_idx - 1].data() + level_rank_offsets[rank_idx - 1] + start_idx, num_vals_to_copy, std::back_inserter(global_mesh_encoding));
                    }

                    /* Update the element offset */
                    current_elem_offset += lvl_rank_num_elems;

                    /* Store the offset of this rank for the next level */
                    level_rank_offsets[rank_idx - 1] += (lvl_rank_bits / (sizeof(uint64_t) * cmc::bits::kCharBit) + (lvl_rank_bits % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0));
                }
            }
        }

        /* Store the global mesh encoding on the root rank */
        this->global_mesh_encoding_ = std::move(global_mesh_encoding);
    }

    /* Store the number of global elements per level */
    this->num_elements_per_level_ = std::move(elems_per_level);
}

//TODO: Make several intra element compression levels possible
template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
CompressionVariableMultiData<T, DIM, N>::GenerateOuputStreams()
{
    cmc_assert(this->output_streams_.empty() == true);

    /* Number of mesh encoding steps */
    const int num_global_mesh_encoding_steps = this->levelwise_encoded_data_.size();

    /* Number of overall encdoing steps */
    const int num_global_encoding_steps = this->levelwise_encoded_data_.size() + (this->HasIntraElementCompression() ? 1 : 0);

    cmc_assert(static_cast<size_t>(num_global_encoding_steps) == this->level_partitioning_info_.size());
    const int num_partition_infos = this->level_partitioning_info_.size();

    /*** Collect the partition table on the root rank ***/
    MPI_Datatype PartitionInfoMPIType;
    CreatePartitionInfoMPIType(&PartitionInfoMPIType);

    /* Allocate memory for the offsets and the encoding lengths */
    std::vector<PartitionInfo> global_partition_info(this->comm_size_ * num_partition_infos);

    /* Allgather the local offsets and encoding lengths */
    const int rv_allgather = MPI_Allgather(this->level_partitioning_info_.data(), num_partition_infos, PartitionInfoMPIType,
                                           global_partition_info.data(), num_partition_infos, PartitionInfoMPIType, this->comm_);
    MPICheckError(rv_allgather);

    /* Free the constructed MPI data type */
    const int rv_type_free = MPI_Type_free(&PartitionInfoMPIType);
    MPICheckError(rv_type_free);

    /* Exchange the mesh level encodings */
    std::vector<std::vector<uint64_t>> offseted_level_mesh_encodings = this->AdjustLevelMeshEncodings(global_partition_info, num_global_mesh_encoding_steps, num_global_encoding_steps);

    /* Collect the global level-wise mesh encoding on the root rank */
    this->CollectMeshEncodingOnTheRootRank(global_partition_info, num_global_encoding_steps, num_global_mesh_encoding_steps, offseted_level_mesh_encodings);

    /* Allocate a vector for the data output */
    this->output_streams_.reserve(num_global_encoding_steps);

    /* Determine the start_offset, e.g. the number of bytes that will be prepended by the root rank to the compressed data */
    const uint64_t start_encoding_offset_ = GetVarHeaderSize(num_global_mesh_encoding_steps, this->HasIntraElementCompression())
                                            + GetPartitionTableSize(num_global_mesh_encoding_steps, this->HasIntraElementCompression(), this->comm_size_)
                                            + this->serialized_entropy_dictionary_.size() * sizeof(uint64_t)
                                            + this->serialized_intra_element_entropy_dictionary_.size() * sizeof(uint64_t)
                                            + ComputeGlobalMeshEncodingLength(this->num_elements_per_level_) * sizeof(uint64_t)
                                            + this->error_mesh_encoding_size_bytes_;
    uint64_t current_byte_offset{start_encoding_offset_};

    /* Define an iterator to the local encoded data */
    auto lvl_iter = levelwise_encoded_data_.begin();

    /* Determine process-locally the offsets writing the data to the file */
    for (int mesh_lvl_idx{0}; mesh_lvl_idx < num_global_mesh_encoding_steps; ++mesh_lvl_idx)
    {
        /* On each mesh encoding level, each process has the opportunity to write data */
        for (int rank_id{0}; rank_id < this->comm_size_; ++rank_id)
        {
            /* Check if it is this processes turn to write data */
            if (rank_id == this->comm_rank_)
            {
                /* If we match the turn of the serialized ouput, we store the current offset */
                output_streams_.emplace_back(static_cast<MPI_Offset>(current_byte_offset), std::move(*lvl_iter));
                
                /* Check if the local length coincides with the length that has been dsitributed to the other ranks */
                cmc_assert(output_streams_.back().data_stream.size() * sizeof(uint64_t) == global_partition_info[this->comm_rank_ * num_global_encoding_steps + mesh_lvl_idx].num_bytes_encoding);

                /* Add the offset for the current data stream */
                current_byte_offset += output_streams_.back().data_stream.size() * sizeof(uint64_t);

                /* Overwrite the vector from which the data has been moved */
                *lvl_iter = std::vector<uint64_t>();

                /* Move to the next local encoded stream */
                ++lvl_iter;
            } else
            {
                /* In case, it is not this processes turn to write data, we just add to the offset */
                const int allgather_access_idx = rank_id * num_global_encoding_steps + mesh_lvl_idx;
                /* Add the length of the corresponding stream */
                current_byte_offset += global_partition_info[allgather_access_idx].num_bytes_encoding;
            }
        }
    }

    /* Potentially, add the intra element compression */
    if (this->HasIntraElementCompression())
    {
        /* Iterate over all ranks and let them store their intra element coding */
        for (int rank_id{0}; rank_id < this->comm_size_; ++rank_id)
        {
            /* Check if it is this processes turn to write data */
            if (rank_id == this->comm_rank_)
            {
                /* If we match the turn of the serialized ouput, we store the current offset */
                this->output_streams_.emplace_back(static_cast<MPI_Offset>(current_byte_offset), std::move(this->intra_element_encoded_data_));
                
                /* Check if the local length coincides with the length that has been dsitributed to the other ranks */
                cmc_assert(this->output_streams_.back().data_stream.size() * sizeof(uint64_t) == global_partition_info[this->comm_rank_ * num_global_encoding_steps + num_global_mesh_encoding_steps].num_bytes_encoding);

                /* Add the offset for the current data stream */
                current_byte_offset += this->output_streams_.back().data_stream.size() * sizeof(uint64_t);
                
                /* Overwrite the vector from which the data has been moved */
                this->intra_element_encoded_data_ = std::vector<uint64_t>();

                /* Move to the next local encoded stream */
                ++lvl_iter;
            } else
            {
                /* In case, it is not this processes turn to write data, we just add to the offset */
                const int allgather_access_idx = rank_id * num_global_encoding_steps + num_global_mesh_encoding_steps;
                /* Add the length of the corresponding stream */
                current_byte_offset += global_partition_info[allgather_access_idx].num_bytes_encoding;
            }
        }
    }

    /* Now after storing the offsets of each rank on each level, we are left with the overall byte count of this compressed variable */
    const SizeType global_byte_count = current_byte_offset;
    this->global_compressed_byte_count_ = global_byte_count;

    /* In case of the root rank, we need to exchange the root level encoding and offset,
     * with an updated stream consisting of the variable header, partition table, huffman codes and level offsets, etc. */
    if (this->comm_rank_ == kRootRank)
    {
        /* Compute the bytes per level */
        const std::vector<SizeType> global_bytes_per_level = ComputeGlobalBytesPerLevel(global_partition_info, num_global_encoding_steps, this->comm_size_);

        /* Allocate a new vector which collects the start stream for the root rank */
        std::vector<uint64_t> root_rank_start_stream;
        root_rank_start_stream.reserve(start_encoding_offset_ + this->output_streams_[0].data_stream.size());

        /* Append the variable header */
        this->AppendVariableHeaderToStream(root_rank_start_stream, global_byte_count, start_encoding_offset_, num_global_mesh_encoding_steps, global_bytes_per_level);

        /* Append the partition table */
        AppendPartitionTableToStream(root_rank_start_stream, global_partition_info, num_global_encoding_steps, this->comm_size_);

        /* Append the Huffman codes from the mesh encoding steps */
        std::copy_n(this->serialized_entropy_dictionary_.begin(), this->serialized_entropy_dictionary_.size(), std::back_inserter(root_rank_start_stream));

        /* Append the Huffman Codes from the mesh encoding steps */
        std::copy_n(this->serialized_intra_element_entropy_dictionary_.begin(), this->serialized_intra_element_entropy_dictionary_.size(), std::back_inserter(root_rank_start_stream));

        /* Append the mesh encoding */
        std::copy_n(this->global_mesh_encoding_.begin(), this->global_mesh_encoding_.size(), std::back_inserter(root_rank_start_stream));

        /* Append the error mesh encoding */
        cmc_assert(this->error_mesh_encoding_.size() * sizeof(uint64_t) == this->error_mesh_encoding_size_bytes_);
        std::copy_n(this->error_mesh_encoding_.begin(), this->error_mesh_encoding_.size(), std::back_inserter(root_rank_start_stream));

        /* Check the anticpiated offset for correctness */
        cmc_assert(start_encoding_offset_ == root_rank_start_stream.size() * sizeof(uint64_t));

        /* Append the root level encoding of the root rank */
        std::copy_n(this->output_streams_[0].data_stream.begin(), this->output_streams_[0].data_stream.size(), std::back_inserter(root_rank_start_stream));

        /* Store the adjusted offset and the new data stream for the root rank */
        this->output_streams_[0].offset = 0;
        this->output_streams_[0].data_stream = std::move(root_rank_start_stream);
    }

    /* Now, all offsets and encoding streams have been added, and the ouput data is fully prepared for writing */
}

inline void 
CheckMPIWriteCorrectness(const MPI_Status* status, const MPI_Datatype datatype, const int expected_num_elems)
{
    int elem_count_{0};
    const int rv_check_write = MPI_Get_count(status, datatype, &elem_count_);
    MPICheckError(rv_check_write);
    if (elem_count_ != expected_num_elems)
    {
        cmc_debug_msg("The expeceted number of bytes could not be written to the file.");
        MPICheckError(MPI_ERR_COUNT);
    }
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
void
CompressionVariableMultiData<T, DIM, N>::WriteCompressedData(const std::string& file_name)
{
    /* Generate the output streams for writing the data to disk */
    this->GenerateOuputStreams();

    cmc_assert(not this->output_streams_.empty());

    /* Check if the output file exists, if so, we delete it */
    const std::filesystem::path output_file_path(file_name);
    if (this->comm_rank_ == kRootRank && std::filesystem::exists(output_file_path))
    {
        const int rv_delete = MPI_File_delete(file_name.c_str(), MPI_INFO_NULL);
        MPICheckError(rv_delete);
    }

    const int rv_barrier_delete = MPI_Barrier(this->comm_);
    MPICheckError(rv_barrier_delete);

    /* Allocate a MPI file handle */
    MPI_File fhandle;

    /* Open the file */
    constexpr int opening_mode = MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL;
    const int rv_open = MPI_File_open(this->comm_, file_name.c_str(), opening_mode, MPI_INFO_NULL, &fhandle);
    MPICheckError(rv_open);

    cmc_debug_msg(this->comm_, "The file ", file_name, " has been opened.");

    /* Set the native data representation */
    const int rv_file_view = MPI_File_set_view(fhandle, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    MPICheckError(rv_file_view);

    /* Pre-Allocate global file storage */
    const uint64_t preallocation_size = this->global_compressed_byte_count_ + (sizeof(uint64_t) - (this->global_compressed_byte_count_ % sizeof(uint64_t)));
    const int rv_file_prealloc = MPI_File_preallocate(fhandle, preallocation_size);
    MPICheckError(rv_file_prealloc);
    cmc_debug_msg("The file has been pre-allocated.");

    /* Allocate a MPI status object */
    MPI_Status status;

    /* Itearte over all data streams and place them in the file */
    for (const auto& stream : this->output_streams_)
    {
        cmc_debug_msg("Compression-Output-Write: File-Offset: ", stream.offset, ", number of bytes: ", stream.data_stream.size() * sizeof(uint64_t));
        /* Move the file handle to the correct position within the file */
        const int rv_pos_fhandle = MPI_File_seek(fhandle, stream.offset, MPI_SEEK_SET);
        MPICheckError(rv_pos_fhandle);

        /* Write the current stream to the file */
        const int rv_write_data_buffer = MPI_File_write(fhandle, stream.data_stream.data(), stream.data_stream.size(), MPI_UINT64_T, &status);
        MPICheckError(rv_write_data_buffer);
        CheckMPIWriteCorrectness(&status, MPI_UINT64_T, static_cast<int>(stream.data_stream.size()));
    }

    /* Close the file */
    const int rv_close = MPI_File_close(&fhandle);
    MPICheckError(rv_close);

    cmc_debug_msg(this->comm_, "The file ", file_name, " has been closed.");
}

}

#endif /* !CMC_AMR_LOSSY_PAR_MULTI_RES_EXTRACTION_HXX */
