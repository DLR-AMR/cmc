#ifndef CMC_AC_COMPRESSION_VARIABLE_HXX
#define CMC_AC_COMPRESSION_VARIABLE_HXX

#include "cmc_config.h"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_variable_utilities.hxx"
#include "utilities/cmc_compression_settings.hxx"
#include "utilities/cmc_iface_compression_adapt_data.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx> 
#include <t8_forest/t8_forest_iterate.h> 
#endif

#include <vector>
#include <functional>

namespace cmc::lossy
{

template <typename T>
class AbstractCompressionVariable;

template<typename T>
using AdaptCreator = std::function<ICompressionAdaptData*(AbstractCompressionVariable<T>*,const CompressionSettings&)>;

template<typename T>
using AdaptDestructor = std::function<void(ICompressionAdaptData*)>;


template <typename T>
class AbstractCompressionVariable
{
public:
    void Compress(const CompressionSettings& settings);

    void WriteVTKFile(const std::string& file_name);

    void SetName(const std::string& name) {name_ = name;};
    const std::string& GetName() const {return name_;};

    size_t Size() const {return data_.size();};
    void PushBack(const T& value) {data_.push_back(value);};
    void PushBack(T&& value) {data_.push_back(std::move(value));};

    AmrMesh& GetAmrMesh() {return mesh_;};
    const AmrMesh& GetAmrMesh() const {return mesh_;};

    VariableUtilities<T>& GetVariableUtilities() {return utilities_;};
    const VariableUtilities<T>& GetVariableUtilities() const {return utilities_;};

    void SetInterpolation(Interpolate<T> interpolation_function) {utilities_.SetInterpolation(interpolation_function);};
    void SetInterpolationSkipMissingValues(InterpolateSkipMissingValues<T> interpolation_function) {utilities_.SetInterpolation(interpolation_function);};

    int EvaluateCoarsening(const int tree_id, const int elem_id, const int num_elements, const std::vector<PermittedError> permitted_errors);
    int EvaluateCoarseningSkipMissingValues(const int tree_id, const int elem_id, const int num_elements, const std::vector<PermittedError> permitted_errors, const T missing_value);

    int LeaveElementUnchanged(const int tree_id, const int elem_id);

    virtual ~AbstractCompressionVariable(){};
protected:
    AbstractCompressionVariable() = default;

    std::string name_; //!< The name of the variable
    std::vector<T> data_; //!< The actual data of the variable
    std::vector<T> data_new_; //!< A helper variable for the adaptation
    AmrMesh mesh_; //!< The mesh on which the variable is defined
    VariableUtilities<T> utilities_; //!< Utilities that are needed in order to track the errors
    AdaptCreator<T> adaptation_creator_; //!< A function pointer which is used to create the wished adaptation structure
    AdaptDestructor<T> adaptation_destructor_; //!< A function pointer which is used to destruct the adaptation structure

private:
    bool IsValidForCompression() const;
    void AllocateCoarseningIteration() {data_new_.reserve(mesh_.GetNumberLocalElements() / (2 << mesh_.GetDimensionality()) + 8);}
    void SwitchToAdaptedData() {data_.swap(data_new_); data_new_.clear();};
    ICompressionAdaptData* CreateAdaptData(const CompressionSettings& settings) {return adaptation_creator_(this, settings);};
    t8_forest_t RepartitionMesh(t8_forest_t adapted_forest);
    void RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest);
};

template <typename T>
inline void
AbstractCompressionVariable<T>::Compress(const CompressionSettings& settings)
{
    cmc_debug_msg("Compression of variable ", this->name_, " (by the means of adaptive coarsening) starts...");

    /* Set up the inaccuracy storage */
    utilities_.SetUpInaccuracyStorage(mesh_.GetNumberLocalElements());

    /* We create the adapt data based on the compression settings, the forest and the variables to consider during the adaptation/coarsening */
    ICompressionAdaptData* adapt_data = this->CreateAdaptData(settings);

    while (adapt_data->IsCompressionProgressing())
    {
        cmc_debug_msg("A coarsening iteration is initialized.");

        /* Initialize/Allocate for a coarsening iteration*/
        this->AllocateCoarseningIteration();
        adapt_data->InitializeCompressionIteration();

        /* Get and indicate to keep the 'previous forest' after the adaptation step */
        t8_forest_t previous_forest = mesh_.GetMesh();
        t8_forest_ref(previous_forest);

        /* Perform a coarsening iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, adapt_data->GetAdaptationFunction(), 0, 0, static_cast<void*>(adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements");

        /* Complete the interpolation step by storing the newly computed adapted data alongside its deviations */
        this->SwitchToAdaptedData();
        utilities_.SwitchDeviations();
        adapt_data->CompleteInterpolation(previous_forest, adapted_forest);

        /* Free the former forest */
        t8_forest_unref(&previous_forest);

        /* Repartition the mesh */
        t8_forest_t partitioned_forest = RepartitionMesh(adapted_forest);

        /* Repartition the data */
        this->RepartitionData(adapted_forest, partitioned_forest);
        adapt_data->RepartitionData(adapted_forest, partitioned_forest);

        cmc_debug_msg("The mesh and the data has been re-partitioned.");

        /* Free the former forest and store the adapted/repartitioned mesh */
        t8_forest_unref(&adapted_forest);
        mesh_.SetMesh(partitioned_forest);

        /* Finalize the comrpession iteration */
        adapt_data->FinalizeCompressionIteration();

        cmc_debug_msg("The coarsening iteration is finished.");
    }

    /* Free the adapt data structure */
    this->adaptation_destructor_(adapt_data);
    cmc_debug_msg("Compression of variable ", this->name_, " is finished.");
}

template <typename T>
inline t8_forest_t
AbstractCompressionVariable<T>::RepartitionMesh(t8_forest_t adapted_forest)
{
    /* Keep the not-partitioned forest */
    t8_forest_ref(adapted_forest);

    /* Allocate a forest */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    const int partition_for_coarsening = 0; //TODO: change to 'one' when partition for coarsening is in t8code
    t8_forest_set_partition(partitioned_forest, adapted_forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    return partitioned_forest;
}

template <typename T>
inline void
AbstractCompressionVariable<T>::RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest)
{
    /* Create an sc_array_t wrapper of the variable's data */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(data_.data()), sizeof(T), data_.size());

    cmc_debug_msg("Number of local data elements before partitioning: ", data_.size());
    cmc_debug_msg("Number of local mesh elements before partitioning: ", t8_forest_get_local_num_leaf_elements(adapted_forest));
    cmc_debug_msg("Size of a single data element: ", in_data->elem_size);

    /* Allocate memory for the partitioned data */
    const t8_locidx_t new_num_elems = t8_forest_get_local_num_leaf_elements(partitioned_forest);
    data_new_ = std::vector<T>(new_num_elems);

    cmc_debug_msg("Number of local data elements after partitioning: ", data_new_.size());
    cmc_debug_msg("Number of local mesh elements after partitioning: ", new_num_elems);

    /* Create a wrapper for the freshly allocated partitioned data */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(data_new_.data()), sizeof(T), data_new_.size());

    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_leaf_elements(adapted_forest)) == data_.size());
    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_leaf_elements(partitioned_forest)) == data_new_.size());

    /* Partition the variables data */
    t8_forest_partition_data(adapted_forest, partitioned_forest, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Set the variable's data to the newly partitioned data */
    SwitchToAdaptedData();
    cmc_debug_msg("Partitioning of mesh and data elements has been finished.");
    
    /* Repartition the inaccuracy tracker */
    utilities_.RepartitionInaccuracyData(adapted_forest, partitioned_forest);
}


template <typename T>
inline int
AbstractCompressionVariable<T>::EvaluateCoarsening(const int tree_id, const int local_elem_id, const int num_elements, const std::vector<PermittedError> permitted_errors)
{
    /* Compute the local element id in the data array */
    const int elem_id = t8_forest_get_tree_element_offset(mesh_.GetMesh(), tree_id) + local_elem_id;
    
    /* Obtain the values to be interpolated */
    const VectorView<T> values(data_.data() + elem_id, num_elements);

    /* Interpolate the values */
    const T interpolated_value = utilities_.Interpolation(values, mesh_.GetMesh(), tree_id, elem_id, num_elements);

    /* Get the deviations of the elements */
    const std::vector<double> previous_deviations = utilities_.GetPreviousDeviations(elem_id, num_elements);

    /* Evaluate whether the interpolation is error-bound-conformal */
    const ErrorCompliance evaluation = utilities_.IsCoarseningErrorCompliant(permitted_errors, values, previous_deviations, interpolated_value);

    if (evaluation.is_error_threshold_satisfied)
    {
        data_new_.push_back(interpolated_value);
        utilities_.StoreInaccuracy(elem_id, evaluation.max_introduced_error);
        return t8::kCoarsenElements;
    } else
    {
        data_new_.push_back(data_[elem_id]);
        utilities_.TransferPreviousDeviation(elem_id);
        return t8::kLeaveElementUnchanged;
    }
}

template <typename T>
inline int
AbstractCompressionVariable<T>::EvaluateCoarseningSkipMissingValues(const int tree_id, const int local_elem_id, const int num_elements, const std::vector<PermittedError> permitted_errors, const T missing_value)
{
    /* Compute the local element id in the data array */
    const int elem_id = t8_forest_get_tree_element_offset(mesh_.GetMesh(), tree_id) + local_elem_id;
    
    /* Obtain the values to be interpolated */
    const VectorView<T> values(data_.data() + elem_id, num_elements);

    /* Interpolate the values */
    const T interpolated_value = utilities_.Interpolation(values, mesh_.GetMesh(), tree_id, elem_id, num_elements, missing_value);

    /* Get the deviations of the elements */
    const std::vector<double> previous_deviations = utilities_.GetPreviousDeviations(elem_id, num_elements);

    /* Evaluate whether the interpolation is error-bound-conformal */
    const ErrorCompliance evaluation = utilities_.IsCoarseningErrorCompliantSkipMissingValues(permitted_errors, values, previous_deviations, interpolated_value, missing_value);

    if (evaluation.is_error_threshold_satisfied)
    {
        data_new_.push_back(interpolated_value);
        utilities_.StoreInaccuracy(elem_id, evaluation.max_introduced_error);
        return t8::kCoarsenElements;
    } else
    {
        data_new_.push_back(data_[elem_id]);
        utilities_.TransferPreviousDeviation(elem_id);
        return t8::kLeaveElementUnchanged;
    }
}

template <typename T>
inline int
AbstractCompressionVariable<T>::LeaveElementUnchanged(const int tree_id, const int local_elem_id) 
{
    /* Compute the local element id in the data array */
    const int elem_id = t8_forest_get_tree_element_offset(mesh_.GetMesh(), tree_id) + local_elem_id;
    
    /* Get the value */
    const T value = data_[elem_id];

    /* Store the unchanged value in the new data vector */
    data_new_.push_back(value);

    /* Transfer the previous inaccuracy */
    utilities_.TransferPreviousDeviation(elem_id);

    return t8::kLeaveElementUnchanged;
}

template <typename T>
inline bool
AbstractCompressionVariable<T>::IsValidForCompression() const 
{
    return true;
}

template <typename T>
inline void
AbstractCompressionVariable<T>::WriteVTKFile(const std::string& file_name)
{
    cmc_debug_msg("The variable ", this->name_, " has been written to the file: ", file_name);


}

}

#endif /* !CMC_AC_COMPRESSION_VARIABLE_HXX */
