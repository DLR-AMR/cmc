#ifndef DECOMPRESSION_CMC_BYTE_DECOMPRESSION_VARIABLE_HXX
#define DECOMPRESSION_CMC_BYTE_DECOMPRESSION_VARIABLE_HXX

#include "cmc_config.h"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_entropy_coder.hxx"
#include "utilities/cmc_byte_compression_arithmetic_encoding.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_compression_schema.hxx"
#include "mesh_compression/cmc_iface_mesh_decoder.hxx"

#ifdef CMC_WITH_T8CODE
#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx> 
#include <t8_forest/t8_forest_iterate.h> 
#include <t8_forest/t8_forest_partition.h>
#endif

#include <vector>
#include <functional>
#include <memory>
#include <algorithm>


namespace cmc::decompression
{

constexpr bool kWriteDecompressionStepToVTK = false;

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

/**
 * @brief A struct holding the data for an extraction process.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template<typename T>
struct RefinementData
{
    RefinementData() = default;
    RefinementData(std::vector<CompressionValue<T>>&& fine_vals)
    : fine_values(std::move(fine_vals)) {};
    RefinementData(const std::vector<CompressionValue<T>>& fine_vals)
    : fine_values(fine_vals) {};
    
    std::vector<CompressionValue<T>> fine_values;
};

/**
 * @brief A struct holding the data for an adapation process that leaves the element unchanged.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template<typename T>
struct UnchangedData
{
    UnchangedData(const CompressionValue<T>& fine_val)
    : fine_value(fine_val) {};
    UnchangedData(CompressionValue<T>&& fine_val)
    : fine_value(std::move(fine_val)) {};
    
    CompressionValue<T> fine_value;
};

/* Forward declarations */
template <typename T>
class AbstractByteDecompressionVariable;
template <typename T>
class IDecompressionAdaptData;

template<typename T>
using AdaptCreator = std::function<IDecompressionAdaptData<T>*(AbstractByteDecompressionVariable<T>*)>;

template<typename T>
using AdaptDestructor = std::function<void(IDecompressionAdaptData<T>*)>;

/**
 * @brief The Interface/Template for a variable that performs lossless compression on the serialized data
 * in a byte-/bit-wise fashion. The compression algorithm is fixed and may be specialized in a derived class
 * with a derived adaptation data (\see ICompressionAdaptData) in order to fit the compression for the 
 * given needs.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template <typename T>
class AbstractByteDecompressionVariable
{
public:
    void Decompress();

    const std::string& GetName() const {return name_;};

    size_t Size() const {return data_.size();};

    const AmrMesh& GetAmrMesh() const {return mesh_;};

    virtual ~AbstractByteDecompressionVariable(){};

    const std::vector<CompressionValue<T>>& GetDecompressedData() const {return data_;};
    
    int GetMaxNumDecompressionIterations() const {return max_num_decompression_iterations_;}

    friend IDecompressionAdaptData<T>;
protected:
    AbstractByteDecompressionVariable() = delete;
    explicit AbstractByteDecompressionVariable(std::vector<uint8_t>&& encoded_data_byte_stream, std::vector<uint8_t>&& encoded_mesh_byte_stream, const int max_num_decompression_iterations)
    : encoded_data_byte_stream_(std::move(encoded_data_byte_stream)), encoded_mesh_byte_stream_(std::move(encoded_mesh_byte_stream)), max_num_decompression_iterations_{max_num_decompression_iterations} {};

    void SetName(const std::string& name) {name_ = name;};
    void SetAmrMesh(const AmrMesh& mesh) {mesh_ = mesh;};
    void SetAmrMesh(AmrMesh&& mesh) {mesh_ = std::move(mesh);};
    void SetData(const std::vector<T>& initial_data);
    void SetData(const std::vector<SerializedCompressionValue<sizeof(T)>>& initial_data);
    void SetData(std::vector<SerializedCompressionValue<sizeof(T)>>&& initial_data);

    const std::vector<uint8_t>& GetEncodedDataStream() const {return encoded_data_byte_stream_;};
    const std::vector<uint8_t>& GetEncodedMeshStream() const {return encoded_mesh_byte_stream_;};

    AdaptCreator<T> adaptation_creator_; //!< A function pointer which is used to create the wished adaptation structure
    AdaptDestructor<T> adaptation_destructor_; //!< A function pointer which is used to destruct the adaptation structure

    std::unique_ptr<mesh_compression::IMeshDecoder> mesh_decoder_{nullptr};

private:
    VectorView<CompressionValue<T>> GetView(const int start_index, const int count) const;
    VectorView<CompressionValue<T>> GetView(const int tree_id, const int lelement_index, const int count) const;
    CompressionValue<T> GetValue(const int tree_id, const int lelement_index) const;

    bool WillNextElementBeRefined() {cmc_assert(mesh_decoder_ != nullptr); return mesh_decoder_->WillNextElementBeRefined();};
    void StoreRefinedValues(const RefinementData<T>& refiend_values);
    void StoreUnchangedElement(const UnchangedData<T>& unchanged_value);

    bool IsValidForDecompression() const;
    void AllocateDecompressionIteration() {data_new_.reserve(mesh_.GetNumberLocalElements() * (2 << mesh_.GetDimensionality()));}
    void SwitchToDecompressedData() {data_.swap(data_new_); data_new_.clear();};
    IDecompressionAdaptData<T>* CreateAdaptData() {return adaptation_creator_(this);};
    t8_forest_t RepartitionMesh(t8_forest_t adapted_forest);
    void RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest);

    std::string name_; //!< The name of the variable
    AmrMesh mesh_; //!< The mesh on which the variable is defined
    std::vector<CompressionValue<T>> data_; //!< The current data of the variable 
    std::vector<CompressionValue<T>> data_new_; //!< A helper variable for the adaptation

    int max_num_decompression_iterations_{0};
    const std::vector<uint8_t> encoded_data_byte_stream_; //!< The encoded byte stream of the variable
    const std::vector<uint8_t> encoded_mesh_byte_stream_; //!< The encoded byte stream of the mesh
};


/**
 * @brief Interface/Template for the adaptation data used within the lossless compression of the variable 
 * 
 * @tparam T The original data type of the underlying data (e.g. float)
 */
template <typename T>
class IDecompressionAdaptData
{
public:
    IDecompressionAdaptData() = delete;
    IDecompressionAdaptData(AbstractByteDecompressionVariable<T>* variable)
    : base_variable_{variable}, encoded_data_byte_stream_{variable->encoded_data_byte_stream_} {};

    virtual bool IsDecompressionProgressing() const = 0;

    virtual std::vector<CompressionValue<T>> DecodeRootLevel(const t8_locidx_t num_local_root_values) = 0;
    virtual void InitializeDecompressionIteration() = 0;
    virtual void FinalizeDecompressionIteration() = 0;
    virtual void CompleteDecompressionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) = 0;
    virtual void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) = 0;

    int ApplyDecompression(const int which_tree, const int lelement_id, const int num_refined_elements);
    int LeaveElementUnchanged(const int which_tree, const int lelement_id);
    bool WillNextElementBeRefined();
    virtual ~IDecompressionAdaptData(){};

    bool IsValidForDeompression() const;
    const AmrMesh& GetAmrMesh() const {return base_variable_->GetAmrMesh();}
    int GetMaxNumDecompressionIterations() const {return base_variable_->GetMaxNumDecompressionIterations();}
    const std::vector<CompressionValue<T>>& GetDecompressedData() const {return base_variable_->data_;};
private:
    AbstractByteDecompressionVariable<T>* const base_variable_{nullptr};

protected:
    virtual RefinementData<T> PerformRefinement(const int which_tree, const int lelement_id, const CompressionValue<T> value, const int num_refined_elements) = 0;
    virtual UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) = 0;

    const std::vector<uint8_t>& encoded_data_byte_stream_;
};

template <typename T>
inline bool
IDecompressionAdaptData<T>::IsValidForDeompression() const
{
    if (base_variable_ == nullptr)
    {
        cmc_err_msg("The pointer to the base variable is not set. Therefore, no decompression can be applied.");
        return false;
    }
    if (encoded_data_byte_stream_.empty())
    {
        cmc_err_msg("The encoded data is emtpy. Therefore, no compression can be applied.");
        return false;
    }

    return true;
}

template <typename T>
inline bool
IDecompressionAdaptData<T>::WillNextElementBeRefined()
{
    return base_variable_->WillNextElementBeRefined();
}

/**
 * @brief This funciton is called during the compression and extracts a value which will be stored on the coarser level
 * (i.e. after coarsening the family of elements). Moreover, the remaining values on the finer level may be adjusted
 * as well. This function calls the custom "PerformRefinement" which needs to be implemented by the derived class.
 * 
 * @tparam T The data type of the underlying data (e.g. float)
 * @param which_tree The lcoal tree id from which the elements are taken
 * @param lelement_id The tree-local start index of the family of elements
 * @param num_elements The number of elements corresponding to this family
 * @return int The return value indicates that this family of elements will be coarsened
 */
template <typename T>
int
IDecompressionAdaptData<T>::ApplyDecompression(const int which_tree, const int lelement_id, const int num_refined_elements)
{
    cmc_assert(which_tree >= 0 && lelement_id >= 0 && num_refined_elements > 1);

    /* Get the corresponding values */
    const CompressionValue<T> value = base_variable_->GetValue(which_tree, lelement_id);

    /* Extract the coarse values, and potentially alter the remaining fine values */
    const RefinementData<T> refined_values = PerformRefinement(which_tree, lelement_id, value, num_refined_elements);

    /* Store the extracted values wihtin the variable */
    base_variable_->StoreRefinedValues(refined_values);

    return cmc::t8::kRefineElement;
}

/**
 * @brief This function is called during the compression and supplies the (potential) altered value for the
 * element (which remains unchanged in the mmesh) after the adaptation which will be stored for the next adaptation
 * iteration. Moreover the left-over value remaining "in the old data vector" can be altered as well, if needed.
 * This function calls the custom "LeaveElementUnchanged" which needs to be implemented by the derived class.
 * 
 * @tparam T The data type of the underlying data
 * @param which_tree The lcoal tree id from which the element is taken
 * @param lelement_id The tree-local index of the element
 * @return int The return value indicates that this element will remain unchanged
 */
template <typename T>
int
IDecompressionAdaptData<T>::LeaveElementUnchanged(const int which_tree, const int lelement_id)
{
    cmc_assert(which_tree >= 0 && lelement_id >= 0);

    /* Get the corresponding value */
    const CompressionValue<T> value = base_variable_->GetValue(which_tree, lelement_id);

    /* Leave the element unchanged */
    const UnchangedData<T> unchanged_value = this->ElementStaysUnchanged(which_tree, lelement_id, value);

    /* Store the unchanged data */
    base_variable_->StoreUnchangedElement(unchanged_value);

    return cmc::t8::kLeaveElementUnchanged;
}

/**
 * @brief The adaptation function which is used for the lossless compression variables.
 * In case a family is passed to this callback, an extraction process is always performed.
 * 
 * @return t8_locidx_t Indicates whether the element stays unchanged or if the family of elements
 * will be coarsened
 */
template<typename T>
inline t8_locidx_t
ByteVariableDecompressionAdaptation (t8_forest_t forest,
                                     t8_forest_t forest_from,
                                     t8_locidx_t which_tree,
                                     const t8_eclass_t tree_class,
                                     t8_locidx_t lelement_id,
                                     const t8_scheme_c * ts,
                                     [[maybe_unused]] const int is_family,
                                     [[maybe_unused]] const int num_elements,
                                     t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    IDecompressionAdaptData<T>* adapt_data = static_cast<IDecompressionAdaptData<T>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Determine whether the element will be refined */
    const bool is_refine_element = adapt_data->WillNextElementBeRefined();

    if (is_refine_element)
    {
        /* If the element will be refined */
        /* Get the number of children elements into which this element will be refined */
        const int num_children = ts->element_get_num_children(tree_class, elements[0]);
        const int ret_val = adapt_data->ApplyDecompression(which_tree, lelement_id, num_children);
        return ret_val;
    } else
    {
        /* If the element stays unchanged */
        const int ret_val = adapt_data->LeaveElementUnchanged(which_tree, lelement_id);
        return ret_val;
    }
}

static int step_count = 0;

template <typename T>
void
WriteData(t8_forest_t forest, const std::vector<CompressionValue<T>>& data)
{
    std::vector<double> double_data;
    double_data.reserve(data.size());

    for (auto val_iter = data.begin(); val_iter != data.end(); ++val_iter)
    {
        double_data.push_back(static_cast<double>(val_iter->template ReinterpretDataAs<T>()));
    }

    t8_vtk_data_field_t vtk_data[1];
    snprintf (vtk_data[0].description, BUFSIZ, "ExampleData");
    vtk_data[0].type = T8_VTK_SCALAR;
    vtk_data[0].data = double_data.data();

    std::string name = std::string("cmc_decompression_hybrid_circlesquare_example_mesh_and_data_step_") + std::to_string(step_count);
    t8_forest_write_vtk_ext (forest, name.c_str(), 0, 0, 0, 0, 0, 0, 0, 1, vtk_data);
    ++step_count;
}

template <typename T>
inline void
AbstractByteDecompressionVariable<T>::Decompress()
{
    cmc_assert(this->IsValidForDecompression());
    cmc_debug_msg("Decompression of variable ", this->name_, " starts...");

    /* Re-Create the base mesh of the variable */
    t8_forest_t base_mesh = mesh_decoder_->DecodeRootMesh();
    mesh_.SetMesh(base_mesh);

    /* We create the adapt data based on the compression settings, the forest and the variables to consider during the adaptation/coarsening */
    IDecompressionAdaptData<T>* adapt_data = this->CreateAdaptData();

    /* Decode the root level values */
    data_ = adapt_data->DecodeRootLevel(mesh_.GetNumberLocalElements());

    if constexpr (kWriteDecompressionStepToVTK)
    {    
        WriteData<T>(mesh_.GetMesh(), data_);
    }

    /* Perform decompression iterations until the mesh is completely reconstructed */
    while (mesh_decoder_->IsDecompressionProgressing())
    {
        cmc_debug_msg("A decompression iteration is initialized.");

        /* Initialize/Allocate for a decompression iteration */
        this->AllocateDecompressionIteration();
        adapt_data->InitializeDecompressionIteration();
        mesh_decoder_->IntializeDecompressionIteration();

        /* Get and indicate to keep the 'previous forest' after the adaptation step */
        t8_forest_t previous_forest = mesh_.GetMesh();
        t8_forest_ref(previous_forest);

        /* Perform a decompression iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, ByteVariableDecompressionAdaptation<T>, 0, 0, static_cast<void*>(adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements");

        /* Complete the interpolation step by storing the newly computed adapted data alongside its deviations */
        adapt_data->CompleteDecompressionIteration(previous_forest, adapted_forest);

        /* Free the former forest */
        t8_forest_unref(&previous_forest);

        /* Switch to the decompressed data */
        this->SwitchToDecompressedData();

        /* Repartition the mesh */
        t8_forest_t partitioned_forest = RepartitionMesh(adapted_forest);

        /* Repartition the data */
        this->RepartitionData(adapted_forest, partitioned_forest);
        adapt_data->RepartitionData(adapted_forest, partitioned_forest);

        cmc_debug_msg("The mesh and the data has been re-partitioned.");

        /* Free the former forest and store the adapted/repartitioned mesh */
        t8_forest_unref(&adapted_forest);
        mesh_.SetMesh(partitioned_forest);

        /* Finalize the decompression iteration */
        adapt_data->FinalizeDecompressionIteration();
        mesh_decoder_->FinalizeDecompressionIteration();

        cmc_debug_msg("The decompression iteration is finished.");
        
        if constexpr (kWriteDecompressionStepToVTK)
        {
            WriteData<T>(mesh_.GetMesh(), data_);
            cmc_debug_msg("The decompression step has been written out in a .vtu file.");
        }
    }

    /* Free the adapt data structure */
    this->adaptation_destructor_(adapt_data);
    cmc_debug_msg("Decompression of variable ", this->name_, " is finished.");
}

template <typename T>
void
AbstractByteDecompressionVariable<T>::SetData(const std::vector<T>& initial_data)
{
    data_.reserve(initial_data.size());

    for (auto val_iter = initial_data.begin(); val_iter != initial_data.end(); ++val_iter)
    {
        data_.emplace_back(CompressionValue<T>(*val_iter));
    }
}

template <typename T>
void
AbstractByteDecompressionVariable<T>::SetData(const std::vector<SerializedCompressionValue<sizeof(T)>>& initial_data)
{
    data_ = initial_data;
}

template <typename T>
void
AbstractByteDecompressionVariable<T>::SetData(std::vector<SerializedCompressionValue<sizeof(T)>>&& initial_data)
{
    data_ = std::move(initial_data);
}

template <typename T>
CompressionValue<T>
AbstractByteDecompressionVariable<T>::GetValue(const int tree_id, const int lelement_index) const
{
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));

    /* Compute the start offset in the local contiguous array */
    const int elem_index = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_index;

    cmc_assert(elem_index >= 0);
    cmc_assert(static_cast<size_t>(elem_index) <= data_.size());

    return data_[elem_index];
}

template <typename T>
VectorView<CompressionValue<T>>
AbstractByteDecompressionVariable<T>::GetView(const int start_index, const int count) const
{
    cmc_assert(start_index >= 0 && count >= 0);
    cmc_assert(static_cast<size_t>(start_index + count) <= data_.size());

    return VectorView(&data_[start_index], count);
}

template <typename T>
VectorView<CompressionValue<T>>
AbstractByteDecompressionVariable<T>::GetView(const int tree_id, const int lelement_index, const int count) const
{
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));

    /* Compute the start offset in the local contiguous array */
    const int start_index = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_index;

    cmc_assert(start_index >= 0 && count >= 0);
    cmc_assert(static_cast<size_t>(start_index + count) <= data_.size());

    return VectorView(&data_[start_index], count);

}

template <typename T>
void
AbstractByteDecompressionVariable<T>::StoreRefinedValues(const RefinementData<T>& refined_values)
{
    /* Store the refined values */
    std::copy_n(refined_values.fine_values.begin(), refined_values.fine_values.size(), std::back_inserter(data_new_));
}

template <typename T>
void
AbstractByteDecompressionVariable<T>::StoreUnchangedElement(const UnchangedData<T>& unchanged_value)
{
    /* Store the element value */
    data_new_.push_back(unchanged_value.fine_value);
}


template <typename T>
inline t8_forest_t
AbstractByteDecompressionVariable<T>::RepartitionMesh(t8_forest_t adapted_forest)
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
AbstractByteDecompressionVariable<T>::RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest)
{
    /* Create an sc_array_t wrapper of the variable's data */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(data_.data()), sizeof(CompressionValue<T>), data_.size());

    cmc_debug_msg("Number of local data elements before partitioning: ", data_.size());
    cmc_debug_msg("Number of local mesh elements before partitioning: ", t8_forest_get_local_num_leaf_elements(adapted_forest));
    cmc_debug_msg("Size of a single data element: ", in_data->elem_size);

    /* Allocate memory for the partitioned data */
    const t8_locidx_t new_num_elems = t8_forest_get_local_num_leaf_elements(partitioned_forest);
    data_new_ = std::vector<CompressionValue<T>>(new_num_elems);

    cmc_debug_msg("Number of local data elements after partitioning: ", data_new_.size());
    cmc_debug_msg("Number of local mesh elements after partitioning: ", new_num_elems);

    /* Create a wrapper for the freshly allocated partitioned data */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(data_new_.data()), sizeof(CompressionValue<T>), data_new_.size());

    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_leaf_elements(adapted_forest)) == data_.size());
    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_leaf_elements(partitioned_forest)) == data_new_.size());

    /* Partition the variables data */
    t8_forest_partition_data(adapted_forest, partitioned_forest, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Set the variable's data to the newly partitioned data */
    SwitchToDecompressedData();
    cmc_debug_msg("Partitioning of mesh and data elements has been finished.");
}

template <typename T>
inline bool
AbstractByteDecompressionVariable<T>::IsValidForDecompression() const 
{
    if (name_.empty())
    {
        cmc_err_msg("The variable needs a name. Therefore, no decompression can be applied.");
        return false;
    }
    if (encoded_data_byte_stream_.empty())
    {
        cmc_err_msg("There is no compressed data attached to the variable. Therefore, no decompression can be applied.");
        return false;
    }
    if (encoded_mesh_byte_stream_.empty())
    {
        cmc_err_msg("There is no compressed mesh attached to the variable. Therefore, no decompression can be applied.");
        return false;
    }
    if (mesh_decoder_ == nullptr)
    {
        cmc_err_msg("The mesh decoder is not set. Therefore, no decompression can be applied.");
        return false;
    }

    return true;
}

}

#endif /* !DECOMPRESSION_CMC_BYTE_DECOMPRESSION_VARIABLE_HXX */
