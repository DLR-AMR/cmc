#ifndef LOSSLESS_CMC_BYTE_COMPRESSION_VARIABLE_HXX
#define LOSSLESS_CMC_BYTE_COMPRESSION_VARIABLE_HXX

#include "cmc_config.h"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_iface_abstract_byte_compression_variable.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_entropy_coder.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_compression_schema.hxx"
#include "mesh_compression/cmc_mesh_encoder.hxx"

#include "mpi/cmc_mpi.hxx"

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


namespace cmc::lossless
{

constexpr bool kCompressStrictlyLevelwise = true;

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

/**
 * @brief A strcut holding the data for an extraction process.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template<typename T>
struct ExtractionData
{
    ExtractionData(const CompressionValue<T>& coarse_val, std::vector<CompressionValue<T>>&& fine_vals)
    : coarse_value(coarse_val), fine_values(std::move(fine_vals)) {};
    ExtractionData(CompressionValue<T>&& coarse_val, std::vector<CompressionValue<T>>&& fine_vals)
    : coarse_value(std::move(coarse_val)), fine_values(std::move(fine_vals)) {};
    
    CompressionValue<T> coarse_value;
    std::vector<CompressionValue<T>> fine_values;
};

/**
 * @brief A strcut holding the data for an adapation process that leaves the element unchanged.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template<typename T>
struct UnchangedData
{
    UnchangedData(const CompressionValue<T>& coarse_val, const CompressionValue<T>& fine_val)
    : coarse_value(coarse_val), fine_value(fine_val) {};
    UnchangedData(CompressionValue<T>&& coarse_val, CompressionValue<T>&& fine_val)
    : coarse_value(std::move(coarse_val)), fine_value(std::move(fine_val)) {};
    
    CompressionValue<T> coarse_value;
    CompressionValue<T> fine_value;
};

/* Forward declarations */
template <typename T>
class AbstractByteCompressionVariable;
template <typename T>
class ICompressionAdaptData;

template<typename T>
using AdaptCreator = std::function<ICompressionAdaptData<T>*(AbstractByteCompressionVariable<T>*)>;

template<typename T>
using AdaptDestructor = std::function<void(ICompressionAdaptData<T>*)>;

/**
 * @brief The Interface/Template for a variable that performs lossless compression on the serialized data
 * in a byte-/bit-wise fashion. The compression algorithm is fixed and may be specialized in a derived class
 * with a derived adaptation data (\see ICompressionAdaptData) in order to fit the compression for the 
 * given needs.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template <typename T>
class AbstractByteCompressionVariable : public IByteCompressionVariable<T>
{
public:
    void Compress() override;

    const std::string& GetName() const override {return name_;};

    size_t Size() const override {return data_.size();};

    const AmrMesh& GetAmrMesh() const override {return mesh_;};

    virtual ~AbstractByteCompressionVariable(){};

    void MoveEncodedEntropyCodesInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_entropy_codes) override
    {
        vec_to_hold_encoded_levelwise_entropy_codes = std::move(buffered_entropy_codes_);
        buffered_entropy_codes_ = std::vector<std::vector<uint8_t>>();
    };

    void MoveEncodedDataInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_data) override
    {
        vec_to_hold_encoded_levelwise_data = std::move(buffered_encoded_data_);
        buffered_encoded_data_ = std::vector<std::vector<uint8_t>>();
    };

    void MoveEncodedMeshInto(std::vector<std::vector<uint8_t>>& vec_to_hold_encoded_levelwise_mesh) override
    {
        vec_to_hold_encoded_levelwise_mesh = std::move(buffered_encoded_mesh_);
        buffered_encoded_mesh_ = std::vector<std::vector<uint8_t>>();
    };

    const std::vector<std::vector<uint8_t>>& GetEncodedEntropyCodes() const override {return buffered_entropy_codes_;}
    const std::vector<std::vector<uint8_t>>& GetEncodedData() const override {return buffered_encoded_data_;};
    const std::vector<std::vector<uint8_t>>& GetEncodedMesh() const override {return buffered_encoded_mesh_;};

    MPI_Comm GetMPIComm() const override {return comm_;};
    int GetInitialMaximumRefinementLevel() const {return max_initial_refinement_level_;}
    virtual CompressionSchema GetCompressionSchema() const = 0;

    friend ICompressionAdaptData<T>;
protected:
    AbstractByteCompressionVariable() = default;

    void SetName(const std::string& name) {name_ = name;};
    void SetAmrMesh(const AmrMesh& mesh) {mesh_ = mesh;};
    void SetAmrMesh(AmrMesh&& mesh) {mesh_ = std::move(mesh);};
    void SetData(const std::vector<T>& initial_data);
    void SetData(const std::vector<SerializedCompressionValue<sizeof(T)>>& initial_data);
    void SetData(std::vector<SerializedCompressionValue<sizeof(T)>>&& initial_data);
    void SetMPIComm(const MPI_Comm comm) {comm_ = comm;};

    AdaptCreator<T> adaptation_creator_; //!< A function pointer which is used to create the wished adaptation structure
    AdaptDestructor<T> adaptation_destructor_; //!< A function pointer which is used to destruct the adaptation structure
    
    std::unique_ptr<mesh_compression::IMeshEncoder> mesh_encoder_{nullptr}; 
private:
    VectorView<CompressionValue<T>> GetView(const int start_index, const int count) const;
    VectorView<CompressionValue<T>> GetView(const int tree_id, const int lelement_index, const int count) const;

    void StoreExtractedValues(const int tree_id, const int lelement_id, const int num_elements, const ExtractionData<T>& extracted_values);
    void StoreUnchangedElement(const int tree_id, const int lelement_id, const UnchangedData<T>& extracted_value);

    void IndicateElementStaysUnchanged() {mesh_encoder_->IndicateElementStaysUnchanged();};
    void IndicateCoarsening() {mesh_encoder_->IndicateCoarsening();};

    bool IsValidForCompression() const;
    void AllocateExtractionIteration() {data_new_.reserve(mesh_.GetNumberLocalElements() / (2 << mesh_.GetDimensionality()) + 8);}
    void SwitchToExtractedData() {data_.swap(data_new_); data_new_.clear();};
    ICompressionAdaptData<T>* CreateAdaptData() {return adaptation_creator_(this);};
    t8_forest_t RepartitionMesh(t8_forest_t adapted_forest);
    void RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest);
    void DetermineInitialMaximumRefinementLevel();

    std::string name_; //!< The name of the variable
    AmrMesh mesh_; //!< The mesh on which the variable is defined
    std::vector<CompressionValue<T>> data_; //!< The current data of the variable 
    std::vector<CompressionValue<T>> data_new_; //!< A helper variable for the adaptation

    std::vector<std::vector<uint8_t>> buffered_entropy_codes_; //!< Level-wise storage of the entropy codes, e.g. LZC or prefix lengths
    std::vector<std::vector<uint8_t>> buffered_encoded_data_; //!< Level-wise storage for the encoded data
    std::vector<std::vector<uint8_t>> buffered_encoded_mesh_; //!< Level-wwise storage for the encoded mesh

    MPI_Comm comm_{MPI_COMM_NULL}; //!< The MPI communicator to use
    int max_initial_refinement_level_{-1};
};


/**
 * @brief Interface/Template for the adaptation data used within the lossless compression of the variable 
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template <typename T>
class ICompressionAdaptData
{
public:
    ICompressionAdaptData() = delete;
    ICompressionAdaptData(AbstractByteCompressionVariable<T>* variable)
    : base_variable_{variable} {};

    bool IsCompressionProgressing() const;

    virtual void InitializeExtractionIteration() = 0;
    virtual void FinalizeExtractionIteration() = 0;
    void FinalizingExtractionIteration();
    virtual void CompleteExtractionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) = 0;
    virtual void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) = 0;

    int ExtractValue(const int which_tree, const int lelement_id, const int num_elements);
    int LeaveElementUnchanged(const int which_tree, const int lelement_id);

    virtual std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeLevelData(const std::vector<CompressionValue<T>>& level_values) const = 0;
    virtual std::pair<std::vector<uint8_t>, std::vector<uint8_t>> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const = 0;
    
    MPI_Comm GetMPIComm() const {return base_variable_->GetMPIComm();};

    virtual ~ICompressionAdaptData(){};

    bool IsValidForCompression() const;

    int GetInitialMaximumRefinementLevel() const {return base_variable_->GetInitialMaximumRefinementLevel();}
    int GetCurrentCompressionStep() const {return compression_step_;}

protected:
    virtual ExtractionData<T> PerformExtraction(const int which_tree, const int lelement_id, const int num_elements, const VectorView<CompressionValue<T>> values) = 0;
    virtual UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) = 0;

    std::unique_ptr<entropy_coding::IByteCompressionEntropyCoder> entropy_coder_{nullptr}; //!< The entropy coder to use in order to encode information
private:
    AbstractByteCompressionVariable<T>* const base_variable_{nullptr};
    int compression_step_{0};
};

template <typename T>
inline bool
ICompressionAdaptData<T>::IsValidForCompression() const
{
    if (base_variable_ == nullptr)
    {
        cmc_err_msg("The pointer to the base variable is not set. Therefore, no compression can be applied.");
        return false;
    }
    if (entropy_coder_ == nullptr)
    {
        cmc_err_msg("The entropy coder is not set. Therefore, no compression can be applied.");
        return false;
    }

    return true;
}

/**
 * @brief This function indicates whether the compression continues or not.
 * The compression is carried out, until only the root elements of the mesh
 * are left over and further coarsening can be applied.
 * 
 * @tparam T The data type of the underlying data (e.g. float)
 * @return true If the compression continues
 * @return false If the compression is finished
 */
template <typename T>
bool
ICompressionAdaptData<T>::IsCompressionProgressing() const
{
    /* Perform the extraction until only the root elements are left */
    return (base_variable_->GetAmrMesh().GetNumberGlobalElements() > base_variable_->GetAmrMesh().GetNumberGlobalTrees());
}

/**
 * @brief This funciton is called during the compression and extracts a value which will be stored on the coarser level
 * (i.e. after coarsening the family of elements). Moreover, the remaining values on the finer level may be adjusted
 * as well. This function calls the custom "PerformExtraction" which needs to be implemented by the derived class.
 * 
 * @tparam T The data type of the underlying data (e.g. float)
 * @param which_tree The lcoal tree id from which the elements are taken
 * @param lelement_id The tree-local start index of the family of elements
 * @param num_elements The number of elements corresponding to this family
 * @return int The return value indicates that this family of elements will be coarsened
 */
template <typename T>
int
ICompressionAdaptData<T>::ExtractValue(const int which_tree, const int lelement_id, const int num_elements)
{
    cmc_assert(which_tree >= 0 && lelement_id >= 0);
    cmc_assert(num_elements > 1);

    /* We indicate to that a coarsening is applied to the family of elements */
    base_variable_->IndicateCoarsening();

    /* Get the corresponding values */
    const VectorView<CompressionValue<T>> values = base_variable_->GetView(which_tree, lelement_id, num_elements);

    /* Extract the coarse values, and potentially alter the remaining fine values */
    const ExtractionData<T> extracted_values = PerformExtraction(which_tree, lelement_id, num_elements, values);

    /* Store the extracted values wihtin the variable */
    base_variable_->StoreExtractedValues(which_tree, lelement_id, num_elements, extracted_values);

    return cmc::t8::kCoarsenElements;
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
ICompressionAdaptData<T>::LeaveElementUnchanged(const int which_tree, const int lelement_id)
{
    cmc_assert(which_tree >= 0 && lelement_id >= 0);

     /* We indicate to that the element stays unchanged */
    base_variable_->IndicateElementStaysUnchanged();

    /* Get the corresponding value */
    const VectorView<CompressionValue<T>> value = base_variable_->GetView(which_tree, lelement_id, 1);

    /* Leave the element unchanged */
    const UnchangedData<T> extracted_values = this->ElementStaysUnchanged(which_tree, lelement_id, value.front());

    /* Store the unchanged data */
    base_variable_->StoreUnchangedElement(which_tree, lelement_id, extracted_values);

    return cmc::t8::kLeaveElementUnchanged;
}

template <typename T>
void
ICompressionAdaptData<T>::FinalizingExtractionIteration()
{
    ++compression_step_;
    this->FinalizeExtractionIteration();
}


inline
bool CheckIfElementIsReadyForCoarsening(const int initial_max_ref_level, const int elem_level, const int compression_step)
{
    if (elem_level == initial_max_ref_level - compression_step)
    {
        return true;
    } else
    {
        return false;
    }
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
LosslessByteCompression (t8_forest_t forest,
                         [[maybe_unused]] t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         [[maybe_unused]] const t8_eclass_t tree_class,
                         t8_locidx_t lelement_id,
                         [[maybe_unused]] const t8_scheme_c * ts,
                         const int is_family,
                         const int num_elements,
                         [[maybe_unused]] t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    ICompressionAdaptData<T>* adapt_data = static_cast<ICompressionAdaptData<T>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Check if a family is supplied to the adaptation function */
    if (is_family == 0 || (kCompressStrictlyLevelwise && not CheckIfElementIsReadyForCoarsening(adapt_data->GetInitialMaximumRefinementLevel(), ts->element_get_level(tree_class, elements[0]), adapt_data->GetCurrentCompressionStep())))
    {
        /* If there is no family, the element stays unchanged */
        const int ret_val = adapt_data->LeaveElementUnchanged(which_tree, lelement_id);
        return ret_val;
    } else
    {
        /* Extract a value of the family and coarsen it */
        const int ret_val = adapt_data->ExtractValue(which_tree, lelement_id, num_elements);
        return ret_val;
    }
}

template <typename T>
inline void
AbstractByteCompressionVariable<T>::Compress()
{
    cmc_assert(this->IsValidForCompression());
    cmc_debug_msg("Lossless compression of variable ", this->name_, " starts...");

    if (kCompressStrictlyLevelwise)
    {
        /* Determine the maximum present refinement level */
        this->DetermineInitialMaximumRefinementLevel();
        cmc_assert(this->GetInitialMaximumRefinementLevel() > 0);
    }

    /* We create the adapt data based on the compression settings, the forest and the variables to consider during the adaptation/coarsening */
    ICompressionAdaptData<T>* adapt_data = this->CreateAdaptData();

    cmc_assert(adapt_data->IsValidForCompression());

    while (adapt_data->IsCompressionProgressing())
    {
        cmc_debug_msg("A coarsening iteration is initialized.");

        /* Initialize/Allocate for a coarsening iteration */
        this->AllocateExtractionIteration();
        adapt_data->InitializeExtractionIteration();
        mesh_encoder_->IntializeCompressionIteration(mesh_.GetNumberLocalElements());

        /* Get and indicate to keep the 'previous forest' after the adaptation step */
        t8_forest_t previous_forest = mesh_.GetMesh();
        t8_forest_ref(previous_forest);

        /* Perform a coarsening iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, LosslessByteCompression<T>, 0, 0, static_cast<void*>(adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements");

        /* Complete the interpolation step by storing the newly computed adapted data alongside its deviations */
        adapt_data->CompleteExtractionIteration(previous_forest, adapted_forest);

        /* Free the former forest */
        t8_forest_unref(&previous_forest);

        /* Encode the data of this level and store it within a buffer */
        auto [encoded_entropy_codes, encoded_data] = adapt_data->EncodeLevelData(data_);
        buffered_entropy_codes_.push_back(std::move(encoded_entropy_codes));
        buffered_encoded_data_.push_back(std::move(encoded_data));

        /* Once the data is buffered, we can overwrite it with the adapted data */
        this->SwitchToExtractedData();

        /* Repartition the mesh */
        t8_forest_t partitioned_forest = RepartitionMesh(adapted_forest);

        /* Get the encoded the mesh adapatations */
        std::vector<uint8_t> encoded_mesh_data = mesh_encoder_->GetPartitionedEncodedLevelData(adapted_forest, partitioned_forest, this->GetMPIComm());
        buffered_encoded_mesh_.push_back(std::move(encoded_mesh_data));

        /* Repartition the data */
        this->RepartitionData(adapted_forest, partitioned_forest);
        adapt_data->RepartitionData(adapted_forest, partitioned_forest);

        cmc_debug_msg("The mesh and the data has been re-partitioned.");

        /* Free the former forest and store the adapted/repartitioned mesh */
        t8_forest_unref(&adapted_forest);
        mesh_.SetMesh(partitioned_forest);

        /* Finalize the comrpession iteration */
        adapt_data->FinalizingExtractionIteration();
        mesh_encoder_->FinalizeCompressionIteration();
        cmc_debug_msg("The coarsening iteration is finished.");
    }

    /* At last, we need to encode the root level data */
    auto [encoded_root_entropy_codes, encoded_root_data] = adapt_data->EncodeRootLevelData(data_);
    buffered_entropy_codes_.push_back(std::move(encoded_root_entropy_codes));
    buffered_encoded_data_.push_back(std::move(encoded_root_data));

    /* At last, we need to encode the root level of the mesh */
    std::vector<uint8_t> encoded_root_mesh = mesh_encoder_->EncodeRootLevelMesh(mesh_.GetMesh());
    buffered_encoded_mesh_.push_back(std::move(encoded_root_mesh));

    /* Free the adapt data structure */
    this->adaptation_destructor_(adapt_data);
    cmc_debug_msg("Compression of variable ", this->name_, " is finished.");
}

template <typename T>
void
AbstractByteCompressionVariable<T>::SetData(const std::vector<T>& initial_data)
{
    data_.reserve(initial_data.size());

    for (auto val_iter = initial_data.begin(); val_iter != initial_data.end(); ++val_iter)
    {
        data_.emplace_back(CompressionValue<T>(*val_iter));
    }
}

template <typename T>
void
AbstractByteCompressionVariable<T>::SetData(const std::vector<SerializedCompressionValue<sizeof(T)>>& initial_data)
{
    data_ = initial_data;
}

template <typename T>
void
AbstractByteCompressionVariable<T>::SetData(std::vector<SerializedCompressionValue<sizeof(T)>>&& initial_data)
{
    data_ = std::move(initial_data);
}

template <typename T>
VectorView<CompressionValue<T>>
AbstractByteCompressionVariable<T>::GetView(const int start_index, const int count) const
{
    cmc_assert(start_index >= 0 && count >= 0);
    cmc_assert(static_cast<size_t>(start_index + count) <= data_.size());

    return VectorView(&data_[start_index], count);
}

template <typename T>
VectorView<CompressionValue<T>>
AbstractByteCompressionVariable<T>::GetView(const int tree_id, const int lelement_index, const int count) const
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
AbstractByteCompressionVariable<T>::StoreExtractedValues(const int tree_id, const int lelement_id, const int num_elements, const ExtractionData<T>& extracted_values)
{
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));
    cmc_assert(static_cast<size_t>(num_elements) == extracted_values.fine_values.size());

    /* Compute the start offset in the local contiguous array */
    const int elem_start_index = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_id;

    cmc_assert(elem_start_index >= 0 && num_elements >= 0);
    cmc_assert(static_cast<size_t>(elem_start_index + num_elements) <= data_.size());

    /* Store the adjusted fine values */
    std::copy_n(extracted_values.fine_values.begin(), num_elements, &data_[elem_start_index]);

    /* Store the extracted coarse value for the next coarsening iteration */
    data_new_.push_back(extracted_values.coarse_value);
}

template <typename T>
void
AbstractByteCompressionVariable<T>::StoreUnchangedElement(const int tree_id, const int lelement_id, const UnchangedData<T>& extracted_value)
{
    cmc_assert(tree_id >= 0 && tree_id < t8_forest_get_num_local_trees(mesh_.GetMesh()));

    /* Get the element id in the contiguous array of all local elements */
    const int elem_id = t8_forest_get_tree_element_offset (mesh_.GetMesh(), tree_id) + lelement_id;
    
    cmc_assert(elem_id >= 0 && static_cast<size_t>(elem_id) < data_.size());

    /* Set the adjusted remaining element value */
    data_[elem_id] = extracted_value.fine_value;

    /* Set the new value for the next coarsening iteration */
    data_new_.push_back(extracted_value.coarse_value);
}


template <typename T>
inline t8_forest_t
AbstractByteCompressionVariable<T>::RepartitionMesh(t8_forest_t adapted_forest)
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
AbstractByteCompressionVariable<T>::RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest)
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
    SwitchToExtractedData();
    cmc_debug_msg("Partitioning of mesh and data elements has been finished.");
}

template <typename T>
inline bool
AbstractByteCompressionVariable<T>::IsValidForCompression() const 
{
    if (name_.empty())
    {
        cmc_err_msg("The variable needs a name in order to store the compressed output. Therefore, no compression can be applied.");
        return false;
    }
    if (data_.empty())
    {
        cmc_err_msg("There is no data attached to the variable. Therefore, no compression can be applied.");
        return false;
    }
    if (not mesh_.IsValid())
    {
        cmc_err_msg("The mesh is not valid. Therefore, no compression can be applied.");
        return false;
    }
    if (mesh_encoder_ == nullptr)
    {
        cmc_err_msg("The mesh encoder is not set. Therefore, no compression can be applied.");
        return false;
    }

#ifdef CMC_ENABLE_MPI
    if (comm_ == MPI_COMM_NULL)
    {
        cmc_err_msg("The MPI-Communicator is not set. Therefore, no compression can be applied.");
        return false;
    }

    MPI_Comm mesh_comm = t8_forest_get_mpicomm(mesh_.GetMesh());
    int are_mpi_comms_equal{0};
    const int ret_val_comm_compare = MPI_Comm_compare(comm_, mesh_comm, &are_mpi_comms_equal);
    MPICheckError(ret_val_comm_compare);
    
    if (are_mpi_comms_equal != MPI_IDENT)
    {
        cmc_err_msg("The MPI-Communicator of the mesh and the variable differs. Therefore, no compression can be applied.");
        return false;
    }

#endif

    return true;
}


template <typename T>
void
AbstractByteCompressionVariable<T>::DetermineInitialMaximumRefinementLevel()
{
    t8_forest_t mesh = this->GetAmrMesh().GetMesh();
    const t8_scheme_c* scheme =  t8_forest_get_scheme(mesh);

    const t8_locidx_t num_local_trees = t8_forest_get_num_local_trees(mesh);
    int val_idx = 0;

    /* Iterate over all elements in all trees */
    for (t8_locidx_t tree_idx = 0; tree_idx < num_local_trees; ++tree_idx)
    {
        const t8_eclass_t tree_class = t8_forest_get_tree_class (mesh, tree_idx);
        const t8_locidx_t  num_elements_in_tree = t8_forest_get_tree_num_leaf_elements (mesh, tree_idx);
        for (t8_locidx_t elem_idx = 0; elem_idx < num_elements_in_tree; ++elem_idx, ++val_idx)
        {
            /* Get the current element */
            const t8_element_t* element = t8_forest_get_leaf_element_in_tree (mesh, tree_idx, elem_idx);

            const int elem_level = scheme->element_get_level(tree_class, element);

            if (max_initial_refinement_level_ < elem_level)
            {
                max_initial_refinement_level_ = elem_level;
            }
        }
    }

    cmc_debug_msg("The maximum present element refinement level is ", max_initial_refinement_level_);
}

}

#endif /* !LOSSLESS_CMC_BYTE_COMPRESSION_VARIABLE_HXX */
