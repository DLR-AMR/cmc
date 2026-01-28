#ifndef LOSSLESS_CMC_BYTE_PAR_DECOMPRESSION_VARIABLE_HXX
#define LOSSLESS_CMC_BYTE_PAR_DECOMPRESSION_VARIABLE_HXX

#include "cmc_config.h"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_huffman_coder.hxx"
#include "utilities/cmc_byte_compression_arithmetic_encoding.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_compression_schema.hxx"
#include "mesh_compression/cmc_iface_mesh_decoder.hxx"

#include <t8.h>
#include <t8_forest/t8_forest.h>
#include <t8_forest/t8_forest_general.h>
#include <t8_schemes/t8_scheme.hxx> 
#include <t8_forest/t8_forest_iterate.h> 
#include <t8_forest/t8_forest_partition.h>

#include <vector>
#include <functional>
#include <memory>
#include <algorithm>
#include <bit>

namespace cmc::decompression::par
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

struct OffsetHint
{
    uint64_t entropy_code_id{0};
    uint64_t global_level_byte_offset{0};
};

struct ProcLevelByteStreamOffsets
{
    std::vector<OffsetHint>> offset_hints;
};

/* Forward declarations */
template <typename T>
class AbstractByteParDecompressionVariable;
template <typename T>
class IParDecompressionAdaptData;

template<typename T>
using AdaptCreator = std::function<IParDecompressionAdaptData<T>*(AbstractByteParDecompressionVariable<T>*)>;

template<typename T>
using AdaptDestructor = std::function<void(IParDecompressionAdaptData<T>*)>;

/**
 * @brief The Interface/Template for a variable that performs lossless compression on the serialized data
 * in a byte-/bit-wise fashion. The compression algorithm is fixed and may be specialized in a derived class
 * with a derived adaptation data (\see ICompressionAdaptData) in order to fit the compression for the 
 * given needs.
 * 
 * @tparam T The origianl data type of the underlying data (e.g. float)
 */
template <typename T>
class AbstractByteParDecompressionVariable
{
public:
    void Decompress(const t8_cmesh_t initial_mesh, const t8_scheme *scheme);
    void DecompressToLevel(const t8_cmesh_t initial_mesh, const t8_scheme *scheme, const int level);

    const std::string& GetName() const {return name_;};

    size_t Size() const {return data_.size();};

    const AmrMesh& GetAmrMesh() const {return mesh_;};

    virtual ~AbstractByteParDecompressionVariable(){};

    const std::vector<CompressionValue<T>>& GetDecompressedData() const {return data_;};
    
    int GetMaxNumDecompressionIterations() const {return max_num_decompression_iterations_;}

    virtual void SetupLevelDecodingStart(const uint8_t* global_var_start, const t8_gloidx_t proc_lvl_elem_offset, const ProcLevelByteStreamOffsets& offset_hints);
    friend IParDecompressionAdaptData<T>;
protected:
    AbstractByteParDecompressionVariable() = delete;
    explicit AbstractByteParDecompressionVariable(std::vector<uint8_t>&& encoded_data_byte_stream, std::vector<uint8_t>&& encoded_mesh_byte_stream, const int max_num_decompression_iterations)
    : encoded_data_byte_stream_(std::move(encoded_data_byte_stream)), encoded_mesh_byte_stream_(std::move(encoded_mesh_byte_stream)), max_num_decompression_iterations_{max_num_decompression_iterations} {};

    explicit AbstractByteParDecompressionVariable(const std::string& name, std::vector<uint8_t>&& global_level_num_elems, std::vector<uint8_t>&& encoded_mesh_stream, 
                                                  std::vector<uint8_t>&& global_level_data_bytes, const uint64_t file_byte_offset_encoded_data, 
                                                  std::vector<ProcLevelByteStreamOffsets>&& level_offset_hints, const std::string& file_name, const MPI_Comm comm)
    : name_(name), global_level_num_elems_(std::move(global_level_num_elems)), encoded_mesh_byte_stream_(std::move(encoded_mesh_stream)),
      global_level_bytes_(std::move(global_level_data_bytes)), file_var_encoded_data_offset_{file_byte_offset_encoded_data}, levelwise_memoffset_hints_(std::move(level_offset_hints)),
      file_name_(file_name), comm_{comm}
    {
        /* Check the amount of levels in the mesh is the same as the amount of data levels*/
        if (global_level_num_elems_.size() != global_level_bytes_.size())
        {
            throw std::invalid_argument("The amount of encoded levels of the mesh and the data encoding do not coincide!");
        }
        /* Check if the compressed file exists */
        if (const std::filesystem::path input_file_path(file_name_); not std::filesystem::exists(input_file_path))
        {
            throw std::invalid_argument("The compressed file does not exist!");
        }
        /* Check if an MPI Communicator is given */
        if (comm_ == MPI_COMM_NULL)
        {
            throw std::invalid_argument("The MPI_Communicator is NULL!");
        }

        max_num_decompression_iterations_ = global_level_num_elems_.size();
    }


    void SetName(const std::string& name) {name_ = name;};
    void SetAmrMesh(const AmrMesh& mesh) {mesh_ = mesh;};
    void SetAmrMesh(AmrMesh&& mesh) {mesh_ = std::move(mesh);};
    void SetData(const std::vector<T>& initial_data);
    void SetData(const std::vector<SerializedCompressionValue<sizeof(T)>>& initial_data);
    void SetData(std::vector<SerializedCompressionValue<sizeof(T)>>&& initial_data);

    const uint8_t* GetEncodedMeshStreamPtr() const {return encoded_mesh_byte_stream_.data();};

    virtual std::vector<T> SetupRootLevelData(const uint8_t* root_level_encoding_ptr, const t8_gloidx_t proc_elem_offset, const t8_gloidx_t proc_elem_count) = 0;

    AdaptCreator<T> adaptation_creator_; //!< A function pointer which is used to create the wished adaptation structure
    AdaptDestructor<T> adaptation_destructor_; //!< A function pointer which is used to destruct the adaptation structure

    std::unique_ptr<mesh_compression::IMeshParDecoder> mesh_decoder_{nullptr};

private:
    VectorView<CompressionValue<T>> GetView(const int start_index, const int count) const;
    VectorView<CompressionValue<T>> GetView(const int tree_id, const int lelement_index, const int count) const;
    CompressionValue<T> GetValue(const int tree_id, const int lelement_index) const;
    t8_forest_t SetupInitialMesh(const t8_cmesh_t cmesh, const t8_scheme *scheme) const;

    bool WillNextElementBeRefined() {cmc_assert(mesh_decoder_ != nullptr); return mesh_decoder_->WillNextElementBeRefined();};
    void StoreRefinedValues(const RefinementData<T>& refiend_values);
    void StoreUnchangedElement(const UnchangedData<T>& unchanged_value);

    bool IsValidForDecompression() const;
    void AllocateDecompressionIteration() {data_new_.reserve(mesh_.GetNumberLocalElements() * (2 << mesh_.GetDimensionality()));}
    void SwitchToDecompressedData() {data_.swap(data_new_); data_new_.clear();};
    IParDecompressionAdaptData<T>* CreateAdaptData() {return adaptation_creator_(this);};
    t8_forest_t RepartitionMesh(t8_forest_t adapted_forest);
    void RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest);

    std::string name_; //!< The name of the variable
    AmrMesh mesh_; //!< The mesh on which the variable is defined

    std::vector<CompressionValue<T>> data_; //!< The current data of the variable 
    std::vector<CompressionValue<T>> data_new_; //!< A helper variable for the adaptation

    std::vector<SizeType> global_level_num_elems_; //!< Indicating the number of global elements per level
    /* The mesh stream will be held completely in memory since it is rather small */
    const std::vector<uint8_t> encoded_mesh_byte_stream_; //!< The encoded byte stream of the mesh
    
    std::vector<SizeType> global_level_bytes_;  //!< Indicating the number of encoded data bytes per level
    const SizeType file_var_encoded_data_offset_; //!< The byte offset in the file to the start of the encoded level data of this variable 

    /* The byte stream indicating hints for the process-local offsets will be held in memory as well since it is rather small */
    const std::vector<ProcLevelByteStreamOffsets> levelwise_memory_hints_;

    std::string file_name_; //!< The name of the file from where the encoded level data will be retrieved 

    MPI_Comm comm_{MPI_COMM_NULL};
    int comm_rank_{-1};
    MPI_Comm shm_comm_{MPI_COMM_NULL};
    int shm_rank_{0}, shm_size_{1};
    MPI_Win lvl_window_;
    int max_num_decompression_iterations_{0};
};


/**
 * @brief Interface/Template for the adaptation data used within the lossless compression of the variable 
 * 
 * @tparam T The original data type of the underlying data (e.g. float)
 */
template <typename T>
class IParDecompressionAdaptData
{
public:
    IParDecompressionAdaptData() = delete;
    IParDecompressionAdaptData(AbstractByteParDecompressionVariable<T>* variable)
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
    virtual ~IParDecompressionAdaptData(){};

    bool IsValidForDeompression() const;
    const AmrMesh& GetAmrMesh() const {return base_variable_->GetAmrMesh();}
    int GetMaxNumDecompressionIterations() const {return base_variable_->GetMaxNumDecompressionIterations();}
    const std::vector<CompressionValue<T>>& GetDecompressedData() const {return base_variable_->data_;};
private:
    AbstractByteParDecompressionVariable<T>* const base_variable_{nullptr};

protected:
    virtual RefinementData<T> PerformRefinement(const int which_tree, const int lelement_id, const CompressionValue<T> value, const int num_refined_elements) = 0;
    virtual UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) = 0;

    const std::vector<uint8_t>& encoded_data_byte_stream_;
};


template <typename T>
inline bool
IParDecompressionAdaptData<T>::IsValidForDeompression() const
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
IParDecompressionAdaptData<T>::WillNextElementBeRefined()
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
IParDecompressionAdaptData<T>::ApplyDecompression(const int which_tree, const int lelement_id, const int num_refined_elements)
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
IParDecompressionAdaptData<T>::LeaveElementUnchanged(const int which_tree, const int lelement_id)
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
ByteVariableParDecompressionAdaptation (t8_forest_t forest,
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
    IParDecompressionAdaptData<T>* adapt_data = static_cast<IParDecompressionAdaptData<T>*>(t8_forest_get_user_data(forest));
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


template <typename T>
void
AbstractByteParDecompressionVariable<T>::DecompressToLevel(const t8_cmesh_t cmesh, const t8_scheme *scheme, const int level)
{
    cmc_assert(level >= 0 && level + 1 <= max_num_decompression_iterations_);

    max_num_decompression_iterations_ = level;

    /* Check if the specified decompression level is possible */
    if (level < 0 || level + 1 > max_num_decompression_iterations_)
    {
        cmc_err_msg("The specified decompression level (", level,") is not in range of all possible decompression levels [0, ", max_num_decompression_iterations_ - 1, "].");
    }

    /* Update the number of compression iterations */
    max_num_decompression_iterations_ = level + 1;

    /* Perform the decompression */
    this->Decompress(cmesh, scheme);
}

template <typename T>
inline t8_forest_t
AbstractByteParDecompressionVariable<T>::SetupInitialMesh(const t8_cmesh_t cmesh, const t8_scheme *scheme) const
{
    /* Create a forest from the cmesh */
    t8_forest_t forest = t8_forest_new_uniform (cmesh, scheme, 0, 1, comm_);
    
    /* ALlocate a new forest */
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
AbstractByteParDecompressionVariable<T>::CreateSharedMemoryComms()
{
    /* Split communicator into shared memory groups */
    const int rv_split_comm = MPI_Comm_split_type(this->comm_, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &(this->shm_comm_));
    MPICheckError(rv_split_comm);

    /* Get the rank within and the size of the shared communicator */
    const int rv_shm_rank_ = MPI_Comm_rank(this->shm_comm_, &(this->shm_rank_));
    MPICheckError(rv_shm_rank_);
    const int rv_shm_size_ = MPI_Comm_size(this->shm_comm_, &(this->shm_size_));
    MPICheckError(rv_shm_size_);
}

/* Read the data from the globally endoded level into a shared memmory array in parallel to store and access it in a shared manner to reduce memory consumption */
template <typename T>
std::tuple<const uint8_t*, t8_gloidx_t, t8_gloidx_t>
AbstractByteParDecompressionVariable<T>::OpenSharedLevelDataWindow(const int level)
{
    /* Get the global offset of the first element */
    const t8_gloidx_t mesh_offset = t8_forest_get_first_local_leaf_element_id(mesh_.GetMesh());

    /* Get the number of local elements */
    const t8_gloidx_t num_local_elems = static_cast<t8_gloidx_t>(t8_forest_get_local_num_leaf_elements(mesh_.GetMesh()));
    
    cmc_assert(global_level_bytes_.size() > static_cast<size_t>(level));

    /* Compute equal distributions for the whole data level */
    const SizeType num_global_lvl_bytes = global_level_bytes_[level];
    const SizeType proc_offset_byte_stream = static_cast<SizeType>(((static_cast<double>(shm_rank_) * static_cast<long double>(num_global_lvl_bytes)) / static_cast<double>(shm_size_)));
    const SizeType next_proc_offset_byte_stream = static_cast<SizeType>(((static_cast<double>(shm_rank_ + 1) * static_cast<long double>(num_global_lvl_bytes)) / static_cast<double>(shm_size_)));
    cmc_assert(proc_offset_byte_stream <= next_proc_offset_byte_stream);
    const int byte_stream_length = static_cast<int>(next_proc_offset_byte_stream - proc_offset_byte_stream);

    /* Store the offset to the start of the contiguous memeory of the global level data */
    const MPI_Aint global_level_start_offset = 0 - proc_offset_byte_stream;

    /* Allocate a window on the shared memory communicators */
    uint8_t* shm_mem{nullptr};
    const int rv_win_alloc = MPI_Win_allocate_shared(static_cast<MPI_Aint>(byte_stream_length), sizeof(uint8_t), MPI_INFO_NULL, this->shm_comm_, &shm_mem, &lvl_window_);
    MPICheckError(rv_win_alloc);

    /* Open the compressed file for reading */
    MPI_File fhandle;
    const int rv_open = MPI_File_open(this->shm_comm_, this->file_name_.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fhandle);
    MPICheckError(rv_open);

    cmc_debug_msg(this->shm_comm_, "The file ", this->file_name_, " has been opened (for comm ", this->shm_comm_, ").");

    /* Move to the correct position in the file for the current level on the current process */
    const MPI_Offset var_level_offset = std::accumulate(this->global_level_bytes_.begin(), std::next(this->global_level_bytes_.begin(), level), 0);
    const MPI_Offset file_offset = this->file_var_encoded_data_offset_ + var_level_offset + proc_offset_byte_stream;
    const int rv_seek_start = MPI_File_seek(fhandle, 0, MPI_SEEK_SET);
    MPICheckError(rv_seek_start);

    /* Read the corresponding data */
    MPI_Status status;
    const int rv_read_data = MPI_File_read(fhandle, shm_mem, byte_stream_length, MPI_UINT8_T, &status);
    MPICheckError(rv_read_data);
    CheckMPIReadCorrectness(&status, MPI_UINT8_T, byte_stream_length);

    cmc_debug_msg(this->shm_comm_, "The process (from comm ", this->shm_comm_, ") has read a part of the data encoding of level ", level, ".");

    /* Close the compressed file */
    const int rv_close = MPI_File_close(&fhandle);
    MPICheckError(rv_close);

    cmc_debug_msg(this->shm_comm_, "The file ", this->file_name_, " has been closed (for comm ", this->shm_comm_, ").");

    /* Next, we synchronize the window */
    const int rv_sync_win = MPI_Win_sync(this->lvl_window_);
    MPICheckError(rv_sync_win);

    /* Impose an additional barrier on the shared communicator */
    const int rv_shm_barrier = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barrier);

    /* Setup the decoding start for this level correctly */
    const uint8_t* global_var_start = shm_mem - global_level_start_offset;

    return std::make_tuple(global_var_start, mesh_offset, num_local_elems);
}

template <typename T>
void
AbstractByteParDecompressionVariable<T>::CloseSharedLevelDataWindow()
{
    /* Impose a barrier to finish all outstanding reads */
    const int rv_shm_barrier = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barrier);

    /* Free the allocated window */
    const int rv_win_free = MPI_Win_free(&(this->lvl_window_));
    MPICheckError(rv_win_free);
}

template <typename T>
SizeType
AbstractByteParDecompressionVariable<T>::DetermineEntropyOffsetsFromMeshOffset(bit_map::BitMapView mesh_lvl_encoding, const size_t proc_elem_offset, [[maybe_unused]] const size_t proc_elem_count)
{
    static_assert(sizeof(uint8_t) == 1);

    #if 0
    //The following code would be helpful, if we opt for the parallel version 
    /* Compute the global start of the local chunk */
    const size_t start_byte = proc_elem_offset / bit_map::kCharBit;
    const size_t start_bit = proc_elem_offset % bit_map::kCharBit;

    /* Set the mask to nullify the bits that belong to the previous process */
    uint8_t UpdateStartByte{0xFF};
    if (start_bit != 0)
    {
        /* In this case the start bit is not the start of a byte but lies somewhere within */
        UpdateStartByte <<= (start_bit - 1);
    }

    /* Compute the global end of the local chunk */
    const size_t end_byte = (proc_elem_offset + proc_elem_count) / bit_map::kCharBit;
    const size_t end_bit = (proc_elem_offset + proc_elem_count) % bit_map::kCharBit;

    /* Set the mask to nullify the bits that belong to the next process */
    uint8_t UpdateEndByte{0xFF};
    if (end_bit != 7)
    {
        UpdateEndByte >>= (bit_map::kCharBit - 1 - end_bit);
    }
    #endif
    
    //Code for the serial iteration through the encoded mesh 
    #if 0
    if (proc_elem_offset == 0)
    {
        /* In case the offset is zero, we are having the start of the encoding, and therefore, there is no entropy offset */
        return 0;
    }

    const size_t prev_proc_elem_offset = (proc_elem_offset != 0 ? proc_elem_offset - 1 : 0);

    /* Compute the global end of the element preceeding the start of this process */
    const size_t prev_proc_end_byte = prev_proc_elem_offset / bit_map::kCharBit;
    const size_t prev_proc_end_bit = prev_proc_elem_offset % bit_map::kCharBit;
    uint8_t UpdateEndByte{0xFF};
    if (prev_proc_end_bit != 7)
    {
        /* Set the mask to nullify the bits that belong to this process */
        UpdateEndByte >>= (bit_map::kCharBit - 1 - prev_proc_end_bit);
    }

    SizeType num_refined_elems{0};

    /* Iterate over all bytes and count the elements that will be refined despite the last byte holding the end bit of the previous process */
    for (size_t byte_iter{0}; byte_iter < prev_proc_end_byte; ++byte_iter)
    {
        num_refined_elems += std::popcount(mesh_lvl_encoding[byte_iter]);
    }

    /* Add up the refined elements from the last byte */
    num_refined_elems += std::popcount(mesh_lvl_encoding[prev_proc_end_byte] & UpdateEndByte);

    /* The number of refined elements coincides with the number of entropy codes */
    return num_refined_elems;
    #endif



    /* In order to determine the number of entropy codes that we need to skip for this process has to be computed in parallel, 
     * because we are not able to retrieve informatin on all global trees, their refinement scheme and the current number of leaf elements per tree. */

    /* Move to the correct bit position within this level */
    mesh_lvl_encoding.MoveToStartBit(proc_elem_offset);

    /* Therefore, we start to iterate locally throught all trees/elemenets and count them */
    uint64_t num_local_entropy_codes{0};

    /* Get the scheme of the mesh */
    const t8_scheme_c* scheme =  t8_forest_get_scheme(this->mesh_.GetMesh());

    /* Get the number of local trees */
    const t8_locidx_t num_trees = t8_forest_get_num_local_trees (this->mesh_.GetMesh());
    for (t8_locidx_t itree{0}; itree < num_trees; ++itree)
    {
        /* Get the tree class from the local tree */
        const t8_eclass_t tree_class = t8_forest_get_tree_class (this->mesh_.GetMesh(), itree);

        /* Check whether the tree refines regularly */
        const bool refines_regular = not scheme->refines_irrgular(tree_class);

        /* Get the number of elements in the tree */
        const t8_locidx_t num_elems_in_tree = t8_forest_get_tree_num_leaf_elements (this->mesh_.GetMesh(), itree);

        /* Check if there are local elements in this tree */
        if (num_elems_in_tree > 0)
        {
            if (refines_regular)
            {
                /* If the tree refines regular, the number of children elements is always the same */
                /* Get the first element in the tree */
                const t8_element_t *element = t8_forest_get_leaf_element_in_tree (this->mesh_.GetMesh(), itree, 0);
                /* Get the number of children elements an element refines to */
                const int num_children = scheme->element_get_num_children(tree_class, elements[0]);

                //TODO: In this case, we can directly iterate through the local bytes of this tree within the encoding and do a population count which can be multiplied afterwards by the num children
                //For now, we just iterate through the elems and check whether they will be refined or not 
                for (t8_locidx_t ielem{0}; ielem < num_elems_in_tree; ++ielem)
                {
                    /* Check whether this element will be refined during this iteration */
                    if (mesh_lvl_encoding.GetNextBit())
                    {
                        /* Add this amount of entropy codes to the counter */
                        num_local_entropy_codes += num_children;
                    }
                }
            } else
            {
                /* If the tree refines irregular, we need to check the number of children elements for each element */
                /* We iterate through all elements and count the number of entropy codes that are stored */
                for (t8_locidx_t ielem{0}; ielem < num_elems_in_tree; ++ielem)
                {
                    /* Check whether this element will be refined during this iteration */
                    if (mesh_lvl_encoding.GetNextBit())
                    {
                        /* Get the element in the tree */
                        const t8_element_t *element = t8_forest_get_leaf_element_in_tree (this->mesh_.GetMesh(), itree, ielem);

                        /* Get the number of elements this element refines to */
                        const int elem_num_children = scheme->element_get_num_children(tree_class, element);

                        /* Add this amount of entropy codes to the counter */
                        num_local_entropy_codes += elem_num_children;
                    }
                }
            }
        }   
    }

    /* After we have counted all local entropy codes, we are going to exchange them with an exclusive scan to determine the entropy offset */
    uint64_t entropy_offset{0};

    /* Peform an exclusive scan to obtain the offset for the entropy codes for each process */
    const int rv_entropy_offset_exscan = MPI_Exscan(&num_local_entropy_codes, &entropy_offset, 1, MPI_UINT64_T, MPI_SUM, this->comm_);
    MPICheckError(rv_entropy_offset_exscan);

    /* The value on the root rank may be undefined, therefore, we explicitly overrite it again */
    if (this->comm_rank_ == kRootRank)
    {
        entropy_offset = 0;
    }

    return entropy_offset;
}

void
SetupLevelDecodingStart(const size_t proc_entropy_offset, const std::vector<ProcLevelByteStreamOffsets>& offset_hints)
{

}


template <typename T>
void
AbstractByteParDecompressionVariable<T>::Decompress(const t8_cmesh_t cmesh, const t8_scheme *scheme)
{
    cmc_assert(this->IsValidForDecompression());
    cmc_debug_msg("Decompression of variable ", this->name_, " starts...");

    /* Split the communicator into shared memory communicators */
    this->CreateSharedMemoryComms();

    /* A decompression step couner (and we start on the root level with the decoding) */
    int decompression_step{0};

    /* Re-Create the base mesh of the variable */
    auto [base_mesh, dimensionality] = mesh_decoder_->DecodeRootLevelMesh(cmesh, scheme);
    mesh_.SetMesh(base_mesh);
    mesh_.SetDimensionality(dimensionality);

    /* Open the root level encoding (the partitioning is based on the partitioning of the root mesh) */
    auto [global_var_start_ptr, proc_mesh_elem_offset, proc_mesh_elem_count] = this->OpenSharedLevelDataWindow(decompression_step);

    /* Decode the root level values */
    data_ = this->SetupRootLevelData(global_var_start_ptr, proc_mesh_elem_offset, proc_mesh_elem_count);

    /* Close the root level window */
    this->CloseSharedLevelDataWindow();

    cmc_debug_msg("The root mesh and the corresponding data has been reconstructed.");

    if constexpr (kWriteDecompressionStepToVTK)
    {
        WriteData<T>(mesh_.GetMesh(), data_);
    }

    /* We create the adapt data based on the compression settings, the forest and the variables to consider during the adaptation/coarsening */
    IParDecompressionAdaptData<T>* adapt_data = this->CreateAdaptData();

    /* Perform decompression iterations until the mesh is completely reconstructed or the specified level is reached */
    while (decompression_step < this->GetMaxNumDecompressionIterations())
    {
        ++decompression_step;
        cmc_debug_msg("A decompression iteration (step ", decompression_step, ") is initialized.");

        /****** Setup of a decompression iteration ******/
        /* We create a view on the encoded global refinement indications on this level */
        bit_map::BitMapView mesh_lvl_encoding = this->mesh_decoder_->GetGlobalMeshEncodingStep(decompression_step);

        /* Open a shared window for the encoded level data */
        const auto [global_var_start_ptr, proc_mesh_elem_offset, proc_mesh_elem_count] = this->OpenSharedLevelDataWindow(decompression_step);

        /* Determine the entropy offset of the process in this level */
        const SizeType entropy_code_offset = DetermineEntropyOffsetsFromMeshOffset(mesh_lvl_encoding, proc_mesh_elem_offset, proc_mesh_elem_count);

        /* Allocate for a decompression iteration */
        this->AllocateDecompressionIteration();

        /* Set the correct view for decoding for the local chunk this process processes */
        this->SetupLevelDecodingStart(global_var_start_ptr, entropy_code_offset, levelwise_memory_hints_[decompression_step]); //TODO: implement in specialization

        /* Initialize the next decompression iteratoion of the adaptation data and the mesh */
        adapt_data->InitializeDecompressionIteration();
        mesh_decoder_->InitializeDecompressionIteration();
        /************************************************/

        /****** Perform the refinement onto the next fienr level as dictated by the mesh encoding ******/
        /* Get and indicate to keep the 'previous forest' after the adaptation step */
        t8_forest_t previous_forest = mesh_.GetMesh();
        t8_forest_ref(previous_forest);

        /* Perform a decompression iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, ByteVariableParDecompressionAdaptation<T>, 0, 0, static_cast<void*>(adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements");

        /* Complete the interpolation step by storing the newly computed adapted data alongside its deviations */
        adapt_data->CompleteDecompressionIteration(previous_forest, adapted_forest);

        /* Free the former forest */
        t8_forest_unref(&previous_forest);

        /* Switch to the decompressed data */
        this->SwitchToDecompressedData();

        /************************************************/

        /****** Equally partition the mesh and the data again ******/
        /* Repartition the mesh */
        t8_forest_t partitioned_forest = RepartitionMesh(adapted_forest);

        /* Repartition the data */
        this->RepartitionData(adapted_forest, partitioned_forest);
        adapt_data->RepartitionData(adapted_forest, partitioned_forest);

        cmc_debug_msg("The mesh and the data has been re-partitioned.");

        /* Free the former forest and store the adapted/repartitioned mesh */
        t8_forest_unref(&adapted_forest);
        mesh_.SetMesh(partitioned_forest);

        /************************************************/

        /****** Clean-Up of the decompression iteration ******/
        /* Finalize the decompression iteration */
        adapt_data->FinalizeDecompressionIteration();
        mesh_decoder_->FinalizeDecompressionIteration();

        cmc_debug_msg("The decompression iteration is finished.");
        
        if constexpr (kWriteDecompressionStepToVTK)
        {
            WriteData<T>(mesh_.GetMesh(), data_);
            cmc_debug_msg("The decompression step has been written out in a .vtu file.");
        }

       /* Close the root level window */
        this->CloseSharedLevelDataWindow();
        /************************************************/
    }

    /* Free the adapt data structure */
    this->adaptation_destructor_(adapt_data);
    cmc_debug_msg("Decompression of variable ", this->name_, " is finished.");
}

template <typename T>
CompressionValue<T>
AbstractByteParDecompressionVariable<T>::GetValue(const int tree_id, const int lelement_index) const
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
AbstractByteParDecompressionVariable<T>::GetView(const int start_index, const int count) const
{
    cmc_assert(start_index >= 0 && count >= 0);
    cmc_assert(static_cast<size_t>(start_index + count) <= data_.size());

    return VectorView(&data_[start_index], count);
}

template <typename T>
VectorView<CompressionValue<T>>
AbstractByteParDecompressionVariable<T>::GetView(const int tree_id, const int lelement_index, const int count) const
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
AbstractByteParDecompressionVariable<T>::StoreRefinedValues(const RefinementData<T>& refined_values)
{
    /* Store the refined values */
    std::copy_n(refined_values.fine_values.begin(), refined_values.fine_values.size(), std::back_inserter(data_new_));
}

template <typename T>
void
AbstractByteParDecompressionVariable<T>::StoreUnchangedElement(const UnchangedData<T>& unchanged_value)
{
    /* Store the element value */
    data_new_.push_back(unchanged_value.fine_value);
}


template <typename T>
inline t8_forest_t
AbstractByteParDecompressionVariable<T>::RepartitionMesh(t8_forest_t adapted_forest)
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
AbstractByteParDecompressionVariable<T>::RepartitionData(t8_forest_t adapted_forest, t8_forest_t partitioned_forest)
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
AbstractByteParDecompressionVariable<T>::IsValidForDecompression() const 
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

#endif /* !LOSSLESS_CMC_BYTE_PAR_DECOMPRESSION_VARIABLE_HXX */
