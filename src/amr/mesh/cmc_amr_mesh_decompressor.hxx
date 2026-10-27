#ifndef AMR_MESH_CMC_AMR_MESH_DECOMPRESSOR_HXX
#define AMR_MESH_CMC_AMR_MESH_DECOMPRESSOR_HXX

#include "cmc.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_bits.hxx"
#include "utilities/cmc_bits_vector.hxx"
#include "utilities/cmc_bits_vector_view.hxx"

#include <t8_forest/t8_forest_general.h>

namespace cmc::mesh
{

constexpr bool kWriteOutIntermediateDecompressionSteps = false;

constexpr uint64_t kT8CodeSupportedMaxDimension = 3;
constexpr uint64_t kT8CodeSupportedMaxRefinementLevel = 65;

class CompressionInfoStruct
{
public:
    CompressionInfoStruct(){} //Null constructor
    CompressionInfoStruct(const uint64_t dimensionality_, const uint64_t init_mpi_comm_size_, const uint64_t mesh_compression_levels_,
                          const std::vector<uint64_t>& global_num_elem_indications_on_levels_, const uint64_t num_bytes_global_mesh_encoding_)
    : dimensionality{dimensionality_}, init_mpi_comm_size{init_mpi_comm_size_}, mesh_compression_levels{mesh_compression_levels_}, global_num_elem_indications_on_levels(global_num_elem_indications_on_levels_),
      num_bytes_global_mesh_encoding{num_bytes_global_mesh_encoding_}
    {}

    cmc::bits::vector_view
    GetGlobalElementIndicationsLevelStart(const int level)
    {
        cmc_assert(level >= 0 && level < static_cast<int>(this->mesh_compression_levels) - 1);
        cmc_assert(shared_data_init_ptr != nullptr);

        int offset{0};

        for (int lvl_idx{0}; lvl_idx < level; ++lvl_idx)
        {
            offset += this->global_num_elem_indications_on_levels[lvl_idx] / (sizeof(uint64_t) * cmc::bits::kCharBit) + (this->global_num_elem_indications_on_levels[lvl_idx] % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0);
        }

        return cmc::bits::vector_view(shared_data_init_ptr + offset);
    }

    int GetMeshRefinementSteps()
    {
        cmc_assert(this->mesh_compression_levels >= 1);
        return static_cast<int>(this->mesh_compression_levels) - 1;
    }

private:
    friend class Decompressor;

    uint64_t dimensionality{0};
    uint64_t init_mpi_comm_size{0};
    uint64_t mesh_compression_levels{0};
    std::vector<uint64_t> global_num_elem_indications_on_levels; 
    uint64_t num_bytes_global_mesh_encoding{0};

    uint64_t* shared_data_init_ptr{nullptr};
    MPI_Win shared_data;
};

constexpr int kMaxDecompressionLevelUndefined = std::numeric_limits<int>::max();

/* Definition of the mesh decompressor class */
class Decompressor
{
public:
    Decompressor() = delete;
    Decompressor(const std::string& file_name, const MPI_Comm comm, const t8_cmesh_t cmesh, const t8_scheme *scheme);

    [[nodiscard]] t8_forest_t ReconstructMesh(const int max_decompression_refinement_level = kMaxDecompressionLevelUndefined);

private:
    void InquireCompressionInfo();
    t8_forest_t RepartitionMesh(t8_forest_t adapted_forest);

    /* File containing the compressed data */
    const std::string file_name_;

    /* The handle to the compressed file */
    MPI_File fhandle_;

    /* Infos extarcted from the preamble */
    CompressionInfoStruct compression_info_;

    /* The current AMR mesh */
    AmrMesh mesh_;

    /* MPI communciators */
    const MPI_Comm comm_{MPI_COMM_NULL};
    int comm_rank_{0}, comm_size_{1};
    MPI_Comm shm_comm_{MPI_COMM_NULL};
    int shm_rank_{0}, shm_size_{1};

    /* Shared memory window */
    MPI_Win lvl_window_;

    /* A flag whether the decompression has been carried out and the file been closed */
    bool is_already_decompressed_{false};
};

inline void 
CheckMPIReadCorrectness(const MPI_Status* status, const MPI_Datatype datatype, const int expected_num_elems)
{
    int elem_count_{0};
    const int rv_check_read = MPI_Get_count(status, datatype, &elem_count_);
    MPICheckError(rv_check_read);
    if (elem_count_ != expected_num_elems)
    {
        cmc_global_msg("The expeceted number of bytes could not be read from the file.");
        MPICheckError(MPI_ERR_COUNT);
    }
}

Decompressor::Decompressor(const std::string& file_name, const MPI_Comm comm, const t8_cmesh_t cmesh, const t8_scheme *scheme)
: file_name_(file_name), comm_{comm}
{
    /* Check if the compressed file exists */
    if (const std::filesystem::path input_file_path(this->file_name_); not std::filesystem::exists(input_file_path))
    {
        throw std::invalid_argument("The compressed file does not exist!");
    }

    /* Check if an MPI Communicator is given */
    if (this->comm_ == MPI_COMM_NULL)
    {
        throw std::invalid_argument("The MPI_Communicator is NULL!");
    }

    /* Gather the MPI rank and size of comm */
    const int rv_rank_ = MPI_Comm_rank(this->comm_, &(this->comm_rank_));
    MPICheckError(rv_rank_);
    const int rv_size_ = MPI_Comm_size(this->comm_, &(this->comm_size_));
    MPICheckError(rv_size_);

    /* Split communicator into shared memory groups */
    const int rv_split_comm = MPI_Comm_split_type(this->comm_, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &(this->shm_comm_));
    MPICheckError(rv_split_comm);

    /* Get the rank within and the size of the shared communicator */
    const int rv_shm_rank_ = MPI_Comm_rank(this->shm_comm_, &(this->shm_rank_));
    MPICheckError(rv_shm_rank_);
    const int rv_shm_size_ = MPI_Comm_size(this->shm_comm_, &(this->shm_size_));
    MPICheckError(rv_shm_size_);

    /* Open the file on the shared memory communicator */
    const int rv_open = MPI_File_open(this->shm_comm_, this->file_name_.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &(this->fhandle_));
    MPICheckError(rv_open);

    /* Set the native data representation */
    const int rv_file_view = MPI_File_set_view(this->fhandle_, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    MPICheckError(rv_file_view);
    
    /* Inquire the compression information */
    this->InquireCompressionInfo();

    /* Create a forest mesh from the given parameters */
    mesh_.SetMesh(t8_forest_new_uniform (cmesh, scheme, 0, 0, this->comm_));
}

void
Decompressor::InquireCompressionInfo()
{
    /* Read the first three uint64_ts from the file */
    std::array<uint64_t, 3> info_vals{0,0,0};

    MPI_Status status;

    /* We read the beginning from the header only on the root rank */
    /* Move to the start of the file */
    const int rv_seek_start = MPI_File_seek(this->fhandle_, 0, MPI_SEEK_SET);
    MPICheckError(rv_seek_start);

    /* Read the first three uint64_ts */
    const int rv_start_bytes = MPI_File_read(this->fhandle_, info_vals.data(), 3, MPI_UINT64_T, &status);
    MPICheckError(rv_start_bytes);
    CheckMPIReadCorrectness(&status, MPI_UINT64_T, 3);
    
    /* Get the mesh dimension */
    const uint64_t mesh_dimension = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(info_vals[0]);

    /* Check for legal input */
    if (mesh_dimension > kT8CodeSupportedMaxDimension) [[unlikely]]
    {
        cmc_err_msg("The dimension of the mesh (DIM: ", mesh_dimension, ") is larger than the maximum supported dimension (MAX_DIM: ", kT8CodeSupportedMaxDimension, ")!");
        throw std::invalid_argument("The mesh dimension is not supported!");
    }

    /* Get the size of the compression MPI Communicator */
    const uint64_t init_mpi_comm_size = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(info_vals[1]);
    
    /* Get the number of mesh compression levels */
    const uint64_t mesh_compression_levels = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(info_vals[2]);

    /* Check for legal input */
    if (mesh_compression_levels > kT8CodeSupportedMaxRefinementLevel) [[unlikely]]
    {
        cmc_err_msg("The compression levels of the mesh (NUM_LVLS: ", mesh_compression_levels, ") extend the maximum supported levels (MAX_NUM_LVLS: ", kT8CodeSupportedMaxRefinementLevel, ")!");
        throw std::invalid_argument("The number of compression levels is too large!");
    }

    
    /* Read in the remaining global number of elements per level */
    const int num_data_to_read = static_cast<int>(mesh_compression_levels + 1);
    std::vector<uint64_t> remaining_vals(num_data_to_read);
    const int rv_remaining_bytes = MPI_File_read(this->fhandle_, remaining_vals.data(), num_data_to_read, MPI_UINT64_T, &status);
    MPICheckError(rv_remaining_bytes);
    CheckMPIReadCorrectness(&status, MPI_UINT64_T, num_data_to_read);

    /* Build up the vector containing the global number of elements on each level (from root to leaf level) */
    std::vector<uint64_t> num_global_elems_per_level;
    num_global_elems_per_level.reserve(mesh_compression_levels);

    uint64_t previous_num_elems{0};

    for (uint64_t lvl_idx{0}; lvl_idx < mesh_compression_levels; ++lvl_idx)
    {
        num_global_elems_per_level.push_back(cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(remaining_vals[lvl_idx]));
        
        if (previous_num_elems > num_global_elems_per_level.back()) [[unlikely]]
        {
            cmc_err_msg("The number of global elements on successive level cannot shrink (PREV_LVL: ", previous_num_elems, "; CURR_LVL: ", num_global_elems_per_level.back(), ")!");
            throw std::invalid_argument("An incorrect table of global elements have been received!"); 
        }

        previous_num_elems = num_global_elems_per_level.back();
    }

    /* Get the length of the global mesh encoding (the last value from the preamble) */
    const uint64_t global_mesh_encoding_bytes = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(remaining_vals.back());

    if (global_mesh_encoding_bytes % sizeof(uint64_t) != 0) [[unlikely]]
    {
        cmc_err_msg("The mesh encoding is of length of a multiple of eight bytes!");
        throw std::invalid_argument("The mesh encoding is of length of a multiple of eight bytes!"); 
    }

    /* Offset to the encoded mesh in the file */
    const uint64_t mesh_encoding_byte_offset = (3 + num_data_to_read) * sizeof(uint64_t);

    if (mesh_encoding_byte_offset % sizeof(uint64_t) != 0) [[unlikely]]
    {
        cmc_err_msg("The mesh encoding needs to start at a multiple of eight bytes from the file start!");
        throw std::invalid_argument("The mesh encoding needs to start at a multiple of eight bytes from the file start!"); 
    }

    /* Create a compression info struct from the obtained information */
    this->compression_info_ = CompressionInfoStruct(mesh_dimension, init_mpi_comm_size, mesh_compression_levels,
                                                    num_global_elems_per_level, global_mesh_encoding_bytes);

    /* We allocate the the shared window from the root rank of the shared memory communicator */
    const MPI_Aint stream_length_in_bytes = (this->shm_rank_ != kRootRank ? 0 : global_mesh_encoding_bytes);
    const MPI_Aint stream_length = stream_length_in_bytes / sizeof(uint64_t);

    /* Allocate a shared window for the compression preamble */
    /* Since the data is relatively small, we allocate it from the root rank only and read it from this process only */
    const int rv_win_alloc = MPI_Win_allocate_shared(stream_length_in_bytes, sizeof(uint64_t), MPI_INFO_NULL, this->shm_comm_, &(this->compression_info_.shared_data_init_ptr), &(this->compression_info_.shared_data));
    MPICheckError(rv_win_alloc);

    if (this->shm_rank_ == kRootRank)
    {
        /* Move to the start of the global mesh encoding in the file */
        const int rv_seek_start_again = MPI_File_seek(this->fhandle_, mesh_encoding_byte_offset, MPI_SEEK_SET);
        MPICheckError(rv_seek_start_again);

        /* Read the complete global mesh encoding from the root rank of the shared communicator */
        const int rv_read_data = MPI_File_read(this->fhandle_, this->compression_info_.shared_data_init_ptr, stream_length, MPI_UINT64_T, &status);
        MPICheckError(rv_read_data);
        CheckMPIReadCorrectness(&status, MPI_UINT64_T, stream_length);
    }

    /* Next, we synchronize the window */
    const int rv_sync_win = MPI_Win_sync(this->compression_info_.shared_data);
    MPICheckError(rv_sync_win);

    /* Impose an additional barrier on the shared communicator */
    const int rv_shm_barrier = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barrier);

    /* Setup the decoding start correctly on all processes */
    this->compression_info_.shared_data_init_ptr -= (this->shm_rank_ != kRootRank ? (global_mesh_encoding_bytes / sizeof(uint64_t))  : 0);
}

t8_forest_t
Decompressor::RepartitionMesh(t8_forest_t adapted_forest)
{
    /* Allocate a forest */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    constexpr int partition_for_coarsening = 0;
    t8_forest_set_partition(partitioned_forest, adapted_forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    return partitioned_forest;
}

struct RefinementIterationData
{
    RefinementIterationData(cmc::bits::vector_view process_local_level_mesh_encoding_)
    : process_local_level_mesh_encoding{process_local_level_mesh_encoding_}
    {}

    bool WillNextElementBeRefined() {return process_local_level_mesh_encoding.GetNextBit();}

private:
    cmc::bits::vector_view process_local_level_mesh_encoding;
};

inline t8_locidx_t
AMRMeshDecompression (t8_forest_t forest,
                      [[maybe_unused]] t8_forest_t forest_from,
                      [[maybe_unused]] t8_locidx_t which_tree,
                      [[maybe_unused]] const t8_eclass_t tree_class,
                      [[maybe_unused]] t8_locidx_t lelement_id,
                      [[maybe_unused]] const t8_scheme_c * ts,
                      [[maybe_unused]] const int is_family,
                      [[maybe_unused]] const int num_elements,
                      [[maybe_unused]] t8_element_t* elements[])
{
    /* Retrieve the adapt_data */
    RefinementIterationData* adapt_data = static_cast<RefinementIterationData*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Check whether we perform a refinement or not */
    if (adapt_data->WillNextElementBeRefined())
    {
        /* If the element will be refined */
        return cmc::t8::kRefineElement;
    } else
    {
        return cmc::t8::kLeaveElementUnchanged;
    }
}

inline t8_forest_t
Decompressor::ReconstructMesh(const int max_decompression_refinement_level)
{
    /* Check if the data has already been decompressed */
    if (this->is_already_decompressed_) [[unlikely]]
    {
        cmc_err_msg("The mesh has already been decompressed! A new decompessor needs to be instantiated!");
    }

    /* Check for legal input */
    if (max_decompression_refinement_level < 0) [[unlikely]]
    {
        cmc_err_msg("The maximum decompression level cannot be negative!");
    }

    // Potentially, write out the root level mesh
    if constexpr (kWriteOutIntermediateDecompressionSteps)
    {
        t8_forest_write_vtk_ext (this->mesh_.GetMesh(), "cmc_mesh_decompression_step_0", 1, 1, 1, 1, 0, 0, 0, 0, NULL);
    }

    /* Determine the number of compression steps */
    const int num_mesh_compression_levels = std::min(this->compression_info_.GetMeshRefinementSteps(), max_decompression_refinement_level);

    /* Reconstruct the mesh level-wise */
    for (int mesh_lvl{0}; mesh_lvl < num_mesh_compression_levels; ++mesh_lvl)
    {
        cmc_debug_msg("The mesh decompression step ", mesh_lvl," starts...");

        /****** Get the correct view on this level for each process ******/

        /* Get the view onto the start of this level's global encosing */
        cmc::bits::vector_view level_mesh_encoding = this->compression_info_.GetGlobalElementIndicationsLevelStart(mesh_lvl);

        /* Get the global offset of the first element */
        const uint64_t mesh_offset = (uint64_t) t8_forest_get_first_local_leaf_element_id(mesh_.GetMesh());

        /* Set the process local mesh encoding start */
        level_mesh_encoding.MoveToOffsetBitInStream(mesh_offset);
    
        /****** Perform the local refinement onto the next finer level as dictated by the mesh encoding ******/
        /* Create the adaptation data */
        RefinementIterationData adapt_data(level_mesh_encoding);

        /* Perform a decompression iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(mesh_.GetMesh(), AMRMeshDecompression, 0, 0, static_cast<void*>(&adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements");

        /****** Partition the mesh ******/
        /* Repartition the mesh */
        t8_forest_t partitioned_forest = this->RepartitionMesh(adapted_forest);

        /* Store the partitioned mesh */
        mesh_.SetMesh(partitioned_forest);

        // Potentially, write out the intermediate refined mesh
        if constexpr (kWriteOutIntermediateDecompressionSteps)
        {
            const std::string vtk_prefix = "cmc_mesh_decompression_step_" + std::to_string(mesh_lvl + 1);
            t8_forest_write_vtk_ext (this->mesh_.GetMesh(), vtk_prefix.c_str(), 1, 1, 1, 1, 0, 0, 0, 0, NULL);
        }

        cmc_debug_msg("The mesh decompression step ", mesh_lvl," has been completed.");
    }

    /* Store the decompressed mesh temporarily */
    t8_forest_t partitioned_decompressed_mesh = this->mesh_.GetMesh();

    /* Impose a barrier to finish all outstanding reads */
    const int rv_shm_barrier = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barrier);

    /* Free the allocated window */
    const int rv_win_free = MPI_Win_free(&(this->compression_info_.shared_data));
    MPICheckError(rv_win_free);

    /* Close the compressed file */
    const int rv_close = MPI_File_close(&this->fhandle_);
    MPICheckError(rv_close);

    /* Set the flag that the data has already been decompressed */
    this->is_already_decompressed_ = true;

    /* Release the pointer to the mesh */
    this->mesh_.SetNullMesh();

    cmc_debug_msg("Decompression of the mesh onto level ", num_mesh_compression_levels, " has been completed.");

    /* Return the decompressed partitioned mesh */
    return partitioned_decompressed_mesh;
}

inline
void
PrintMeshCompressionInfo(const std::string& file_name)
{
    /* The handle to the compressed file */
    MPI_File fhandle;

    /* Open the file on each processs */
    const int rv_open = MPI_File_open(MPI_COMM_SELF, file_name.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fhandle);
    MPICheckError(rv_open);

    /* Set the native data representation */
    const int rv_file_view = MPI_File_set_view(fhandle, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
    MPICheckError(rv_file_view);

    /* Read the first three uint64_ts from the file */
    std::array<uint64_t, 3> info_vals{};

    MPI_Status status;

    /* We read the beginning from the header only on the root rank */
    /* Move to the start of the file */
    const int rv_seek_start = MPI_File_seek(fhandle, 0, MPI_SEEK_SET);
    MPICheckError(rv_seek_start);

    /* Read the first three uint64_ts */
    const int rv_start_bytes = MPI_File_read(fhandle, info_vals.data(), 3, MPI_UINT64_T, &status);
    MPICheckError(rv_start_bytes);
    CheckMPIReadCorrectness(&status, MPI_UINT64_T, 3);
    
    /* Get the mesh dimension */
    const uint64_t mesh_dimension = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(info_vals[0]);

    /* Check for legal input */
    if (mesh_dimension > kT8CodeSupportedMaxDimension) [[unlikely]]
    {
        cmc_err_msg("The dimension of the mesh (DIM: ", mesh_dimension, ") is larger than the maximum supported dimension (MAX_DIM: ", kT8CodeSupportedMaxDimension, ")!");
        throw std::invalid_argument("The mesh dimension is not supported!");
    }

    /* Get the size of the compression MPI Communicator */
    const uint64_t init_mpi_comm_size = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(info_vals[1]);
    
    /* Get the number of mesh compression levels */
    const uint64_t mesh_compression_levels = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(info_vals[2]);

    /* Check for legal input */
    if (mesh_compression_levels > kT8CodeSupportedMaxRefinementLevel) [[unlikely]]
    {
        cmc_err_msg("The compression levels of the mesh (NUM_LVLS: ", mesh_compression_levels, ") extend the maximum supported levels (MAX_NUM_LVLS: ", kT8CodeSupportedMaxRefinementLevel, ")!");
        throw std::invalid_argument("The number of compression levels is too large!");
    }

    /* Read in the remaining global number of elements per level */
    const int num_data_to_read = static_cast<int>(mesh_compression_levels + 1);
    std::vector<uint64_t> remaining_vals(num_data_to_read);
    const int rv_remaining_bytes = MPI_File_read(fhandle, remaining_vals.data(), num_data_to_read, MPI_UINT64_T, &status);
    MPICheckError(rv_remaining_bytes);
    CheckMPIReadCorrectness(&status, MPI_UINT64_T, num_data_to_read);

    cmc_global_msg("mesh_compression_levels: ", mesh_compression_levels);

    /* Build up the vector containing the global number of elements on each level (from root to leaf level) */
    std::vector<uint64_t> num_global_elems_per_level;
    num_global_elems_per_level.reserve(mesh_compression_levels);

    uint64_t previous_num_elems{0};

    for (uint64_t lvl_idx{0}; lvl_idx < mesh_compression_levels; ++lvl_idx)
    {
        num_global_elems_per_level.push_back(cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(remaining_vals[lvl_idx]));
        
        cmc_global_msg("Num global elems per level: ", num_global_elems_per_level.back());

        if (previous_num_elems > num_global_elems_per_level.back()) [[unlikely]]
        {
            cmc_err_msg("The number of global elements on successive level cannot shrink (PREV_LVL: ", previous_num_elems, "; CURR_LVL: ", num_global_elems_per_level.back(), ")!");
            throw std::invalid_argument("An incorrect table of global elements have been received!"); 
        }

        previous_num_elems = num_global_elems_per_level.back();
    }

    /* Get the length of the global mesh encoding (the last value from the preamble) */
    const uint64_t global_mesh_encoding_bytes = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(remaining_vals.back());

    if (global_mesh_encoding_bytes % sizeof(uint64_t) != 0) [[unlikely]]
    {
        cmc_err_msg("The mesh encoding is of length of a multiple of eight bytes!");
        throw std::invalid_argument("The mesh encoding is of length of a multiple of eight bytes!"); 
    }

    /* Offset to the encoded mesh in the file */
    const uint64_t mesh_encoding_byte_offset = (3 + num_data_to_read) * sizeof(uint64_t);

    if (mesh_encoding_byte_offset % sizeof(uint64_t) != 0) [[unlikely]]
    {
        cmc_err_msg("The mesh encoding needs to start at a multiple of eight bytes from the file start!");
        throw std::invalid_argument("The mesh encoding needs to start at a multiple of eight bytes from the file start!"); 
    }

    /* Print the gathered information */
    cmc_msg("Compression Information retrieved from: ", file_name);
    cmc_msg("\t Global Overall Mesh Encoding Size: ", global_mesh_encoding_bytes, " bytes");
    cmc_msg("\t Mesh Dimensionality: ", mesh_dimension, "D");
    cmc_msg("\t Compression MPI Communicator Size: ", init_mpi_comm_size);
    cmc_msg("\t Mesh Compression Levels: ", mesh_compression_levels);
    cmc_msg("\t Mesh Compression Level - Global ELement Count: ");
    for (uint64_t lvl_idx{0}; lvl_idx < mesh_compression_levels; ++lvl_idx)
    {
        cmc_msg("\t\t Level ", lvl_idx, ": ", num_global_elems_per_level[lvl_idx], " elements");
    }

    /* Close the compressed file */
    const int rv_close = MPI_File_close(&fhandle);
    MPICheckError(rv_close);
}

}

#endif /* !AMR_MESH_CMC_AMR_MESH_DECOMPRESSOR_HXX */
