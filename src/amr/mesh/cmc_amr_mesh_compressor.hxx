#ifndef AMR_MESH_CMC_AMR_MESH_COMPRESSOR_HXX
#define AMR_MESH_CMC_AMR_MESH_COMPRESSOR_HXX

#include "cmc.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_log_functions.hxx"
#include <t8_forest/t8_forest_partition.h>
#include "utilities/cmc_bits.hxx"
#include "utilities/cmc_bits_vector.hxx"

#include <string>
#include <cmath>
#include <filesystem>

namespace cmc::mesh
{

constexpr bool kWriteOutIntermediateCompressionSteps = false;

constexpr int kRootRank = 0;
constexpr int kTagMeshEncoding = 1000;

inline constexpr int32_t kMaxPresentElementLevelUnknown = -1;
inline constexpr int32_t kMeshDimensionUnknown = -1;

class Compressor
{
public:
    Compressor() = delete;
    Compressor(t8_forest_t forest)
    : mesh_(forest) {
        /* Get the MPI communicator from the mesh */
        this->comm_ = t8_forest_get_mpicomm(forest);
        /* Gather the initial number of global leaf elements */
        this->num_init_global_elements_ = t8_forest_get_global_num_leaf_elements(forest);
        /* Get the rank and the size of the communicator */
        const int rv_size = sc_MPI_Comm_size(this->comm_, &(this->comm_size_));
        MPICheckError(rv_size);
        const int rv_rank = sc_MPI_Comm_rank(this->comm_, &(this->comm_rank_));
        MPICheckError(rv_rank);
        /* Increment the reference count, since we do not have ownership of the forest mesh */
        t8_forest_ref(forest);
        /* Get the cmesh */
        t8_cmesh_t cmesh = t8_forest_get_cmesh (this->mesh_.GetMesh());
        /* Check whether the dimensionality matches */
        this->mesh_dimension_ = t8_cmesh_get_dimension (cmesh);
    }

    void CompressMeshStructure(const bool enfore_strictly_level_wise_compression = false);
    void WriteCompressedMeshStructure(const std::string& file_name, const bool block_mpi_processes_until_file_is_written = true);
    
    
    void SetMaximumInitialElementLevel(const int max_init_elem_level);
    void SetIsMeshAlreadyPartitionedForCoarsening(const bool is_partitioned_for_coarsening);

private:
    void DetermineMaxInitElementLevel();
    bool IsCompressionProgressing() const;
    void RepartitionMeshForCoarsening(t8_forest_t adapted_mesh);
    bool IsAlreadyPartitionedForCoarsening() const;
    void GenerateOuputStreams();
    void AppendVariableHeaderToStream(std::vector<uint64_t>& stream, const int mesh_compression_levels);
    std::vector<std::vector<uint64_t>> AdjustLevelMeshEncodings(const std::vector<uint64_t>& partition_info, const int num_global_encoding_steps);
    void CollectMeshEncodingOnTheRootRank(const std::vector<uint64_t>& global_partition_info, const uint64_t num_global_encoding_steps, const std::vector<std::vector<uint64_t>>& offseted_level_mesh_encodings);
    void InsertFileHeader(std::vector<uint64_t>& stream, const int mesh_compression_levels);

    /* The input mesh for the compression */
    AmrMesh mesh_;

    uint64_t num_init_global_elements_{0};
    int32_t max_init_elem_level_{kMaxPresentElementLevelUnknown};
    int32_t mesh_dimension_{kMeshDimensionUnknown};
    bool is_initially_partitioned_for_coarsening_{false};

    /* The MPI communicator to use (which is extracted from the mesh) */
    MPI_Comm comm_{MPI_COMM_NULL};
    int comm_size_{1}, comm_rank_{0};

    /* Members filled by the mesh coarsening compression */
    std::vector<cmc::bits::vector> coarsening_indications_;
    std::vector<uint64_t> level_partitioning_info_;
    
    /* Mesh Encoding after the compression */
    std::vector<uint64_t> num_elements_per_level_;
    std::vector<uint64_t> global_mesh_encoding_; //Only filled by the root rank

    /* Data to be written out generated after compression */
    uint64_t global_compressed_file_size_in_bytes_{0};
    std::vector<uint64_t> root_rank_compression_stream_;
};

inline void
Compressor::SetMaximumInitialElementLevel(const int max_init_elem_level)
{
    max_init_elem_level_ = max_init_elem_level;
}

inline void
Compressor::SetIsMeshAlreadyPartitionedForCoarsening(const bool is_partitioned_for_coarsening)
{
    is_initially_partitioned_for_coarsening_ = is_partitioned_for_coarsening;
}

inline bool
Compressor::IsAlreadyPartitionedForCoarsening() const
{
    return is_initially_partitioned_for_coarsening_;
}

inline bool
Compressor::IsCompressionProgressing() const
{
    return (mesh_.GetNumberGlobalElements() > mesh_.GetNumberGlobalTrees());
}

struct CoarseningIterationData
{
    CoarseningIterationData() = delete;
    CoarseningIterationData(const int32_t init_max_level, const int32_t coarsening_step, const uint64_t num_global_elems, const int32_t mesh_dimensionality)
    : current_coarsening_level{init_max_level - coarsening_step}
    {
        coarsening_indications.Reserve(num_global_elems / std::pow<int32_t>(int32_t{2}, mesh_dimensionality) + 1);
    }

    void LeaveElementUnchanged()
    {
        /* Indicate that no coarsening has been performed */
        this->coarsening_indications.AppendUnsetBit();
    }

    void PerformExtraction()
    {
        /* Indicate that corsening has been performed */
        this->coarsening_indications.AppendSetBit();
    }

    const int64_t current_coarsening_level;
    cmc::bits::vector coarsening_indications;
};

inline void
Compressor::RepartitionMeshForCoarsening(t8_forest_t adapted_mesh)
{
    /** Partition the mesh **/

    /* Allocate a forest */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    constexpr int partition_for_coarsening = 1;
    t8_forest_set_partition(partitioned_forest, adapted_mesh, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    /* Free the former forest and store the adapted/repartitioned mesh */
    mesh_.SetMesh(partitioned_forest);
}

inline void
Compressor::DetermineMaxInitElementLevel()
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

constexpr bool
CheckIfElementIsEligibleForCoarsening(const int current_coarsening_level, const int element_level)
{
    return (current_coarsening_level == element_level);
}

inline t8_locidx_t
CoarsenMeshLevelWise (t8_forest_t forest,
                      [[maybe_unused]] t8_forest_t forest_from,
                      [[maybe_unused]] t8_locidx_t which_tree,
                      const t8_eclass_t tree_class,
                      [[maybe_unused]] t8_locidx_t lelement_id,
                      const t8_scheme_c * ts,
                      const int is_family,
                      [[maybe_unused]] const int num_elements,
                      t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    CoarseningIterationData* adapt_data = static_cast<CoarseningIterationData*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Check if a family is supplied to the adaptation function */
    if (is_family == 0 || (not CheckIfElementIsEligibleForCoarsening(adapt_data->current_coarsening_level, ts->element_get_level(tree_class, elements[0]))))
    {
        /* If there is no family, the element stays unchanged */
        adapt_data->LeaveElementUnchanged();
        return cmc::t8::kLeaveElementUnchanged;
    } else
    {
        /* Extract a value of the family and coarsen it */
        adapt_data->PerformExtraction();
        return cmc::t8::kCoarsenElements;
    }
}

inline t8_locidx_t
CoarsenMesh (t8_forest_t forest,
             [[maybe_unused]] t8_forest_t forest_from,
             [[maybe_unused]] t8_locidx_t which_tree,
             [[maybe_unused]] const t8_eclass_t tree_class,
             [[maybe_unused]] t8_locidx_t lelement_id,
             [[maybe_unused]] const t8_scheme_c* ts,
             const int is_family,
             [[maybe_unused]] const int num_elements,
             [[maybe_unused]] t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    CoarseningIterationData* adapt_data = static_cast<CoarseningIterationData*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Check if a family is supplied to the adaptation function */
    if (is_family == 0)
    {
        /* If there is no family, the element stays unchanged */
        adapt_data->LeaveElementUnchanged();
        return cmc::t8::kLeaveElementUnchanged;
    } else
    {
        /* Extract a value of the family and coarsen it */
        adapt_data->PerformExtraction();
        return cmc::t8::kCoarsenElements;
    }
}

inline void
Compressor::CompressMeshStructure(const bool enfore_strictly_level_wise_compression)
{
    cmc_global_msg("The global AMR mesh compression begins.");
    if (enfore_strictly_level_wise_compression) [[unlikely]]
    {
        cmc_debug_msg("The mesh compression is performed strictly level-wise.");
        /* Potentially, gather the maximum present element level if has not been set already */
        this->DetermineMaxInitElementLevel();
        /* Allocate for the expected levelwise data-streams */
        this->coarsening_indications_.reserve(this->max_init_elem_level_ + 1);
    } else
    {
        cmc_debug_msg("The mesh compression is performed as fast as possible disregarding a strictly level-wise compression.");
    }

    // Potentially, write out the leaf level mesh
    if constexpr (kWriteOutIntermediateCompressionSteps)
    {
        t8_forest_write_vtk_ext (this->mesh_.GetMesh(), "cmc_mesh_compression_step_0", 1, 1, 1, 1, 0, 0, 0, 0, NULL);
    }

    /* Potentially, perform partition for coarsening */
    if (not this->IsAlreadyPartitionedForCoarsening())
    {
        this->RepartitionMeshForCoarsening(this->mesh_.GetMesh());
    }

    /* We also store a partition info for the leaf level */
    this->level_partitioning_info_.push_back(static_cast<uint64_t>(t8_forest_get_local_num_leaf_elements(mesh_.GetMesh())));

    int32_t compression_step{0};

    // Compress until the root level
    while(this->IsCompressionProgressing())
    {
        cmc_debug_msg("A coarsening iteration is initialized.");
        
        /* Get the number of global elements */
        const uint64_t num_global_elems = static_cast<uint64_t>(t8_forest_get_global_num_leaf_elements(mesh_.GetMesh()));

        // 1) Allocate an extraction iteration 
        CoarseningIterationData adapt_data(this->max_init_elem_level_, compression_step, num_global_elems, this->mesh_dimension_);

        // 2)  Adapt
        /* Get the correct adaptation function pointer */
        cmc::t8::AdaptationFn adapt_fn = (enfore_strictly_level_wise_compression ? CoarsenMeshLevelWise : CoarsenMesh);

        /* Perform a coarsening iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(mesh_.GetMesh(), adapt_fn, 0, 0, static_cast<void*>(&adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements");

        // 3) Store data to be encoded later
        this->coarsening_indications_.push_back(std::move(adapt_data.coarsening_indications));

        // We need to store the number of local elements after the adaptation for the partition info (before repartioning)
        this->level_partitioning_info_.push_back(static_cast<uint64_t>(t8_forest_get_local_num_leaf_elements(adapted_forest)));

        // 4) Repartition the mesh
        this->RepartitionMeshForCoarsening(adapted_forest);

        // Potentially, write out the intermediate caorsened mesh
        if constexpr (kWriteOutIntermediateCompressionSteps)
        {
            const std::string vtk_prefix = "cmc_mesh_compression_step_" + std::to_string(compression_step + 1);
            t8_forest_write_vtk_ext (this->mesh_.GetMesh(), vtk_prefix.c_str(), 1, 1, 1, 1, 0, 0, 0, 0, NULL);
        }

        cmc_debug_msg("This coarsening iteration is finished.");
        ++compression_step;
    }

    cmc_global_msg("The global AMR mesh has been compressed.");
}

inline uint64_t
ComputeGlobalMeshEncodingNumBytes(const std::vector<uint64_t>& global_elems_per_level)
{
    uint64_t num_vals{0};
    for (auto lvl_iter = global_elems_per_level.begin(); lvl_iter != global_elems_per_level.end(); ++lvl_iter)
    {
        num_vals += (*lvl_iter) / (sizeof(uint64_t) * cmc::bits::kCharBit) + ((*lvl_iter) % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0);
    }
    
    return num_vals * sizeof(uint64_t);
}

inline uint64_t
GetVarHeaderSize(const int mesh_compression_levels)
{
    cmc_assert(mesh_compression_levels > 0);
    return (sizeof(uint64_t)   //Dimensionality
            + sizeof(uint64_t) //Size of MPI communicator
            + sizeof(uint64_t) //Num Mesh Compression Levels
            + ((mesh_compression_levels - 1) * sizeof(uint64_t)) //Global Element Indications per level
            + sizeof(uint64_t) //Global number of initial leaf elements
            + sizeof(uint64_t) //Num Bytes global mesh encoding
           );
}

[[maybe_unused]] inline void
Compressor::AppendVariableHeaderToStream(std::vector<uint64_t>& stream, const int mesh_compression_levels)
{
    /* Store the dimensionality */
    stream.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(this->mesh_dimension_)));

    /* Store the size of the MPI Communicator */
    stream.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(this->comm_size_)));

    /* Store the number of mesh compresion levels */
    stream.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(mesh_compression_levels)));

    /* Append the global num elements per level */
    for (auto lvl_elem_iter = this->num_elements_per_level_.begin(); lvl_elem_iter != this->num_elements_per_level_.end(); ++lvl_elem_iter)
    {
        stream.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(*lvl_elem_iter)));
    }

    /* Push back the number of initial global mesh elements */
    stream.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(this->num_init_global_elements_));

    /* Store the length of the global mesh encoding */
    cmc_assert(not this->global_mesh_encoding_.empty());
    stream.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(this->global_mesh_encoding_.size() * sizeof(uint64_t))));
}

inline void
Compressor::InsertFileHeader(std::vector<uint64_t>& stream, const int mesh_compression_levels)
{
    cmc_assert(stream.size() >= GetVarHeaderSize(mesh_compression_levels) / sizeof(uint64_t));

    int offset{0};

    /* Store the dimensionality */
    stream[offset] = cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(this->mesh_dimension_));
    ++offset;

    /* Store the size of the MPI Communicator */
    stream[offset] = cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(this->comm_size_));
    ++offset;

    /* Store the number of mesh compresion levels */
    stream[offset] = cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(mesh_compression_levels));
    ++offset;

    /* Append the global num elements per level */
    for (auto lvl_elem_iter = this->num_elements_per_level_.begin(); lvl_elem_iter != this->num_elements_per_level_.end(); ++lvl_elem_iter)
    {
        stream[offset] = cmc::bits::ConvertToBigEndian<uint64_t>(static_cast<uint64_t>(*lvl_elem_iter));
        ++offset;
    }

    /* Push back the number of initial global mesh elements */
    stream[offset] = cmc::bits::ConvertToBigEndian<uint64_t>(this->num_init_global_elements_);
    ++offset;

    /* Store the length of the global mesh encoding */
    const uint64_t num_global_mesh_encoding_bytes = ComputeGlobalMeshEncodingNumBytes(this->num_elements_per_level_);
    stream[offset] = cmc::bits::ConvertToBigEndian<uint64_t>(num_global_mesh_encoding_bytes);
}

inline std::vector<std::vector<uint64_t>>
Compressor::AdjustLevelMeshEncodings(const std::vector<uint64_t>& partition_info, const int num_global_encoding_steps)
{
    const int num_coarsening_encodings = static_cast<int>(this->coarsening_indications_.size());
    cmc_assert(num_global_encoding_steps - 1 == num_coarsening_encodings);

    std::vector<std::vector<uint64_t>> offseted_level_mesh_encodings;
    offseted_level_mesh_encodings.reserve(num_coarsening_encodings);

    for (int lvl_idx{0}; lvl_idx < num_coarsening_encodings; ++lvl_idx)
    {
        uint64_t lvl_elem_offset{0};

        /* Count this levels offset */
        for (int rank_idx{0}; rank_idx < this->comm_rank_; ++rank_idx)
        {
            const int access_idx = rank_idx * num_global_encoding_steps + lvl_idx;
            lvl_elem_offset += partition_info[access_idx];
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

inline uint64_t
ComputeGlobalMeshEncodingLength(const std::vector<uint64_t> global_elems_per_level)
{
    uint64_t num_bytes{0};
    for (auto lvl_iter = global_elems_per_level.begin(); lvl_iter != global_elems_per_level.end(); ++lvl_iter)
    {
        num_bytes += (*lvl_iter) / (sizeof(uint64_t) * cmc::bits::kCharBit) + ((*lvl_iter) % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0);
    }

    cmc_debug_msg("Computed number of bytes for the global mesh encoding: ", num_bytes);
    
    return num_bytes;
}

inline void
Compressor::CollectMeshEncodingOnTheRootRank(const std::vector<uint64_t>& global_partition_info, const uint64_t num_global_encoding_steps,  const std::vector<std::vector<uint64_t>>& offseted_level_mesh_encodings)
{
    /* Compute the number of mesh encoding levels */
    const int num_levels = this->coarsening_indications_.size();
    cmc_assert(num_levels == static_cast<int>(num_global_encoding_steps) - 1);

    /* Number of mesh encoding steps */
    const int num_partition_infos = static_cast<int>(this->coarsening_indications_.size() + 1);

    /* Allocate an output vector */
    this->num_elements_per_level_.clear();
    this->num_elements_per_level_.reserve(this->coarsening_indications_.size());

    cmc_assert(num_levels >= 1);

    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        uint64_t global_elems_level{0};
        /* Count this levels offset */
        for (int rank_idx{0}; rank_idx < this->comm_size_; ++rank_idx)
        {
            const int access_idx = rank_idx * num_global_encoding_steps + lvl_idx;
            global_elems_level += global_partition_info[access_idx];
        }

        this->num_elements_per_level_.push_back(global_elems_level);
        cmc_debug_msg("On level: ", lvl_idx, ", number of global elements: ", global_elems_level);
    }

    if (this->comm_rank_ != kRootRank)
    {
        /* Count overall length of accumulated message */
        uint64_t num_bytes_msg{0};
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
        const int rv_send = MPI_Ssend(message.data(), message.size(), MPI_UINT64_T, kRootRank, kTagMeshEncoding, this->comm_);
        MPICheckError(rv_send);
    } else
    {
        /* Compute the number of bytes for the global level-wise mesh encoding */
        uint64_t mesh_encoding_vals{0};
        for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
        {
            mesh_encoding_vals += (this->num_elements_per_level_[lvl_idx] / (sizeof(uint64_t) * cmc::bits::kCharBit) + (this->num_elements_per_level_[lvl_idx] % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0));
        }

        /* Compute the size of the file header */
        const uint64_t var_header_size_vals = GetVarHeaderSize(num_partition_infos) / sizeof(uint64_t);

        /* Allocate the encoded output stream on the root rank */
        std::vector<uint64_t> global_encoding(mesh_encoding_vals + var_header_size_vals, uint64_t{0});

        /* Write the variable header to the stream */
        this->InsertFileHeader(global_encoding, num_partition_infos);

        /* The mesh encoding starts directly after the variable header */
        uint64_t root_rank_encoding_offset = var_header_size_vals;

        /* Write the root-local mesh encodings per level to the stream */
        for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
        {
            /* Copy the root-local part into the global mesh encoding */
            std::copy_n(offseted_level_mesh_encodings[lvl_idx].begin(), offseted_level_mesh_encodings[lvl_idx].size(), std::next(global_encoding.begin(), root_rank_encoding_offset));

            /* Update the write offset */
            root_rank_encoding_offset += this->num_elements_per_level_[lvl_idx] / (sizeof(uint64_t) * cmc::bits::kCharBit) + (this->num_elements_per_level_[lvl_idx] % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0);
        }

        /* Compute the number of expected messages */
        const int num_expected_msgs = this->comm_size_ - 1;

        /* Declare a vector to receive the messages */
        std::vector<uint64_t> recv_msg;

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
            recv_msg = std::vector<uint64_t>(num_elements);

            /* Receive the message */
            const int rv_recv = MPI_Recv(recv_msg.data(), num_elements, MPI_UINT64_T, source_rank, kTagMeshEncoding, this->comm_, MPI_STATUS_IGNORE);
            MPICheckError(rv_recv);

            uint64_t global_start_offset{var_header_size_vals};
            uint64_t intra_msg_offset{0};

            /* Processs the message */
            for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
            {
                uint64_t current_elem_offset{0};
                /* Iterate over the level from all other ranks and append their encodings */
                for (int rank_idx{0}; rank_idx < source_rank; ++rank_idx)
                {
                    /* Compute the access index for this level on this rank */
                    const int access_idx = rank_idx * num_global_encoding_steps + lvl_idx;

                    /* Get the number of bits/elements on this level for this rank */
                    const uint64_t lvl_other_rank_num_elems = global_partition_info[access_idx];

                    /* Update the element offset */
                    current_elem_offset += lvl_other_rank_num_elems;
                }

                /* Now, we have reached the correct offset for this rank's message on this level */
                const int rank_access_idx = source_rank * num_global_encoding_steps + lvl_idx;
                const uint64_t lvl_rank_num_elems = global_partition_info[rank_access_idx];

                if (lvl_rank_num_elems > 0) [[likely]]
                {
                    const uint64_t lvl_rank_bits = lvl_rank_num_elems + (current_elem_offset % (sizeof(uint64_t) * cmc::bits::kCharBit));
                    const uint64_t num_vals_to_copy = lvl_rank_bits / (sizeof(uint64_t) * cmc::bits::kCharBit) + (lvl_rank_bits % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0);

                    const uint64_t global_encoding_start_offset = global_start_offset + current_elem_offset / (sizeof(uint64_t) * cmc::bits::kCharBit);

                    /** To ensure correctness, we always need to interleave the first and the last value in order to not overwrite the contents from the other ranks **/
                    /* Interleave the first value */
                    global_encoding[global_encoding_start_offset] |= recv_msg[intra_msg_offset];

                    /* Copy the intermediate values between start and end */
                    if (num_vals_to_copy >= 3)
                    {
                        std::copy_n(&(recv_msg[intra_msg_offset + 1]), num_vals_to_copy - 2, std::next(global_encoding.begin(), global_encoding_start_offset + 1));
                    }

                    /* Interleave the last value */
                    if (num_vals_to_copy >= 2)
                    {
                        global_encoding[global_encoding_start_offset + num_vals_to_copy - 1] |= recv_msg[intra_msg_offset + num_vals_to_copy - 1];
                    }

                    /* Update the element offset */
                    current_elem_offset += lvl_rank_num_elems;

                    /* Update the intra msg offset for the next level */
                    intra_msg_offset += (lvl_rank_bits / (sizeof(uint64_t) * cmc::bits::kCharBit) + (lvl_rank_bits % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0));
                }

                /* Update the offset to the start of the next level */
                global_start_offset += this->num_elements_per_level_[lvl_idx] / (sizeof(uint64_t) * cmc::bits::kCharBit) + (this->num_elements_per_level_[lvl_idx] % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0);
            }
        }
        
        /* After all expeceted messages have been received and processed, the encoding is completed */
        this->root_rank_compression_stream_ = std::move(global_encoding);
    }
}

inline void
Compressor::GenerateOuputStreams()
{
    /* Number of mesh encoding steps */
    const int num_partition_infos = static_cast<int>(this->coarsening_indications_.size() + 1);

    /* Allocate memory for the offsets and the encoding lengths */
    std::vector<uint64_t> global_partition_info(this->comm_size_ * num_partition_infos);

    /* Reverse */
    std::reverse(this->level_partitioning_info_.begin(), this->level_partitioning_info_.end());

    /* Allgather the local offsets and encoding lengths */
    const int rv_allgather = MPI_Allgather(this->level_partitioning_info_.data(), num_partition_infos, MPI_UINT64_T,
                                           global_partition_info.data(), num_partition_infos, MPI_UINT64_T, this->comm_);
    MPICheckError(rv_allgather);

    /* Exchange the mesh level encodings */
    std::vector<std::vector<uint64_t>> offseted_level_mesh_encodings = this->AdjustLevelMeshEncodings(global_partition_info, num_partition_infos);

    /* Collect the global level-wise mesh encoding on the root rank */
    this->CollectMeshEncodingOnTheRootRank(global_partition_info, num_partition_infos, offseted_level_mesh_encodings);

    /* Compute the global file sizes on all ranks */
    this->global_compressed_file_size_in_bytes_ = ComputeGlobalMeshEncodingNumBytes(this->num_elements_per_level_) + GetVarHeaderSize(num_partition_infos);

    if (this->comm_rank_ == kRootRank)
    {
        cmc_assert(this->global_compressed_file_size_in_bytes_ == this->root_rank_compression_stream_.size() * sizeof(uint64_t));
    }
}

inline void 
CheckMPIWriteCorrectness(const MPI_Status* status, const MPI_Datatype datatype, const int expected_num_elems)
{
    int elem_count_{0};
    const int rv_check_write = MPI_Get_count(status, datatype, &elem_count_);
    MPICheckError(rv_check_write);
    if (elem_count_ != expected_num_elems)
    {
        cmc_msg("The expeceted number of values could not be written to the file.");
        MPICheckError(MPI_ERR_COUNT);
    }
}

inline void
Compressor::WriteCompressedMeshStructure(const std::string& file_name, const bool block_mpi_processes_until_file_is_written)
{
    /* Generate the output streams for writing the data to disk */
    this->GenerateOuputStreams();

    /* We write the file from the root rank only */
    if (this->comm_rank_ == kRootRank)
    {
        /* Check if the output file exists, if so, we delete it, since the preallcoation fails if the file already exists */
        const std::filesystem::path output_file_path(file_name);
        if (std::filesystem::exists(output_file_path))
        {
            const int rv_delete = MPI_File_delete(file_name.c_str(), MPI_INFO_NULL);
            MPICheckError(rv_delete);
        }

        /* Allocate a MPI file handle */
        MPI_File fhandle;

        /* Open the file */
        constexpr int opening_mode = MPI_MODE_CREATE | MPI_MODE_WRONLY | MPI_MODE_EXCL;
        const int rv_open = MPI_File_open(MPI_COMM_SELF, file_name.c_str(), opening_mode, MPI_INFO_NULL, &fhandle);
        MPICheckError(rv_open);

        cmc_debug_msg(MPI_COMM_SELF, "The file ", file_name, " has been opened.");

        /* Set the native data representation */
        const int rv_file_view = MPI_File_set_view(fhandle, 0, MPI_BYTE, MPI_BYTE, "native", MPI_INFO_NULL);
        MPICheckError(rv_file_view);
        
        /* Pre-Allocate global file storage */
        const uint64_t preallocation_size = this->global_compressed_file_size_in_bytes_;
        cmc_debug_msg("File pre-allocation size: ", preallocation_size);
        const int rv_file_prealloc = MPI_File_preallocate(fhandle, preallocation_size);
        MPICheckError(rv_file_prealloc);
        cmc_debug_msg("The file has been pre-allocated.");

        /* Allocate an MPI status object */
        MPI_Status status;

        cmc_assert(preallocation_size == this->root_rank_compression_stream_.size() * sizeof(uint64_t));
        /* Write out the single stream making up the whole compressed file */
        const int rv_pos_fhandle = MPI_File_seek(fhandle, 0, MPI_SEEK_SET);
        MPICheckError(rv_pos_fhandle);

        /* Write the current stream to the file */
        const int rv_write_data_buffer = MPI_File_write(fhandle, this->root_rank_compression_stream_.data(), root_rank_compression_stream_.size(), MPI_UINT64_T, &status);
        MPICheckError(rv_write_data_buffer);
        CheckMPIWriteCorrectness(&status, MPI_UINT64_T, static_cast<int>(this->root_rank_compression_stream_.size()));
        cmc_debug_msg("The file has been written.");

        /* Close the file */
        const int rv_close = MPI_File_close(&fhandle);
        MPICheckError(rv_close);
        cmc_debug_msg("The file ", file_name, " has been closed.");
    }

    if (block_mpi_processes_until_file_is_written)
    {
        const int rv_barrier_wait_file_completion = MPI_Barrier(this->comm_);
        MPICheckError(rv_barrier_wait_file_completion);
        cmc_debug_msg(this->comm_, "The file ", file_name, " output has been completed.");
    }

    cmc_global_msg("The compressed global AMR mesh has been written to: ", file_name);
}

}

#endif /* !AMR_MESH_CMC_AMR_MESH_COMPRESSOR_HXX */
