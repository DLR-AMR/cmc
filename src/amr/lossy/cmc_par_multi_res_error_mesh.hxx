#ifndef CMC_AMR_LOSSY_PAR_MULTI_RES_ERROR_MESH_HXX
#define CMC_AMR_LOSSY_PAR_MULTI_RES_ERROR_MESH_HXX

#include "cmc.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"
#include "utilities/cmc_bits_vector.hxx"
#include "utilities/cmc_bits_vector_view.hxx"
#include "utilities/cmc_error_domain.hxx"

#include <t8_forest/t8_forest_partition.h>

#include <vector>
#include <cstdint>
#include <limits>

namespace cmc::par::lossy::multi_res
{

using bfloat16_t = uint16_t;

constexpr bfloat16_t kBfloat16Zero{0x0000};

inline const MPI_Datatype MPI_CMC_BFLOAT16_T = MPI_UINT16_T;

constexpr inline float
GetFloat(const bfloat16_t bfloat_value)
{
    static_assert(sizeof(uint32_t) == sizeof(float));
    const uint32_t fvalue = static_cast<uint32_t>(bfloat_value) << 16;
    return std::bit_cast<float>(fvalue);
}

constexpr inline bfloat16_t
GetBfloat16(const float fvalue)
{
    static_assert(sizeof(uint32_t) == sizeof(float));
    const uint32_t uvalue = std::bit_cast<uint32_t>(fvalue);
    const uint16_t bfloat_value = static_cast<uint16_t>(uvalue >> 16);
    return bfloat_value;
}

class ErrorMesh
{
public:
    template <typename T> ErrorMesh(t8_forest_t init_mesh, const std::vector<ErrorDomain>& error_domains, const std::vector<T>& init_data); //Compression Constructor 
    ErrorMesh(t8_forest_t root_level_mesh, const uint64_t* encoded_error_mesh, const MPI_Comm shm_comm, const int shm_rank, const int shm_size); //Decompression Constructor

    float GetPermittedAbsError(const t8_scheme_c* scheme, const int global_tree_id, const t8_eclass_t tree_class, const t8_element_t* element) const;

    int GetUniformErrorRefinementLevel(const int global_tree_id) const;

    std::vector<uint64_t> GetSerializedErrorMeshCopy() const {return error_mesh_encoding_;}
    const std::vector<uint64_t>& GetSerializedErrorMesh() const {return error_mesh_encoding_;}
    uint64_t GetSerializedErrorMeshEncodingSizeBytes() const {return encoded_error_mesh_size_;}

    void DestructErrorMeshCollectively();

    constexpr static float kCoarseningRelTolerance = 0.1;
    constexpr static int kNumObligatoryCoarseningIterations = 3;

private:
    template <typename T> void CreateErrorMesh(t8_forest_t init_mesh, const std::vector<ErrorDomain>& error_domains, const std::vector<T>& init_data);
    template <typename T> std::pair<t8_forest_t, std::vector<bfloat16_t>> GetCoarseErrors(t8_forest_t init_mesh, const std::vector<ErrorDomain>& error_domains, const std::vector<T>& init_data);
    void GatherTreeElementsOffsetAndErrors(t8_forest_t mesh, const uint64_t* permitted_abs_error_ptr);
    void ReconstructErrorMesh(t8_forest_t root_level_mesh, const uint64_t* encoded_error_mesh_start_ptr);
    t8_forest_t ReconstructCoarseMesh(t8_forest_t root_level_mesh, const std::vector<uint64_t>& global_elems_per_level, cmc::bits::vector_view encoded_mesh_view);

    MPI_Comm comm_{MPI_COMM_NULL};
    int comm_rank_{0}, comm_size_{0};

    MPI_Comm shm_comm_{MPI_COMM_NULL};
    int shm_rank_{0}, shm_size_{0};

    int num_global_trees_{0};

    MPI_Win window_uniform_permitted_abs_errors_;
    bfloat16_t* uniform_permitted_abs_errors_{nullptr};

    MPI_Win window_tree_elem_offsets_;
    uint64_t* tree_elem_offsets_{nullptr};

    MPI_Win window_uniform_tree_levels_;
    uint8_t* uniform_tree_levels_{nullptr};

    uint64_t encoded_error_mesh_size_;
    std::vector<uint64_t> error_mesh_encoding_; //Only the root rank in comm_ hold the serialization
};

template <typename T> 
ErrorMesh::ErrorMesh(t8_forest_t init_mesh, const std::vector<ErrorDomain>& error_domains, const std::vector<T>& init_data)
{
    /* Construct the error mesh and serialize it */
    this->CreateErrorMesh(init_mesh, error_domains, init_data);
}

ErrorMesh::ErrorMesh(t8_forest_t root_level_mesh, const uint64_t* encoded_error_mesh, const MPI_Comm shm_comm, const int shm_rank, const int shm_size)
: shm_comm_{shm_comm}, shm_rank_{shm_rank}, shm_size_{shm_size}
{
    this->ReconstructErrorMesh(root_level_mesh, encoded_error_mesh);
}

inline float
ErrorMesh::GetPermittedAbsError(const t8_scheme_c* scheme, const int global_tree_id, const t8_eclass_t tree_class, const t8_element_t* element) const
{
    cmc_assert(global_tree_id < this->num_global_trees_);

    /* Get the previous tree_offset */
    const uint64_t prior_tree_offset = this->tree_elem_offsets_[global_tree_id];

    /* Get the uniform refinement level of the corresponding tree */
    const int uniform_tree_level = this->uniform_tree_levels_[global_tree_id];

    /* Get the linear id of this element in the uniform refinement of this tree */
    const int lin_idx = scheme->element_get_linear_id(tree_class, element, uniform_tree_level);

    /* Compute the accessor for the permitted absolute errors */
    const uint64_t error_offset = prior_tree_offset + lin_idx;

    /* Return the float of the permitted absolute error */
    return GetFloat(this->uniform_permitted_abs_errors_[error_offset]);
}

inline int
ErrorMesh::GetUniformErrorRefinementLevel(const int global_tree_id) const
{
    cmc_assert(global_tree_id < this->num_global_trees_);
    return this->uniform_tree_levels_[global_tree_id];
}

inline bool IsFloatZero(const bfloat16_t value)
{
    /* There might be positive and negative float zero, therefore, we zero the (potential) sign bit and comapre afterwards */
    return ((value & 0x7FFF) == 0);
}

struct CoarseningErrors
{
    std::vector<bfloat16_t> current_permitted_errors;
    std::vector<bfloat16_t> coarse_permitted_errors;
};

template<typename T>
inline t8_locidx_t
CollectCoarseErrors (t8_forest_t forest,
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
    CoarseningErrors* adapt_data = static_cast<CoarseningErrors*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Compute the start offset in the local contiguous array of the data */
    const int local_start_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

    /* Check if a family is supplied to the adaptation function */
    if (is_family)
    {
        /* Compute the relative deviation from the currently permitted errors in order to assess the coarsening */
        bfloat16_t fam_min_bf16 = std::numeric_limits<bfloat16_t>::max();

        /* Since the absolute error thresholds are all positive, we can directly compare the bfloat16_t instead converting them to floats first */
        for (int elem_idx{0}; elem_idx < num_elements; ++elem_idx)
        {
            if (fam_min_bf16 > adapt_data->current_permitted_errors[local_start_index + elem_idx])
            {
                fam_min_bf16 = adapt_data->current_permitted_errors[local_start_index + elem_idx];
            }
        }
        
        /* If the minimum is zero */
        if (IsFloatZero(fam_min_bf16)) [[unlikely]]
        {
            /* Check if the other element values are also zero */
            for (int elem_idx{0}; elem_idx < num_elements; ++elem_idx)
            {
                if (IsFloatZero(adapt_data->current_permitted_errors[local_start_index + elem_idx]))
                {
                    goto LeaveThisElementUnchanged;
                }
            }

            /* If all element errors are zero, we are able to coarsen it */
            adapt_data->coarse_permitted_errors.push_back(kBfloat16Zero);
            return cmc::t8::kCoarsenElements;
        }

        const float fam_min = GetFloat(fam_min_bf16);

        /* Compute the relative deviation and whether we are able to coarsen the permitted error locally.
         * There are no zero values left in the error, otherwise we would not reach this section */
        for (int elem_idx{0}; elem_idx < num_elements; ++elem_idx)
        {
            const float elem_float_err = GetFloat(adapt_data->current_permitted_errors[local_start_index + elem_idx]);
            /* Compute the relative deviation */
            const float rel_deviation_elem = std::abs(elem_float_err - fam_min) / elem_float_err;    
            if (rel_deviation_elem > ErrorMesh::kCoarseningRelTolerance)
            {
                /* In case the deviation is "too large", we cannot coarsen the permitted errors */
               goto LeaveThisElementUnchanged; 
            }
        }

        /* Extract a value of the family and coarsen it */
        adapt_data->coarse_permitted_errors.push_back(fam_min_bf16);
        return cmc::t8::kCoarsenElements;
    }

LeaveThisElementUnchanged:
    /* The element stays unchanged */
    adapt_data->coarse_permitted_errors.push_back(adapt_data->current_permitted_errors[local_start_index]);
    return cmc::t8::kLeaveElementUnchanged;
}

inline bool
IsCoarseningErrorsProgressing(const t8_gloidx_t previous_num_elems, const t8_gloidx_t current_num_elems, const t8_gloidx_t num_global_trees_coarse_error_mesh)
{
    return (previous_num_elems > current_num_elems && current_num_elems > num_global_trees_coarse_error_mesh);
}

std::vector<uint8_t>
GatherMaximumRefinementLevelsPerTree(t8_forest_t mesh)
{
    /* Get the communicator of the mesh */
    const MPI_Comm comm = t8_forest_get_mpicomm(mesh);

    /* Get the number of global trees */
    const int num_global_trees = t8_forest_get_num_global_trees(mesh);

    /* Allocate a vector for exchanging */
    std::vector<uint8_t> max_lvl_per_tree(num_global_trees, uint8_t{0});

    /* Get the number of local trees */
    const int num_local_trees = t8_forest_get_num_local_trees(mesh);

    /* Get the scheme of the mesh */
    const t8_scheme *scheme = t8_forest_get_scheme (mesh);

    for (int tree_idx{0}; tree_idx < num_local_trees; ++tree_idx)
    {
        /* Convert the local to a global tree id */
        const int global_tree_idx = t8_forest_global_tree_id (mesh, tree_idx);

        /* Get the number of elements in this tree */
        const int num_tree_local_elems = t8_forest_get_tree_num_leaf_elements(mesh, tree_idx);

        /* Get the tree class */
        const t8_eclass_t tree_class = t8_forest_get_eclass(mesh, tree_idx);

        for (int elem_idx{0}; elem_idx < num_tree_local_elems; ++elem_idx)
        {
            /* Get the element in the tree */
            const t8_element_t* elem = t8_forest_get_leaf_element_in_tree(mesh, tree_idx, elem_idx);

            /* Get the level of this element */
            const uint8_t elem_level = static_cast<uint8_t>(scheme->element_get_level (tree_class, elem));

            /* Check whether it is greater than the currently found level */
            max_lvl_per_tree[global_tree_idx] = (max_lvl_per_tree[global_tree_idx] < elem_level ? elem_level : max_lvl_per_tree[global_tree_idx]);
        }
    }

    /* Allocate a vector for exchanged maximum refinement levels */
    std::vector<uint8_t> exchanged_max_lvl_per_tree(num_global_trees, uint8_t{0});

    /* After all local elements have been checked, we gather the global maxima */
    const int rv_allredc = MPI_Allreduce(max_lvl_per_tree.data(), exchanged_max_lvl_per_tree.data(), num_global_trees, MPI_UINT8_T, MPI_MAX, comm);
    MPICheckError(rv_allredc);

    return exchanged_max_lvl_per_tree;
}

struct RecursiveRefinementData
{
    RecursiveRefinementData() = delete;
    RecursiveRefinementData(const std::vector<uint8_t>& global_tree_levels_)
    : global_tree_levels{global_tree_levels_} {}

    const std::vector<uint8_t>& global_tree_levels;
};

struct SharedMemMessages
{
    std::vector<bfloat16_t> permitted_errors;
    int tag;
};

inline void
AppendUniformRefinementLevelsPerTree(std::vector<uint64_t>& error_mesh_serialization, const std::vector<uint8_t>& tree_levels)
{
    constexpr size_t levels_per_uint64_t = 8;
    const int num_complete_values = tree_levels.size() / levels_per_uint64_t;
    const int num_levels_in_incomplete_value = tree_levels.size() % levels_per_uint64_t;

    for (int idx{0}; idx < num_complete_values; ++idx)
    {
        /* Compute the value holding eight tree levels */
        const uint64_t val = (static_cast<uint64_t>(tree_levels[levels_per_uint64_t * idx]) << 56) |
                             (static_cast<uint64_t>(tree_levels[levels_per_uint64_t * idx + 1]) << 48) |
                             (static_cast<uint64_t>(tree_levels[levels_per_uint64_t * idx + 2]) << 40) |
                             (static_cast<uint64_t>(tree_levels[levels_per_uint64_t * idx + 3]) << 32) |
                             (static_cast<uint64_t>(tree_levels[levels_per_uint64_t * idx + 4]) << 24) |
                             (static_cast<uint64_t>(tree_levels[levels_per_uint64_t * idx + 5]) << 16) |
                             (static_cast<uint64_t>(tree_levels[levels_per_uint64_t * idx + 6]) << 8) |
                             (static_cast<uint64_t>(tree_levels[levels_per_uint64_t * idx + 7]));

        /* Store the value containing the levels */
        error_mesh_serialization.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(val));
    }

    /* Potentially, add an incomplete value with the remaining tree levels */
    if (num_levels_in_incomplete_value != 0)
    {
        uint64_t val{0};
        int shift_param{56};
        for (int idx{0}; idx < num_levels_in_incomplete_value; ++idx)
        {
            val |= (static_cast<uint64_t>(tree_levels[levels_per_uint64_t * num_complete_values + idx]) << shift_param);
            shift_param -= 8;
        }

        error_mesh_serialization.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(val));
    }
}

inline void
AppendCoarsePermittedAbsErrors(std::vector<uint64_t>& error_mesh_serialization, const uint64_t num_global_elems_error_mesh, const std::vector<bfloat16_t>& permitted_abs_errors)
{
    cmc_assert(num_global_elems_error_mesh == permitted_abs_errors.size());

    constexpr size_t errors_per_uint64_t = 4;
    const int num_complete_values = num_global_elems_error_mesh / errors_per_uint64_t;
    const int num_errors_in_incomplete_value = num_global_elems_error_mesh % errors_per_uint64_t;

    for (int idx{0}; idx < num_complete_values; ++idx)
    {
        /* Compute the value holding four permitted abs errors */
        const uint64_t val = (static_cast<uint64_t>(cmc::bits::ConvertToBigEndian<bfloat16_t>(permitted_abs_errors[errors_per_uint64_t * idx])) << 48) |
                             (static_cast<uint64_t>(cmc::bits::ConvertToBigEndian<bfloat16_t>(permitted_abs_errors[errors_per_uint64_t * idx + 1])) << 32) |
                             (static_cast<uint64_t>(cmc::bits::ConvertToBigEndian<bfloat16_t>(permitted_abs_errors[errors_per_uint64_t * idx + 2])) << 16) |
                             (static_cast<uint64_t>(cmc::bits::ConvertToBigEndian<bfloat16_t>(permitted_abs_errors[errors_per_uint64_t * idx + 3])));

        /* Store the value containing the levels */
        error_mesh_serialization.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(val));
    }

    /* Potentially, add an incomplete value with the remaining permitted errors */
    if (num_errors_in_incomplete_value != 0)
    {
        uint64_t val{0};
        int shift_param{48};
        for (int idx{0}; idx < num_errors_in_incomplete_value; ++idx)
        {
            val |= (static_cast<uint64_t>(cmc::bits::ConvertToBigEndian<bfloat16_t>(permitted_abs_errors[errors_per_uint64_t * num_complete_values + idx])) << shift_param);
            shift_param -= 16;
        }
        error_mesh_serialization.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(val));
    }
}

template <typename T>
struct FirstErrorCoarseningData
{
    FirstErrorCoarseningData(const int num_local_elems, const std::vector<ErrorDomain>& error_domains_, const std::vector<T>& initial_data)
    : error_domains{error_domains_}, init_data{initial_data}
    {
        permitted_abs_errors.reserve(num_local_elems);
    }

    std::vector<bfloat16_t> permitted_abs_errors;
    const std::vector<ErrorDomain>& error_domains;
    const std::vector<T>& init_data;
};

template<typename T>
inline t8_locidx_t
CoarsenAndEvaluatePermimttedErrors (t8_forest_t forest,
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
    FirstErrorCoarseningData<T>* adapt_data = static_cast<FirstErrorCoarseningData<T>*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    float min_abs_err{std::numeric_limits<float>::max()};

    /* Check all error domains */
    for (auto err_iter = adapt_data->error_domains.begin(); err_iter != adapt_data->error_domains.end(); ++err_iter)
    {
        /* Check if any of the considered elements are within the error domain */
        if (err_iter->IsAnyElementWithinDomain(forest_from, which_tree, tree_class, lelement_id, ts, num_elements, elements))
        {
            /* Get the permitted error */
            const PermittedError error_criterion = err_iter->GetPermittedError();

            /* Check if it is an absolute or relative error criterion */
            if (error_criterion.criterion == CompressionCriterion::AbsoluteErrorThreshold)
            {
                /* In case of an absolute error criterion, we can just compare the prescribed error with the current minumum */
                if (min_abs_err > error_criterion.error)
                {
                    min_abs_err = error_criterion.error;
                }
            } else
            {
                cmc_assert(error_criterion.criterion == CompressionCriterion::RelativeErrorThreshold);

                /* Get the local offset */
                const int local_start_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

                /* In case of a relative error, we need to compute all resulting absolute deviations */
                for (int elem_idx{0}; elem_idx < num_elements; ++elem_idx)
                {
                    if constexpr (std::is_signed_v<T>)
                    {
                        /* In case the data type is signed */
                        const float abs_err = std::fabs(static_cast<float>(adapt_data->init_data[local_start_index + elem_idx]) * error_criterion.error);
                        
                        /* Check if the deviation is smaller than the currently permitted abs error */
                        if (min_abs_err > abs_err)
                        {
                            min_abs_err = abs_err;
                        }
                    } else
                    {
                        /* In case the data type is unsigned */
                        const float abs_err = static_cast<float>(adapt_data->init_data[local_start_index + elem_idx]) * std::fabs(error_criterion.error);
                        
                        /* Check if the deviation is smaller than the currently permitted abs error */
                        if (min_abs_err > abs_err)
                        {
                            min_abs_err = abs_err;
                        }
                    }
                }
            }
        }
    }

    /* Store the computed absolute minimum for this element/family */
    adapt_data->permitted_abs_errors.push_back(GetBfloat16(min_abs_err));

    /* If it is a family, it will be coarsened, otherwise the element stays unchanged */
    if (is_family)
    {
        return cmc::t8::kCoarsenElements;
    } else
    {
        return cmc::t8::kLeaveElementUnchanged;
    }
}


template<typename T>
inline std::pair<t8_forest_t, std::vector<T>>
RepartitionErrorData(t8_forest_t adapted_mesh, std::vector<T>& adapted_data)
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

    return std::make_pair(partitioned_forest, std::move(partitioned_data));
}

struct ErrorCoarseningAdaptData
{
   std::vector<bfloat16_t> current_errors;
   std::vector<bfloat16_t> coarse_errors;
};

template<typename T>
inline t8_locidx_t
CoarsenPermimttedErrors (t8_forest_t forest,
                         t8_forest_t forest_from,
                         t8_locidx_t which_tree,
                         [[maybe_unused]] const t8_eclass_t tree_class,
                         t8_locidx_t lelement_id,
                         [[maybe_unused]] const t8_scheme_c * ts,
                         const int is_family,
                         const int num_elements,
                         [[maybe_unused]] t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    ErrorCoarseningAdaptData* adapt_data = static_cast<ErrorCoarseningAdaptData*>(t8_forest_get_user_data(forest));
    cmc_assert(adapt_data != nullptr);

    /* Get the local offset */
    const int local_start_index = t8_forest_get_tree_element_offset (forest_from, which_tree) + lelement_id;

    if (is_family)
    {
        /* Gather the minumum of all elements and coarsen the family */
        bfloat16_t min_abs_err{std::numeric_limits<bfloat16_t>::max()};

        /* SInce the absolute errot thresholds are all positive, we can directly compare the bfloat16_t instead converting them to floats first */
        for (int elem_idx{0}; elem_idx < num_elements; ++elem_idx)
        {
            if (min_abs_err > adapt_data->current_errors[local_start_index + elem_idx])
            {
                min_abs_err = adapt_data->current_errors[local_start_index + elem_idx];
            }
        }

        /* And we store the minumum for the coarser level */
        adapt_data->coarse_errors.push_back(min_abs_err);
        
        return cmc::t8::kCoarsenElements;
    } else
    {
        /* In case it is not a family, we leave the element unchanged */
        adapt_data->coarse_errors.push_back(adapt_data->current_errors[local_start_index ]);
        return cmc::t8::kLeaveElementUnchanged;
    }
}

template <typename T>
std::pair<t8_forest_t, std::vector<bfloat16_t>>
ErrorMesh::GetCoarseErrors(t8_forest_t init_mesh, const std::vector<ErrorDomain>& error_domains, const std::vector<T>& init_data)
{
    /* Reference the init mesh, since we do not own it */
    t8_forest_ref(init_mesh);

    /* Create the adapt data for the initial error coarsening */
    FirstErrorCoarseningData init_adapt_data(t8_forest_get_local_num_leaf_elements(init_mesh), error_domains, init_data);

    /* Perform the first coarsening iteration while evaluating the permitted errors */
    t8_forest_t coarsened_mesh = t8_forest_new_adapt(init_mesh, CoarsenAndEvaluatePermimttedErrors<T>, 0, 0, &init_adapt_data);

    /* Repartition the mesh and the data */
    auto [partitioned_mesh, partitoned_data] = RepartitionErrorData<bfloat16_t>(coarsened_mesh, init_adapt_data.permitted_abs_errors);

    /* Next, we perform some more obligatory coarsening iterations */
    constexpr int num_coarsening_iterations = ErrorMesh::kNumObligatoryCoarseningIterations - 1;

    t8_forest_t mesh = partitioned_mesh;
    std::vector<bfloat16_t> permitted_abs_errors = std::move(partitoned_data);

    /* Perform the obligatory coarsening steps */
    for (int coarsening_step{0}; coarsening_step < num_coarsening_iterations; ++coarsening_step)
    {
        ErrorCoarseningAdaptData adapt_data;
        adapt_data.current_errors = std::move(permitted_abs_errors);
        adapt_data.coarse_errors.reserve(t8_forest_get_local_num_leaf_elements(mesh));

        /* Perform the consecutive obligatory coarsening iterations */
        t8_forest_t coarsened_error_mesh = t8_forest_new_adapt(mesh, CoarsenPermimttedErrors<T>, 0, 0, &adapt_data);

        /* Repartition the mesh and the data */
        auto [coarse_partitioned_mesh, coarse_partitoned_data] = RepartitionErrorData<bfloat16_t>(coarsened_error_mesh, adapt_data.coarse_errors);

        mesh = coarse_partitioned_mesh;
        permitted_abs_errors = std::move(coarse_partitoned_data);
        cmc_debug_msg("Error coarsening iteration finished, step: ", coarsening_step);
    }

    /* After the obligatory coarsening steps have been performed, we try to further coarsen the permitted errors */    
    int64_t previous_num_elements = std::numeric_limits<int64_t>::max();

    t8_gloidx_t num_global_elems_coarse_error_mesh = t8_forest_get_global_num_leaf_elements(mesh);
    t8_gloidx_t num_global_trees_coarse_error_mesh = t8_forest_get_num_global_trees(mesh);

    /* Get the coarsest possible error mesh without restraining the permitted abs error too much */
    while(IsCoarseningErrorsProgressing(previous_num_elements, num_global_elems_coarse_error_mesh, num_global_trees_coarse_error_mesh))
    {
        previous_num_elements = num_global_elems_coarse_error_mesh;

        CoarseningErrors adapt_data;
        adapt_data.current_permitted_errors = std::move(permitted_abs_errors);
        adapt_data.coarse_permitted_errors.reserve(num_global_elems_coarse_error_mesh);
        
        /* Coarsen the permitted errors if possible */
        t8_forest_t coarsened_error_mesh = t8_forest_new_adapt(mesh, CollectCoarseErrors<T>, 0, 0, &adapt_data);

        /* Repartition the mesh and the data */
        auto [coarse_partitioned_mesh, coarse_partitoned_data] = RepartitionErrorData<bfloat16_t>(coarsened_error_mesh, adapt_data.coarse_permitted_errors);

        /* Store the coarsened errors and mesh */
        mesh = coarse_partitioned_mesh;
        permitted_abs_errors = std::move(coarse_partitoned_data);

        /* Update the global number of elements */
        num_global_elems_coarse_error_mesh = t8_forest_get_global_num_leaf_elements(mesh);

        cmc_debug_msg("Additional coarsening iteration is finished");
    }

    #if 1
    /* Write out the error mesh */
    std::vector<double> abs_errors_double;
    abs_errors_double.reserve(permitted_abs_errors.size());
    for(auto iter=permitted_abs_errors.begin(); iter != permitted_abs_errors.end(); ++iter)
    {
        abs_errors_double.push_back(static_cast<double>(GetFloat(*iter)));
    }
    t8_vtk_data_field_t vtk_data[1];
    snprintf (vtk_data[0].description, BUFSIZ, "PermittedAbsError");
    vtk_data[0].type = T8_VTK_SCALAR;
    vtk_data[0].data = abs_errors_double.data();

    t8_forest_write_vtk_ext (mesh, "cmc_test_error_mesh", 1, 1, 1, 1, 0, 0, 0, 1, vtk_data);
    #endif

    /* Finally, we received the coarse error mesh alognside the permitted absolute errors */
    return std::make_pair(mesh, std::move(permitted_abs_errors));
}

inline bool
IsErrorMeshCoarseningProgressing(t8_forest_t mesh)
{
    return (t8_forest_get_global_num_leaf_elements(mesh) > t8_forest_get_num_global_trees(mesh));
}

inline t8_locidx_t
CoarsenAll(t8_forest_t forest,
                         [[maybe_unused]]t8_forest_t forest_from,
                         [[maybe_unused]] t8_locidx_t which_tree,
                         [[maybe_unused]] const t8_eclass_t tree_class,
                         [[maybe_unused]] t8_locidx_t lelement_id,
                         [[maybe_unused]] const t8_scheme_c * ts,
                         const int is_family,
                         [[maybe_unused]] const int num_elements,
                         [[maybe_unused]] t8_element_t * elements[])
{
    cmc::bits::vector* coarsening_indications = static_cast<cmc::bits::vector*>(t8_forest_get_user_data(forest));

    if (is_family)
    {
        coarsening_indications->AppendSetBit();
        return cmc::t8::kCoarsenElements;
    } else
    {
        coarsening_indications->AppendUnsetBit();
        return cmc::t8::kLeaveElementUnchanged;
    }
}


inline std::pair<std::vector<uint64_t>, std::vector<std::vector<uint64_t>>>
AdjustLevelMeshEncodings(const std::vector<uint64_t>& num_local_elems_per_level, const std::vector<cmc::bits::vector>& level_indications, const MPI_Comm comm, const int comm_rank, const int comm_size)
{
    /* Exchange the local encoding lengths */
    const int num_levels = num_local_elems_per_level.size();
    std::vector<uint64_t> level_partition(comm_size * num_levels);

    /* Allgather the local offsets */
    const int rv_allgather = MPI_Allgather(num_local_elems_per_level.data(), num_levels, MPI_UINT64_T,
                                           level_partition.data(), num_levels, MPI_UINT64_T, comm);
    MPICheckError(rv_allgather);

    std::vector<std::vector<uint64_t>> offseted_level_mesh_encodings;
    offseted_level_mesh_encodings.reserve(num_levels);

    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        uint64_t lvl_elem_offset{0};

        /* Count this levels offset */
        for (int rank_idx{0}; rank_idx < comm_rank; ++rank_idx)
        {
            const int access_idx = rank_idx * num_levels + lvl_idx;
            lvl_elem_offset += level_partition[access_idx];
        }

        /* We only store the mesh level encoding if there are actual signficant bits, otherwise we store an empty vector */
        if (level_indications[lvl_idx].size() == 0) [[unlikely]]
        {
           offseted_level_mesh_encodings.emplace_back(std::vector<uint64_t>());
           continue;
        }

        /* Compute the local offset */
        const int shift = lvl_elem_offset % (sizeof(uint64_t) * cmc::bits::kCharBit);

        /* Get the offseted byte stream */
        offseted_level_mesh_encodings.emplace_back(level_indications[lvl_idx].GetSerializedOffsetByteStreamBE(shift));
    }

    return std::make_pair(std::move(level_partition), std::move(offseted_level_mesh_encodings));
}

constexpr int kTagErrorMeshEncoding = 1001;

inline std::pair<std::vector<uint64_t>, std::vector<uint64_t>>
CollectMeshEncodingOnTheRootRank(const std::vector<uint64_t>& num_local_elems_per_level, const std::vector<std::vector<uint64_t>>& offseted_level_mesh_encodings, const MPI_Comm comm, const int comm_rank, const int comm_size)
{
    constexpr int kRootRank = 0;

    /* Compute the number of mesh encoding levels */
    const int num_levels = offseted_level_mesh_encodings.size();

    /* Allocate an output vector */
    std::vector<uint64_t> elems_per_level;
    elems_per_level.reserve(offseted_level_mesh_encodings.size());

    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        uint64_t global_elems_level{0};
        /* Count this levels offset */
        for (int rank_idx{0}; rank_idx < comm_size; ++rank_idx)
        {
            const int access_idx = rank_idx * num_levels + lvl_idx;
            global_elems_level += num_local_elems_per_level[access_idx];
        }

        elems_per_level.push_back(global_elems_level);
    }

    if (comm_rank != kRootRank)
    {
        /* Count overall length of accumulated message */
        uint64_t num_vals_msg{0};
        for (const auto& lvl_mesh_encoding : offseted_level_mesh_encodings)
        {
            num_vals_msg += lvl_mesh_encoding.size();
        }

        std::vector<uint64_t> message;
        message.reserve(num_vals_msg);

        /* Setup the message accordingly */
        for (const auto& lvl_mesh_encoding : offseted_level_mesh_encodings)
        {
            std::copy_n(lvl_mesh_encoding.begin(), lvl_mesh_encoding.size(), std::back_inserter(message));
        }

        /* Send the message to the root rank */
        const int rv_send = MPI_Send(message.data(), message.size(), MPI_UINT64_T, kRootRank, kTagErrorMeshEncoding, comm);
        MPICheckError(rv_send);
    } else
    {
        /* Allocate a vector collecting the messages from the other ranks */
        const int num_expected_msgs = comm_size - 1;

        std::vector<std::vector<uint64_t>> recv_messages(num_expected_msgs);
        recv_messages.reserve(num_expected_msgs);

        /* Receive all messages from the other ranks collecting their process-local mesh encoding for all levels */
        for (int msg_idx{0}; msg_idx < num_expected_msgs; ++msg_idx)
        {
            /* Wait for a message to be received  */
            MPI_Status probe_status;
            const int rv_probe = MPI_Probe(MPI_ANY_SOURCE, kTagErrorMeshEncoding, comm, &probe_status);
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
            const int rv_recv = MPI_Recv(recv_messages[source_rank - 1].data(), num_elements, MPI_UINT64_T, source_rank, kTagErrorMeshEncoding, comm, MPI_STATUS_IGNORE);
            MPICheckError(rv_recv);
        }
        
        /* After all expeceted messages have been received, we process them and create the global level-wise mesh encoding */

        /* Compute the number of bytes for the global level-wise mesh encoding */
        uint64_t mesh_encoding_vals{0};
        for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
        {
            mesh_encoding_vals += (elems_per_level[lvl_idx] / (sizeof(uint64_t) * cmc::bits::kCharBit) + (elems_per_level[lvl_idx] % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0));
        }

        /* Allocate mesh encoding */
        std::vector<uint64_t> global_mesh_encoding;
        global_mesh_encoding.reserve(mesh_encoding_vals);

        std::vector<uint64_t> level_rank_offsets(comm_size - 1, 0);

        for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
        {
            uint64_t current_elem_offset{0};

            /* Copy the root-local part into the global mesh encoding */
            std::copy_n(offseted_level_mesh_encodings[lvl_idx].begin(), offseted_level_mesh_encodings[lvl_idx].size(), std::back_inserter(global_mesh_encoding));

            cmc_assert(not global_mesh_encoding.empty());
            if (global_mesh_encoding.empty()) [[unlikely]] {global_mesh_encoding.push_back(uint64_t{0});}
            
            /* Update the offset by the root local elements */
            current_elem_offset += num_local_elems_per_level[lvl_idx];

            /* Iterate over the root level from all other ranks and append their encodings */
            for (int rank_idx{1}; rank_idx < comm_size; ++rank_idx)
            {
                const int access_idx = rank_idx * num_levels + lvl_idx;

                const uint64_t lvl_rank_num_elems = num_local_elems_per_level[access_idx];

                if (lvl_rank_num_elems > 0) [[likely]]
                {
                    const uint64_t lvl_rank_bits = lvl_rank_num_elems + (current_elem_offset % (sizeof(uint64_t) * cmc::bits::kCharBit));
                    uint64_t start_idx{0};
                    uint64_t num_vals_to_copy = lvl_rank_bits / (sizeof(uint64_t) * cmc::bits::kCharBit) + (lvl_rank_bits % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0);

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
                        std::copy_n(&(recv_messages[rank_idx - 1].operator[](level_rank_offsets[rank_idx - 1] + start_idx)), num_vals_to_copy, std::back_inserter(global_mesh_encoding));
                    }

                    /* Update the element offset */
                    current_elem_offset += lvl_rank_num_elems;

                    /* Store the offset of this rank for the next level */
                    level_rank_offsets[rank_idx - 1] += (lvl_rank_bits / (sizeof(uint64_t) * cmc::bits::kCharBit) + (lvl_rank_bits % (sizeof(uint64_t) * cmc::bits::kCharBit) != 0 ? 1 : 0));
                }
            }
        }

        return std::make_pair(std::move(global_mesh_encoding), std::move(elems_per_level));
    }

    return std::make_pair(std::vector<uint64_t>(), std::move(elems_per_level));
}


inline std::tuple<t8_forest_t, std::vector<uint64_t>, std::vector<uint64_t>, int>
SerializeErrorMeshOnTheRootRank(t8_forest_t mesh, const MPI_Comm comm, const int comm_rank, const int comm_size)
{
    std::vector<cmc::bits::vector> coarsening_indications;
    std::vector<uint64_t> num_local_level_elems;

    t8_forest_t coarse_mesh = mesh;

    int adaptation_step_count{0};

    while (IsErrorMeshCoarseningProgressing(coarse_mesh))
    {   
        /* Allocate a bit-field */
        cmc::bits::vector level_indications;
        level_indications.Reserve(t8_forest_get_local_num_leaf_elements(coarse_mesh));

        /* Perform the coarsening to serialize this level's refinement structure */
        t8_forest_t coarsened_mesh = t8_forest_new_adapt(coarse_mesh, CoarsenAll, 0, 0, &level_indications);

        /* Get number of local elements */
        const int num_local_elems = t8_forest_get_local_num_leaf_elements(coarsened_mesh); 
        num_local_level_elems.push_back(num_local_elems);

        /* Store the bit-field as the coarsening indications */
        coarsening_indications.push_back(std::move(level_indications));

        /* Allocate a forest */
        t8_forest_t partitioned_forest;
        t8_forest_init(&partitioned_forest);

        /* Partition the forest */
        constexpr int partition_for_coarsening = 1;
        t8_forest_set_partition(partitioned_forest, coarsened_mesh, partition_for_coarsening);
        t8_forest_commit(partitioned_forest);

        /* Store the mesh for the next caorsening iteration */
        coarse_mesh = partitioned_forest;

        /* Update the step count */
        ++adaptation_step_count;
    }

    /* Reverse the partition table and the indications such that they run from the coarse to the fine level */
    std::reverse(num_local_level_elems.begin(), num_local_level_elems.end());
    std::reverse(coarsening_indications.begin(), coarsening_indications.end());

    /* Offset the mesh encodings correctly and gather the global level-wise partitioning */
    auto [level_partition_table, offseted_level_mesh_encodings] = AdjustLevelMeshEncodings(num_local_level_elems, coarsening_indications, comm, comm_rank, comm_size);

    /* Gather the global encoding on the root rank and retrieve the information about the number of global elements per level */
    auto [mesh_encdoing, global_elements_per_level] = CollectMeshEncodingOnTheRootRank(level_partition_table, offseted_level_mesh_encodings, comm, comm_rank, comm_size);

    return std::make_tuple(coarse_mesh, std::move(mesh_encdoing), std::move(global_elements_per_level), adaptation_step_count);
}


inline t8_locidx_t
RefineErrorMesh(t8_forest_t forest,
                [[maybe_unused]]t8_forest_t forest_from,
                [[maybe_unused]] t8_locidx_t which_tree,
                [[maybe_unused]] const t8_eclass_t tree_class,
                [[maybe_unused]] t8_locidx_t lelement_id,
                [[maybe_unused]] const t8_scheme_c * ts,
                [[maybe_unused]] const int is_family,
                [[maybe_unused]] const int num_elements,
                [[maybe_unused]] t8_element_t * elements[])
{
    cmc::bits::vector_view* refinement_indications = static_cast<cmc::bits::vector_view*>(t8_forest_get_user_data(forest));

    if (refinement_indications->GetNextBit())
    {
        /* Refine the element */
        return cmc::t8::kRefineElement;
    } else
    {
        /* Leave the element unchanged */
        return cmc::t8::kLeaveElementUnchanged;
    }
}

t8_forest_t
ErrorMesh::ReconstructCoarseMesh(t8_forest_t root_level_mesh, const std::vector<uint64_t>& global_elems_per_level, cmc::bits::vector_view encoded_mesh_view)
{
    /** Reconstruct the coarse error mesh **/
    const int num_levels = global_elems_per_level.size();
    uint64_t lvl_bit_offset{0};

    t8_forest_t mesh = root_level_mesh;

    for (int lvl_idx{0}; lvl_idx < num_levels; ++lvl_idx)
    {
        /* Get the global number of elements on this level */
        const uint64_t lvl_global_elems = global_elems_per_level[lvl_idx];

        /* Determine the global element offsets */
        const t8_gloidx_t elem_offset = t8_forest_get_first_local_leaf_element_id(mesh);

        /* Move to the correct bit position for the current level and partitioning */
        encoded_mesh_view.MoveToOffsetBitInStream(lvl_bit_offset + elem_offset);

        /* Adapt the mesh */
        t8_forest_t refined_mesh = t8_forest_new_adapt(mesh, RefineErrorMesh, 0, 0, &encoded_mesh_view);

        /* Repartition the mesh */
        t8_forest_t refined_partitioned_mesh;
        t8_forest_init(&refined_partitioned_mesh);
        t8_forest_set_partition(refined_partitioned_mesh, refined_mesh, 0);
        /* Create the new forest */
        t8_forest_commit(refined_partitioned_mesh);

        /* Set the mesh for the enxt refinement iteration */
        mesh = refined_partitioned_mesh;

        /* Compute the next level's offset */
        lvl_bit_offset += ((lvl_global_elems / 64) + (lvl_global_elems % 64 != 0 ? 1 : 0)) * 64;
    }

    return mesh;
}

void
ErrorMesh::GatherTreeElementsOffsetAndErrors(t8_forest_t mesh, const uint64_t* permitted_abs_error_ptr)
{
    /** Allocate and determine the global tree offsets **/
    const uint64_t num_coarse_error_mesh_elems = t8_forest_get_local_num_leaf_elements(mesh); 
    const uint64_t num_global_coarse_mesh_elems = t8_forest_get_global_num_leaf_elements(mesh);

    std::vector<uint32_t> num_elem_copies(num_coarse_error_mesh_elems, 0);

    std::vector<uint64_t> global_tree_offsets(this->num_global_trees_, 0);

    /* We iterate through the mesh and insert the necessary errors */
    const int num_current_local_trees = t8_forest_get_num_local_trees(mesh);

    /* Get the scheme of the mesh */
    const t8_scheme* scheme = t8_forest_get_scheme(mesh);

    uint64_t proc_offset{0};
    for (int tree_idx{0}, local_offset{0}; tree_idx < num_current_local_trees; ++tree_idx)
    {   
        /* Get the global tree id */
        const t8_gloidx_t gtree_id = t8_forest_global_tree_id (mesh, tree_idx);

        /* Get the prescribed element level */
        const int uniform_level = static_cast<int>(this->uniform_tree_levels_[gtree_id]);

        /* Get the number of elements in the tree */
        const int num_current_local_elems = t8_forest_get_tree_num_leaf_elements(mesh, tree_idx);
        
        uint64_t num_uniform_elems{0};

        /* Get the tree class */
        const t8_eclass_t tree_class = t8_forest_get_eclass(mesh, tree_idx);

        for (int elem_idx{0}; elem_idx < num_current_local_elems; ++elem_idx, ++local_offset)
        {
            /* Get the current element from the tree */
            const t8_element_t* element = t8_forest_get_leaf_element_in_tree (mesh, tree_idx, elem_idx);

            /* Count the leaf element that will be inserted within the uniform refinement */
            const int num_leaves = scheme->element_count_leaves(tree_class, element, uniform_level);

            /* Store the amount of copies that needs to be inserted */
            num_elem_copies[local_offset] = static_cast<uint32_t>(num_leaves);

            /* Count the inserted elements */
            num_uniform_elems += num_leaves;
        }

        /* Store the number of elements on the uniform level */
        global_tree_offsets[gtree_id] = num_uniform_elems;
        proc_offset += num_uniform_elems;
    }

    /* Now, we have vector indicating the local amount of copies for the permitted absolute errors */

    /* Exchange the global tree_offsets */
    std::vector<uint64_t> exchanged_global_tree_offsets(this->num_global_trees_, 0);
    const int rv_allgather_tree_offsets = MPI_Allreduce(global_tree_offsets.data(), exchanged_global_tree_offsets.data(), this->num_global_trees_, MPI_UINT64_T, MPI_SUM, this->comm_);
    MPICheckError(rv_allgather_tree_offsets);

    /* Next, we need to perform an exclusive scan, in order to obtain the global offsets */
    uint64_t acc_tree_offset{0};
    for (uint64_t celem_idx{0}; celem_idx < this->num_global_trees_; ++celem_idx)
    {
        const uint64_t num_uniform_elems_in_tree = exchanged_global_tree_offsets[celem_idx];
        exchanged_global_tree_offsets[celem_idx] = acc_tree_offset;
        acc_tree_offset += num_uniform_elems_in_tree;
    }

    /* Exchange the process local element lengths */
    std::vector<int> num_proc_local_elems(this->comm_size_, 0);
    const int num_local_elems = num_coarse_error_mesh_elems;
    const int rv_allgather_lengths = MPI_Allgather(&num_local_elems, 1, MPI_INT, num_proc_local_elems.data(), 1, MPI_INT, this->comm_);
    MPICheckError(rv_allgather_lengths);

    cmc_assert(std::reduce(num_proc_local_elems.begin(), num_proc_local_elems.end()) == num_global_coarse_mesh_elems);

    /* Exchange the amount of copies globally */
    std::vector<uint32_t> num_global_elem_copies(num_global_coarse_mesh_elems, 0);

    std::vector<int> gather_v_displ;
    gather_v_displ.reserve(this->comm_size_);
    int displ_offset{0};
    for (int idx{0}; idx < this->comm_size_; ++idx)
    {
        gather_v_displ.push_back(displ_offset);
        displ_offset += num_proc_local_elems[idx];
    }

    /* Exchange the amount of copies per element */
    const int rv_allgatherv = MPI_Allgatherv(num_elem_copies.data(), num_local_elems, MPI_UINT32_T, num_global_elem_copies.data(), num_proc_local_elems.data(), gather_v_displ.data(), MPI_UINT32_T, this->comm_);
    MPICheckError(rv_allgatherv);

    /* Compute the intervals for the ranks in shm_comm to process the levels */
    const uint64_t tree_start_offset = static_cast<uint64_t>(((static_cast<double>(this->shm_rank_) * static_cast<long double>(this->num_global_trees_)) / static_cast<double>(this->shm_size_)));
    const uint64_t tree_end_offset = static_cast<uint64_t>(((static_cast<double>(this->shm_rank_ + 1) * static_cast<long double>(this->num_global_trees_)) / static_cast<double>(this->shm_size_)));
    const uint64_t num_vals_tree_levels = tree_end_offset - tree_start_offset;

    /* Copy the tree levels into the shared memeory region */
    std::copy_n(exchanged_global_tree_offsets.data() + tree_start_offset, num_vals_tree_levels, this->tree_elem_offsets_ + tree_start_offset);

    /* Compute the shared comm intervals for inserting the copies the per element in the coarse error mesh */
    static_assert(sizeof(bfloat16_t) * 4 == sizeof(uint64_t));

    const uint64_t num_abs_error_vals = num_global_coarse_mesh_elems / 4;
    const uint64_t val_start_offset = static_cast<uint64_t>(((static_cast<double>(this->shm_rank_) * static_cast<long double>(num_abs_error_vals)) / static_cast<double>(this->shm_size_)));
    const uint64_t val_end_offset = static_cast<uint64_t>(((static_cast<double>(this->shm_rank_ + 1) * static_cast<long double>(num_abs_error_vals)) / static_cast<double>(this->shm_size_)));
    const uint64_t num_vals_errors = val_end_offset - val_start_offset;

    uint64_t elem_idx = val_start_offset * 4;

    /* Insert the permitted abs errors into the shared memory region */
    /* A mask to zero all bits except the relevant bits for the permitted abs errors */
    const uint64_t error_mask{0x000000000000FFFF};

    /* Compute the starting position within the shared memory region to insert the copies of the element's permitted absolute errors */
    uint32_t start_err_offset = std::reduce(num_global_elem_copies.data(), num_global_elem_copies.data() + val_start_offset * 4, uint32_t{0});
    for (uint64_t val_id{val_start_offset}; val_id < val_end_offset; ++val_id)
    {
        const uint64_t current_val = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(permitted_abs_error_ptr + val_id));

        /* Extarct the four permitted absolute error bounds from the value */
        std::array<bfloat16_t, 4> errors{
            cmc::bits::ConvertBigEndianToNativeEndianness<bfloat16_t>(static_cast<bfloat16_t>(current_val >> 48)),
            cmc::bits::ConvertBigEndianToNativeEndianness<bfloat16_t>(static_cast<bfloat16_t>((current_val >> 32) & error_mask)),
            cmc::bits::ConvertBigEndianToNativeEndianness<bfloat16_t>(static_cast<bfloat16_t>((current_val >> 16) & error_mask)),
            cmc::bits::ConvertBigEndianToNativeEndianness<bfloat16_t>(static_cast<bfloat16_t>(current_val & error_mask))
        };

        /* Insert the copies of the permitted errors */
        const uint32_t num_copies_1 = num_global_elem_copies[elem_idx];
        std::fill_n(this->uniform_permitted_abs_errors_ + start_err_offset, num_copies_1, errors[0]);
        start_err_offset += num_copies_1;

        const uint32_t num_copies_2 = num_global_elem_copies[elem_idx + 1];
        std::fill_n(this->uniform_permitted_abs_errors_ + start_err_offset, num_copies_2, errors[1]);
        start_err_offset += num_copies_2;

        const uint32_t num_copies_3 = num_global_elem_copies[elem_idx + 2];
        std::fill_n(this->uniform_permitted_abs_errors_ + start_err_offset, num_copies_3, errors[2]);
        start_err_offset += num_copies_3;

        const uint32_t num_copies_4 = num_global_elem_copies[elem_idx + 3];
        std::fill_n(this->uniform_permitted_abs_errors_ + start_err_offset, num_copies_4, errors[3]);
        start_err_offset += num_copies_4;

        /* Update the element offset to extract the next errors from the next value */
        elem_idx += 4;
    }

    /* The last rank needs to put the permitted errors from the potential incomplete value into the shared memory region */
    if (this->shm_rank_ == this->shm_size_ - 1)
    {
        uint32_t end_err_offset = std::reduce(num_global_elem_copies.data(), num_global_elem_copies.data() + num_abs_error_vals * 4, uint32_t{0});

        const int num_errs_incomplete_val = num_global_coarse_mesh_elems % 4;
        if (num_errs_incomplete_val != 0) [[likely]]
        {
            const uint64_t current_val_iv = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(permitted_abs_error_ptr + num_abs_error_vals));

            int shift{48};
            for (int err_idx{0}; err_idx < num_errs_incomplete_val; ++err_idx)
            {
                const uint32_t num_copies_iv = num_global_elem_copies[num_abs_error_vals * 4 + err_idx];
                const bfloat16_t err_val = cmc::bits::ConvertBigEndianToNativeEndianness<bfloat16_t>(static_cast<bfloat16_t>((current_val_iv >> shift) & error_mask));

                std::fill_n(this->uniform_permitted_abs_errors_ + end_err_offset, num_copies_iv, err_val);
                end_err_offset += num_copies_iv;
                shift -= 16;
            }
        }
    }

    /* Now, the tree offsets and the uniform permitted absolute errors are copied to the shared memory region */
    const int rv_shm_barr = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barr);
}

void
ErrorMesh::ReconstructErrorMesh(t8_forest_t root_level_mesh, const uint64_t* encoded_error_mesh_start_ptr)
{
    /* Reference the root level mesh, since we do not own it */
    t8_forest_ref(root_level_mesh);

    /* Get the MPI communicator of the mesh */
    const MPI_Comm comm = t8_forest_get_mpicomm(root_level_mesh);
    int rank, comm_size;
    /* Get the rank within and the size of the initial communicator */
    const int rv_rank = MPI_Comm_rank(comm, &rank);
    MPICheckError(rv_rank);
    const int rv_comm_size = MPI_Comm_size(comm, &comm_size);
    MPICheckError(rv_comm_size);

    /* Store the communicators */
    this->comm_ = comm;
    this->comm_rank_ = rank;
    this->comm_size_ = comm_size;

    const int kShmCommRootRank = 0;

    int offset{0};

    /* Get the global number of trees */
    const uint64_t num_global_trees = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_error_mesh_start_ptr + offset));
    ++offset;
    this->num_global_trees_ = num_global_trees;

    /* Get the global number of elements in the coarse error mesh */
    const uint64_t num_coarse_error_mesh_global_elems = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_error_mesh_start_ptr + offset));
    ++offset;

    /* Get the global number of uniform elements */
    const uint64_t num_global_elems = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_error_mesh_start_ptr + offset));
    ++offset;

    /* Get the number of bytes for the mesh encoding */
    const uint64_t num_mesh_encoding_bytes = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_error_mesh_start_ptr + offset));
    ++offset;
    cmc_assert(num_mesh_encoding_bytes % sizeof(uint64_t) == 0);

    /* Allocate the window for the permitted absolute errors */
    const uint64_t shared_win_proc_byte_length_permitted_abs_err = (this->shm_rank_ != kShmCommRootRank ? 0 : num_global_elems * sizeof(bfloat16_t));
    /* Create a shared memory window for the permitted absolute errors */
    bfloat16_t* shm_mem{nullptr};
    const int rv_win_alloc = MPI_Win_allocate_shared(static_cast<MPI_Aint>(shared_win_proc_byte_length_permitted_abs_err), sizeof(bfloat16_t), MPI_INFO_NULL, this->shm_comm_, &shm_mem, &(this->window_uniform_permitted_abs_errors_));
    MPICheckError(rv_win_alloc);
    /* Global start_ptr to the shared window data */
    bfloat16_t* shared_permitted_abs_error_ptr = (this->shm_rank_ != kShmCommRootRank ? shm_mem - num_global_elems : shm_mem);
    this->uniform_permitted_abs_errors_ = shared_permitted_abs_error_ptr;

    /* Create a shared window for the tree element offsets */
    const uint64_t shared_win_proc_byte_length_tree_elem_count = (this->shm_rank_ != kShmCommRootRank ? 0 : (num_global_trees + 1) * sizeof(uint64_t));
    uint64_t* shm_mem_tree_elem_count{nullptr};
    const int rv_win_alloc_tree_elem_count = MPI_Win_allocate_shared(static_cast<MPI_Aint>(shared_win_proc_byte_length_tree_elem_count), sizeof(uint64_t), MPI_INFO_NULL, this->shm_comm_, &shm_mem_tree_elem_count, &(this->window_tree_elem_offsets_));
    MPICheckError(rv_win_alloc_tree_elem_count);
    /* Global start_ptr to the shared window data */
    uint64_t* shared_permitted_tree_counts_ptr = (this->shm_rank_ != kShmCommRootRank ? shm_mem_tree_elem_count - (num_global_trees + 1) : shm_mem_tree_elem_count);
    this->tree_elem_offsets_ = shared_permitted_tree_counts_ptr;

    /* Allocate the shared memory window for the uniform tree levels */
    const uint64_t shared_win_proc_byte_length_tree_levels = (this->shm_rank_ != kShmCommRootRank ? 0 : num_global_trees * sizeof(uint8_t));
    uint8_t* shm_mem_tree_levels{nullptr};
    const int rv_win_alloc_tree_levels = MPI_Win_allocate_shared(static_cast<MPI_Aint>(shared_win_proc_byte_length_tree_levels), sizeof(uint8_t), MPI_INFO_NULL, this->shm_comm_, &shm_mem_tree_levels, &(this->window_uniform_tree_levels_));
    MPICheckError(rv_win_alloc_tree_levels);
    /* Global start_ptr to the shared window data */
    uint8_t* shared_permitted_tree_levels_ptr = (this->shm_rank_ != kShmCommRootRank ? shm_mem_tree_levels - num_global_trees : shm_mem_tree_levels);
    this->uniform_tree_levels_ = shared_permitted_tree_levels_ptr;

    /* Compute the intervals for the ranks in shm_comm to process the levels */
    const uint64_t num_values_tree_levels = num_global_trees / 8; //There are eight levels within one value
    const int start_offset = static_cast<uint64_t>(((static_cast<double>(this->shm_rank_) * static_cast<long double>(num_values_tree_levels)) / static_cast<double>(this->shm_size_)));
    const int end_offset =  static_cast<uint64_t>(((static_cast<double>(this->shm_rank_ + 1) * static_cast<long double>(num_values_tree_levels)) / static_cast<double>(this->shm_size_)));
    const int num_vals_to_read = end_offset - start_offset;

    /* Set the offset correctly */
    offset += start_offset;

    /* Create an offset for tree index */
    int tree_offset = start_offset * 8;

    /* A mask to zero all bits except the relevant bits for the tree-level */
    const uint64_t tree_lvl_mask{0x00000000000000FF};

    /* Get the uniform refinement levels per tree from the values */
    for (int idx{0}; idx < num_vals_to_read; ++idx)
    {
        /* Get the value */
        const uint64_t val = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_error_mesh_start_ptr + offset));

        /* Extract the eight encoded tree levels */
        shared_permitted_tree_levels_ptr[tree_offset] = static_cast<uint8_t>((val >> 56));
        ++tree_offset;
        shared_permitted_tree_levels_ptr[tree_offset] = static_cast<uint8_t>((val >> 48) & tree_lvl_mask);
        ++tree_offset;
        shared_permitted_tree_levels_ptr[tree_offset] = static_cast<uint8_t>((val >> 40) & tree_lvl_mask);
        ++tree_offset;
        shared_permitted_tree_levels_ptr[tree_offset] = static_cast<uint8_t>((val >> 32) & tree_lvl_mask);
        ++tree_offset;
        shared_permitted_tree_levels_ptr[tree_offset] = static_cast<uint8_t>((val >> 24) & tree_lvl_mask);
        ++tree_offset;
        shared_permitted_tree_levels_ptr[tree_offset] = static_cast<uint8_t>((val >> 16) & tree_lvl_mask);
        ++tree_offset;
        shared_permitted_tree_levels_ptr[tree_offset] = static_cast<uint8_t>((val >> 8) & tree_lvl_mask);
        ++tree_offset;
        shared_permitted_tree_levels_ptr[tree_offset] = static_cast<uint8_t>(val & tree_lvl_mask);
        ++tree_offset;

        /* Update the value accessor */
        ++offset;
    }

    /* The last rank may potentially add the tree-levels from a not full value */
    if (this->shm_rank_ == this->shm_size_ - 1)
    {
        const int num_levels_incomplete_val = num_global_trees % 8;
        if (num_levels_incomplete_val != 0) [[likely]]
        {
            const uint64_t val = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_error_mesh_start_ptr + 4 + num_values_tree_levels));

            tree_offset = num_values_tree_levels * 8;
            int shift{56};
            for (int lvl_idx{0}; lvl_idx < num_levels_incomplete_val; ++lvl_idx)
            {
                shared_permitted_tree_levels_ptr[tree_offset] =  static_cast<uint8_t>((val >> shift) & tree_lvl_mask);
                ++tree_offset;
                shift -= 8;
            }
        }
    }

    /* Reset the offset such that each ranks starts with the next encdoing section */
    const int mesh_encoding_offset = 4 + num_values_tree_levels + (num_global_trees % 8 != 0 ? 1 : 0);

    /* Get the number of adaptation level */
    const uint64_t num_adaptation_levels = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_error_mesh_start_ptr + mesh_encoding_offset));

    std::vector<uint64_t> global_elems_per_level;
    global_elems_per_level.reserve(num_adaptation_levels);
    for (int lvl_idx{0}; lvl_idx < num_adaptation_levels; ++lvl_idx)
    {
        const uint64_t lvl_elems = cmc::bits::ConvertBigEndianToNativeEndianness<uint64_t>(*(encoded_error_mesh_start_ptr + mesh_encoding_offset + 1 + lvl_idx));
        global_elems_per_level.push_back(lvl_elems);
    }

    /* Next within the encoding is the mesh serialization */
    cmc::bits::vector_view mesh_encoding(encoded_error_mesh_start_ptr + mesh_encoding_offset + 1 + num_adaptation_levels);

    /* Reconstruct the coarse error mesh and generate the tree element offset table */
    t8_forest_t coarse_error_mesh = this->ReconstructCoarseMesh(root_level_mesh, global_elems_per_level, mesh_encoding);

    cmc_assert(num_mesh_encoding_bytes % sizeof(uint64_t) == 0);

    /* Update the offset for accessing the permitted abs errors */
    const int mesh_errors_offset = mesh_encoding_offset + 1 + num_adaptation_levels + (num_mesh_encoding_bytes / sizeof(uint64_t));

    /* Gather the tree offesets and the uniform errors */
    this->GatherTreeElementsOffsetAndErrors(coarse_error_mesh, encoded_error_mesh_start_ptr + mesh_errors_offset);

    /* Afterwards, we can delete the mesh */
    t8_forest_unref(&coarse_error_mesh);
}

template <typename T>
void
ErrorMesh::CreateErrorMesh(t8_forest_t init_mesh, const std::vector<ErrorDomain>& error_domains, const std::vector<T>& init_data)
{
    constexpr int kRootRank = 0;

    /* Get the MPI communicator of the mesh */
    const MPI_Comm comm = t8_forest_get_mpicomm(init_mesh);
    int rank, comm_size;
    /* Get the rank within and the size of the initial communicator */
    const int rv_rank = MPI_Comm_rank(comm, &rank);
    MPICheckError(rv_rank);
    const int rv_comm_size = MPI_Comm_size(comm, &comm_size);
    MPICheckError(rv_comm_size);

    /* Split communicator into shared memory groups */
    int shm_rank, shm_size;
    MPI_Comm shm_comm;
    const int rv_split_comm = MPI_Comm_split_type(comm, MPI_COMM_TYPE_SHARED, 0, MPI_INFO_NULL, &shm_comm);
    MPICheckError(rv_split_comm);

    /* Get the rank within and the size of the shared communicator */
    const int rv_shm_rank = MPI_Comm_rank(shm_comm, &shm_rank);
    MPICheckError(rv_shm_rank);
    const int rv_shm_size = MPI_Comm_size(shm_comm, &shm_size);
    MPICheckError(rv_shm_size);

    /* Store the communicators */
    this->comm_ = comm;
    this->comm_rank_ = rank;
    this->comm_size_ = comm_size;
    this->shm_comm_ = shm_comm;
    this->shm_rank_ = shm_rank;
    this->shm_size_ = shm_size;

    /* Define the root ranks for the communicators */
    const int kCommRootRank = 0;
    const int kShmCommRootRank = 0;

    /* Create an inter-communicator of the root ranks within the shared memory communicators */
    const bool kParticipateInInterCommunication = (shm_rank == kShmCommRootRank);
    const int is_root_shm_comm_rank = (shm_rank == kShmCommRootRank ? 1 : MPI_UNDEFINED);
    MPI_Comm root_shm_comm;
    const int rv_split_root_shm_comm = MPI_Comm_split(comm, is_root_shm_comm_rank, 0, &root_shm_comm);
    MPICheckError(rv_split_root_shm_comm);

    /* Get the size of the communicator */
    int root_shm_comm_rank{0}, root_shm_comm_size{0};
    if (kParticipateInInterCommunication)
    {
        const int rv_root_shm_rank = MPI_Comm_rank(root_shm_comm, &root_shm_comm_rank);
        MPICheckError(rv_root_shm_rank);
        const int rv_root_shm_size = MPI_Comm_size(root_shm_comm, &root_shm_comm_size);
        MPICheckError(rv_root_shm_size);
    }

    /* Get the coarsest possible mesh alongside the errors */
    auto [mesh, local_permitted_abs_errors] = this->GetCoarseErrors(init_mesh, error_domains, init_data);
    cmc_assert(static_cast<size_t>(t8_forest_get_local_num_leaf_elements(mesh)) == local_permitted_abs_errors.size());

    /* Get the number of global elements in the coarse error mesh */
    const uint64_t num_coarse_error_mesh_global_elems = t8_forest_get_global_num_leaf_elements(mesh);

    /* We need to gather and exchange the maximum present refinement levels per tree */
    const std::vector<uint8_t> max_present_refinement_level_per_tree = GatherMaximumRefinementLevelsPerTree(mesh);

    /* Get the number of global trees */
    const t8_gloidx_t num_global_trees = t8_forest_get_num_global_trees(mesh);

    /* Store the number of global trees */
    this->num_global_trees_ = num_global_trees;

    std::vector<uint64_t> global_tree_offsets(num_global_trees, 0);

    /* Create the refinement adapt data */
    std::vector<bfloat16_t> recursively_refined_errors;

    /* Get the scheme of the mesh */
    const t8_scheme* scheme = t8_forest_get_scheme (mesh);

    /* We iterate through the mesh and insert the necessary errors */
    const int num_current_local_trees = t8_forest_get_num_local_trees(mesh);

    uint64_t proc_offset{0};
    for (int tree_idx{0}, local_offset{0}; tree_idx < num_current_local_trees; ++tree_idx)
    {   
        /* Get the global tree id */
        const t8_gloidx_t gtree_id = t8_forest_global_tree_id (mesh, tree_idx);

        /* Get the prescribed element level */
        const int uniform_level = static_cast<int>(max_present_refinement_level_per_tree[gtree_id]);

        /* Get the number of elements in the tree */
        const int num_current_local_elems = t8_forest_get_tree_num_leaf_elements(mesh, tree_idx);
        
        uint64_t num_uniform_elems{0};

        /* Get the tree class */
        const t8_eclass_t tree_class = t8_forest_get_eclass(mesh, tree_idx);

        for (int elem_idx{0}; elem_idx < num_current_local_elems; ++elem_idx, ++local_offset)
        {
            /* Get the current element from the tree */
            const t8_element_t* element = t8_forest_get_leaf_element_in_tree (mesh, tree_idx, elem_idx);

            /* Count the leaf element that will be inserted within the uniform refinement */
            const int num_leaves = scheme->element_count_leaves(tree_class, element, uniform_level);

            /* Get the permitted error for this element */
            const bfloat16_t permitted_abs_error = local_permitted_abs_errors[local_offset];

            /* Copy the permitted error num_leaves times to the error vector */
            std::fill_n(std::back_inserter(recursively_refined_errors), num_leaves, permitted_abs_error);

            /* Count the inserted elements */
            num_uniform_elems += num_leaves;
        }

        /* Store the number of elements on the uniform level */
        global_tree_offsets[gtree_id] = num_uniform_elems;
        proc_offset += num_uniform_elems;
    }

    /* Exchange the global tree_offsets */
    std::vector<uint64_t> exchanged_global_tree_offsets(num_global_trees, 0);
    const int rv_allgather_tree_offsets = MPI_Allreduce(global_tree_offsets.data(), exchanged_global_tree_offsets.data(), num_global_trees, MPI_UINT64_T, MPI_SUM, comm);
    MPICheckError(rv_allgather_tree_offsets);

    /* Perform an exclusive scan to determine the global offsets for the contiguous permitted abs error sequences */
    uint64_t global_proc_offset{0};
    const int rv_exscan_proc_offset = MPI_Exscan(&proc_offset, &global_proc_offset, 1, MPI_UINT64_T, MPI_SUM, comm);
    MPICheckError(rv_exscan_proc_offset);
    
    /* Since the value for the root rank is theoretically undefined after the exclusive scan, we explicitly set it again */
    if (rank == kRootRank)
    {
        global_proc_offset = 0;
    }

    /* Gather the not recursively refiend error mesh data on the root rank */
    const int num_local_elems_coarse_error_mesh = t8_forest_get_local_num_leaf_elements(mesh);
    std::vector<int> coarse_error_mesh_byte_offsets(comm_size, 0);
    const int rv_gather_offsets = MPI_Gather(&num_local_elems_coarse_error_mesh, 1, MPI_INT,
                                             coarse_error_mesh_byte_offsets.data(), 1, MPI_INT,  
                                             kCommRootRank, comm);
    MPICheckError(rv_gather_offsets);
    
    /* Compute the global number of elements */
    uint64_t num_global_elems{0};
    for (int tree_idx{0}; tree_idx < num_global_trees; ++tree_idx)
    {
        num_global_elems += exchanged_global_tree_offsets[tree_idx];
    }

    /** ALlocate shared memory windows for all neccessary information regarding the error retrieval */
    /* We allocate the window from the shared comm root rank */
    const uint64_t shared_win_proc_byte_length_permitted_abs_err = (shm_rank != kShmCommRootRank ? 0 : num_global_elems * sizeof(bfloat16_t));
    /* Create a shared memory window for the permitted absolute errors */
    bfloat16_t* shm_mem{nullptr};
    const int rv_win_alloc = MPI_Win_allocate_shared(static_cast<MPI_Aint>(shared_win_proc_byte_length_permitted_abs_err), sizeof(bfloat16_t), MPI_INFO_NULL, shm_comm, &shm_mem, &(this->window_uniform_permitted_abs_errors_));
    MPICheckError(rv_win_alloc);
    /* Global start_ptr to the shared window data */
    bfloat16_t* shared_permitted_abs_error_ptr = (shm_rank != kShmCommRootRank ? shm_mem - num_global_elems : shm_mem);
    this->uniform_permitted_abs_errors_ = shared_permitted_abs_error_ptr;

    /* Create a shared window for the tree element counts */
    const uint64_t shared_win_proc_byte_length_tree_elem_count = (shm_rank != kShmCommRootRank ? 0 : (num_global_trees + 1) * sizeof(uint64_t));
    uint64_t* shm_mem_tree_elem_count{nullptr};
    const int rv_win_alloc_tree_elem_count = MPI_Win_allocate_shared(static_cast<MPI_Aint>(shared_win_proc_byte_length_tree_elem_count), sizeof(uint64_t), MPI_INFO_NULL, shm_comm, &shm_mem_tree_elem_count, &(this->window_tree_elem_offsets_));
    MPICheckError(rv_win_alloc_tree_elem_count);
    /* Global start_ptr to the shared window data */
    uint64_t* shared_permitted_tree_counts_ptr = (shm_rank != kShmCommRootRank ? shm_mem_tree_elem_count - (num_global_trees + 1) : shm_mem_tree_elem_count);
    this->tree_elem_offsets_ = shared_permitted_tree_counts_ptr;

    /* Create a shared window for the uniform refinement levels */
    const uint64_t shared_win_proc_byte_length_tree_levels = (shm_rank != kShmCommRootRank ? 0 : num_global_trees * sizeof(uint8_t));
    uint8_t* shm_mem_tree_levels{nullptr};
    const int rv_win_alloc_tree_levels = MPI_Win_allocate_shared(static_cast<MPI_Aint>(shared_win_proc_byte_length_tree_levels), sizeof(uint8_t), MPI_INFO_NULL, shm_comm, &shm_mem_tree_levels, &(this->window_uniform_tree_levels_));
    MPICheckError(rv_win_alloc_tree_levels);
    /* Global start_ptr to the shared window data */
    uint8_t* shared_permitted_tree_levels_ptr = (shm_rank != kShmCommRootRank ? shm_mem_tree_levels - num_global_trees : shm_mem_tree_levels);
    this->uniform_tree_levels_ = shared_permitted_tree_levels_ptr;
    /** End of shared memory allocations */
    
    /* Declare a vector for the messages to be sent via the inter-communciator */
    std::vector<SharedMemMessages> inter_msgs;

    /* Gather the shared data from the intra-communicators on the corresponding root rank */
    if (shm_rank != kShmCommRootRank)
    {
        if (global_proc_offset > std::numeric_limits<int>::max())
        {
            cmc_err_msg("The error mesh is too large to be communicated in the implemented fashion.");
        }
        
        /* Set the global offset as a tag */
        const int kSendOffsetAsTag = static_cast<int>(global_proc_offset);

        /* In this case, the process sends the data to the root rank of the intra-communicator */
        const int rv_send_shm = MPI_Send(recursively_refined_errors.data(), proc_offset, MPI_CMC_BFLOAT16_T, kShmCommRootRank, kSendOffsetAsTag, shm_comm);
        MPICheckError(rv_send_shm);
    } else
    {
        inter_msgs.reserve(shm_size);

        /* In case we are the root rank, we collect the data from the intra-communicator */
        const int num_messages = shm_size - 1;

        /* Receive all messages from the intra-communicator */
        for (int msg_idx{0}; msg_idx < num_messages; ++msg_idx)
        {
            MPI_Status probe_status;

            /* Probe for a message */
            const int rv_probe_msg = MPI_Probe(MPI_ANY_SOURCE, MPI_ANY_TAG, shm_comm, &probe_status);
            MPICheckError(rv_probe_msg);

            /* Get the sending rank from the shared commnicator and the tag */
            const int source_rank = probe_status.MPI_SOURCE;
            const int msg_tag = probe_status.MPI_TAG;

            /* Get the count */
            int num_elements{0};
            const int rv_count = MPI_Get_count(&probe_status, MPI_CMC_BFLOAT16_T, &num_elements);
            MPICheckError(rv_count);

            /* Allocate the vector */
            inter_msgs.emplace_back();
            inter_msgs.back().permitted_errors = std::vector<bfloat16_t>(num_elements);
            inter_msgs.back().tag = msg_tag;

            /* Receive the message */
            const int rv_recv = MPI_Recv(inter_msgs.back().permitted_errors.data(), num_elements, MPI_CMC_BFLOAT16_T, source_rank, msg_tag, shm_comm, MPI_STATUS_IGNORE);
            MPICheckError(rv_recv);
        }

        /* Afterwards, we append the message from this root rank to be sent as well */
        if (global_proc_offset > std::numeric_limits<int>::max())
        {
            cmc_err_msg("The error mesh is too large to be communicated in the implemented fashion.");
        }
        inter_msgs.emplace_back();
        inter_msgs.back().permitted_errors = recursively_refined_errors;
        inter_msgs.back().tag = static_cast<int>(global_proc_offset);
    }

    /* The processes can put their data of the local permitted absolute errrors into the shared memory window already */
    std::copy_n(recursively_refined_errors.data(), proc_offset, shared_permitted_abs_error_ptr + global_proc_offset);


    #if 1

    /* The root rank from the shared memory communicators is busy receiving the data,
     * therefore, the other ranks put the data regarding the tree element counts and the levels into the windows */
    if (shm_size > 1 && shm_rank != kShmCommRootRank)
    {
        /* Determine the intervals the the processes need to write to */
        const int num_procs_to_write = shm_size - 1;

        const int start_offset = static_cast<int>(((static_cast<double>(shm_rank - 1) * static_cast<long double>(num_global_trees)) / static_cast<double>(num_procs_to_write)));
        const int end_offset =  static_cast<int>(((static_cast<double>(shm_rank) * static_cast<long double>(num_global_trees)) / static_cast<double>(num_procs_to_write)));

        const int num_vals_to_write = end_offset - start_offset;
        cmc_assert(num_vals_to_write >= 0);

        if (num_vals_to_write > 0)
        {
            /* Copy the refinement level data into the shared memory window */
            std::copy_n(max_present_refinement_level_per_tree.data() + start_offset, num_vals_to_write, this->uniform_tree_levels_ + start_offset);
        }

        /* Compute the initial offset to the start index */
        uint64_t tree_offset_count{0};
        for (int tree_idx{0}; tree_idx < start_offset; ++tree_idx)
        {
            tree_offset_count += exchanged_global_tree_offsets[tree_idx];
        }

        /* Accumulate and put the data into the shared memory window */
        for (int tree_idx{0}; tree_idx < num_vals_to_write; ++tree_idx)
        {
            this->tree_elem_offsets_[start_offset + tree_idx] = tree_offset_count;
            tree_offset_count += exchanged_global_tree_offsets[start_offset + tree_idx];
        }

        /* The last rank adds the overall count at last */
        if (shm_rank == shm_size - 1)
        {
            this->tree_elem_offsets_[num_global_trees] = tree_offset_count;
        }
    } else if (shm_size == 1)
    {
        cmc_assert(static_cast<size_t>(num_global_trees) == exchanged_global_tree_offsets.size());
        cmc_assert(static_cast<size_t>(num_global_trees) == max_present_refinement_level_per_tree.size());

        /* Only the root rank is available and puts the data into the windows */
        std::copy_n(max_present_refinement_level_per_tree.data(), num_global_trees, this->uniform_tree_levels_);

        /* Compute and store the tree offset table */
        uint64_t tree_offset_count{0};
        for (int tree_idx{0}; tree_idx < num_global_trees; ++tree_idx)
        {
            this->tree_elem_offsets_[tree_idx] = tree_offset_count;
            tree_offset_count += exchanged_global_tree_offsets[tree_idx];
        }
        /* Store the overall element count at last */
        this->tree_elem_offsets_[num_global_trees] = tree_offset_count;
    }

    #else
    //Let only the root rank write the data in any case
    if (shm_rank == kShmCommRootRank)
    {
        cmc_assert(static_cast<size_t>(num_global_trees) == exchanged_global_tree_offsets.size());
        cmc_assert(static_cast<size_t>(num_global_trees) == max_present_refinement_level_per_tree.size());

        /* Only the root rank is available and puts the data into the windows */
        std::copy_n(max_present_refinement_level_per_tree.data(), num_global_trees, this->uniform_tree_levels_);

        /* Compute and store the tree offset table */
        uint64_t tree_offset_count{0};
        for (int tree_idx{0}; tree_idx < num_global_trees; ++tree_idx)
        {
            this->tree_elem_offsets_[tree_idx] = tree_offset_count;
            tree_offset_count += exchanged_global_tree_offsets[tree_idx];
        }
        /* Store the overall element count at last */
        this->tree_elem_offsets_[num_global_trees] = tree_offset_count;
    }
    
    #endif

    /* Stop before moving onto sending the data via the inter-communicator */
    const int rv_barrier_send = MPI_Barrier(comm);
    MPICheckError(rv_barrier_send);

    if (kParticipateInInterCommunication)
    {
        std::vector<MPI_Request> requests_inter_msgs;
        requests_inter_msgs.reserve(root_shm_comm_size - 1);

        int msg_idx{0};

        /* Share the data on the root ranks of the intra communicators on write it into the shared memory window */
        for (auto inter_msg_iter = inter_msgs.begin(); inter_msg_iter != inter_msgs.end(); ++inter_msg_iter)
        {
            /* Send the message to all other root processes of the intra-communicators */
            for (int root_shm_rank_idx{0}; root_shm_rank_idx < root_shm_comm_size; ++root_shm_rank_idx)
            {
                if (root_shm_rank_idx != root_shm_comm_rank) [[likely]]
                {
                    /* Create a new request */
                    requests_inter_msgs.emplace_back();

                    /* Send the message non-blocking */
                    const int rv_isend_inter_msg = MPI_Isend(inter_msg_iter->permitted_errors.data(), inter_msg_iter->permitted_errors.size(), MPI_CMC_BFLOAT16_T, root_shm_rank_idx, inter_msg_iter->tag, root_shm_comm, &requests_inter_msgs[msg_idx]);
                    MPICheckError(rv_isend_inter_msg);

                    /* Update the message count */
                    ++msg_idx;
                }
            }
        }

        /* Impose to a barrier to be able to receive all messages afterwards */
        const int rv_barrier_stage_inter_msgs = MPI_Barrier(root_shm_comm);
        MPICheckError(rv_barrier_stage_inter_msgs);

        /* Receive all messages and copy them into the shared memory window */
        bool continue_receiving = true;

        /* Start receiving the messages into the shared memory window */
        while(continue_receiving)
        {
            /* Check if there is a message to receive */
            int is_message_waiting{0};
            MPI_Status status;
            const int rv_probe_ret_val = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, root_shm_comm, &is_message_waiting, &status);
            MPICheckError(rv_probe_ret_val);

            /* If there is a message, we receive it */
            if (is_message_waiting)
            {
                /* Get the count permitted abs errors */
                int count{0};
                const int rv_count_inter_msg = MPI_Get_count(&status, MPI_CMC_BFLOAT16_T, &count);
                MPICheckError(rv_count_inter_msg);

                /* The global offset is given as the tag */
                const int proc_offset = status.MPI_TAG;
                const int source_rank = status.MPI_SOURCE;

                /* Define the pointer to the shared memory window to the correct position */
                bfloat16_t* recv_data_ptr = shared_permitted_abs_error_ptr + proc_offset;

                /* Receive the actual message */
                const int rv_recv_inter_msg = MPI_Recv(recv_data_ptr, count, MPI_CMC_BFLOAT16_T, source_rank, proc_offset, root_shm_comm, MPI_STATUS_IGNORE);
                MPICheckError(rv_recv_inter_msg);
            }

            /* Update the loop flag, we continue until no more messages are queued */
            continue_receiving = is_message_waiting;
        }
    }

    /* Next, we synchronize the windows */
    const int rv_sync_win_abs_errors = MPI_Win_sync(this->window_uniform_permitted_abs_errors_);
    MPICheckError(rv_sync_win_abs_errors);

    const int rv_sync_win_offsets = MPI_Win_sync(this->window_tree_elem_offsets_);
    MPICheckError(rv_sync_win_offsets);

    const int rv_sync_win_levels = MPI_Win_sync(this->window_uniform_tree_levels_);
    MPICheckError(rv_sync_win_levels);

    /* Impose a global barrier after the global permitted abs errors have been reproduced in shared memeory windows */
    const int rv_barrier_replicated = MPI_Barrier(comm);
    MPICheckError(rv_barrier_replicated);

    /** At this moment, each shared memory communicator has a window with the global absolute permitted errors,
     * the tree levels and the tree element counts replicated **/

    /* Gather the serialized caorse error mesh on the root rank and retireve some additional information needed for the encoding */
    auto [root_mesh, mesh_encoding, elems_per_level, num_encoding_levels] = SerializeErrorMeshOnTheRootRank(mesh, comm, rank, comm_size);
    
    /* After the serialization of the mesh, we can get rid of it */
    t8_forest_unref(&root_mesh);

    /* Computed mesh serialization */
    uint64_t mesh_serialization_count{0};
    for (const auto num_elems_lvl : elems_per_level)
    {
        mesh_serialization_count += ((num_elems_lvl / 64) + (num_elems_lvl % 64 != 0 ? 1 : 0)) * sizeof(uint64_t);
    }

    /* Compute the size of the error mesh encoding */
    const uint64_t computed_mesh_encoding_size_bytes = sizeof(uint64_t) * 4 +
                                                       ((num_global_trees / 8) + (num_global_trees % 8 != 0 ? 1 : 0)) * sizeof(uint64_t) +
                                                       sizeof(uint64_t) + sizeof(uint64_t) * elems_per_level.size() +
                                                       mesh_serialization_count +
                                                       ((num_coarse_error_mesh_global_elems / 4) + (num_coarse_error_mesh_global_elems % 4 != 0 ? 1 : 0)) * sizeof(uint64_t);

    /* Store the error mesh encoding size */
    this->encoded_error_mesh_size_ = computed_mesh_encoding_size_bytes;

    /* Compute gather displacements */
    std::vector<int> gather_v_displ;
    if (rank == kCommRootRank)
    {
        gather_v_displ.reserve(comm_size);

        int displ_offset{0};
        for (int idx{0}; idx < comm_size; ++idx)
        {
            gather_v_displ.push_back(displ_offset);
            displ_offset += coarse_error_mesh_byte_offsets[idx];
        }
    }
    
    /* Gather the caorse error data */
    std::vector<bfloat16_t> gathered_coarse_error_data;
    if (rank == kCommRootRank)
    {
        gathered_coarse_error_data = std::vector<bfloat16_t>(num_coarse_error_mesh_global_elems);
    }
    const int rv_coarse_error_mesh_data = MPI_Gatherv(local_permitted_abs_errors.data(), num_local_elems_coarse_error_mesh, MPI_CMC_BFLOAT16_T,
                                                      gathered_coarse_error_data.data(), coarse_error_mesh_byte_offsets.data(), gather_v_displ.data(), MPI_CMC_BFLOAT16_T,
                                                      kCommRootRank, comm);
    MPICheckError(rv_coarse_error_mesh_data);

    /* In case of the root rank, we encode the data */
    if (rank == kCommRootRank)
    {
        std::vector<uint64_t> error_mesh_serialization;

        /* Store the global number of trees */
        error_mesh_serialization.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(num_global_trees));

        /* Store the global number of elements in the coarse error mesh */
        error_mesh_serialization.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(num_coarse_error_mesh_global_elems));
        
        /* Store the global numer of uniform elements */
        error_mesh_serialization.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(num_global_elems));

        /* Compute the length of the mesh encoding */
        const uint64_t num_mesh_encoding_bytes = mesh_encoding.size() * sizeof(uint64_t);

        /* Store the number of bytes for the mesh encoding */
        error_mesh_serialization.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(num_mesh_encoding_bytes));
        //cmc_err_msg("Stop here");
        /* Store the uniform refinement level per tree */
        cmc_assert(num_global_trees == max_present_refinement_level_per_tree.size());

        /* Append the uniform tree levels to the serialization */
        AppendUniformRefinementLevelsPerTree(error_mesh_serialization, max_present_refinement_level_per_tree);

        /* Store the number of adaptation levels */
        const uint64_t num_adaptation_levels = elems_per_level.size();
        error_mesh_serialization.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(num_adaptation_levels));

        /* Store the number of elements per level */
        for (uint64_t lvl_idx{0}; lvl_idx < num_adaptation_levels; ++lvl_idx)
        {
            error_mesh_serialization.push_back(cmc::bits::ConvertToBigEndian<uint64_t>(elems_per_level[lvl_idx]));
        }

        /* Append the coarse error mesh encoding */
        std::copy_n(mesh_encoding.data(), mesh_encoding.size(), std::back_inserter(error_mesh_serialization));

        /* Append the permitted abs errors on the coarse mesh */
        AppendCoarsePermittedAbsErrors(error_mesh_serialization, num_coarse_error_mesh_global_elems, gathered_coarse_error_data);

        /* Check whether the computation of the size equals the actual length */
        if (computed_mesh_encoding_size_bytes != error_mesh_serialization.size() * sizeof(uint64_t))
        {
            cmc_err_msg("The size of the error mesh encoding is unexpecetd and, therefore, corrupts the file layout!");
        }

        /* Store the serialized data */
        this->error_mesh_encoding_ = std::move(error_mesh_serialization);
    }

    /* Free the shared root communciator */
    if (kParticipateInInterCommunication)
    {
        const int rv_free_root_shm_comm = MPI_Comm_free(&root_shm_comm);
        MPICheckError(rv_free_root_shm_comm);
    }
}


inline void
ErrorMesh::DestructErrorMeshCollectively()
{
    /** Free the MPI Windows **/
    
    /* Impose a barrier to finish all outstanding reads */
    const int rv_shm_barrier = MPI_Barrier(this->shm_comm_);
    MPICheckError(rv_shm_barrier);

    /* Free the allocated window */
    const int rv_win_free_errors = MPI_Win_free(&(this->window_uniform_permitted_abs_errors_));
    MPICheckError(rv_win_free_errors);
    this->uniform_permitted_abs_errors_ = nullptr;

    const int rv_win_free_tree_elem_counts = MPI_Win_free(&(this->window_tree_elem_offsets_));
    MPICheckError(rv_win_free_tree_elem_counts);
    this->tree_elem_offsets_ = nullptr;

    const int rv_win_free_tree_levels = MPI_Win_free(&(this->window_uniform_tree_levels_));
    MPICheckError(rv_win_free_tree_levels);
    this->uniform_tree_levels_ = nullptr;

    /* Free the shared memory communicator */
    const int rv_free_shm_comm = MPI_Comm_free(&(this->shm_comm_));
    MPICheckError(rv_free_shm_comm);
}


}


#endif /* !CMC_AMR_LOSSY_PAR_MULTI_RES_ERROR_MESH_HXX */
