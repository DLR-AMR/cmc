#ifndef CMC_PAR_MULTI_RES_EXTRACTION_HXX
#define CMC_PAR_MULTI_RES_EXTRACTION_HXX

#include "cmc.hxx"
#include "amr/lossless/cmc_par_multi_res_extraction_util.hxx"
#include "utilities/cmc_huffman_coder.hxx"
#include "mpi/cmc_mpi.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#include "t8code/cmc_t8_adaptation_callbacks.hxx"

#include <string>
#include <span>

namespace cmc::par::lossless::multi_res
{

template<ArithmeticType T>
class CompressionVariable
{
public:
    CompressionVariable() = delete;
    CompressionVariable(std::string name, t8_forest_t forest, std::span<T> data, const int num_data_per_element);

    void Compress();

private:
    void DetermineMaxInitElementLevel();
    void PerformIntraElementCompression();
    bool IsCompressionProgressing() const;
    void Repartition(t8_forest_t& adapted_mesh, std::vector<T>& adapted_data);
    bool HasIntraElementCompression() const;
    
    std::string name_;
    AmrMesh mesh_;

    const std::span<T> init_data_;
    const int32_t num_data_per_element_{1};
    int32_t max_init_elem_level_{kMaxPresentElementLevelUnknown};

    MPI_Comm comm_{MPI_COMM_NULL};

    std::vector<T> data_;
    std::vector<cmc::bits::vector> coarsening_indications_;
    std::vector<std::vector<LevelEncodingData<T>>> residual_encodings_;

    std::vector<std::vector<uint8_t>> levelwise_encoded_data_;
};

template<ArithmeticType T>
inline bool
CompressionVariable<T>::IsCompressionProgressing() const
{
    return (mesh_.GetNumberGlobalElements() > mesh_.GetNumberGlobalTrees());
}

template<ArithmeticType T>
inline bool
CompressionVariable<T>::HasIntraElementCompression() const
{
    return (num_data_per_element_ > 1);
}

template<ArithmeticType T>
struct CoarseningIterationData
{
    CoarseningIterationData(const std::span<T> current_data, const int32_t init_max_level, const int32_t coarsening_step)
    : data(current_data), current_coarsening_level{init_max_level - coarsening_step}
    {
        coarse_level_data.reserve(current_data_size / 4 + 1);
        coarsening_indications.Reserve(current_data_size / 4 + 1);
        residual_encodings.reserve(current_data_size / 4 + 1);
    }

    void LeaveElementUnchanged(const int local_idx);
    void PerformExtraction(const int local_idx, const int num_elements);

    const std::span<T> data;
    const int64_t current_coarsening_level;
    std::vector<T> coarse_level_data;
    cmc::bits::vector coarsening_indications;
    std::vector<LevelEncodingData<T>> residual_encodings;
};

template<ArithmeticType T>
inline void
CoarseningIterationData<T>::LeaveElementUnchanged(const int local_idx)
{
    /* Get the value from the current level */
    coarse_level_data.push_back(this->data[local_idx]);

    /* Indicate that no corsening has been performed */
    this->coarsening_indications.AppendUnsetBit();
}

template<ArithmeticType T>
inline void
CoarseningIterationData<T>::PerformExtraction(const int local_idx, const int num_elements)
{
    /* Indicate that corsening has been performed */
    this->coarsening_indications.AppendSetBit();

    /* Perform the multi-resolution extraction */
    std::array<T, num_elements + 2> predictors;
    std::copy_n(&(this->data[local_idx]), num_elements, predictors.begin())

    /* Append the arithmetic mean as a predictor */
    predictors[num_elements] = ComputeArithmeticMean<T>(predictors.data(), num_elements);

    /* Append the mid-range as a predictor */
    predictors[num_elements + 1] = ComputeMidRange<T>(predictors.data(), num_elements);

    int current_lzc{-1};
    T lzc_maximizing_predictor{};

    /* For all predictors, we evaluate the one that gives us the overall maximum number of leading zeros */
    for (int pred_idx{0}; pred_idx < num_elements + 2; ++pred_idx)
    {
        int cumulative_lzc{0};
        /* Check the predictor for all element values */
        for (int elem_idx{0}; elem_idx < num_elements; ++elem_idx)
        {
            /* Compute the residual */
            const auto [_, residual] = ComputeIntegerResidual<T>(predictors[pred_idx], predictors[elem_idx]);
            /* Add the LZC */
            cumulative_lzc += cmc::bits::GetLZC(residual);
        }

        /* If the predictor maximizes the LZC, we store it */
        if (current_lzc < cumulative_lzc)
        {
            current_lzc = cumulative_lzc;
            lzc_maximizing_predictor = predictors[pred_idx];
        }
    }

    /* We store the lzc_maximizing_predictor for the next coarse level */
    this->coarse_level_data.push_back(lzc_maximizing_predictor);

    /* Allocate the coarsening data for the coarsening of this family of elements */
    this->residual_encodings.emplace_back();
    /* Define a reference for freshly allocated LevelEncdoingData struct for the ease of notation */
    LevelEncodingData<T>& coarsening_data = this->residual_encodings.back();

    /* Store the number of elements */
    coarsening_data.num_elements = num_elements;

    /* After we have found the "best" predictor, we compute the residuals and the entropy symbols */
    for (int elem_idx{0}; elem_idx < num_elements; ++elem_idx)
    {
        /* Compute residual */
        const auto [is_pred_greater, residual] = ComputeIntegerResidual<T>(lzc_maximizing_predictor, predictors[elem_idx]);

        /* Determine the entropy code */
        coarsening_data.entropy_symbols[elem_idx] = CreateEntropySymbol(is_pred_greater, residual);

        /* Store the residual */
        coarsening_data.residuals[elem_idx] = residual;
    }
}

template<ArithmeticType T>
inline void
CompressionVariable<T>::Repartition(t8_forest_t& adapted_mesh, std::vector<T>& adapted_data)
{
    /** Partition the mesh **/
    /* Keep the not-partitioned forest */
    t8_forest_ref(adapted_forest);

    /* Allocate a forest */
    t8_forest_t partitioned_forest;
    t8_forest_init(&partitioned_forest);

    /* Partition the forest */
    const int partition_for_coarsening = 0; //TODO: change to 'one' when partition for coarsening is in t8code
    t8_forest_set_partition(partitioned_forest, adapted_forest, partition_for_coarsening);
    t8_forest_commit(partitioned_forest);

    /** Partition the adapted data **/
    /* Create an sc_array_t wrapper of the variable's data */
    sc_array_t* in_data = sc_array_new_data (static_cast<void*>(adapted_data.data()), sizeof(T), adapted_data.size());

    /* Allocate an output vector for the partitioned data */
    this->data_ = std::vector<T>(t8_forest_get_local_num_leaf_elements(partitioned_forest));

    /* Create a wrapper for the freshly allocated partitioned data */
    sc_array_t* out_data = sc_array_new_data (static_cast<void*>(this->data_.data()), sizeof(CompressionValue<T>), this->data_.size());

    /* Partition the variables data */
    t8_forest_partition_data(adapted_forest, partitioned_forest, in_data, out_data);

    /* Destroy the array wrappers */
    sc_array_destroy(in_data);
    sc_array_destroy(out_data);

    /* Free the former forest and store the adapted/repartitioned mesh */
    t8_forest_unref(&adapted_forest);
    mesh_.SetMesh(partitioned_forest);
}

template<ArithmeticType T>
void
CompressionVariable<T>::DetermineMaxInitElementLevel()
{
    /* If the maximum initial element level is not known, we need to gather it */
    if (this->max_init_elem_level_ != kMaxPresentElementLevelUnknown)
    {
        t8_forest_t mesh = this->GetAmrMesh().GetMesh();
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
void
CompressionVariable<T>::PerformIntraElementCompression()
{
    /* In case more than a single value is given on each element */
    if (num_data_per_element_ > 1)
    {
        //TODO
    }
}

constexpr bool
CheckIfElementIsEligibleForCoarsening(const int current_coarsening_level, const int element_level)
{
    return (current_coarsening_level == element_level);
}

template<typename T>
inline t8_locidx_t
LosslessMultiResCompression (t8_forest_t forest,
                             t8_forest_t forest_from,
                             t8_locidx_t which_tree,
                             const t8_eclass_t tree_class,
                             t8_locidx_t lelement_id,
                             const t8_scheme_c * ts,
                             const int is_family,
                             const int num_elements,
                             [[maybe_unused]] t8_element_t * elements[])
{
    /* Retrieve the adapt_data */
    CoarseningIterationData<T>* adapt_data = static_cast<CoarseningIterationData<T>*>(t8_forest_get_user_data(forest));
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
        adapt_data->PerformExtraction(local_start_index, num_elements);
        return cmc::t8::kCoarsenElements;
    }
}

template<ArithmeticType T>
void
CompressionVariable<T>::Compress()
{
    cmc_debug_msg("The lossless multi-resolution compression on variable ", this->name_, " is performed.");
    /* Potentially, gather the maximum present element level */
    this->DetermineMaxInitElementLevel();
    
    /* Allocate for the expected levelwise data-streams */
    this->coarsening_indications_.reserve(this->max_init_elem_level_ + 1);
    this->residual_encodings_.reserve(this->max_init_elem_level_ + 1);

    /* Potentially, perform intra-element compression, such that we obtain one data point per element */
    this->PerformIntraElementCompression();

    /* Potentially, perform partition for coarsening */
    this->Repartition(...);

    int32_t compression_step{0};

    // Compress until the root level
    while(this->IsCompressionProgressing())
    {
        cmc_debug_msg("A coarsening iteration is initialized.");

        // 1) Allocate an extraction iteration 
        CoarseningIterationData<T> adapt_data(std::span(this->data_), this->max_init_elem_level_, compression_step);

        // 2)  Adapt
        /* Perform a coarsening iteration */
        t8_forest_t adapted_forest = t8_forest_new_adapt(previous_forest, LosslessMultiResCompression<T>, 0, 0, static_cast<void*>(adapt_data));
        cmc_debug_msg("The mesh adaptation step is finished; resulting in ", t8_forest_get_global_num_leaf_elements(adapted_forest), " global elements");

        // 3) Store data to be encoded later
        this->coarsening_indications_.push_back(std::move(adapt_data.coarsening_indications));
        this->residual_encodings_.push_back(std::move(adapt_data.residual_encodings));

        // 4) Repartition the mesh and the data for the next iteration
        this->Repartition(adapted_forest, adapt_data.coarse_level_data);
        cmc_debug_msg("The mesh and the data has been re-partitioned.");

        cmc_debug_msg("This coarsening iteration is finished.");
        ++compression_step;
    }

    /* Encode the data that has been collected */
    this->EncodeData();

    cmc_debug_msg("The lossless multi-resolution compression of variable ", this->name_, " has been completed.");
}

template<ArithmeticType T>
std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>>
ExchangeEntropySymbols(const std::vector<std::vector<LevelEncodingData<T>>>& levelwise_residuals, const MPI_Comm comm)
{
    /* Get the number of all possible entropy symbols */
    constexpr int num_entropy_symbols = GetNumEntropySymbols<T>();

    /* Set the array and zero intialiaze the frequencies */
    std::array<uint64_t, num_entropy_symbols> entropy_symbol_frequencies{};

    /* Iterate through all entropy codes and accumulate their frequencies */
    for (size_t lvl_idx{0}; lvl_idx < levelwise_residuals.size(); ++lvl_idx)
    {
        /* Iterate through all coarsening data on this level */
        for (size_t coarsening_idx{0}; coarsening_idx < levelwise_residuals[lvl_idx].size(); ++coarsening_idx)
        {
            /* Iterate over all entropy codes from this coarsening data */
            for (int entropy_sym_idx{0}; entropy_sym_idx < levelwise_residuals[lvl_idx][coarsening_idx].num_elements; ++entropy_sym_idx)
            {
                /* Convert the symbol to the corresponding array index */
                const int array_idx = MapEntropySymbolToArrayIndex<T>(levelwise_residuals[lvl_idx][coarsening_idx].entropy_symbols[entropy_sym_idx]);
                /* Update the frequency */
                ++entropy_symbol_frequencies[array_idx];
            }
        }
    }

    /* Contract the redundant full LZC symobls and add the process end symbol */
    AddProcessEndSymbol<T>(entropy_symbol_frequencies, levelwise_residuals.size());

    /* After all entropy symbol frequencies have been coolected, we exchange them */
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

template<ArithmeticType T>
void
CompressionVariable<T>::EncodeData()
{
    /* Collect and exchange all entropy symbols */
    const std::vector<cmc::entropy_coding::huffman::EntropySymbol<SymbolType>> entropy_symbols = ExchangeEntropySymbols<T>(this->residual_encodings);
    
    /* Create a Huffman encoder */
    cmc::entropy_coding::huffman::HuffmanCoder<SymbolType> entropy_coder(entropy_symbols);

    /* Number of overall encdoing steps */
    const int num_encoding_steps = coarsening_indications_.size() + 1 + (this->HasIntraElementCompression() ? 1 : 0);

    /* Allocate an output vector for the levelwise encoding */
    std::vector<std::vector<uint8_t>> encoded_levelwise_data;
    encoded_levelwise_data.reserve(num_encoding_steps);

    /* We need to encode the data resididng on the root level */
    //TODO: Encode root level 
    //...

    /* We encode the data from the root level to the leaf level */
    lvl_iter = coarsening_indications_.rbegin();
    res_iter = residual_encodings_.rbegin();

    for (int step_idx{0}; step_idx < coarsening_indications_.size(); ++step_idx, ++lvl_iter, ++res_iter)
    {
        /* Allocate a bits::vector to store this level's encoded data */
        cmc::bits::vector lvl_data;
        lvl_data.reserve(res_iter->size() * sizeof(LevelEncodingData<T>) * cmc::bits::kCharBit);

        /* We define a view onto refinement indications */
        cmc::bits::vector_view lvl_view(*lvl_iter);

        /* Define a reference on the coarsening data for the ease of notation */
        const std::vector<LevelEncodingData<T>>& lvl_coarsening_data = *res_iter;

        /* Number of refinement indication bits on this level */
        const size_t num_elems = lvl_iter->size();
        
        int coarsening_data_idx{0};

        /* Iterate over this level's refinement indications */
        for (size_t elem_idx{0}; elem_idx < num_elems; ++elem_idx)
        {
            if (lvl_view.GetNextBit() == true)
            {
                /* In case a coarsening step has been performed */
                /* Set the bit for refinement */
                lvl_data.AppendSetBit();

                /* Define a refernce for the ease of notation */
                const LevelEncodingData<T>& coarse_data = lvl_coarsening_data[coarsening_data_idx];

                /* Store the entropy codes for this family of elements */
                for (int child_elem_idx{0}; child_elem_idx < coarse_data.num_elements; ++child_elem_idx)
                {
                    /* Encode the entropy symbol */
                    const cmc::entropy_coding::huffman::HuffmanCode code = entropy_coder.EncodeSymbol(coarse_data.entropy_symbols[child_elem_idx]);

                    /* Serialize the encoded entropy symbol */
                    lvl_data.AppendBits(code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - code.code_length), 0);

                    /* Retrieve the LZC from the entropy symbol */
                    const int lzc = GetLZCFromEntropySymbol(coarse_data.entropy_symbols[child_elem_idx]);

                    /* We do not need to encode the implicit given one-bit following the LZC */
                    if (lzc + 1 < sizeof(T) * cmc::bits::kCharBit) [[likely]]
                    {
                        /* Append the significant reisdual bits */
                        lvl_data.AppendBits(coarse_data.residuals[child_elem_idx], lzc + 1, 0);
                    }
                }

                /* Update the coarsening data accessing index */
                ++coarsening_data_idx;
            } else
            {
                /* In case the element remained unchanged */
                /* Set the bit that the element remains unchanged */
                lvl_data.AppendUnsetBit();
            }
        }

        /* At the end of the local encoding of the level, we append the process-end symbol */
        const cmc::entropy_coding::huffman::HuffmanCode process_lvl_end_code = entropy_coder.EncodeSymbol(kProcessEndSymbol<T>);
        /* Serialize the encoded process end symbol */
        lvl_data.AppendBits(process_lvl_end_code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - process_lvl_end_code.code_length), 0);

        /* We store this level's encoding in the variable's buffer */
        levelwise_encoded_data_.push_back(lvl_data.GetSerializedByteStream());
    }
}

}

#endif /* !CMC_PAR_MULTI_RES_EXTRACTION_HXX */
