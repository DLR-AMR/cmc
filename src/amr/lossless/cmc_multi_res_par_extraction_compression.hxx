#ifndef CMC_MULTI_RES_PAR_EXTRACTION_COMPRESSION_HXX
#define CMC_MULTI_RES_PAR_EXTRACTION_COMPRESSION_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_compression_values.hxx"
#include "utilities/cmc_interpolation_fn.hxx"
#include "utilities/cmc_serialization.hxx"
#include "utilities/cmc_huffman_coder.hxx"
#include "utilities/cmc_multi_res_extraction_util.hxx"
#include "utilities/cmc_lossless_multi_res_extraction_residual_computation.hxx"

#include "amr/lossless/cmc_byte_par_compression_variable.hxx"
#include "amr/lossless/cmc_multi_res_par_extraction_util.hxx"

#include "mpi/cmc_mpi_data.hxx"

#include <utility>
#include <vector>
#include <array>
#include <algorithm>

namespace cmc::lossless::par::multi_res
{

/* Indicate whether unchanged elements have an encoded zero reisdual or no residual at all (in this case information from the mesh refinement bits is needed)
 * Moreover, this according setting needs to set in the decompression routine */
constexpr bool kEncodeLessResiduals = false;

template <typename T>
class MultiResAdaptData : public cmc::lossless::par::ICompressionAdaptData<T>
{
public:
    MultiResAdaptData() = delete;
    MultiResAdaptData(AbstractByteCompressionVariable<T>* variable)
    : ICompressionAdaptData<T>(variable) {};

    void InitializeExtractionIteration() override;
    void FinalizeExtractionIteration() override;
    void CompleteExtractionIteration(const t8_forest_t previous_forest, const t8_forest_t adapted_forest) override;
    void RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest) override;

    std::vector<uint8_t> EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const override;
    std::vector<uint8_t> EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const override;

    std::vector<uint8_t> StorePartitionTableOnTheRootRank(const t8_forest_t coarsened_forest) const override;
    std::vector<uint8_t> StoreRootLevelPartitionTableOnTheRootRank(const t8_forest_t coarsened_forest) const override;

protected:
    ExtractionData<T> PerformExtraction(const int which_tree, const int lelement_id, const int num_elements, const VectorView<CompressionValue<T>> values) override;
    UnchangedData<T> ElementStaysUnchanged(const int which_tree, const int lelement_id, const CompressionValue<T>& value) override;

private:
    std::vector<cmc::entropy_coding::huffman::HuffmanSymbol<uint32_t>> CollectGlobalSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values, bit_map::BitMapView is_residual_present, bit_map::BitMapView residual_indications) const;

    bit_map::BitMap resdiual_order_indications_;
    bit_map::BitMap residual_presence_indications_;
    int count_adaptation_step_{0};
    uint32_t num_levelwise_entropy_codes_{0};
    uint32_t num_local_entropy_bytes_encoded_level_data_{0};
};

template <typename T>
void
MultiResAdaptData<T>::InitializeExtractionIteration()
{
    resdiual_order_indications_ = bit_map::BitMap();
    residual_presence_indications_ = bit_map::BitMap();
    num_levelwise_entropy_codes_ = 0;
    num_local_entropy_bytes_encoded_level_data_ = 0;
}

template <typename T>
void
MultiResAdaptData<T>::FinalizeExtractionIteration()
{
    ++count_adaptation_step_;
}

template <typename T>
void
MultiResAdaptData<T>::CompleteExtractionIteration([[maybe_unused]] const t8_forest_t previous_forest, [[maybe_unused]] const t8_forest_t adapted_forest)
{
    //Nothing to be done here!
}

template <typename T>
void
MultiResAdaptData<T>::RepartitionData(const t8_forest_t adapted_forest, const t8_forest_t partitioned_forest)
{
    //Currently, nothing to be done here!
}

template<typename T>
T
GetCoarseApproximationMaximizingResidualsLZC(const VectorView<CompressionValue<T>> values)
{
    cmc_assert(values.size() >= 1);

    T current_best_predictor{T()};
    int current_max_lzc{-1};

    /* Try each value from the view as a predictor */
    for (auto pred_iter = values.begin(); pred_iter != values.end(); ++pred_iter)
    {
        /* Convert the compression value back to its origianl type */
        const T predictor = pred_iter->template ReinterpretDataAs<T>();

        /* Compute the LZC in the residuals for this approximation */
        const int pred_lzc = cmc::lossless::multi_res::GetCumulativeResidualsLZC<T>(predictor, values);

        /* Potentially, update the current predictor */
        if (pred_lzc > current_max_lzc)
        {
            current_best_predictor = predictor;
            current_max_lzc = pred_lzc;
        }
    }

    /* Convert the view to actual values of the underlying data type */
    const std::vector<T> converted_vals = ConvertCompressionValues<T>(values);

    /* Try the mid-range as an predictor */
    const T mid_range = InterpolateToMidRange<T>(converted_vals);

    /* Compute the LZC in the residuals for this approximation */
    const int mid_range_lzc = cmc::lossless::multi_res::GetCumulativeResidualsLZC<T>(mid_range, values);

    /* Potentially, update the current predictor */
    if (mid_range_lzc > current_max_lzc)
    {
        current_best_predictor = mid_range;
        current_max_lzc = mid_range_lzc;
    }

    /* Try the arithmetic mean as an predictor */
    const T mean = InterpolateToArithmeticMean<T>(converted_vals);

    /* Compute the LZC in the residuals for this approximation */
    const int mean_lzc = cmc::lossless::multi_res::GetCumulativeResidualsLZC<T>(mean, values);

    /* Potentially, update the current predictor */
    if (mean_lzc > current_max_lzc)
    {
        current_best_predictor = mean;
        current_max_lzc = mean_lzc;
    }

    return current_best_predictor;
}

template <typename T>
ExtractionData<T>
MultiResAdaptData<T>::PerformExtraction([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, const int num_elements, const VectorView<CompressionValue<T>> values)
{
    /* Get the coarse approximation for these values which maximizes the cumulative residual LZC */
    const T coarse_approximation = GetCoarseApproximationMaximizingResidualsLZC<T>(values);

    std::vector<CompressionValue<T>> fine_values;
    fine_values.reserve(values.size());

    /* Compute the residuals for the given predictor */
    for (auto val_iter = values.begin(); val_iter != values.end(); ++val_iter)
    {
        /* Compute the residual between approximation and actual value */
        auto [is_approximation_greater, residual] = cmc::lossless::multi_res::ComputeResidual<T>(coarse_approximation, *val_iter);

        /* Store whether the approximation is greater than the real value */
        resdiual_order_indications_.AppendBit(is_approximation_greater);

        /* Store the residual */
        fine_values.push_back(residual);

        residual_presence_indications_.AppendSetBit();
    }

    /* Store the number of levelwise entropy codes */
    num_levelwise_entropy_codes_ += num_elements;

    return ExtractionData<T>(CompressionValue<T>(coarse_approximation), std::move(fine_values));
}

template <typename T>
UnchangedData<T>
MultiResAdaptData<T>::ElementStaysUnchanged([[maybe_unused]] const int which_tree, [[maybe_unused]] const int lelement_id, const CompressionValue<T>& value)
{
    residual_presence_indications_.AppendUnsetBit();
    return UnchangedData<T>(value, CompressionValue<T>());
}

template <typename T>
std::vector<cmc::entropy_coding::huffman::HuffmanSymbol<uint32_t>>
MultiResAdaptData<T>::CollectGlobalSymbolFrequenciesForEntropyCoding(const std::vector<CompressionValue<T>>& level_byte_values, bit_map::BitMapView is_residual_present, bit_map::BitMapView residual_indications) const
{
    /* Get the MPI communicator */
    const MPI_Comm comm = this->GetMPIComm();
    int rank{0}, size{0};

    /* Get the rank of the process in the communicator */
    const int rv_rank = MPI_Comm_rank(comm, &rank);
    MPICheckError(rv_rank);

    /* Get the size of the communicator */
    const int rv_size = MPI_Comm_size(comm, &size);
    MPICheckError(rv_size);

    /******** Compute the local symbol frequencies ********/
    /* Vector that stores the frequencies */
    std::vector<entropy_coding::huffman::FrequencyType> frequencies(GetNumEntropySymbols<T>(), 0);

    /* Iterate over all values and encode them */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /* If there is no family that could be coarsened, we do not need to encode a residual for the element since the value remains unchanged */
        if (is_residual_present.GetNextBit() == false)
        {
            continue;
        }

        /* Get the current value */
        CompressionValue<T> val = *val_iter;

        /* Get the LZC */
        const uint32_t first_one_bit = val.GetNumberLeadingZeros();

        /* Get the info whether the capproximation was greater or smaller */
        const bool next_residual_indication = residual_indications.GetNextBit();

        /* Check whether the byte value has been fully extracted */
        if (val.GetNumberLeadingZeros() == sizeof(T) * bit_map::kCharBit)
        {
            /* Update the count for the fully extracted residuals */
            ++frequencies[GetFullyExtractedSymbol<T>()];
        } else
        {
            /* Update the frequency for the corresponding symbol */
            ++frequencies[ConvertToSymbolInFrequencyTable<T>(next_residual_indication, first_one_bit)];
        }
    }

    /* Add a one frequency for the process end (on this local process) */
    ++frequencies[GetProcessEndSymbol<T>()];
    /******** END of Compute the local symbol frequencies ********/

    /******** Exchange the frequencies globally ********/
    /* Allocate a receive buffer for the frequencies */
    std::vector<entropy_coding::huffman::FrequencyType> exchanged_frequencies(GetNumEntropySymbols<T>(), 0);
    cmc_assert(exchanged_frequencies.size() == frequencies.size());

    /* Get the corresponding MPI data type */
    MPI_Datatype mpi_freq_type = ConvertToMPIType<entropy_coding::huffman::FrequencyType>();

    /* Reduce the frequencies globally */
    const int rv_allreduce = MPI_Allreduce(frequencies.data(), exchanged_frequencies.data(), frequencies.size(), mpi_freq_type, MPI_SUM, comm);
    MPICheckError(rv_allreduce);
    /******** END of Exchange the frequencies globally ********/


    /* Allocate the local symbol frequency table locally on all processes */
    std::vector<cmc::entropy_coding::huffman::HuffmanSymbol<uint32_t>> sym_freq_table;
    sym_freq_table.reserve(GetNumEntropySymbols<T>());

    /* Fill the frequency table */
    for (uint32_t iter{0}; iter < GetNumEntropySymbols<T>(); ++iter)
    {
        if (exchanged_frequencies[iter] > 0)
        {
            sym_freq_table.emplace_back(iter, exchanged_frequencies[iter]);
        }
    }

    return sym_freq_table;
}

/**
 * @brief  We use an arithmetic encoder to encode the position of the first "one-bit" in the compression value
 * 
 * @tparam T The original data type of the underlying data (e.g. float)
 * @param level_byte_values The remaining "fine compression" values after an extraction iteration
 * @return std::vector<uint8_t> The encoded data stream
 */
template <typename T>
std::vector<uint8_t>
MultiResAdaptData<T>::EncodeLevelData(const std::vector<CompressionValue<T>>& level_byte_values) const
{
    cmc_debug_msg("The encoding of the CompressionValues after the parallel multi-resolution extraction iteration starts...");
    
    /* Get the rank of the mpi process within the communicator */
    const MPI_Comm comm = this->GetMPIComm();
    int rank{0};
    int ret_val = MPI_Comm_rank(comm, &rank);
    MPICheckError(ret_val);

    /* Define the root rank */
    const int root_rank = 0;

    /* Get a view on whether there is a residual for this element */
    bit_map::BitMapView is_residual_present(residual_presence_indications_);

    /* Get a view on the residual indciation flags */
    bit_map::BitMapView residual_indications(resdiual_order_indications_);

    /* Collect the global symbol frequency table */
    const auto sym_freq_table = CollectGlobalSymbolFrequenciesForEntropyCoding(level_byte_values, is_residual_present, residual_indications);

    /* Setup the Huffman coder */
    entropy_coding::huffman::HuffmanCoder<uint32_t> huff_coder(sym_freq_table);

    /* The encoded data will be stored in a BitVector */
    cmc::bit_vector::BitVector encoding;
    encoding.Reserve(3 * level_byte_values.size());

    /****** Encode the entropy codes interleaved with the encoded significant bits ******/

    /* Iterate over all values and encode them */
    for (auto val_iter = level_byte_values.begin(); val_iter != level_byte_values.end(); ++val_iter)
    {
        /* If there is no family that could be coarsened, we do not need to encode a residual for the element since the value remains unchanged */
        if (is_residual_present.GetNextBit() == false)
        {
            continue;
        }

        /* Get the current value */
        CompressionValue<T> val = *val_iter;

        /* Get the LZC */
        const uint32_t first_one_bit = val.GetNumberLeadingZeros();

        /* Get the info whether the capproximation was greater or smaller */
        bool next_residual_indication = residual_indications.GetNextBit();

        /* Check whether the byte value has been fully extracted */
        if (val.GetNumberLeadingZeros() == sizeof(T) * bit_map::kCharBit)
        {
            /* We only store a single zero residual for plus and minus residual operations */
            next_residual_indication = false;
        }

        /*** Enctropy Coding ***/
        /* We encode the entropy code with the Huffman coder */
        const auto symbol = ConvertToSymbolInFrequencyTable<T>(next_residual_indication, first_one_bit);

        /* Encode the symbol */
        const auto [huff_code, code_num_bits] = huff_coder.EncodeSymbol(symbol);

        /* Append the entropy code to the encoded stream */
        encoding.AppendBits(huff_code, code_num_bits);
        /******/

        /*** Significant Bits Coding ***/
        /* The first one is implicitly given by the leading zero count; therefore we set the "front bit" in order to discard the leading zeros and the following one */
        val.SetFrontBit(first_one_bit + 1);
        
        /* If there are remaining bits in the residual, append them to the encoded residuals of this level */
        if (not val.IsEmpty())
        {
            encoding.AppendBits(val.GetSignificantBitsInBigEndianOrdering(), val.GetCountOfSignificantBits());
        }
        /******/
    }

    /* At last, we encode the process boundary */
    const auto [proc_end_code, proc_end_code_num_bits] = huff_coder.EncodeSymbol(GetProcessEndSymbol<T>());
    encoding.AppendBits(proc_end_code, proc_end_code_num_bits);

    /* Set up the BitVector holding the encoded data for further use */
    encoding.TrimToContent();

    /* Count the additional process boundary symbol */
    ++num_levelwise_entropy_codes_;
    /****** END of Encode the entropy codes interleaved with the encoded significant bits ******/

    /* Store the encoded data */
    num_local_entropy_bytes_encoded_level_data_ = static_cast<uint32_t>(encoding.size());

    /* The global count of bytes for this level needs to be gathered */
    const uint64_t local_bytes = encoding.size();
    uint64_t global_bytes{0};
    const int rv_reduce_bytes = MPI_Reduce(&local_bytes, &global_bytes, 1, MPI_UINT64_T, MPI_SUM, root_rank, comm);
    MPICheckError(rv_reduce_bytes);

    /* Allocate an output vector */
    std::vector<uint8_t> encoded_data;

    /* The root rank writes a level header which stores for example the huffmann tree */
    if (rank == root_rank)
    {
        /* Serialize the Huffmann Coder Symbol Frequency Table */
        std::vector<uint8_t> sym_freq_table_serialized = huff_coder.GetSerializedSymbolFrequencyTable();
        const uint64_t sym_freq_table_num_bytes = sym_freq_table_serialized.size();
        
        /* Allcoate memory for the root rank's encoded stream */
        encoded_data.reserve(local_bytes + 3 * sizeof(uint64_t) + sym_freq_table_num_bytes);

        /* Compute the global umber of bytes describing this level */
        const uint64_t level_global_num_bytes = global_bytes + sym_freq_table_num_bytes + 3 * sizeof(uint64_t);

        /* Store the global level bytes count */
        PushBackValueToByteStream<uint64_t>(encoded_data, level_global_num_bytes);

        /* Store the global encoding number of bytes */
        PushBackValueToByteStream<uint64_t>(encoded_data, global_bytes);

        /* Store the amount of bytes for the serialized Huffmann symbol frequency table */
        PushBackValueToByteStream<uint64_t>(encoded_data, sym_freq_table_num_bytes);

        /* Next we store the serialized entropy table */
        std::copy_n(sym_freq_table_serialized.begin(), sym_freq_table_num_bytes, std::back_inserter(encoded_data));

        /* Finally, we store the root rank's local encoding */
        std::copy_n(encoding.begin(), local_bytes, std::back_inserter(encoded_data));
    } else
    {
        /* In case the process is not the root rank, we only store the local encoding without additional information */
        encoding.MoveDataInto(encoded_data);
    }

    cmc_debug_msg("The encoding of this level of the parallel multi-resolution extraction compression has been completed.");

    return encoded_data;
}

/**
 * We store a serialized partition table of the mesh and the entropy codes on the root rank, the other ranks obtain an empty vector
 */
template <typename T>
std::vector<uint8_t>
MultiResAdaptData<T>::StorePartitionTableOnTheRootRank(const t8_forest_t coarsened_forest) const
{
    /* Get the rank and size of the mpi process within the communicator */
    const MPI_Comm comm = this->GetMPIComm();
    int rank{0};
    int rv_rank = MPI_Comm_rank(comm, &rank);
    MPICheckError(rv_rank);
    int size{1};
    int rv_size = MPI_Comm_size(comm, &size);
    MPICheckError(rv_size);

    /* Get the global offset of the processes */
    std::array<uint32_t, 3> offset_values;

    /* Store the number of local elements */
    offset_values[0] = static_cast<uint32_t>(t8_forest_get_local_num_leaf_elements(coarsened_forest));

    /* Store the entropy offset */
    offset_values[1] = num_levelwise_entropy_codes_;

    /* Store the amount of bytes */
    offset_values[2] = num_local_entropy_bytes_encoded_level_data_;

    /* Declare an output vector for the gathering */
    std::vector<uint32_t> gathered_values;

    /* Allocate the vector on the root rank */
    if (rank == kRootRank)
    {
        gathered_values = std::vector<uint32_t>(3 * size);
    }

    /* Gather all the process local data */
    const int rv_gather = MPI_Gather(offset_values.data(), 3, MPI_UINT32_T, gathered_values.data(), 3, MPI_UINT32_T, kRootRank, comm);
    MPICheckError(rv_gather);

    std::vector<uint8_t> serialized_partition_table;
    if (rank == kRootRank)
    {
        serialized_partition_table.reserve(3 * size * sizeof(uint32_t));

        /* We need to exclusively scan the data to indicate the serialized positions */
        uint32_t num_elems{0};
        uint32_t num_entropy_codes{0};
        uint32_t num_encoded_bytes{0};

        /* Store the default start offset */
        PushBackValueToByteStream<uint32_t>(serialized_partition_table, num_elems);
        PushBackValueToByteStream<uint32_t>(serialized_partition_table, num_entropy_codes);
        PushBackValueToByteStream<uint32_t>(serialized_partition_table, num_encoded_bytes);

        /* Sum up the offsets for the next processes */
        for (int rank_id{0}; rank_id < size - 1; ++rank_id)
        {
            /* Compute the next offsets */
            num_elems += gathered_values[rank_id * 3];
            num_entropy_codes += gathered_values[rank_id * 3 + 1];
            num_encoded_bytes += gathered_values[rank_id * 3 + 2];

            /* Store the offsets */
            PushBackValueToByteStream<uint32_t>(serialized_partition_table, num_elems);
            PushBackValueToByteStream<uint32_t>(serialized_partition_table, num_entropy_codes);
            PushBackValueToByteStream<uint32_t>(serialized_partition_table, num_encoded_bytes);
        }
    }

    return serialized_partition_table;
}

template <typename T>
std::vector<uint8_t>
MultiResAdaptData<T>::StoreRootLevelPartitionTableOnTheRootRank(const t8_forest_t coarsened_forest) const
{
/* Get the rank and size of the mpi process within the communicator */
    const MPI_Comm comm = this->GetMPIComm();
    int rank{0};
    int rv_rank = MPI_Comm_rank(comm, &rank);
    MPICheckError(rv_rank);
    int size{1};
    int rv_size = MPI_Comm_size(comm, &size);
    MPICheckError(rv_size);

    /* Store the number of local elements */
    const uint32_t local_elem_count = static_cast<uint32_t>(t8_forest_get_local_num_leaf_elements(coarsened_forest));

    /* Declare an output vector for the gathering */
    std::vector<uint32_t> gathered_values;

    /* Allocate the vector on the root rank */
    if (rank == kRootRank)
    {
        gathered_values = std::vector<uint32_t>(size);
    }

    /* Gather all the process local data */
    const int rv_gather = MPI_Gather(&local_elem_count, 1, MPI_UINT32_T, gathered_values.data(), 1, MPI_UINT32_T, kRootRank, comm);
    MPICheckError(rv_gather);

    std::vector<uint8_t> serialized_partition_table;
    if (rank == kRootRank)
    {
        serialized_partition_table.reserve(3 * size * sizeof(uint32_t));

        /* We need to exclusively scan the data to indicate the serialized positions */
        uint32_t num_elems{0};

        /* Store the default start offset */
        PushBackValueToByteStream<uint32_t>(serialized_partition_table, num_elems);

        /* Sum up the offsets for the next processes */
        for (int rank_id{0}; rank_id < size - 1; ++rank_id)
        {
            /* Compute the next offsets */
            num_elems += gathered_values[rank_id];

            /* Store the offsets */
            PushBackValueToByteStream<uint32_t>(serialized_partition_table, num_elems);
        }
    }

    return serialized_partition_table;
}

template <typename T>
std::vector<uint8_t>
MultiResAdaptData<T>::EncodeRootLevelData(const std::vector<CompressionValue<T>>& root_level_values) const
{
    cmc_debug_msg("The encoding of the root level values of the parallel multi-resolution compression starts.");

    /* Get the rank of the mpi process within the communicator */
    const MPI_Comm comm = this->GetMPIComm();
    int rank{0};
    int ret_val = MPI_Comm_rank(comm, &rank);
    MPICheckError(ret_val);

    /* Define the root rank */
    const int root_rank = 0;

    /* Get the number of global root level values */
    const uint64_t num_local_root_values = root_level_values.size();
    uint64_t num_global_root_values{0};

    /* Get the global number of root values */
    const int rv_reduce_num_vals = MPI_Reduce(&num_local_root_values, &num_global_root_values, 1, MPI_UINT64_T, MPI_SUM, root_rank, comm);
    MPICheckError(rv_reduce_num_vals);

    /* Declare an output vector */
    std::vector<uint8_t> encoded_stream;

    if (rank == root_rank)
    {
        /* Allocate an additional uint64_t for storing the number of root values */
        encoded_stream.reserve(num_local_root_values * sizeof(T) + sizeof(uint64_t));
        PushBackValueToByteStream(encoded_stream, num_global_root_values);
    } else
    {
        /* Allocate memory for all local root values */
        encoded_stream.reserve(num_local_root_values * sizeof(T));
    }

    /* Store the root values in a not encoded fashion */
    for (auto val_iter = root_level_values.begin(); val_iter != root_level_values.end(); ++val_iter)
    {
        const T val = val_iter->template ReinterpretDataAs<T>();
        PushBackValueToByteStream(encoded_stream, val);
    }

    cmc_debug_msg("The encoding of the root level of the parllel multi-resolution compression has been finished.");

    return encoded_stream;
}


template <typename T>
inline ICompressionAdaptData<T>*
CreateMultiResExtractionAdaptationClass(AbstractByteCompressionVariable<T>* abstract_var)
{
    return new MultiResAdaptData<T>(abstract_var);
}

template <typename T>
inline void
DestroyMultiResExtractionAdaptationClass(ICompressionAdaptData<T>* iadapt_data)
{
    delete iadapt_data;
}

template<class T>
class CompressionVariable : public AbstractByteCompressionVariable<T>
{
public:
    CompressionVariable() = delete;

    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, const std::vector<T>& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_leaf_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(variable_data);
        StoreMeshMPIComm();
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreateMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::MeshEncoder>();  
    };

    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, const std::vector<CompressionValue<T>>& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_leaf_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(variable_data);
        StoreMeshMPIComm();
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreateMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::MeshEncoder>();  
    };

    CompressionVariable(const std::string& name, t8_forest_t initial_mesh, std::vector<CompressionValue<T>>&& variable_data)
    : AbstractByteCompressionVariable<T>()
    {
        if (static_cast<size_t>(t8_forest_get_local_num_leaf_elements(initial_mesh)) != variable_data.size())
        {
            throw std::invalid_argument("The number of local mesh elements does not match the amount of data elements.");
        }

        this->SetName(name);
        this->SetAmrMesh(AmrMesh(initial_mesh));
        this->SetData(std::move(variable_data));
        StoreMeshMPIComm();
        AbstractByteCompressionVariable<T>::adaptation_creator_ = CreateMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::adaptation_destructor_ = DestroyMultiResExtractionAdaptationClass<T>;
        AbstractByteCompressionVariable<T>::mesh_encoder_ = std::make_unique<mesh_compression::MeshEncoder>();  
    };

    CompressionSchema GetCompressionSchema() const override
    {
        return CompressionSchema::ParallelMultiResExtraction;
    }

private:
    void StoreMeshMPIComm(){this->SetMPIComm(t8_forest_get_mpicomm(this->GetAmrMesh().GetMesh()));};

};

}

#endif /* !CMC_MULTI_RES_PAR_EXTRACTION_COMPRESSION_HXX */
