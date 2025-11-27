#ifndef CMC_IFACE_ABSTRACT_MESH_ENCODER_HXX
#define CMC_IFACE_ABSTRACT_MESH_ENCODER_HXX

#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_serialization.hxx"
#include "t8code/cmc_t8_mesh.hxx"
#ifdef CMC_ENABLE_MPI
#include "mpi/cmc_mpi.hxx"
#endif

#include <vector>
#include <utility>

namespace cmc::mesh_compression
{

/* Forward declaration */
struct ExchangeData;

class IAbstractMeshEncoder
{
public:

    void IntializeCompressionIteration(const t8_locidx_t num_local_elems);
    void FinalizeCompressionIteration();
    void IndicateCoarsening();
    void IndicateElementStaysUnchanged();
    int GetNumberOfCompressionIterations() const {return compression_iteration_count_;};

    std::vector<uint8_t> GetEncodedLevelData();

#ifdef CMC_ENABLE_MPI
    std::vector<uint8_t> GetPartitionedEncodedLevelData(t8_forest_t adapted_mesh, t8_forest_t partitioned_mesh, MPI_Comm comm);
#endif

    virtual ~IAbstractMeshEncoder(){};

protected:
    IAbstractMeshEncoder() = default;
    
private:
#ifdef CMC_ENABLE_MPI
    std::vector<ExchangeData> DetermineSendMessages(const uint64_t num_global_bits, const uint64_t global_bit_offset, const int comm_size);
    std::pair<uint64_t, uint64_t> DetermineSendingRange(const uint64_t current_local_bit_idx, const uint64_t current_global_bit_idx, const uint64_t num_global_byte_packages, const int comm_size);
#endif

    bit_map::BitMap encoded_level_data_;
    int compression_iteration_count_{0};
};

inline void
IAbstractMeshEncoder::IntializeCompressionIteration(const t8_locidx_t num_local_elems)
{
    encoded_level_data_ = bit_map::BitMap();
    encoded_level_data_.Reserve(num_local_elems);
}

inline void
IAbstractMeshEncoder::IndicateCoarsening()
{
    encoded_level_data_.AppendSetBit();
}

inline void
IAbstractMeshEncoder::IndicateElementStaysUnchanged()
{
    encoded_level_data_.AppendUnsetBit();
}

inline void
IAbstractMeshEncoder::FinalizeCompressionIteration()
{
    ++compression_iteration_count_;
}


inline std::vector<uint8_t>
IAbstractMeshEncoder::GetEncodedLevelData()
{
    /* Get the vectorized data of the bitmap and return it */
    std::vector<uint8_t> encoded_mesh;

    encoded_level_data_.MoveDataInto(encoded_mesh);

    return encoded_mesh;
}

#ifdef CMC_ENABLE_MPI

struct ExchangeData
{
    ExchangeData(const int rank_, const int tag_, const bit_map::BitMap& data_)
    : rank{rank_}, tag{tag_}, data{data_} {};

    ExchangeData(const int rank_, const int tag_, bit_map::BitMap&& data_)
    : rank{rank_}, tag{tag_}, data{std::move(data_)} {};
    
    ExchangeData(const int rank_, const int tag_, const size_t num_bits)
    : rank{rank_}, tag{tag_}, data(num_bits) {};

    int rank{-1};
    int tag{0};
    bit_map::BitMap data;
};

inline int
GetReceivingRank(const uint64_t current_global_bit_idx, const uint64_t num_global_byte_packages, const int comm_size)
{
    /* Determine global byte package holding the bit */
    const uint64_t global_byte_package_id = current_global_bit_idx / bit_map::kCharBit;

    /* Determine the rank which will hold this byte package */
    const int rank = static_cast<int>(((double) global_byte_package_id * (double) comm_size) / (long double) num_global_byte_packages);

    return rank;
}

inline std::pair<uint64_t, uint64_t>
IAbstractMeshEncoder::DetermineSendingRange(const uint64_t current_local_bit_idx, const uint64_t current_global_bit_idx, const uint64_t num_global_byte_packages, const int comm_size)
{
    /* Get the rank which will receive the current bit */
    const int rank = GetReceivingRank(current_global_bit_idx, num_global_byte_packages, comm_size);
    cmc_debug_msg("Receiving rank is: ", rank);

    /* Compute the first global bit of the next rank */
    const uint64_t computed_end_byte = static_cast<uint64_t>(((double) (rank + 1) * (long double) num_global_byte_packages) / (double) comm_size);
    const uint64_t end_byte = (computed_end_byte > 0 ? computed_end_byte : 1);
    const uint64_t rank_end_bit = static_cast<uint64_t>(bit_map::kCharBit) * end_byte;

    cmc_debug_msg("rank_end_bit: ", rank_end_bit, ", current_global_bit_idx: ", current_global_bit_idx);
    cmc_assert(rank_end_bit >= current_global_bit_idx);

    /* Determine the number of contiguous bits to send to this process alongside */
    const uint64_t range = (rank_end_bit - current_global_bit_idx) <= (encoded_level_data_.size() - current_local_bit_idx) ? (rank_end_bit - current_global_bit_idx) : (encoded_level_data_.size() - current_local_bit_idx);
    cmc_debug_msg("Computed range is: ", range);
    return std::make_pair(rank, range);
}

std::vector<ExchangeData>
IAbstractMeshEncoder::DetermineSendMessages(const uint64_t num_global_bits, const uint64_t global_bit_offset, const int comm_size)
{
    /* Determine the number of global byte packages */
    const uint64_t num_global_byte_packages = (num_global_bits / bit_map::kCharBit) + (num_global_bits % bit_map::kCharBit != 0 ? 1 : 0);

    /* Get the number of encoded bits */
    const uint64_t num_local_encoded_bits = static_cast<uint64_t>(encoded_level_data_.size());

    /* Get a view on the local encoding bits */
    bit_map::BitMapView indication_bits_view(encoded_level_data_);

    std::vector<ExchangeData> send_data;

    uint64_t processed_bits{0};
    uint64_t& current_local_bit_idx = processed_bits;

    /* Iterate over all local bits and determine their rank-affiliation */
    while (processed_bits < num_local_encoded_bits)
    {
        cmc_debug_msg("Calling: DetermineSendingRange(): with: ", current_local_bit_idx,", ",global_bit_offset + current_local_bit_idx,", ", num_global_byte_packages,", ", comm_size);
        /* Get the rank the current bit will be sent to and the amount of contiguous bits to send */
        auto [recv_rank, num_bits] = DetermineSendingRange(current_local_bit_idx, global_bit_offset + current_local_bit_idx, num_global_byte_packages, comm_size);

        /* Get this mnumber of bits */
        std::vector<uint8_t> bits = indication_bits_view.GetNextNumberOfBits(num_bits);

        /* Emplace a send data struct with the computed properties */
        send_data.emplace_back(recv_rank, num_bits, bit_map::BitMap(bits, num_bits));

        processed_bits += num_bits;
    }

    return send_data;
}


inline bit_map::BitMap
ProcessMessages(std::vector<ExchangeData>& recv_messages)
{
    /* After all messages have been transferred and received, we sort the messages based on the source rank */
    std::sort(recv_messages.begin(), recv_messages.end(), [](const ExchangeData& a, const ExchangeData& b){
        return a.rank < b.rank;
    });

    /* Count the bits */
    size_t bit_count{0};
    for (auto msg_iter = recv_messages.begin(); msg_iter != recv_messages.end(); ++msg_iter)
    {
        bit_count += static_cast<size_t>(msg_iter->tag);
    }

    /* We stack the bits from the different ranks in a contiguous BitMap */
    bit_map::BitMap local_indication_bits;
    local_indication_bits.Reserve(bit_count);

    /* Add the bits contiguously */
    for (auto msg_iter = recv_messages.begin(); msg_iter != recv_messages.end(); ++msg_iter)
    {
        local_indication_bits.AppendBits(msg_iter->data);
    }

    return local_indication_bits;
}


inline std::vector<uint8_t>
IAbstractMeshEncoder::GetPartitionedEncodedLevelData(t8_forest_t adapted_mesh, t8_forest_t partitioned_mesh, MPI_Comm comm)
{
    /* Get the rank of the mpi process within the communicator */
    int comm_size{0};
    int ret_val = MPI_Comm_size(comm, &comm_size);
    MPICheckError(ret_val);

    int rank;
    ret_val = MPI_Comm_rank(comm, &rank);
    MPICheckError(ret_val);

    /* Define the root rank */
    const int root_rank = 0;

    /* Get the number of global bits/elements */
    const uint64_t num_global_bits = t8_forest_get_global_num_leaf_elements (adapted_mesh);

    /* The number of local bits should coincide the number of local elements in the adapted forest */
    cmc_assert(static_cast<uint64_t>(encoded_level_data_.size()) == static_cast<uint64_t>(t8_forest_get_local_num_leaf_elements(adapted_mesh)));
    cmc_assert(t8_forest_get_global_num_leaf_elements (adapted_mesh) == t8_forest_get_global_num_leaf_elements (partitioned_mesh));

    /* Get the process-local global element offset */
    const uint64_t global_bit_offset = static_cast<uint64_t>(t8_forest_get_first_local_leaf_element_id (adapted_mesh)); 

    /* Gather all messages that need to be sent */
    const std::vector<ExchangeData> send_messages = DetermineSendMessages(num_global_bits, global_bit_offset, comm_size);

    /* Allocate MPI Requests for the messages */
    std::vector<MPI_Request> send_requests(send_messages.size());

    /* Stage all messages */
    int send_idx{0};
    for (auto msg_iter = send_messages.begin(); msg_iter != send_messages.end(); ++msg_iter, ++send_idx)
    {
        /* We only send a message, if it contains bits */
        if (msg_iter->tag != 0)
        {
            /* Get the number of bytes to be send within this message */
            const int count_bytes = static_cast<int>(msg_iter->data.size_bytes());

            /* Actually send the message */
            ret_val = MPI_Isend(msg_iter->data.data(), count_bytes, MPI_UINT8_T, msg_iter->rank, msg_iter->tag, comm, &send_requests[send_idx]);
            MPICheckError(ret_val);
        }
    }

    /* Wait until all messages are staged */
    MPI_Request barrier_req;
    ret_val = MPI_Ibarrier(comm, &barrier_req);
    MPICheckError(ret_val);

    std::vector<ExchangeData> recv_messages;
    recv_messages.reserve(send_messages.size());

    bool continue_receiving = true;

    /* Start receiving the messages */
    while(continue_receiving)
    {
        /* Check if the barrier has been reached by all processes */
        int is_barrier_complete;
        const int ret_val_barrier_check = MPI_Test(&barrier_req, &is_barrier_complete, MPI_STATUS_IGNORE);
        MPICheckError(ret_val_barrier_check);

        /* Check if there is a message to receive */
        int is_message_waiting{0};
        MPI_Status status;
        const int probe_ret_val = MPI_Iprobe(MPI_ANY_SOURCE, MPI_ANY_TAG, comm, &is_message_waiting, &status);
        MPICheckError(probe_ret_val);

        /* If there is a message, we receive it */
        if (is_message_waiting)
        {
            /* Get the count of symbol-frequencies */
            int count{0};
            const int ret_val_count = MPI_Get_count(&status, MPI_UINT8_T, &count);
            MPICheckError(ret_val_count);

            /* Get the number of bits */
            const int num_bits = status.MPI_TAG;
            const int source_rank = status.MPI_SOURCE;

            cmc_assert(static_cast<size_t>(count) == (num_bits / bit_map::kCharBit) + (num_bits % bit_map::kCharBit != 0 ? 1 : 0));

            /* Allocate a struct for receiving the message */
            recv_messages.emplace_back(source_rank, num_bits, static_cast<size_t>(num_bits));

            /* Receive the actual message */
            const int ret_val_recv = MPI_Recv(recv_messages.back().data.data_overwrite(), count, MPI_UINT8_T, source_rank, num_bits, comm, MPI_STATUS_IGNORE);
            MPICheckError(ret_val_recv);
        }

        /* Update the loop flag, we continue until no more messages are queued */
        continue_receiving = is_message_waiting || !(is_barrier_complete);
    }

    /* Process the messages to build a contiguous BitMap */
    bit_map::BitMap contiguous_bits = ProcessMessages(recv_messages);

    if (rank != root_rank)
    {
        /* Get the vectorized data of the bitmap and return it */
        std::vector<uint8_t> encoded_mesh;
        if (not contiguous_bits.IsEmpty())
        {
            contiguous_bits.MoveDataInto(encoded_mesh);
        }

        return encoded_mesh;
    } else
    {
        /* In case of the root rank, we add a small preamble indicating the number of bits on this level */
        std::vector<uint8_t> encoded_stream;
        encoded_stream.reserve(contiguous_bits.size_bytes() + sizeof(uint64_t));
        cmc_debug_msg("contiguous_bits.size_bytes(): ", contiguous_bits.size_bytes());
        /* Store the number of global bits */
        PushBackValueToByteStream<uint64_t>(encoded_stream, num_global_bits);

        /* Afterwards the encoded mesh data is appended */
        if (not contiguous_bits.IsEmpty())
        {
            std::copy_n(contiguous_bits.begin_bytes(), contiguous_bits.size_bytes(), std::back_inserter(encoded_stream));
        }
        
        return encoded_stream;
    }
}

#endif

}


#endif /* !CMC_IFACE_ABSTRACT_MESH_ENCODER_HXX */
