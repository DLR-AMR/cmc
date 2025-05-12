#ifndef CMC_MPI_DATA_HXX
#define CMC_MPI_DATA_HXX

#include "utilities/cmc_utilities.hxx"
#include "mpi/cmc_mpi.hxx"
#include "utilities/cmc_log_functions.hxx"

namespace cmc
{

constexpr int kLowestPermittedID = 0;
constexpr int kHighestPermittedID = 1023;
constexpr int kLowestNotPermittedID = kHighestPermittedID + 1;
constexpr int kSubandIDOffset = 1024;
constexpr int kSubbandIDFactor = 1024;
constexpr int kMortonIDOffset = 1024 * 1024 + kSubandIDOffset;

constexpr inline
int
CreateSubbandID(const int id, const int subband_id)
{
    cmc_assert(kLowestPermittedID <= id && id < kLowestNotPermittedID);
    cmc_assert(kLowestPermittedID <= subband_id && subband_id < kLowestNotPermittedID);

    return id * kSubbandIDFactor + subband_id;
}

constexpr inline
int
CreateMortonIndicesTag(const int id)
{
    return id + kMortonIDOffset;
}

constexpr inline
int
CreateDataTag(const int id)
{
    return id;
}

constexpr inline 
int
GetIDFromTag(const int msg_tag)
{
    if (msg_tag >= kMortonIDOffset)
    {
        return msg_tag - kMortonIDOffset;
    } else
    {
        return msg_tag;
    }
}

constexpr inline 
int
GetActualVariableIDFromTag(const int msg_tag)
{
    /* Check if the message corresponded to a message containing Morton Indices */
    if (msg_tag >= kMortonIDOffset)
    {
        const int var_comm_id = msg_tag - kMortonIDOffset;

        /* Check if it is a subband variable */
        if (var_comm_id >= kLowestNotPermittedID)
        {
            /* Thsi integer division erases the subband information which was added for the communication and we obtain the initial id */
            return static_cast<int>(var_comm_id / kSubbandIDFactor);
        } else
        {
            return var_comm_id;
        }
    } else
    {
        /* Check if it is a subband variable */
        if (msg_tag >= kLowestNotPermittedID)
        {
            /* Thsi integer division erases the subband information which was added for the communication and we obtain the initial id */
            return static_cast<int>(msg_tag / kSubbandIDFactor);
        } else
        {
            return msg_tag;
        }
    }
}

constexpr inline
bool IsMessageADataMessage(const int msg_tag)
{
    return (msg_tag < kMortonIDOffset);
}

constexpr inline
bool IsMessageAMortonIndicesMessage(const int msg_tag)
{
    return (msg_tag >= kMortonIDOffset);
}

constexpr int MortonIndicesTagOffset = 1024;

MPI_Datatype ConvertCmcTypeToMPIType(const CmcType type);

template<typename T>
constexpr MPI_Datatype ConvertToMPIType();


class DataOffsets
{
    using iterator = std::vector<MortonIndex>::iterator;
    using const_iterator = std::vector<MortonIndex>::const_iterator;
public:

    DataOffsets() = delete;
    DataOffsets(const size_t num_elements)
    : offsets_(num_elements){}

    DataOffsets(const DataOffsets& other) = default;
    DataOffsets& operator=(const DataOffsets& other) = default;
    DataOffsets(DataOffsets&& other) = default;
    DataOffsets& operator=(DataOffsets&& other) = default;

    template<typename U> auto operator[](U index) const -> std::enable_if_t<std::is_integral_v<U>, MortonIndex> {return offsets_[index];}
    template<typename U> auto operator[](U index) -> std::enable_if_t<std::is_integral_v<U>, MortonIndex&> {return offsets_[index];}

    iterator Begin() { return offsets_.begin(); }
    iterator End() { return offsets_.end(); }
    const_iterator Begin() const { return offsets_.begin(); }
    const_iterator End() const { return offsets_.end(); }
    const_iterator CBegin() const { return offsets_.cbegin(); }
    const_iterator CEnd() const { return offsets_.cend(); }

    size_t size() const {return offsets_.size();}

    MortonIndex* data()
    {
        return offsets_.data();
    }

    void SetMPIComm(const MPI_Comm comm) {comm_ = comm;}

private:
    std::vector<MortonIndex> offsets_;
    MPI_Comm comm_{MPI_COMM_WORLD};
};


constexpr int kVarMessageMPIErr = -1;

template<class T>
class VariableMessage
{
public:
    VariableMessage() = default;
    VariableMessage(const int rank, const int variable_id)
    : rank_{rank}, variable_id_{variable_id}{}
    VariableMessage(const int rank, const int variable_id, const int num_elements)
    : rank_{rank}, variable_id_{variable_id}, data_(num_elements), morton_indices_(num_elements){}

    VariableMessage(const VariableMessage& other) = default;
    VariableMessage& operator=(const VariableMessage& other) = default;
    VariableMessage(VariableMessage&& other) = default;
    VariableMessage& operator=(VariableMessage&& other) = default;

    ~VariableMessage(){
        cmc_debug_msg("VariableMessage of rank ", rank_, " for var ", variable_id_, " with data size = ", data_.size(), " and morton indices size ", morton_indices_.size(), " is destroyed!");
    }
    void* GetInitialDataPtr() {return static_cast<void*>(data_.data());}
    void* GetInitialMortonIndicesPtr() {return static_cast<void*>(morton_indices_.data());}

    int GetRank() const {return rank_;}
    int GetVariableID() const {return variable_id_;}

    std::vector<T>::iterator DataBegin() {return data_.begin();}
    std::vector<T>::const_iterator DataBegin() const {return data_.begin();}
    std::vector<T>::iterator DataEnd() {return data_.end();}
    std::vector<T>::const_iterator DataEnd() const {return data_.end();}

    std::vector<MortonIndex>::iterator MortonIndicesBegin() {return morton_indices_.begin();}
    std::vector<MortonIndex>::const_iterator MortonIndicesBegin() const {return morton_indices_.begin();}
    std::vector<MortonIndex>::iterator MortonIndicesEnd() {return morton_indices_.end();}
    std::vector<MortonIndex>::const_iterator MortonIndicesEnd() const {return morton_indices_.end();}

    int rank_{kVarMessageMPIErr};
    int variable_id_{kVarMessageMPIErr};
    std::vector<T> data_;
    std::vector<MortonIndex> morton_indices_;
};

#if 0
inline
int
CreateMortonIndicesTag(const int variable_id)
{
    return variable_id + MortonIndicesTagOffset;
}

inline
int
CreateVariableDataTag(const int variable_id) 
{
    return variable_id;
}

inline 
int
GetVariableIDFromTag(const int msg_tag)
{
    return (msg_tag >= MortonIndicesTagOffset ? msg_tag - MortonIndicesTagOffset : msg_tag);
}

inline
bool IsMessageADataMessage(const int msg_tag)
{
    return (msg_tag < MortonIndicesTagOffset);
}
#endif

template<typename T>
static 
void
DetachVarMessage(VariableMessage<T>& source, VariableMessage<T>& destination)
{
    destination.rank_ = source.rank_;
    destination.variable_id_ = source.variable_id_;
    std::swap(source.data_, destination.data_);
    std::swap(source.morton_indices_, destination.morton_indices_);
}

using VarMessage = std::variant<VariableMessage<int8_t>, VariableMessage<char>, VariableMessage<int16_t>, VariableMessage<int32_t>, VariableMessage<float>, VariableMessage<double>,
                                 VariableMessage<uint8_t>, VariableMessage<uint16_t>, VariableMessage<uint32_t>, VariableMessage<int64_t>, VariableMessage<uint64_t>>;

VarMessage
CreateVariableMessage(const int rank, const int variable_id, const CmcType type, const int num_elements);

class VariableSendMessage
{
public:
    template<typename T> VariableSendMessage(VariableMessage<T>&& var_message)
    : type_{ConvertToCmcType<T>()}, message_{std::move(var_message)}{}

    VariableSendMessage(const VariableSendMessage& other) = default;
    VariableSendMessage& operator=(const VariableSendMessage& other) = default;
    VariableSendMessage(VariableSendMessage&& other) = default;
    VariableSendMessage& operator=(VariableSendMessage&& other) = default;

    std::pair<MPI_Request, MPI_Request> Send(const MPI_Comm comm);
    VarMessage& GetInternalVariant() {return message_;}
    const VarMessage& GetInternalVariant() const {return message_;}

private:
    CmcType type_;
    VarMessage message_;
};

class
VariableRecvMessage
{
public:
    template<typename T> VariableRecvMessage(VariableMessage<T>&& var_message)
    : type_{ConvertToCmcType<T>()}, message_{std::move(var_message)}{}

    VariableRecvMessage(const CmcType type, VarMessage&& message)
    : type_{type}, message_{std::move(message)} {};

    VariableRecvMessage(const VariableRecvMessage& other) = default;
    VariableRecvMessage& operator=(const VariableRecvMessage& other) = default;
    VariableRecvMessage(VariableRecvMessage&& other) = default;
    VariableRecvMessage& operator=(VariableRecvMessage&& other) = default;

    void* GetInitialDataPtr();
    void* GetInitialMortonIndicesPtr();

    int GetSendingRank() const;
    int GetVariableID() const;

    VarMessage& GetInternalVariant() {return message_;}
    const VarMessage& GetInternalVariant() const {return message_;}
private:
    CmcType type_;
    VarMessage message_;
};

template<typename T>
constexpr MPI_Datatype
ConvertToMPIType()
{
    if constexpr (std::is_same_v<T, int8_t>)
    {
        return MPI_INT8_T;
    } else if constexpr (std::is_same_v<T, char>)
    {
        return MPI_CHAR;
    } else if constexpr (std::is_same_v<T, int16_t>)
    {   
        return MPI_INT16_T;
    } else if constexpr (std::is_same_v<T, int32_t>)
    {
        return MPI_INT32_T;
    } else if constexpr (std::is_same_v<T, float>)
    {
        return MPI_FLOAT;
    } else if constexpr (std::is_same_v<T, double>)
    {
        return MPI_DOUBLE;
    } else if constexpr (std::is_same_v<T, uint8_t>)
    {
        return MPI_UINT8_T;
    } else if constexpr (std::is_same_v<T, uint16_t>)
    {
        return MPI_UINT16_T;
    } else if constexpr (std::is_same_v<T, uint32_t>)
    {
        return MPI_UINT32_T;
    } else if constexpr (std::is_same_v<T, int64_t>)
    {
        return MPI_INT64_T;
    } else if constexpr (std::is_same_v<T, uint64_t>)
    {
        return MPI_UINT64_T;
    } else
    {
        cmc_err_msg("There is no supported conversion to an MPI_Datatype.");
        return MPI_DATATYPE_NULL;
    }
}

}

#endif /* !CMC_MPI_DATA_HXX */
