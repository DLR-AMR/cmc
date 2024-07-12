#include "mpi/cmc_mpi_data.hxx"

namespace cmc
{

MPI_Datatype
ConvertCmcTypeToMPIType(const CmcType type)
{
    switch (type)
    {
        case CmcType::Int8_t:
            return MPI_INT8_T;
        break;
        case CmcType::Char:
            return MPI_CHAR;
        break;
        case CmcType::Int16_t:
            return MPI_INT16_T;
        break;
        case CmcType::Int32_t:
            return MPI_INT32_T;
        break;
        case CmcType::Float:
            return MPI_FLOAT;
        break;
        case CmcType::Double:
            return MPI_DOUBLE;
        break;
        case CmcType::Uint8_t:
            return MPI_UINT8_T;
        break;
        case CmcType::Uint16_t:
            return MPI_UINT16_T;
        break;
        case CmcType::Uint32_t:
            return MPI_UINT32_T;
        break;
        case CmcType::Int64_t:
            return MPI_INT64_T;
        break;
        case CmcType::Uint64_t:
            return MPI_UINT64_T;
        break;
        default:
            return MPI_DATATYPE_NULL;
        
    }
}

VarMessage
CreateVariableMessage(const int rank, const int variable_id, const CmcType type, const int num_elements)
{
 switch (type)
    {
        case CmcType::Int8_t:
            return VariableMessage<int8_t>(rank, variable_id, num_elements);
        break;
        case CmcType::Char:
            return VariableMessage<char>(rank, variable_id, num_elements);;
        break;
        case CmcType::Int16_t:
            return VariableMessage<int16_t>(rank, variable_id, num_elements);;
        break;
        case CmcType::Int32_t:
            return VariableMessage<int32_t>(rank, variable_id, num_elements);;
        break;
        case CmcType::Float:
            return VariableMessage<float>(rank, variable_id, num_elements);;
        break;
        case CmcType::Double:
            return VariableMessage<double>(rank, variable_id, num_elements);;
        break;
        case CmcType::Uint8_t:
            return VariableMessage<uint8_t>(rank, variable_id, num_elements);;
        break;
        case CmcType::Uint16_t:
            return VariableMessage<uint16_t>(rank, variable_id, num_elements);;
        break;
        case CmcType::Uint32_t:
            return VariableMessage<uint32_t>(rank, variable_id, num_elements);;
        break;
        case CmcType::Int64_t:
            return VariableMessage<int64_t>(rank, variable_id, num_elements);;
        break;
        case CmcType::Uint64_t:
            return VariableMessage<uint64_t>(rank, variable_id, num_elements);;
        break;
        default:
            cmc_err_msg("A variable of the given type cannot be constructed.");
            return VariableMessage<int8_t>(-1, -1, 1);
        
    }   
}

struct SendMessage
{
public:
    SendMessage() = delete;
    SendMessage(const MPI_Comm comm)
    : comm_{comm} {};

    std::pair<MPI_Request, MPI_Request> operator()(const VariableMessage<int8_t>& msg) {
        MPI_Request morton_req, data_req;
        /* Send the Morton indices */
        int err = MPI_Isend(msg.morton_indices_.data(), msg.morton_indices_.size(), MPI_MORTON_INDEX_T, msg.rank_, CreateMortonIndicesTag(msg.variable_id_), comm_, &morton_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.morton_indices_.size(), " Morton indices to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        /* Send the actual data */
        err = MPI_Isend(msg.data_.data(), msg.data_.size(), ConvertToMPIType<int8_t>(), msg.rank_, CreateDataTag(msg.variable_id_), comm_, &data_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.data_.size(), " data values (of type: int8_t) to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        return std::make_pair(morton_req, data_req);
    }
    std::pair<MPI_Request, MPI_Request> operator()(const VariableMessage<char>& msg) {
        MPI_Request morton_req, data_req;
        /* Send the Morton indices */
        int err = MPI_Isend(msg.morton_indices_.data(), msg.morton_indices_.size(), MPI_MORTON_INDEX_T, msg.rank_, CreateMortonIndicesTag(msg.variable_id_), comm_, &morton_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.morton_indices_.size(), " Morton indices to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        /* Send the actual data */
        err = MPI_Isend(msg.data_.data(), msg.data_.size(), ConvertToMPIType<char>(), msg.rank_, CreateDataTag(msg.variable_id_), comm_, &data_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.data_.size(), " data values (of type: char) to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        return std::make_pair(morton_req, data_req);
    }
    std::pair<MPI_Request, MPI_Request> operator()(const VariableMessage<int16_t>& msg) {
        MPI_Request morton_req, data_req;
        /* Send the Morton indices */
        //cmc_debug_msg("Morton indices that going to be sent");
        //for (auto miter = msg.morton_indices_.begin(); miter != msg.morton_indices_.end(); ++miter)
        //{
        //    std::cout << *miter << ", ";
        //}
        int err = MPI_Isend(msg.morton_indices_.data(), msg.morton_indices_.size(), MPI_MORTON_INDEX_T, msg.rank_, CreateMortonIndicesTag(msg.variable_id_), comm_, &morton_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.morton_indices_.size(), " Morton indices to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        //cmc_debug_msg("Data that is going to be sent");
        //for (auto miter = msg.data_.begin(); miter != msg.data_.end(); ++miter)
        //{
        //    std::cout << *miter << ", ";
        //}
        /* Send the actual data */
        err = MPI_Isend(msg.data_.data(), msg.data_.size(), ConvertToMPIType<int16_t>(), msg.rank_, CreateDataTag(msg.variable_id_), comm_, &data_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.data_.size(), " data values (of type: int16_t) to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        return std::make_pair(morton_req, data_req);
    }
    std::pair<MPI_Request, MPI_Request> operator()(const VariableMessage<int32_t>& msg) {
        MPI_Request morton_req, data_req;
        /* Send the Morton indices */
        int err = MPI_Isend(msg.morton_indices_.data(), msg.morton_indices_.size(), MPI_MORTON_INDEX_T, msg.rank_, CreateMortonIndicesTag(msg.variable_id_), comm_, &morton_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.morton_indices_.size(), " Morton indices to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        /* Send the actual data */
        err = MPI_Isend(msg.data_.data(), msg.data_.size(), ConvertToMPIType<int32_t>(), msg.rank_, CreateDataTag(msg.variable_id_), comm_, &data_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.data_.size(), " data values (of type: int32_t) to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        return std::make_pair(morton_req, data_req);
    }
    std::pair<MPI_Request, MPI_Request> operator()(const VariableMessage<float>& msg) {
        MPI_Request morton_req, data_req;
        /* Send the Morton indices */
        int err = MPI_Isend(msg.morton_indices_.data(), msg.morton_indices_.size(), MPI_MORTON_INDEX_T, msg.rank_, CreateMortonIndicesTag(msg.variable_id_), comm_, &morton_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.morton_indices_.size(), " Morton indices to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        /* Send the actual data */
        err = MPI_Isend(msg.data_.data(), msg.data_.size(), ConvertToMPIType<float>(), msg.rank_, CreateDataTag(msg.variable_id_), comm_, &data_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.data_.size(), " data values (of type: float) to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        return std::make_pair(morton_req, data_req);
    }
    std::pair<MPI_Request, MPI_Request> operator()(const VariableMessage<double>& msg) {
        MPI_Request morton_req, data_req;
        /* Send the Morton indices */
        int err = MPI_Isend(msg.morton_indices_.data(), msg.morton_indices_.size(), MPI_MORTON_INDEX_T, msg.rank_, CreateMortonIndicesTag(msg.variable_id_), comm_, &morton_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.morton_indices_.size(), " Morton indices to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        /* Send the actual data */
        err = MPI_Isend(msg.data_.data(), msg.data_.size(), ConvertToMPIType<double>(), msg.rank_, CreateDataTag(msg.variable_id_), comm_, &data_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.data_.size(), " data values (of type: double) to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        return std::make_pair(morton_req, data_req);
    }
    std::pair<MPI_Request, MPI_Request> operator()(const VariableMessage<uint8_t>& msg) {
        MPI_Request morton_req, data_req;
        /* Send the Morton indices */
        int err = MPI_Isend(msg.morton_indices_.data(), msg.morton_indices_.size(), MPI_MORTON_INDEX_T, msg.rank_, CreateMortonIndicesTag(msg.variable_id_), comm_, &morton_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.morton_indices_.size(), " Morton indices to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        /* Send the actual data */
        err = MPI_Isend(msg.data_.data(), msg.data_.size(), ConvertToMPIType<uint8_t>(), msg.rank_, CreateDataTag(msg.variable_id_), comm_, &data_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.data_.size(), " data values (of type: uint8_t) to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        return std::make_pair(morton_req, data_req);
    }
    std::pair<MPI_Request, MPI_Request> operator()(const VariableMessage<uint16_t>& msg) {
        MPI_Request morton_req, data_req;
        /* Send the Morton indices */
        int err = MPI_Isend(msg.morton_indices_.data(), msg.morton_indices_.size(), MPI_MORTON_INDEX_T, msg.rank_, CreateMortonIndicesTag(msg.variable_id_), comm_, &morton_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.morton_indices_.size(), " Morton indices to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        /* Send the actual data */
        err = MPI_Isend(msg.data_.data(), msg.data_.size(), ConvertToMPIType<uint16_t>(), msg.rank_, CreateDataTag(msg.variable_id_), comm_, &data_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.data_.size(), " data values (of type: uint16_t) to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        return std::make_pair(morton_req, data_req);
    }
    std::pair<MPI_Request, MPI_Request> operator()(const VariableMessage<uint32_t>& msg) {
        MPI_Request morton_req, data_req;
        /* Send the Morton indices */
        int err = MPI_Isend(msg.morton_indices_.data(), msg.morton_indices_.size(), MPI_MORTON_INDEX_T, msg.rank_, CreateMortonIndicesTag(msg.variable_id_), comm_, &morton_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.morton_indices_.size(), " Morton indices to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        /* Send the actual data */
        err = MPI_Isend(msg.data_.data(), msg.data_.size(), ConvertToMPIType<uint32_t>(), msg.rank_, CreateDataTag(msg.variable_id_), comm_, &data_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.data_.size(), " data values (of type: uint32_t) to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        return std::make_pair(morton_req, data_req);
    }
    std::pair<MPI_Request, MPI_Request> operator()(const VariableMessage<int64_t>& msg) {
        MPI_Request morton_req, data_req;
        /* Send the Morton indices */
        int err = MPI_Isend(msg.morton_indices_.data(), msg.morton_indices_.size(), MPI_MORTON_INDEX_T, msg.rank_, CreateMortonIndicesTag(msg.variable_id_), comm_, &morton_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.morton_indices_.size(), " Morton indices to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        /* Send the actual data */
        err = MPI_Isend(msg.data_.data(), msg.data_.size(), ConvertToMPIType<int64_t>(), msg.rank_, CreateDataTag(msg.variable_id_), comm_, &data_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.data_.size(), " data values (of type: int64_t) to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        return std::make_pair(morton_req, data_req);
    }
    std::pair<MPI_Request, MPI_Request> operator()(const VariableMessage<uint64_t>& msg) {
        MPI_Request morton_req, data_req;
        /* Send the Morton indices */
        int err = MPI_Isend(msg.morton_indices_.data(), msg.morton_indices_.size(), MPI_MORTON_INDEX_T, msg.rank_, CreateMortonIndicesTag(msg.variable_id_), comm_, &morton_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.morton_indices_.size(), " Morton indices to rank ", msg.rank_, " for variable ", msg.variable_id_);

        /* Send the actual data */
        err = MPI_Isend(msg.data_.data(), msg.data_.size(), ConvertToMPIType<uint64_t>(), msg.rank_, CreateDataTag(msg.variable_id_), comm_, &data_req);
        MPICheckError(err);
        cmc_debug_msg("Send ", msg.data_.size(), " data values (of type: uint64_t) to rank ", msg.rank_, " for variable ", msg.variable_id_);
        
        return std::make_pair(morton_req, data_req);
    }

private:
    const MPI_Comm comm_;
};

std::pair<MPI_Request, MPI_Request>  
VariableSendMessage::Send(const MPI_Comm comm)
{
    return std::visit(SendMessage(comm), message_);
}

void*
VariableRecvMessage::GetInitialDataPtr()
{
    return std::visit([](auto& msg) -> void* {
        return msg.GetInitialDataPtr();
    }, message_);
}

void*
VariableRecvMessage::GetInitialMortonIndicesPtr()
{
    return std::visit([](auto& msg) -> void* {
        return msg.GetInitialMortonIndicesPtr();
    }, message_);
}

int
VariableRecvMessage::GetSendingRank() const
{
    return std::visit([](auto& msg) -> int {
        return msg.GetRank();
    }, message_);
}

int
VariableRecvMessage::GetVariableID() const
{
    return std::visit([](auto& msg) -> int {
        return msg.GetVariableID();
    }, message_);
}

}
