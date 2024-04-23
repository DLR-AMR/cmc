#ifndef CMC_MPI_DATA_HXX
#define CMC_MPI_DATA_HXX

#include "utilities/cmc_utilities.hxx"
#include "mpi/cmc_mpi.hxx"

namespace cmc
{

constexpr int kVarMessageMPIErr = -1;

template<class T>
class VariableMessage
{
public:
    VariableMessage() = default;
    VariableMessage(int receiving_rank, const int variable_id)
    : receiving_rank_{receiving_rank}, variable_id_{variable_id}{};

    int receiving_rank_{kVarMessageMPIErr};
    int variable_id_{kVarMessageMPIErr};
    std::vector<T> data_;
    std::vector<MortonIndex> morton_indices_;
};

template<typename T>
static 
void
DetachVarMessage(VariableMessage<T>& source, VariableMessage<T>& destination)
{
    destination.receiving_rank_ = source.receiving_rank_;
    destination.variable_id_ = source.variable_id_;
    std::swap(source.data_, destination.data_);
    std::swap(source.morton_indices_, destination.morton_indices_);
}

using VarMessage = std::variant<VariableMessage<int8_t>, VariableMessage<char>, VariableMessage<int16_t>, VariableMessage<int32_t>, VariableMessage<float>, VariableMessage<double>,
                                 VariableMessage<uint8_t>, VariableMessage<uint16_t>, VariableMessage<uint32_t>, VariableMessage<int64_t>, VariableMessage<uint64_t>>;

class VariableSendMessage
{
public:
    //TODO: type needs to be set corretcly
    template<typename T> VariableSendMessage(VariableMessage<T>&& var_message)
    : type_{TypeUndefined}, message_{std::move(var_message)} {};


    //MPI_Status* SendMessage();
    
private:
    CmcType type_;
    VarMessage message_;
};

}

#endif /* !CMC_MPI_DATA_HXX */
