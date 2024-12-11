#ifndef CMC_INPUT_VARIABLE_FORWARD_HXX
#define CMC_INPUT_VARIABLE_FORWARD_HXX

#include "utilities/cmc_utilities.hxx"
#include "mpi/cmc_mpi_data.hxx"
#include "t8code/cmc_t8_mpi.hxx"

#include <map>

namespace  cmc
{

template<typename T>
using ReceiverMap = std::map<int, VariableMessage<T>>;

/* Forward declarations */
template<class T> class InputVariable;
class InputVar;

using CmcInputVariable = std::variant<InputVariable<int8_t>, InputVariable<char>, InputVariable<int16_t>, InputVariable<int32_t>, InputVariable<float>, InputVariable<double>,
                                      InputVariable<uint8_t>, InputVariable<uint16_t>, InputVariable<uint32_t>, InputVariable<int64_t>, InputVariable<uint64_t>>;

using CmcDefaultDataType = float;
inline CmcType GetDefaultCmcType() { return CmcType::Float; };

//template<class T> ReceiverMap<T> GatherDataToBeDistributed(const InputVariable<T>& variable, const DataOffsets& offsets);

template<class T> std::vector<InputVariable<T>> ExtractSubVariables (const InputVariable<T>& variable, const Dimension split_dimension);
    
template<class T> InputVariable<T> HollowCopy (const InputVariable<T>& variable);

InputVar MetaCopy(const InputVar& variable);

class AccessKey;

struct IndexReduction;

class UpdateLinearIndices;

}

#endif /* !CMC_INPUT_VARIABLE_FORWARD_HXX */
