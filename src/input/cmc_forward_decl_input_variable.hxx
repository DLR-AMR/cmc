#ifndef INPUT_CMC_FORWARD_DECL_INPUT_VARIABLE
#define INPUT_CMC_FORWARD_DECL_INPUT_VARIABLE

#include "utilities/cmc_utilities.hxx"

#ifdef CMC_ENABLE_MPI
#include "mpi/cmc_mpi_data.hxx"
#endif

#include <map>
#include <variant>

namespace cmc::input
{

/** Forward declarations **/
#ifdef CMC_ENABLE_MPI
template<typename T>
using ReceiverMap = std::map<int, VariableMessage<T>>;
#endif

template<class T> class Variable;
class Var;

using GeneralVariable = std::variant<Variable<int8_t>, Variable<char>, Variable<int16_t>, Variable<int32_t>, Variable<float>, Variable<double>,
                                      Variable<uint8_t>, Variable<uint16_t>, Variable<uint32_t>, Variable<int64_t>, Variable<uint64_t>>;

using CmcDefaultDataType = float;
inline CmcType GetDefaultCmcType() { return CmcType::Float; }

template<class T> std::vector<Variable<T>> ExtractSubVariables (const Variable<T>& variable, const Dimension split_dimension);
    
template<class T> Variable<T> HollowCopy (const Variable<T>& variable);

Var MetaCopy(const Var& variable);

class AccessKey;

struct IndexReduction;

class UpdateLinearIndices;


}




#endif /* !INPUT_CMC_FORWARD_DECL_INPUT_VARIABLE */
