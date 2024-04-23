#ifndef CMC_T8_DATA_VARIABLES_FORWARD_HXX
#define CMC_T8_DATA_VARIABLES_FORWARD_HXX

#include "utilities/cmc_utilities.hxx"

namespace cmc
{

template<class T>
class Variable;

using CmcVariable = std::variant<Variable<int8_t>, Variable<char>, Variable<int16_t>, Variable<int32_t>, Variable<float>, Variable<double>,
                                 Variable<uint8_t>, Variable<uint16_t>, Variable<uint32_t>, Variable<int64_t>, Variable<uint64_t>>;


class Var;

}

#endif /* !CMC_T8_DATA_VARIABLES_FORWARD_HXX */
