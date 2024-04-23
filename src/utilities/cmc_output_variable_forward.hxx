#ifndef CMC_OUTPUT_VARIABLE_FORWARD_HXX
#define CMC_OUTPUT_VARIABLE_FORWARD_HXX

#include "utilities/cmc_utilities.hxx"

namespace  cmc
{

/* Forward declarations */
template<class T> class OutputVariable;
class OutputVar;

using CmcOutputVariable = std::variant<OutputVariable<int8_t>, OutputVariable<char>, OutputVariable<int16_t>, OutputVariable<int32_t>, OutputVariable<float>, OutputVariable<double>,
                                      OutputVariable<uint8_t>, OutputVariable<uint16_t>, OutputVariable<uint32_t>, OutputVariable<int64_t>, OutputVariable<uint64_t>>;


}

#endif /* !CMC_OUTPUT_VARIABLE_FORWARD_HXX */
