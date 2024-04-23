#include "t8code/cmc_t8_data_variables.hxx"
#include "utilities/cmc_log_functions.h"

namespace cmc
{

/* Explicit template instantiations */
template class Variable<int8_t>;
template class Variable<char>;
template class Variable<int16_t>;
template class Variable<int32_t>;
template class Variable<float>;
template class Variable<double>;
template class Variable<uint8_t>;
template class Variable<uint16_t>;
template class Variable<uint32_t>;
template class Variable<int64_t>;
template class Variable<uint64_t>;

}
