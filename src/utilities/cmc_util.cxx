#include "cmc_util.h"

[[noreturn]] void cmc_exit(const char* _err_msg, const char* _location)
{
    std::cout << "CMC_EXIT is invoked..." << std::endl << _err_msg << std::endl << "Error Occurence:"  << _location << std::endl;
    std::exit(EXIT_FAILURE);
}

[[noreturn]] void cmc_abort(const char* _err_msg, const char* _location)
{
    std::cout << "CMC_ABORT is invoked..." << std::endl << _err_msg << std::endl << "Error Occurence:"  << _location << std::endl;
    std::abort();
}

#if 0
std::string
cmc_convert_from_fortran_string(const char* f_string, int length)
{
    assert(f_string != nullptr);
    return std::string(f_string, length);
}
#endif
