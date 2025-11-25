#ifndef CMC_SERIAL_IO_UTIL_HXX
#define CMC_SERIAL_IO_UTIL_HXX

#include <string>
#include <cstdint>

namespace cmc::compression_io::serial
{

constexpr size_t kNumHeaderInfo = 6;

using attr_type_t = int32_t;

inline std::string
CreateFileName(const std::string& file_prefix, const std::string& variable_name)
{
    return file_prefix + "_" + variable_name + ".cmc";
}

}

#endif /* !CMC_SERIAL_IO_UTIL_HXX */
