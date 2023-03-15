#ifndef CMC_UTIL_H
#define CMC_UTIL_H

/* Standard Headers */
#include <iostream>
#include <string>
#include <cmath>
#include <memory>
#include <array>
#include <vector>
#include <cassert>
#include <algorithm>
#include <variant>
#include <limits>
#include <cstring>
#include <functional>
#include <cstddef>

#include "cmc_constants_definitions.h"

/* Define a typedef of variant able to hold severeal data types */
typedef std::variant<std::byte, int8_t, char, int16_t, int32_t, float, double, uint8_t, uint16_t, uint32_t, int64_t, uint64_t> cmc_universal_type_t;

[[noreturn]] void cmc_exit(const char* _err_msg, const char* _location);

[[noreturn]] void cmc_abort(const char* _err_msg, const char* _location);

template<typename T>
inline auto
cmc_approx_cmp(T num1, T num2, T epsilon = 2*std::numeric_limits<T>::epsilon())
    -> std::enable_if_t<std::is_signed_v<T>, bool>
{
    if (std::abs(num1 - num2) <= epsilon)
    {
        return true;
    } else
    {
        return false;
    }
}

template<typename T>
inline auto
cmc_approx_cmp(T num1, T num2, T epsilon = 2*std::numeric_limits<T>::epsilon())
    -> std::enable_if_t<std::is_unsigned_v<T>, bool>
{
    const T _diff = (num1 > num2) ? num1 - num2 : num2 - num1;
    if (_diff <= epsilon)
    {
        return true;
    } else
    {
        return false;
    }
}

template<typename T>
auto cmc_get_universal_data(const cmc_universal_type_t& _universal_data)
    -> std::enable_if_t<std::is_arithmetic_v<T>,T>
{
    return std::visit([](auto&& value) -> T {
        return static_cast<T>(value);
    }, _universal_data);
}

#if 0
std::string
cmc_convert_from_fortran_string(const char* f_string, int length);
#endif

#endif /* CMC_UTIL_H */