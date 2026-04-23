#ifndef CMC_UTILITIES_HXX
#define CMC_UTILITIES_HXX

#include "cmc_config.h"

#include <iostream>
#include <cstddef>
#include <cstdint>
#include <variant>
#include <vector>
#include <array>
#include <cassert>
#include <limits>
#include <string>
#include <cstring>

#define CMC_MACRO_EXPANSION(x) #x
#define CMC_MACRO_EXPANSION2(x) CMC_MACRO_EXPANSION(x)
#define CMC_FILE_LOCATION __FILE__ ": " CMC_MACRO_EXPANSION2(__LINE__)

namespace cmc
{

#ifdef CMC_ENABLE_DEBUG
#define cmc_assert(condition) assert(condition)
#define cmc_static_assert(condition) static_assert(condition)
#else
#define cmc_assert(condition) ((void)0)
#define cmc_static_assert(condition) ((void)0)
#endif

[[noreturn]] void cmc_exit(const char* _err_msg, const char* _location);

[[noreturn]] void cmc_abort(const char* _err_msg, const char* _location);

enum CmcType {TypeUndefined = -1, Int8_t, Char, Int16_t, Int32_t, Float, Double, Uint8_t, Uint16_t, Uint32_t, Int64_t, Uint64_t, NumTypes};

/* A type for holding 'arbitrary' data, in particular any possible CmcType */
using CmcUniversalType = std::variant<int8_t, char, int16_t, int32_t, float, double, uint8_t, uint16_t, uint32_t, int64_t, uint64_t>;

/* Helper function for the variant */
template<class...>
constexpr bool always_false_v = false;

template<class... Ts>
struct overloaded : Ts...
{
  using Ts::operator()...;

  /* Prevent implicit type conversions */
  template<typename T>
  constexpr void operator()(T) const
  {
    static_assert(always_false_v<T>, "Unsupported type");
  }
};

template<class... Ts>
overloaded(Ts...) -> overloaded<Ts...>;

inline constexpr
size_t CmcTypeToBytes(const CmcType type)
{
    cmc_assert(type > CmcType::TypeUndefined && type < CmcType::NumTypes);
    constexpr std::array<size_t, CmcType::NumTypes> type_to_bytes = {sizeof(int8_t), sizeof(char), sizeof(int16_t), sizeof(int32_t), sizeof(float), sizeof(double), sizeof(uint8_t), sizeof(uint16_t), sizeof(uint32_t), sizeof(int64_t), sizeof(uint64_t)};
    return type_to_bytes[type];
}

inline constexpr
bool IsSameType(const CmcType& type1, const CmcType& type2)
{
    return (type1 == type2);
}

template<typename T>
auto ApproxCompare(const T& value1, const T& value2, T&& epsilon = 4 * std::numeric_limits<T>::epsilon())
 -> std::enable_if_t<std::is_floating_point_v<T>, bool>
{
    return (std::abs(value1 - value2) < epsilon ? true : false);
}

template<typename T>
auto ApproxCompare(const T& value1, const T& value2)
 -> std::enable_if_t<std::is_integral_v<T>, bool>
{
    return (value1 == value2);
}

template<typename T>
auto GetUniversalDataAs(const CmcUniversalType& universal_data)
    -> std::enable_if_t<std::is_arithmetic_v<T>, T>
{
    return std::visit([](auto&& value) -> T {
        return static_cast<T>(value);
    }, universal_data);
}

template<typename T>
constexpr CmcType
ConvertToCmcType()
{
    if constexpr (std::is_same_v<T, int8_t>)
    {
        return CmcType::Int8_t;
    } else if constexpr (std::is_same_v<T, char>)
    {
        return CmcType::Char;
    } else if constexpr (std::is_same_v<T, int16_t>)
    {   
        return CmcType::Int16_t;
    } else if constexpr (std::is_same_v<T, int32_t>)
    {
        return CmcType::Int32_t;
    } else if constexpr (std::is_same_v<T, float>)
    {
        return CmcType::Float;
    } else if constexpr (std::is_same_v<T, double>)
    {
        return CmcType::Double;
    } else if constexpr (std::is_same_v<T, uint8_t>)
    {
        return CmcType::Uint8_t;
    } else if constexpr (std::is_same_v<T, uint16_t>)
    {
        return CmcType::Uint16_t;
    } else if constexpr (std::is_same_v<T, uint32_t>)
    {
        return CmcType::Uint32_t;
    } else if constexpr (std::is_same_v<T, int64_t>)
    {
        return CmcType::Int64_t;
    } else if constexpr (std::is_same_v<T, uint64_t>)
    {
        return CmcType::Uint64_t;
    } else
    {
        std::cout << "[cmc] ERROR: The template parameter could not be converted to a CmcType." << std::endl;
        std::exit(EXIT_FAILURE);
        return CmcType::TypeUndefined;
    }
}

[[noreturn]] inline void
cmc_exit(const char* _err_msg, const char* _location)
{
    std::cout << "CMC_EXIT is invoked..." << std::endl << _err_msg << std::endl << "Error Occurence:"  << _location << std::endl;
    std::exit(EXIT_FAILURE);
}

[[noreturn]] inline void
cmc_abort(const char* _err_msg, const char* _location)
{
    std::cout << "CMC_ABORT is invoked..." << std::endl << _err_msg << std::endl << "Error Occurence:"  << _location << std::endl;
    std::abort();
}

}

#endif /* !CMC_UTILITIES_HXX */
