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

#define CMC_MACRO_EXPANSION(x) #x
#define CMC_MACRO_EXPANSION2(x) CMC_MACRO_EXPANSION(x)
#define CMC_FILE_LOCATION __FILE__ ": " CMC_MACRO_EXPANSION2(__LINE__)

namespace cmc {

#if CMC_ENABLE_DEBUG
#define cmc_assert(condition) assert(condition)
#else
#define cmc_assert(condition) ((void)0)
#endif

[[noreturn]] void cmc_exit(const char* _err_msg, const char* _location);

[[noreturn]] void cmc_abort(const char* _err_msg, const char* _location);

enum Dimension {DimensionUndefined = -1, Lon = 0, Lat = 1, Lev = 2, Time = 3, NumCoordinates};

enum CmcType {TypeUndefined = -1, Int8_t, Char, Int16_t, Int32_t, Float, Double, Uint8_t, Uint16_t, Uint32_t, Int64_t, Uint64_t, NumTypes};

//TODO: Maybe make a nicer design of the layout, such that there is a pattern between layouts 
enum DataLayout {LayoutUndefined, Lat_Lon, Lon_Lat, Lat_Lev, Lev_Lat, Lon_Lev, Lev_Lon, _InternEnd2DLayouts,
                  Lat_Lon_Lev, Lat_Lev_Lon, Lev_Lat_Lon, Lev_Lon_Lat, Lon_Lev_Lat, Lon_Lat_Lev, _InternEnd3DLayouts};

enum DataRepresentation {RepresentationUndefined, SpaceFillingCurve, CartesianCoordinates, HyperslabCoordinates};

/* A type for holding 'arbitrary' data, in particular any possible CmcType */
using CmcUniversalType = std::variant<int8_t, char, int16_t, int32_t, float, double, uint8_t, uint16_t, uint32_t, int64_t, uint64_t>;

typedef int64_t MortonIndex;
typedef int64_t LinearIndex;
typedef int64_t DomainIndex;

enum DataFormat {FormatUndefined, LinearFormat, CartesianFormat, HyperslabFormat};

/* Helper function for the variant */
template<class... Ts>
struct overloaded : Ts... { using Ts::operator()...; };
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

DataLayout
GetDataLayoutAfterRemoval(const DataLayout initial_layout, const Dimension removed_dimension);

std::vector<Dimension>
GetDimensionVectorFromLayout(const DataLayout layout);

int
GetDimensionalityOfDataLayout(const DataLayout layout);

}

#endif /* !CMC_UTILITIES_HXX */
