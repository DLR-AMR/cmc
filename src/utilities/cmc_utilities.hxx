#ifndef CMC_UTILITIES_HXX
#define CMC_UTILITIES_HXX

#include <cstddef>
#include <variant>

namespace cmc {

enum Dimension {DimensionUndefined = -1, Lon = 0, Lat = 1, Lev = 2, Time = 3, NumCoordinates};

enum DataLayout {LayoutUndefined, Lat_Lon, Lon_Lat, Lat_Lev, Lev_Lat, Lon_Lev, Lev_Lon, _InternEnd2DLayouts,
                  Lat_Lon_Lev, Lat_Lev_Lon, Lev_Lat_Lon, Lev_Lon_Lat, Lon_Lev_Lat, Lon_Lat_Lev};

enum CmcType {TypeUndefined = -1, Byte, Int8_t, Char, Int16_t, Int32_t, Float, Double, Uint8_t, Uint16_t, Uint32_t, Int64_t, Uint64_t, NumTypes};

/* A type for holding 'arbitrary' data, in particular any possible CmcType */
typedef std::variant<std::byte, int8_t, char, int16_t, int32_t, float, double, uint8_t, uint16_t, uint32_t, int64_t, uint64_t> CmcUniversalType;

inline constexpr
bool isSameType(const CmcType& type1, const CmcType& type2)
{
    return (type1.index() == type2.index());
}

inline constexpr
size_t CmcTypeToBytes(const CmcType type)
{
    constexpr std::array<size_t, CmcType::NumTypes> type_to_bytes = {sizeof(std::byte), sizeof(int8_t), sizeof(char), sizeof(int16_t), sizeof(int32_t), sizeof(float), sizeof(double), sizeof(uint8_t), sizeof(uint16_t), sizeof(uint32_t), sizeof(int64_t), sizeof(uint64_t)};
    return type_to_bytes[type];
}

struct DimensionInterval
{
    constexpr DimensionInterval(Dimension dimension, int start_idx, int end_idx)
    : dim{dimension}, start_index{start_idx}, end_index{end_idx}{};
    Dimension dim{Dimension::DimensionUndefined};
    int start_index{0};
    int end_index{0};
};

class GeoDomain
{
public:
    /* 1D Constructor */
    constexpr GeoDomain(const Dimension dimension1, const int start_idx_dim1, const int end_idx_dim1);
    /* 2D Constructor */
    constexpr GeoDomain(const Dimension dimension1, const int start_idx_dim1, const int end_idx_dim1,
              const Dimension dimension2, const int start_idx_dim2, const int end_idx_dim2);
    /* 3D Constructor */
    constexpr GeoDomain(const Dimension dimension1, const int start_idx_dim1, const int end_idx_dim1,
              const Dimension dimension2, const int start_idx_dim2, const int end_idx_dim2,
              const Dimension dimension3, const int start_idx_dim3, const int end_idx_dim3);
    /* 4D constructor */
    constexpr GeoDomain(const Dimension dimension1, const int start_idx_dim1, const int end_idx_dim1,
              const Dimension dimension2, const int start_idx_dim2, const int end_idx_dim2,
              const Dimension dimension3, const int start_idx_dim3, const int end_idx_dim3,
              const Dimension dimension4, const int start_idx_dim4, const int end_idx_dim4);

    GeoDomain(const GeoDomain& other) = default;
    GeoDomain& operator=(const GeoDomain& other) = default;
    GeoDomain(GeoDomain&& other) = default;
    GeoDomain& operator=(GeoDomain&& other) = default;

    ~GeoDomain() = default;

    constexpr int GetDimensionality() const;
    void AddDimension(const Dimension dimension, const int start_idx, const int end_idx);
private:
    std::vector<DimensionInterval> domain_;
};

template <typename T>
class CoordinateArray
{
public:
    T& operator[](size_t idx){
        cmc_assert(idx < Dimension::NumCoordinates);
        return coordinates[idx];
    }
    const T& operator[](size_t idx) const {
        cmc_assert(idx < Dimension::NumCoordinates);
        return coordinates[idx];
    }

    auto begin(){
        return coordinates.begin();
    }
    auto end()
    {
        return coordinates.end();
    }

    std::array<T, Dimension::NumCoordinates> coordinates{0, 0, 0, 0};
};

struct CartesianCoordinate
{
    CartesianCoordinate(const int64_t x_coordinate, const int64_t y_coordinate = 0, const int64_t z_coordinate = 0)
    : x{x_coordinate}, y{y_coordinate}, z{z_coordinate} {};

    int64_t x{0}, y{0}, z{0};
};

struct Hyperslab
{
    std::vector<int64_t> start_of_boxed_domain;
    std::vector<int64_t> count_of_boxed_domain;
};

struct GlobalDomain
{
    CoordinateArray<int64_t> start_coordinates;
    CoordinateArray<int64_t> end_coordinates;
    int dimensionality{0};
};

}
#endif /* !CMC_UTILITIES_HXX */
