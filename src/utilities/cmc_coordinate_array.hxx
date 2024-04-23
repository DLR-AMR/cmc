#ifndef CMC_COORDINATE_ARRAY_HXX
#define CMC_COORDINATE_ARRAY_HXX

#include "utilities/cmc_utilities.hxx"

#include <array>

namespace cmc
{

template <typename T>
class CoordinateArray
{
public:
    CoordinateArray() = default;
    CoordinateArray(const T fill_value)
    : coordinates{fill_value, fill_value, fill_value, fill_value}{};

    T& operator[](int idx){
        cmc_assert(idx < Dimension::NumCoordinates);
        return coordinates[idx];
    }
    const T& operator[](int idx) const {
        cmc_assert(idx < Dimension::NumCoordinates);
        return coordinates[idx];
    }

    std::array<T, Dimension::NumCoordinates>::iterator begin()
    {
        return coordinates.begin();
    }
    std::array<T, Dimension::NumCoordinates>::iterator end()
    {
        return coordinates.end();
    }

    std::array<T, Dimension::NumCoordinates>::const_iterator begin() const
    {
        return coordinates.begin();
    }
    std::array<T, Dimension::NumCoordinates>::const_iterator end() const
    {
        return coordinates.end();
    }
    
    std::array<T, Dimension::NumCoordinates> coordinates{0, 0, 0, 0};
};

}

#endif /* !CMC_COORDINATE_ARRAY_HXX */
