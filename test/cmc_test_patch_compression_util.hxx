#include <vector>
#include <cmath>
#include <cstdint>

namespace cmc::test
{

/* Forward declaration of the data generator */
template <uint32_t DIM_X, uint32_t DIM_Y, uint32_t DIM_Z, uint32_t DIM_T>
std::vector<float>
GenerateExampleData_MD();

/* Typedef for 1D data */
template <uint32_t DIM_X>
auto const GenerateExampleData_1D = GenerateExampleData_MD<DIM_X, 1, 1, 1>;

/* Typedef for 2D data */
template <uint32_t DIM_X, uint32_t DIM_Y>
auto const GenerateExampleData_2D = GenerateExampleData_MD<DIM_X, DIM_Y, 1, 1>;

/* Typedef for 3D data */
template <uint32_t DIM_X, uint32_t DIM_Y, uint32_t DIM_Z>
auto const GenerateExampleData_3D = GenerateExampleData_MD<DIM_X, DIM_Y, DIM_Z, 1>;

/* Typedef for 4D data */
template <uint32_t DIM_X, uint32_t DIM_Y, uint32_t DIM_Z, uint32_t DIM_T>
auto const GenerateExampleData_4D = GenerateExampleData_MD<DIM_X, DIM_Y, DIM_Z, DIM_T>;

template <uint32_t DIM_X, uint32_t DIM_Y, uint32_t DIM_Z, uint32_t DIM_T>
std::vector<float>
GenerateExampleData_MD()
{
    static_assert(DIM_X >= 1 && DIM_Y >= 1 && DIM_Z >= 1 && DIM_T >= 1);

    std::vector<float> data(DIM_X * DIM_Y * DIM_Z * DIM_T);

    constexpr float t_scaling = (DIM_T > 1 ? (2 * std::numbers::pi_v<float>) / static_cast<float>(DIM_T - 1) : 0.0);
    constexpr float z_scaling = (DIM_Z > 1 ? (2 * std::numbers::pi_v<float>) / static_cast<float>(DIM_Z - 1) : 0.0);
    constexpr float y_scaling = (DIM_Y > 1 ? (2 * std::numbers::pi_v<float>) / static_cast<float>(DIM_Y - 1) : 0.0);
    constexpr float x_scaling = (DIM_X > 1 ? (2 * std::numbers::pi_v<float>) / static_cast<float>(DIM_X - 1) : 0.0);

    for (uint32_t t{0}; t < DIM_T; ++t)
    {
        const uint32_t t_offset = t * DIM_Z * DIM_Y * DIM_X;
        for (uint32_t z{0}; z < DIM_Z; ++z)
        {
            const uint32_t z_offset = z * DIM_Y * DIM_X;
            for (uint32_t y{0}; y < DIM_Y; ++y)
            {
                const uint32_t y_offset = y * DIM_X;
                for (uint32_t x{0}; x < DIM_X; ++x)
                {
                    const uint32_t idx = t_offset + z_offset + y_offset + x;
                    data[idx] = std::sin(x_scaling * x) + std::sin(y_scaling * y) + std::sin(z_scaling * z) + std::sin(t_scaling * t);
                }   
            }   
        }
    }

    return data;
}



}