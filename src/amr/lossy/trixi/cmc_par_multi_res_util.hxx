#ifndef CMC_LOSSY_TRIXI_PAR_MULTI_RES_UTIL_HXX
#define CMC_LOSSY_TRIXI_PAR_MULTI_RES_UTIL_HXX

#include "cmc.hxx"

#include <array>
#include <cstdint>
#include <vector>

namespace cmc::par::lossy::multi_res::idw::trixi
{

template<int32_t DIM>
requires Dimension<DIM>
constexpr
inline int
GetNumIntraElemControlValues()
{
    if constexpr (DIM == 1)
    {
        return 2;
    } else if constexpr (DIM == 2)
    {
        return 4;
    } else if constexpr (DIM == 3)
    {
        return 8;
    } else if constexpr (DIM == 4)
    {
        return 16;
    } else
    {
        return 1;
    }
}

template <uint32_t P>
constexpr int32_t
integral_pow(const int32_t base)
{
    if constexpr (P == 0)
    {
        return 1;
    }
    if constexpr (P == 1)
    {
        return base;
    }

    const int temp = integral_pow<P / 2>(base);

    if constexpr (P % 2 == 0)
    {
        return temp * temp;
    } else
    {
        return base * temp * temp;
    }
}

template <int32_t DIM>
requires Dimension<DIM>
constexpr int kNumIntraElemControlValues = GetNumIntraElemControlValues<DIM>();

template <int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr std::array<int, kNumIntraElemControlValues<DIM>>
GetControlValueIndices()
{
    static_assert(DIM >= 1 || DIM < 4, "Currently, the supplied dimension is not supported.");

    std::array<int, kNumIntraElemControlValues<DIM>> indices{};

    if constexpr (DIM >= 1)
    {
        indices[0] = 0;
        indices[1] = (N - 1);
    }
    
    if constexpr (DIM >= 2)
    {
        indices[2] = (N - 1) * N;
        indices[3] = (N * N - 1);
    } 

    if constexpr (DIM >= 3)
    {
        indices[4] = (N * N * (N - 1));
        indices[5] = (N * N * (N - 1) + N - 1);
        indices[6] = (N * N * (N - 1) + (N - 1) * N);
        indices[7] = (N * N * (N - 1) + (N * N - 1));
    }

    return indices;
}

template <int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr int kGetNumPointsPerElement()
{
    return integral_pow<DIM>(N);
}

template <int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr int kNumInitDataPredictions = std::max(kGetNumPointsPerElement<DIM, N>() - GetNumIntraElemControlValues<DIM>(), 0);


template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr
inline bool IsControlValueIndex(const int check_idx)
{
    constexpr std::array<int, kNumIntraElemControlValues<DIM>> control_value_indices = GetControlValueIndices<DIM, N>();
    for (int idx{0}; idx < kNumIntraElemControlValues<DIM>; ++idx)
    {
        if (check_idx == control_value_indices[idx])
        {
            return true;
        }
    } 
    return false;
}

template<int32_t DIM, int32_t N>
requires Dimension<DIM>
constexpr
inline std::array<int, kNumInitDataPredictions<DIM, N>>
GetPredictionReferenceCoordinates()
{
    std::array<int, kNumInitDataPredictions<DIM, N>> coords_indices{};

    /* Get the control value indices */
    constexpr std::array<int, kNumIntraElemControlValues<DIM>> control_value_indices = GetControlValueIndices<DIM, N>();

    /* Iterate over all element indices */
    constexpr int kNumElemPoints = kGetNumPointsPerElement<DIM, N>();
    int access_idx{0};
    for (int data_idx{0}; data_idx < kNumElemPoints; ++data_idx)
    {
        if (not IsControlValueIndex<DIM, N>(data_idx))
        {
            coords_indices[access_idx] = data_idx;
            ++access_idx;
        }
    }

    return coords_indices;
}

template<ArithmeticType T, int32_t DIM, int32_t N>
requires Dimension<DIM>
struct IntraElemEncodingData
{
    std::array<SymbolType, kNumIntraElemControlValues<DIM>> control_quantization_bins{};
    std::array<SymbolType, kNumInitDataPredictions<DIM, N>> init_quantization_bins{};
    std::vector<T> unpredictable_values{};
};

}

#endif /* !CMC_LOSSY_TRIXI_PAR_MULTI_RES_UTIL_HXX */
