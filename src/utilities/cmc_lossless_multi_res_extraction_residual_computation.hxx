#ifndef CMC_LOSSLESS_MULTI_RES_EXTRACTION_RESIDUAL_COMPUTATION_HXX
#define CMC_LOSSLESS_MULTI_RES_EXTRACTION_RESIDUAL_COMPUTATION_HXX

#include "utilities/cmc_vector_view.hxx"
#include "utilities/cmc_byte_compression_values.hxx"

#include <type_traits>
#include <cstring>
#include <cstddef>
#include <utility>

namespace cmc::lossless::multi_res
{

/* A typedef for the sake of brevity */
template <typename T>
using CompressionValue = SerializedCompressionValue<sizeof(T)>;

template <typename T>
auto
GetResidualsLZC(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint32_t), int>
{
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint32_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint32_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    /* Compute the unsigned residual and determine the LZC */
    const uint32_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);
    return diff_compression_val.GetNumberLeadingZeros();
}

template <typename T>
auto
GetResidualsLZC(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint16_t), int>
{
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint16_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint16_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    /* Compute the unsigned residual and determine the LZC */
    const uint16_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);
    return diff_compression_val.GetNumberLeadingZeros();
}

template <typename T>
auto
GetResidualsLZC(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint8_t), int>
{
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint8_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint8_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    /* Compute the unsigned residual and determine the LZC */
    const uint8_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);
    return diff_compression_val.GetNumberLeadingZeros();
}

template <typename T>
auto
GetResidualsLZC(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint64_t), int>
{
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint64_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint64_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    /* Compute the unsigned residual and determine the LZC */
    const uint64_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

    /* Convert the difference to a CompressionValue and determine the LZC */
    const CompressionValue<T> diff_compression_val(diff);
    return diff_compression_val.GetNumberLeadingZeros();
}

template <typename T>
auto
GetCumulativeResidualsLZC(const T& approximation, const VectorView<CompressionValue<T>>& real_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint32_t), int>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Get an unsigned representation of the approximation value */
    uint32_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = real_values.begin(); val_iter != real_values.end(); ++val_iter)
    {
        /* Get an unsigned representation of the real value */
        uint32_t unsigned_real_value;
        std::memcpy(&unsigned_real_value, val_iter->GetMemoryForReading().data(), N);

        /* Compute the unsigned residual and determine the LZC */
        const uint32_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();
    }

    return cumulative_lzc;
}

template <typename T>
auto
GetCumulativeResidualsLZC(const T& approximation, const VectorView<CompressionValue<T>>& real_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint16_t), int>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Get an unsigned representation of the approximation value */
    uint16_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = real_values.begin(); val_iter != real_values.end(); ++val_iter)
    {
        /* Get an unsigned representation of the real value */
        uint16_t unsigned_real_value;
        std::memcpy(&unsigned_real_value, val_iter->GetMemoryForReading().data(), N);

        /* Compute the unsigned residual and determine the LZC */
        const uint16_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();
    }

    return cumulative_lzc;
}

template <typename T>
auto
GetCumulativeResidualsLZC(const T& approximation, const VectorView<CompressionValue<T>>& real_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint8_t), int>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Get an unsigned representation of the approximation value */
    uint8_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = real_values.begin(); val_iter != real_values.end(); ++val_iter)
    {
        /* Get an unsigned representation of the real value */
        uint8_t unsigned_real_value;
        std::memcpy(&unsigned_real_value, val_iter->GetMemoryForReading().data(), N);

        /* Compute the unsigned residual and determine the LZC */
        const uint8_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();
    }

    return cumulative_lzc;
}

template <typename T>
auto
GetCumulativeResidualsLZC(const T& approximation, const VectorView<CompressionValue<T>>& real_values)
 -> std::enable_if_t<sizeof(T) == sizeof(uint64_t), int>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    int cumulative_lzc = 0;

    /* Get an unsigned representation of the approximation value */
    uint64_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Iterate over all values, compute the unsigned residuals and determine the LZC*/
    for (auto val_iter = real_values.begin(); val_iter != real_values.end(); ++val_iter)
    {
        /* Get an unsigned representation of the real value */
        uint64_t unsigned_real_value;
        std::memcpy(&unsigned_real_value, val_iter->GetMemoryForReading().data(), N);

        /* Compute the unsigned residual and determine the LZC */
        const uint64_t diff = (unsigned_real_value >= unsigned_approximation ? unsigned_real_value - unsigned_approximation : unsigned_approximation - unsigned_real_value);

        /* Convert the difference to a CompressionValue and determine the LZC */
        const CompressionValue<T> diff_compression_val(diff);

        cumulative_lzc += diff_compression_val.GetNumberLeadingZeros();
    }

    return cumulative_lzc;
}


template <typename T>
auto
ComputeResidual(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint32_t), std::pair<bool, CompressionValue<T>>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint32_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint32_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    const bool is_approximation_greater = (unsigned_approximation >= unsigned_real_value);

    /* Compute the residual */
    const uint32_t diff = (unsigned_approximation >= unsigned_real_value ? unsigned_approximation - unsigned_real_value : unsigned_real_value - unsigned_approximation);

    /* Convert the difference to a CompressionValue */
    const CompressionValue<T> diff_compression_val(diff);

    /* Return the flag and the residual */
    return std::make_pair(is_approximation_greater, diff_compression_val);
}


template <typename T>
auto
ComputeResidual(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint8_t), std::pair<bool, CompressionValue<T>>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint8_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint8_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    const bool is_approximation_greater = (unsigned_approximation >= unsigned_real_value);

    /* Compute the residual */
    const uint8_t diff = (unsigned_approximation >= unsigned_real_value ? unsigned_approximation - unsigned_real_value : unsigned_real_value - unsigned_approximation);

    /* Convert the difference to a CompressionValue */
    const CompressionValue<T> diff_compression_val(diff);

    /* Return the flag and the residual */
    return std::make_pair(is_approximation_greater, diff_compression_val);
}

template <typename T>
auto
ComputeResidual(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint16_t), std::pair<bool, CompressionValue<T>>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint16_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint16_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    const bool is_approximation_greater = (unsigned_approximation >= unsigned_real_value);

    /* Compute the residual */
    const uint16_t diff = (unsigned_approximation >= unsigned_real_value ? unsigned_approximation - unsigned_real_value : unsigned_real_value - unsigned_approximation);

    /* Convert the difference to a CompressionValue */
    const CompressionValue<T> diff_compression_val(diff);

    /* Return the flag and the residual */
    return std::make_pair(is_approximation_greater, diff_compression_val);
}

template <typename T>
auto
ComputeResidual(const T& approximation, const CompressionValue<T>& real_value)
 -> std::enable_if_t<sizeof(T) == sizeof(uint64_t), std::pair<bool, CompressionValue<T>>>
{
    /* Byte size of the data type */
    constexpr size_t N = sizeof(T);

    /* Get an unsigned representation of the approximation value */
    uint64_t unsigned_approximation;
    std::memcpy(&unsigned_approximation, &approximation, N);

    /* Get an unsigned representation of the real value */
    uint64_t unsigned_real_value;
    std::memcpy(&unsigned_real_value, real_value.GetMemoryForReading().data(), N);

    const bool is_approximation_greater = (unsigned_approximation >= unsigned_real_value);

    /* Compute the residual */
    const uint64_t diff = (unsigned_approximation >= unsigned_real_value ? unsigned_approximation - unsigned_real_value : unsigned_real_value - unsigned_approximation);

    /* Convert the difference to a CompressionValue */
    const CompressionValue<T> diff_compression_val(diff);

    /* Return the flag and the residual */
    return std::make_pair(is_approximation_greater, diff_compression_val);
}


}


#endif /* !CMC_LOSSLESS_MULTI_RES_EXTRACTION_RESIDUAL_COMPUTATION_HXX */
