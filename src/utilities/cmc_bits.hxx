#ifndef CMC_BITS_HXX
#define CMC_BITS_HXX

#include "cmc.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <bit>
#include <cstdint>
#include <concepts>
#include <cstddef>
#include <span>
#include <vector>
#include <algorithm>
#include <type_traits>
#include <array>
#include <cstring>

namespace cmc::bits
{

constexpr int32_t kCharBit = 8;

//Extend below with std__is__arithmetic
template<typename T>
concept OneByteType = (sizeof(T) == 1);

template<typename T>
concept TwoByteType = (sizeof(T) == 2);

template<typename T>
concept FourByteType = (sizeof(T) == 4);

template<typename T>
concept EightByteType = (sizeof(T) == 8);

template<typename T>
concept IntegerType = (std::is_integral_v<T>);

template<typename T>
concept UnsignedIntegerType = (std::is_unsigned_v<T> && std::is_integral_v<T>);

template<typename T>
concept ArithmeticType = (std::is_arithmetic_v<T> && std::is_fundamental_v<T>);

template<typename T>
concept PointerType = (std::is_pointer_v<T> && not std::is_member_function_pointer_v<T>);

template<IntegerType T>
constexpr T
ConvertToBigEndian(const T value)
{
    static_assert(std::endian::native == std::endian::big || std::endian::native == std::endian::little,
                  "Only little-endian and big-endian systems are suppoprted!");

    if constexpr (std::endian::native == std::endian::big)
    {
        /* The value is already big-endian */
        return value;
    } else
    {
        /* In case the values are little-endian */
        return std::byteswap(value);
    }
}

template<IntegerType T>
constexpr T
ConvertToLittleEndian(const T value)
{
    static_assert(std::endian::native == std::endian::big || std::endian::native == std::endian::little,
                  "Only little-endian and big-endian systems are suppoprted!");

    if constexpr (std::endian::native == std::endian::little)
    {
        /* The value is already little-endian */
        return value;
    } else
    {
        /* In case the values are big-endian */
        return std::byteswap(value);
    }
}

template<IntegerType T>
constexpr T
ConvertBigEndianToNativeEndianness(const T value)
{
    static_assert(std::endian::native == std::endian::big || std::endian::native == std::endian::little,
                  "Only little-endian and big-endian systems are suppoprted!");

    if constexpr (std::endian::native == std::endian::big)
    {
        /* The value is already big-endian */
        return value;
    } else
    {
        /* In case the values are little-endian */
        return std::byteswap(value);
    }
}

template<IntegerType T>
inline const T
ConvertBigEndianBytePositionToNativeEndiannessBytePosition(const int bigendian_byte_idx)
{
    static_assert(std::endian::native == std::endian::big || std::endian::native == std::endian::little,
                  "Only little-endian and big-endian systems are suppoprted!");
    if (sizeof(T) <= bigendian_byte_idx) [[unlikely]]
    {
        cmc_err_msg("The specified byte index is larger than the overall amount of bytes in type T!");
    }

    if constexpr (std::endian::native == std::endian::big)
    {
        return bigendian_byte_idx;
    } else
    {
        return sizeof(T) - (1 + bigendian_byte_idx);
    }
}

template<typename T>
class CRValue;

template<OneByteType T>
class CRValue<T>
{
public:

private:
    uint8_t bytes_;
    uint8_t front_bit_;
    uint8_t end_bit_;
};

template<TwoByteType T>
class CRValue<T>
{
public:
    
private:
    uint16_t bytes_;
    uint8_t front_bit_;
    uint8_t end_bit_;
};

template<FourByteType T>
class CRValue<T>
{
public:
    
private:
    uint32_t bytes_;
    uint16_t front_bit_;
    uint16_t end_bit_;
};

template<EightByteType T>
class CRValue<T>
{
public:
    
private:
    uint64_t bytes_;
    uint32_t front_bit_;
    uint32_t end_bit_;
};

struct
SignificantBitsIndicator
{
    uint8_t front_bit_{0};
    uint8_t end_bit_{0};
};

#if 0
template<typename T>
T
GetCommonPrefix(std::span<T> values)
{
    /* Get the maximum possible prefix length */
    const auto max_prefix_length_iter = std::max_element(values.begin(), indicators.end(), [](const SignificantBitsIndicator val1, const SignificantBitsIndicator val2)
                                                        {
                                                            return val1.end_bit_ < val2.end_bit_;
                                                        });
    cmc_assert(max_prefix_length_iter !=  indicators.end());
    const int max_prefix_length = *max_prefix_length_iter;

}


template<typename T>
struct SRValue;

template<OneByteType T>
struct SRValue<T>
{
    uint8_t value(0);
};

template<TwoByteType T>
struct SRValue<T>
{
    uint16_t value(0);
};

template<FourByteType T>
struct SRValue<T>
{
    uint32_t value(0);
};

template<EightByteType T>
struct SRValue<T>
{
    uint64_t value(0);
};

//TODO: Adapt for different types
template<FourByteType T>
T
GetCommonPrefix(const std::span<T> values)
{
    if (values.size() < 2) [[unlikely]]
    {
        return ...;
    }

    /* Generate header full */
    const uint32_t prefix_cmp{std::bitcast<uint32_t>(values.front())};

    std::vector<int> zcounts(values.size(), kCharBit * sizeof(T));

    /* Peform bitwise XOR over all values */
    for (size_t i=1; i < values.size(); ++i)
    {
        /* Check the common prefix */
        const uint32_t xor_prefix = prefix_cmp ^ std::bitcast<uint32_t>(values[i]);

        /* Check the leading zero count */
        const int zero_count = std::countl_zero(xor_prefix);

        /* Store the value */
        zcounts[i] = zero_count;
    }

    /* Check the minimum applicable prefix length */
    const auto min_zero_count_iter = std::min_element(zcounts.begin(), zcounts.end());

    cmc_assert(min_zero_count_iter != zcounts.end());

    /* Check the prefix length */
    const int prefix_length = *min_zero_count_iter;

    cmc_assert(prefix_length >= 0 && prefix_length <= sizeof(T) * kCharBit);

    /* Build the prefix */
    const T prefix_value = std::bitcast<T>(prefix_cmp & ((~uint32_t{0}) << (sizeof(T) * kCharBit - prefix_length)));

    /* Return the compression value */
    ....
}
#endif

template<OneByteType T, OneByteType U>
constexpr uint8_t
IntegerAddition(const T summand1, const U summand2)
{
    return std::bit_cast<uint8_t>(summand1) + std::bit_cast<uint8_t>(summand2);
}

template<TwoByteType T, TwoByteType U>
constexpr uint16_t
IntegerAddition(const T summand1, const U summand2)
{
    return std::bit_cast<uint16_t>(summand1) + std::bit_cast<uint16_t>(summand2);
}

template<FourByteType T, FourByteType U>
constexpr uint32_t
IntegerAddition(const T summand1, const U summand2)
{
    return std::bit_cast<uint32_t>(summand1) + std::bit_cast<uint32_t>(summand2);
}

template<EightByteType T, EightByteType U>
constexpr uint64_t
IntegerAddition(const T summand1, const U summand2)
{
    return std::bit_cast<uint64_t>(summand1) + std::bit_cast<uint64_t>(summand2);
}

template<OneByteType T, OneByteType U>
constexpr uint8_t
IntegerSubtraction(const T minuend, const U subtrahend)
{
    return std::bit_cast<uint8_t>(minuend) - std::bit_cast<uint8_t>(subtrahend);
}

template<TwoByteType T, TwoByteType U>
constexpr uint16_t
IntegerSubtraction(const T minuend, const U subtrahend)
{
    return std::bit_cast<uint16_t>(minuend) - std::bit_cast<uint16_t>(subtrahend);
}

template<FourByteType T, FourByteType U>
constexpr uint32_t
IntegerSubtraction(const T minuend, const U subtrahend)
{
    return std::bit_cast<uint32_t>(minuend) - std::bit_cast<uint32_t>(subtrahend);
}

template<EightByteType T, EightByteType U>
constexpr uint64_t
IntegerSubtraction(const T minuend, const U subtrahend)
{
    return std::bit_cast<uint64_t>(minuend) - std::bit_cast<uint64_t>(subtrahend);
}

template<OneByteType T, OneByteType U>
constexpr bool
BinGreaterThan(const T value, const U cmp)
{
    return (std::bit_cast<uint8_t>(value) > std::bit_cast<uint8_t>(cmp));
}

template<TwoByteType T, TwoByteType U>
constexpr bool
BinGreaterThan(const T value, const U cmp)
{
    return (std::bit_cast<uint16_t>(value) > std::bit_cast<uint16_t>(cmp));
}

template<FourByteType T, FourByteType U>
constexpr bool
BinGreaterThan(const T value, const U cmp)
{
    return (std::bit_cast<uint32_t>(value) > std::bit_cast<uint32_t>(cmp));
}

template<EightByteType T, EightByteType U>
constexpr bool
BinGreaterThan(const T value, const U cmp)
{
    return (std::bit_cast<uint64_t>(value) > std::bit_cast<uint64_t>(cmp));
}

template<OneByteType T, OneByteType U>
constexpr bool
BinLessThan(const T value, const U cmp)
{
    return (std::bit_cast<uint8_t>(value) < std::bit_cast<uint8_t>(cmp));
}

template<TwoByteType T, TwoByteType U>
constexpr bool
BinLessThan(const T value, const U cmp)
{
    return (std::bit_cast<uint16_t>(value) < std::bit_cast<uint16_t>(cmp));
}

template<FourByteType T, FourByteType U>
constexpr bool
BinLessThan(const T value, const U cmp)
{
    return (std::bit_cast<uint32_t>(value) < std::bit_cast<uint32_t>(cmp));
}

template<EightByteType T, EightByteType U>
constexpr bool
BinLessThan(const T value, const U cmp)
{
    return (std::bit_cast<uint64_t>(value) < std::bit_cast<uint64_t>(cmp));
}

template<OneByteType T, OneByteType U>
constexpr std::pair<bool, uint8_t>
ComputeIntegerResidual(const T approximation, const U real_value)
{
    /* Check if the approximaiton is greater than the real value */
    const bool is_approx_greater = BinGreaterThan(approximation, real_value);

    /* Compute the difference in inetegr arithmetic */
    const uint8_t diff = is_approx_greater ? IntegerSubtraction(approximation, real_value) : IntegerSubtraction(real_value, approximation);

    /* Return the residual */
    return std::make_pair(is_approx_greater, diff);
}

template<TwoByteType T, TwoByteType U>
constexpr std::pair<bool, uint16_t>
ComputeIntegerResidual(const T approximation, const U real_value)
{
    /* Check if the approximaiton is greater than the real value */
    const bool is_approx_greater = BinGreaterThan(approximation, real_value);

    /* Compute the difference in inetegr arithmetic */
    const uint16_t diff = is_approx_greater ? IntegerSubtraction(approximation, real_value) : IntegerSubtraction(real_value, approximation);

    /* Return the residual */
    return std::make_pair(is_approx_greater, diff);
}

template<FourByteType T, FourByteType U>
constexpr std::pair<bool, uint32_t>
ComputeIntegerResidual(const T approximation, const U real_value)
{
    /* Check if the approximaiton is greater than the real value */
    const bool is_approx_greater = BinGreaterThan(approximation, real_value);

    /* Compute the difference in inetegr arithmetic */
    const uint32_t diff = is_approx_greater ? IntegerSubtraction(approximation, real_value) : IntegerSubtraction(real_value, approximation);

    /* Return the residual */
    return std::make_pair(is_approx_greater, diff);
}

template<EightByteType T, EightByteType U>
constexpr std::pair<bool, uint64_t>
ComputeIntegerResidual(const T approximation, const U real_value)
{
    /* Check if the approximaiton is greater than the real value */
    const bool is_approx_greater = BinGreaterThan(approximation, real_value);

    /* Compute the difference in inetegr arithmetic */
    const uint64_t diff = is_approx_greater ? IntegerSubtraction(approximation, real_value) : IntegerSubtraction(real_value, approximation);

    /* Return the residual */
    return std::make_pair(is_approx_greater, diff);
}

template<typename T>
constexpr int
GetLZC(const T residual)
{
    return std::countl_zero(residual);
}


template<OneByteType T>
inline std::array<uint8_t, 1>
SerializeValueBE(const T value)
{
    /* Convert the value */
    const uint8_t converted_value = ConvertToBigEndian(std::bit_cast<uint8_t>(value));
    /* Create a byte pointer to the data */
    const uint8_t* byte_ptr = reinterpret_cast<const uint8_t*>(&converted_value);
    /* Copy the bytes into the array */
    const std::array<uint8_t, 1> serialized{*byte_ptr};
    return serialized;
}

template<TwoByteType T>
inline std::array<uint8_t, 2>
SerializeValueBE(const T value)
{
    /* Convert the value */
    const uint16_t converted_value = ConvertToBigEndian(std::bit_cast<uint16_t>(value));
    /* Create a byte pointer to the data */
    const uint8_t* byte_ptr = reinterpret_cast<const uint8_t*>(&converted_value);
    /* Copy the bytes into the array */
    const std::array<uint8_t, 2> serialized{*byte_ptr, *(byte_ptr+1)};
    return serialized;
}

template<FourByteType T>
inline std::array<uint8_t, 4>
SerializeValueBE(const T value)
{
    /* Convert the value */
    const uint32_t converted_value = ConvertToBigEndian(std::bit_cast<uint32_t>(value));
    /* Create a byte pointer to the data */
    const uint8_t* byte_ptr = reinterpret_cast<const uint8_t*>(&converted_value);
    /* Copy the bytes into the array */
    const std::array<uint8_t, 4> serialized{*byte_ptr, *(byte_ptr+1), *(byte_ptr+2), *(byte_ptr+3)};
    return serialized;
}

template<EightByteType T>
inline std::array<uint8_t, 8>
SerializeValueBE(const T value)
{
    /* Convert the value */
    const uint64_t converted_value = ConvertToBigEndian(std::bit_cast<uint64_t>(value));
    /* Create a byte pointer to the data */
    const uint8_t* byte_ptr = reinterpret_cast<const uint8_t*>(&converted_value);
    /* Copy the bytes into the array */
    const std::array<uint8_t, 8> serialized{*byte_ptr, *(byte_ptr+1), *(byte_ptr+2), *(byte_ptr+3), *(byte_ptr+4), *(byte_ptr+5), *(byte_ptr+6), *(byte_ptr+7)};
    return serialized;
}

template <OneByteType T, PointerType Iter>
inline T
DeserializeValueBE(Iter pos)
{
    /* Create a byte pointer from the iterator */
    const uint8_t* pos_ = reinterpret_cast<const uint8_t*>(&(*pos));
    /* Declare input value */
    uint8_t value;
    /* Copy the bytes over to the type */
    std::memcpy(&value, pos_, sizeof(T));
    /* Convert the value to the native endiannes */
    const T native_value = std::bit_cast<T>(ConvertBigEndianToNativeEndianness(value));
    return native_value;
}

template <TwoByteType T, PointerType Iter>
inline T
DeserializeValueBE(Iter pos)
{
    /* Create a byte pointer from the iterator */
    const uint8_t* pos_ = reinterpret_cast<const uint8_t*>(&(*pos));
    /* Declare input value */
    uint16_t value;
    /* Copy the bytes over to the type */
    std::memcpy(&value, pos_, sizeof(T));
    /* Convert the value to the native endiannes */
    const T native_value = std::bit_cast<T>(ConvertBigEndianToNativeEndianness(value));
    return native_value;
}

template <FourByteType T, PointerType Iter>
inline T
DeserializeValueBE(Iter pos)
{
    /* Create a byte pointer from the iterator */
    const uint8_t* pos_ = reinterpret_cast<const uint8_t*>(&(*pos));
    /* Declare input value */
    uint32_t value;
    /* Copy the bytes over to the type */
    std::memcpy(&value, pos_, sizeof(T));
    /* Convert the value to the native endiannes */
    const T native_value = std::bit_cast<T>(ConvertBigEndianToNativeEndianness(value));
    return native_value;
}

template <EightByteType T, PointerType Iter>
inline T
DeserializeValueBE(Iter pos)
{
    /* Create a byte pointer from the iterator */
    const uint8_t* pos_ = reinterpret_cast<const uint8_t*>(&(*pos));
    /* Declare input value */
    uint64_t value;
    /* Copy the bytes over to the type */
    std::memcpy(&value, pos_, sizeof(T));
    /* Convert the value to the native endiannes */
    const T native_value = std::bit_cast<T>(ConvertBigEndianToNativeEndianness(value));
    return native_value;
}

template <OneByteType T>
inline T
DeserializeValueBE(const uint8_t* pos)
{
    /* Declare input value */
    uint8_t value;
    /* Copy the bytes over to the type */
    std::memcpy(&value, pos, sizeof(T));
    /* Convert the value to the native endiannes */
    const T native_value = std::bit_cast<T>(ConvertBigEndianToNativeEndianness(value));
    return native_value;
}

template <TwoByteType T>
inline T
DeserializeValueBE(const uint8_t* pos)
{
    /* Declare input value */
    uint16_t value;
    /* Copy the bytes over to the type */
    std::memcpy(&value, pos, sizeof(T));
    /* Convert the value to the native endiannes */
    const T native_value = std::bit_cast<T>(ConvertBigEndianToNativeEndianness(value));
    return native_value;
}

template <FourByteType T>
inline T
DeserializeValueBE(const uint8_t* pos)
{
    /* Declare input value */
    uint32_t value;
    /* Copy the bytes over to the type */
    std::memcpy(&value, pos, sizeof(T));
    /* Convert the value to the native endiannes */
    const T native_value = std::bit_cast<T>(ConvertBigEndianToNativeEndianness(value));
    return native_value;
}

template <EightByteType T>
inline T
DeserializeValueBE(const uint8_t* pos)
{
    /* Declare input value */
    uint64_t value;
    /* Copy the bytes over to the type */
    std::memcpy(&value, pos, sizeof(T));
    /* Convert the value to the native endiannes */
    const T native_value = std::bit_cast<T>(ConvertBigEndianToNativeEndianness(value));
    return native_value;
}

template <OneByteType T>
inline void
SerializeBEToByteStream(std::vector<uint8_t>& byte_stream, const T value)
{
    /* Serialize the value correctly in big-endian */
    const uint8_t converted_value = ConvertToBigEndian(std::bit_cast<uint8_t>(value));
    /* Create a byte pointer to the data */
    const uint8_t* byte_ptr = reinterpret_cast<const uint8_t*>(&converted_value);
    /* Copy the data to the byte stream */
    std::copy_n(byte_ptr, sizeof(T), std::back_inserter(byte_stream));
}

template <TwoByteType T>
inline void
SerializeBEToByteStream(std::vector<uint8_t>& byte_stream, const T value)
{
    /* Serialize the value correctly in big-endian */
    const uint16_t converted_value = ConvertToBigEndian(std::bit_cast<uint16_t>(value));
    /* Create a byte pointer to the data */
    const uint8_t* byte_ptr = reinterpret_cast<const uint8_t*>(&converted_value);
    /* Copy the data to the byte stream */
    std::copy_n(byte_ptr, sizeof(T), std::back_inserter(byte_stream));
}

template <FourByteType T>
inline void
SerializeBEToByteStream(std::vector<uint8_t>& byte_stream, const T value)
{
    /* Serialize the value correctly in big-endian */
    const uint32_t converted_value = ConvertToBigEndian(std::bit_cast<uint32_t>(value));
    /* Create a byte pointer to the data */
    const uint8_t* byte_ptr = reinterpret_cast<const uint8_t*>(&converted_value);
    /* Copy the data to the byte stream */
    std::copy_n(byte_ptr, sizeof(T), std::back_inserter(byte_stream));
}

template <EightByteType T>
inline void
SerializeBEToByteStream(std::vector<uint8_t>& byte_stream, const T value)
{
    /* Serialize the value correctly in big-endian */
    const uint64_t converted_value = ConvertToBigEndian(std::bit_cast<uint64_t>(value));
    /* Create a byte pointer to the data */
    const uint8_t* byte_ptr = reinterpret_cast<const uint8_t*>(&converted_value);
    /* Copy the data to the byte stream */
    std::copy_n(byte_ptr, sizeof(T), std::back_inserter(byte_stream));
}

}

#endif /* !CMC_BITS_HXX */
