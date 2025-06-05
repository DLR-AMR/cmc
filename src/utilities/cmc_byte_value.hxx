#ifndef CMC_BYTE_VALUE_HXX
#define CMC_BYTE_VALUE_HXX

#include "utilities/cmc_endian.hxx"
#include "utilities/cmc_byte_iteration.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <array>
#include <cstddef>
#include <cstring>
#include <climits>

namespace cmc
{

constexpr int kFloatBytes = static_cast<int>(sizeof(float));
constexpr int kDoubleBytes = static_cast<int>(sizeof(double));


constexpr uint8_t kZeroByte = 0x00;
constexpr uint8_t kOneByte = 0xFF;
constexpr uint8_t kLowBit = 0x01;
constexpr uint8_t kHighBit = 0x80;

constexpr int kNumBitsPerByte = CHAR_BIT;

constexpr std::array<uint8_t, CHAR_BIT + 1> LowBitMask{0xFF, 0x7F, 0x3F, 0x1F, 0x0F, 0x07, 0x03, 0x01, 0x00};

constexpr int kByteError = -1;

/* This class stores the bytes of a given value in the native order */
template<int N>
class Serialized
{
public:
    Serialized(){
        bytes_.fill((uint8_t)0U);
    };

    template<typename T> Serialized(const T& value);
    Serialized(std::array<uint8_t, N>&& serialized_value);

    ~Serialized() = default;

    Serialized(const Serialized& other) = default;
    Serialized& operator=(const Serialized& other) = default;
    Serialized(Serialized&& other) = default;
    Serialized& operator=(Serialized&& other) = default;

    uint8_t& operator[](const int index)
    {
        cmc_assert(index < N);
        return bytes_[index];
    }
    const uint8_t& operator[](const int index) const
    {
        cmc_assert(index < N);
        return bytes_[index];
    }

    friend Serialized operator^(const Serialized& value1, const Serialized& value2)
    {
        std::array<uint8_t, N> val;
        for (int i = 0; i < N; ++i)
        {
            val[i] = value1[i] ^ value2[i];
        }

        return val;
    }

    Serialized& operator^=(const Serialized& rhs)
    {
        for (int i = 0; i < N; ++i)
        {
            bytes_[i] ^= rhs.bytes_[i];
        }
        return *this;
    }

    int GetNumberLeadingZeros() const;
    
    int GetNumberTrailingZeros() const;

    const std::array<uint8_t, N>& GetMemoryForReading() const
    {
        return bytes_;
    }


private:
    std::array<uint8_t, N> bytes_;
};

template<int N>
template<typename T>
Serialized<N>::Serialized(const T& value)
{
    cmc_assert(N == static_cast<int>(sizeof(value)));
    std::memcpy(this->bytes_.data(), &value, sizeof(value));
};

template<int N>
Serialized<N>::Serialized(std::array<uint8_t, N>&& serialized_value)
: bytes_{std::move(serialized_value)} {};


template<int N>
int
Serialized<N>::GetNumberLeadingZeros() const
{
    int num_leading_zeros = 0;
    for (int byte_id = GetMSBByteStart<N>(*this); MSBContinueIteration<N>(byte_id, *this); MSBByteIncrement(byte_id))
    {
        /* Check whether the full byte is zero */
        if (bytes_[byte_id] == kZeroByte)
        {
            num_leading_zeros += CHAR_BIT;
        } else
        {
            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if (bytes_[byte_id] & (kHighBit >> bit_index))
                {
                    /* If true, the value does hold not a zero at this position */
                    return num_leading_zeros;
                } else
                {
                    ++num_leading_zeros;
                }
            }
        }
    }
    return num_leading_zeros;
}

template<int N>
int
Serialized<N>::GetNumberTrailingZeros() const
{
    int num_trailing_zeros = 0;
    for (int byte_id = GetLSBByteStart<N>(*this); LSBContinueIteration<N>(byte_id, *this); LSBByteIncrement(byte_id))
    {
        /* Check whether the full byte is zero */
        if (bytes_[byte_id] == kZeroByte)
        {
            num_trailing_zeros += CHAR_BIT;
        } else
        {
            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if (bytes_[byte_id] & (kLowBit << bit_index))
                {
                    /* If true, the value does hold a zero at this position */
                    return num_trailing_zeros;
                } else
                {
                    ++num_trailing_zeros;
                }
            }
        }
    }
    return num_trailing_zeros;
}

inline uint8_t
GetBitMaskToClearBit(const int index)
{
    return ~(kLowBit << index);
}

inline bool
CheckIfBitIsSet(const uint8_t byte, const int bit_index)
{
    return (byte & (kLowBit << bit_index) ? true : false);
}

template <typename T>
std::array<uint8_t, sizeof(T)>
SerializeValue(const T& value, const Endian desired_endianness = Endian::Big)
{
    std::array<uint8_t, sizeof(T)> serialized_value;
    Serialized<sizeof(T)> value_(value);

    /* Assign the bytes in the desired order */
    switch (desired_endianness)
    {
        case Endian::Big:
            for (int byte_id = GetMSBByteStart<sizeof(T)>(), index = 0; MSBContinueIteration<sizeof(T)>(byte_id); MSBByteIncrement(byte_id), ++index)
            {
                serialized_value[index] = value_[byte_id];
            }
        break;
        case Endian::Little:
            for (int byte_id = GetLSBByteStart<sizeof(T)>(), index = 0; LSBContinueIteration<sizeof(T)>(byte_id); LSBByteIncrement(byte_id), ++index)
            {
                serialized_value[index] = value_[byte_id];
            }
        break;
        default:
            cmc_err_msg("The given endianness is not recognized.");
    }

    return serialized_value;
}

template <typename T, typename Iter>
std::array<uint8_t, sizeof(T)>
DeserializeValue(Iter pos, const Endian endianness_of_bytes)
{
    std::array<uint8_t, sizeof(T)> serialized_value;

    /* Assign the bytes compliant to the native order */
    if (Endian::Native == endianness_of_bytes)
    {
        /* If the endianness_of_bytes of the serialized value equals the native format, we can just copy the values over */
        std::copy_n(pos, sizeof(T), serialized_value.begin());
    } else
    {
        /* If the endianness_of_bytes does not equal the native endianness, the byte sequence needs to be reversed */
        const auto end_iter = pos + sizeof(T);
        int index = sizeof(T) - 1;
        for (auto iter = pos; iter != end_iter; ++iter, --index)
        {
            serialized_value[index] = *iter;
        }
    }
    return serialized_value;
}

template <typename T, typename Iter>
T
ReconstructSerializedValue(Iter pos, const Endian endianness_of_bytes)
{
    /* Deserialize the bytes in the native order */
    const std::array<uint8_t, sizeof(T)> deserialized_memory = DeserializeValue<T>(pos, endianness_of_bytes);

    /* Copy the memory to an actual type */
    T reconstructed_value;
    std::memcpy(&reconstructed_value, deserialized_memory.data(), sizeof(T));

    return reconstructed_value;
}

}

#endif /* !CMC_BYTE_VALUE_HXX */
