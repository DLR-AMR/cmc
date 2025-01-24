#ifndef CMC_BYTES_HXX
#define CMC_BYTES_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_serialized_compression_value_forward.hxx"

#include <array>
#include <cstddef>
#include <cstring>
#include <climits>

namespace cmc
{

constexpr int NIBBLE_BIT = 4;

enum Endian 
{
    /* These macros are GCC macros. Other compilers might not work. */
    Little = __ORDER_LITTLE_ENDIAN__,
    Big = __ORDER_BIG_ENDIAN__,
    Native = __BYTE_ORDER__
};

constexpr bool IsLittleEndian = (Endian::Native == Endian::Little);
constexpr bool IsBigEndian = (Endian::Native == Endian::Big);
constexpr bool IsUnsupportedEndianness = (IsLittleEndian || IsBigEndian);

constexpr int kFloat = static_cast<int>(sizeof(float));
constexpr int kDouble = static_cast<int>(sizeof(double));


constexpr uint8_t kZeroByte = 0x00;
constexpr uint8_t kLowBit = 0x01;
constexpr uint8_t kHighBit = 0x80;

constexpr int kNumBits = CHAR_BIT;

constexpr std::array<uint8_t, CHAR_BIT + 1> LowBitMask{0xFF, 0x7F, 0x3F, 0x1F, 0x0F, 0x07, 0x03, 0x01, 0x00};

constexpr int kPrefixError = -1;

/* This class stores the bytes of a given value in the native order */
template<int N>
class Prefix
{
public:
    Prefix(){
        prefix_mem_.fill((uint8_t)0U);
    };

    template<typename T> Prefix(const T& value);
    Prefix(std::array<uint8_t, N>&& serialized_value);

    ~Prefix() = default;

    Prefix(const Prefix& other) = default;
    Prefix& operator=(const Prefix& other) = default;
    Prefix(Prefix&& other) = default;
    Prefix& operator=(Prefix&& other) = default;

    uint8_t& operator[](const int index)
    {
        cmc_assert(index < N);
        return prefix_mem_[index];
    }
    const uint8_t& operator[](const int index) const
    {
        cmc_assert(index < N);
        return prefix_mem_[index];
    }

    friend Prefix operator^(const Prefix& prefix1, const Prefix& prefix2)
    {
        std::array<uint8_t, N> prefix;
        for (int i = 0; i < N; ++i)
        {
            prefix[i] = prefix1[i] ^ prefix2[i];
        }

        return prefix;
    }

    Prefix& operator^=(const Prefix& rhs)
    {
        for (int i = 0; i < N; ++i)
        {
            prefix_mem_[i] ^= rhs.prefix_mem_[i];
        }
        return *this;
    }

    int GetNumberLeadingZeros() const;
    
    int GetNumberTrailingZeros() const;

    const std::array<uint8_t, N>& GetMemoryForReading() const
    {
        return prefix_mem_;
    }

    friend class CompressionValue<N>;
private:
    std::array<uint8_t, N> prefix_mem_;
};

template<int N>
template<typename T>
Prefix<N>::Prefix(const T& value)
{
    cmc_assert(N == static_cast<int>(sizeof(value)));
    std::memcpy(this->prefix_mem_.data(), &value, sizeof(value));
};

template<int N>
Prefix<N>::Prefix(std::array<uint8_t, N>&& serialized_value)
: prefix_mem_{std::move(serialized_value)} {};

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


template<typename T>
static inline
int GetMSBytePosition(const T& value)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return sizeof(value) - 1;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return 0;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kPrefixError;
    }
}

template<typename T>
static inline
int GetLSBytePosition(const T& value)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return 0;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return sizeof(value) - 1;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kPrefixError;
    }
}



template<int N>
inline
constexpr int GetMSBByteStart([[maybe_unused]] const Prefix<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return N - 1;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return 0;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kPrefixError;
    }
}

template<int N>
inline
constexpr int GetMSBByteEnd([[maybe_unused]] const Prefix<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return 0;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return N - 1;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kPrefixError;
    }
}

static inline
void MSBByteIncrement(int& iterator)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        --iterator;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        ++iterator;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
    }
}

template<int N>
inline
bool MSBContinueIteration(const int iterator, [[maybe_unused]] const Prefix<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return iterator >= GetMSBByteEnd(prefix);
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return iterator <= GetMSBByteEnd(prefix);
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return false;
    }
}


template<int N>
inline
constexpr int GetLSBByteStart([[maybe_unused]] const Prefix<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return 0;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return N - 1;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kPrefixError;
    }
}

template<int N>
inline
constexpr int GetLSBByteEnd([[maybe_unused]] const Prefix<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return N - 1;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return 0;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return kPrefixError;
    }
}

static inline
void LSBByteIncrement(int& iterator)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        ++iterator;
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        --iterator;
    } else
    {
        cmc_err_msg("The native endianness is not supported");
    }
}

template<int N>
inline
bool LSBContinueIteration(const int iterator, [[maybe_unused]] const Prefix<N>& prefix)
{
    if constexpr (IsLittleEndian)
    {
        /* If little endian */
        return iterator <= GetLSBByteEnd(prefix);
    } else if constexpr (IsBigEndian)
    {
        /* If big endian */
        return iterator >= GetLSBByteEnd(prefix);
    } else
    {
        cmc_err_msg("The native endianness is not supported");
        return false;
    }
}


template<int N>
int
Prefix<N>::GetNumberLeadingZeros() const
{
    int num_leading_zeros = 0;
    for (int byte_id = GetMSBByteStart(*this); MSBContinueIteration(byte_id, *this); MSBByteIncrement(byte_id))
    {
        /* Check whether the full byte is zero */
        if (prefix_mem_[byte_id] == kZeroByte)
        {
            num_leading_zeros += CHAR_BIT;
        } else
        {
            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if (prefix_mem_[byte_id] & (kHighBit >> bit_index))
                {
                    /* If true, the prefix does hold not a zero at this position */
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
Prefix<N>::GetNumberTrailingZeros() const
{
    int num_trailing_zeros = 0;
    for (int byte_id = GetLSBByteStart(*this); LSBContinueIteration(byte_id, *this); LSBByteIncrement(byte_id))
    {
        /* Check whether the full byte is zero */
        if (prefix_mem_[byte_id] == kZeroByte)
        {
            num_trailing_zeros += CHAR_BIT;
        } else
        {
            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if (prefix_mem_[byte_id] & (kLowBit << bit_index))
                {
                    /* If true, the prefix does hold a zero at this position */
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

template <typename T>
std::array<uint8_t, sizeof(T)>
SerializeValue(const T& value, const Endian desired_endianness = Endian::Big)
{
    //cmc_debug_msg("data to be serialized: ", value);
    std::array<uint8_t, sizeof(T)> serialized_value;
    Prefix<sizeof(T)> value_(value);

    /* Assign the bytes in the desired order */
    switch (desired_endianness)
    {
        case Endian::Big:
            for (int byte_id = GetMSBByteStart(value_), index = 0; MSBContinueIteration(byte_id, value_); MSBByteIncrement(byte_id), ++index)
            {
                serialized_value[index] = value_[byte_id];
            }
        break;
        case Endian::Little:
            for (int byte_id = GetLSBByteStart(value_), index = 0; LSBContinueIteration(byte_id, value_); LSBByteIncrement(byte_id), ++index)
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

#endif /* !CMC_BYTES_HXX */
