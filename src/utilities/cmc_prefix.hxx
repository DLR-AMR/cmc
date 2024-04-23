#ifndef CMC_PREFIX_HXX
#define CMC_PREFIX_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.h"

#include <array>
#include <climits>
#include <cstring>
#include <bitset>

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

constexpr int PREFIX_ERROR = -1;

constexpr std::array<uint8_t, CHAR_BIT> LowBitMask{0xFF, 0x7F, 0x3F, 0x1F, 0x0F, 0x07, 0x03, 0x01};

template<int N>
class CompressionValue;

template<int N>
class Prefix;

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
        return PREFIX_ERROR;
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
        return PREFIX_ERROR;
    }
}

template<int N>
template<typename T>
Prefix<N>::Prefix(const T& value)
{
    if (N != static_cast<int>(sizeof(value)))
    {
        cmc_debug_msg("Hier gehts nicht.");// N ist = ", N, " und value ist: ", value, " mit groesse : ", sizeof(T));
    }
    cmc_assert(N == static_cast<int>(sizeof(value)));
    std::memcpy(this->prefix_mem_.data(), &value, sizeof(value));
};

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
        return PREFIX_ERROR;
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
        return PREFIX_ERROR;
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
        return PREFIX_ERROR;
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
        return PREFIX_ERROR;
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
Prefix<N>::Prefix(std::array<uint8_t, N>&& serialized_value)
: prefix_mem_{std::move(serialized_value)} {};

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
                    /* If true, the prefix does hold a zero at this position */
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

template<int N>
struct CommonPrefix 
{
    Prefix<N> prefix;
    int bit_length{0};
};

template<int N>
inline CommonPrefix<N>
GetCommonPrefixOnly(const Prefix<N>& prefix1, const Prefix<N>& prefix2)
{
    CommonPrefix<N> cprefix;

    for (int byte_id = GetMSBByteStart(prefix1); MSBContinueIteration(byte_id, prefix1); MSBByteIncrement(byte_id))
    {
        if (prefix1[byte_id] == prefix2[byte_id])
        {
            /* If a whole byte coincides */
            cprefix.prefix[byte_id] = prefix2[byte_id];
            cprefix.bit_length += CHAR_BIT;
        } else
        {
            /* Only copy the equal part of the byte */

            const uint8_t xor_bits = prefix1[byte_id] ^ prefix2[byte_id];

            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if (xor_bits & (kHighBit >> bit_index))
                {
                    /* The first unequal bit has been found */
                    cprefix.prefix[byte_id] = prefix1[byte_id] & ~(xor_bits | LowBitMask[bit_index]);
                    return cprefix;
                } else
                {
                    ++cprefix.bit_length;
                }
            }
        }
    }

    return cprefix;
}

constexpr int kNoPrefixIndicationBit = -1;

template<int N>
CompressionValue<N>
GetCommonPrefix(const CompressionValue<N>& value1, const CompressionValue<N>& value2);

template<int N>
class CompressionValue
{
public:
    CompressionValue()
    : prefix_(Prefix<N>()), front_bit_{0}, trail_bit_{N * CHAR_BIT}{};
    template<typename T> CompressionValue(const T& value);
    CompressionValue(std::array<uint8_t, N>&& serialized_value);
    CompressionValue(std::array<uint8_t, N>&& serialized_value, const int trail_bit);
    CompressionValue(Prefix<N>&& prefix, const int front_bit, const int tail_bit)
    : prefix_{std::move(prefix)}, front_bit_{front_bit}, trail_bit_{tail_bit} {};
    CompressionValue(const std::vector<uint8_t>& serialized_prefix, const int num_bits);

    ~CompressionValue() = default;

    CompressionValue(const CompressionValue& other) = default;
    CompressionValue& operator=(const CompressionValue& other) = default;
    CompressionValue(CompressionValue&& other) = default;
    CompressionValue& operator=(CompressionValue&& other) = default;

    void ToggleTailUntilNextUnsetBit();

    void ClearNextSetBitFromTail();

    int GetNumberTrailingZeros() const;
    int GetNumberLeadingZeros() const;
    void UpdateTrailBitCount();
    void UpdateFrontBitCount();
    void SetFrontBit(const int front_bit);
    int GetTrailBit() const;
    int GetFrontBit() const;
    int GetCountOfSignificantBits() const;
    bool IsEmpty() const;
    std::vector<uint8_t> GetSignificantBitsInBigEndianOrdering(const int start_offset = 0) const;
    std::vector<uint8_t> GetSignificantBitsInBigEndianOrdering_Offset_4_Stride_12() const;
    void ApplyPrefix(const std::vector<uint8_t>& serialized_prefix, const int num_bits);

    template<typename T> 
    auto ReinterpretDataAs() const
        -> std::enable_if_t<sizeof(T) == N, T>
    {
        T value;
        std::memcpy(&value, prefix_.prefix_mem_.data(), sizeof(T));
        return value;
    }

    friend CompressionValue
    GetCommonPrefix <> (const CompressionValue<N>& value1, const CompressionValue<N>& value2);

    uint8_t& operator[](const int index)
    {
        cmc_assert(index < N);
        return prefix_.prefix_mem_[index];
    }
    const uint8_t& operator[](const int index) const
    {
        cmc_assert(index < N);
        return prefix_.prefix_mem_[index];
    }

    const std::array<uint8_t, N>& GetMemoryForReading() const
    {
        return prefix_.GetMemoryForReading();
    }

private:
    Prefix<N> prefix_;
    int front_bit_{0};
    int trail_bit_{N * CHAR_BIT};
};

template<int N>
CompressionValue<N>::CompressionValue(const std::vector<uint8_t>& serialized_prefix, const int num_bits)
{
    cmc_assert(N >= serialized_prefix.size());

    int bits_set = 0;
    for (int byte_id = GetMSBByteStart(prefix_), index = 0; MSBContinueIteration(byte_id, prefix_); MSBByteIncrement(byte_id), ++index)
    {
        prefix_[byte_id] = serialized_prefix[index];
        bits_set += CHAR_BIT;
        if (bits_set >= num_bits) {break;}
    }

    trail_bit_ = N * CHAR_BIT - num_bits;
};


/* This functions needs the serialized prefix to be aligned at the high bits ( e.g. four bit prefix: 0x(p1 p2 p3 p4 0 0 0 0)b */
template<int N>
void
CompressionValue<N>::ApplyPrefix(const std::vector<uint8_t>& serialized_prefix, const int num_bits)
{
    if (num_bits <= 0) {return;}

    int bits_written = 0;

    const int offset_shifting_ = (trail_bit_  - static_cast<int>(trail_bit_ / CHAR_BIT) * CHAR_BIT);
    const int offset_shifting = (offset_shifting_ > 0 ? CHAR_BIT - offset_shifting_ : 0);

    for (int byte_id = GetMSBByteStart(prefix_), iter = 0, vec_index = 0;
         MSBContinueIteration(byte_id, prefix_);
         MSBByteIncrement(byte_id), ++iter)
    {
        /* Iterate until we are in the byte holding the bit after the current tail bit */
        if (trail_bit_ <= (N - iter - 1) * CHAR_BIT) {continue;}

        if (offset_shifting > 0 && bits_written >= CHAR_BIT - offset_shifting)
        {
            prefix_[byte_id] |= (serialized_prefix[vec_index - 1] << (CHAR_BIT - offset_shifting));
            bits_written += offset_shifting;
            if (bits_written >= num_bits) {break;}
        }

        prefix_[byte_id] |= (serialized_prefix[vec_index] >> offset_shifting);
        bits_written += CHAR_BIT - offset_shifting;
        if (bits_written >= num_bits) {break;}
        ++vec_index;
    }

    trail_bit_ -= num_bits;
}

template<int N>
CompressionValue<N>
GetCommonPrefix(const CompressionValue<N>& value1, const CompressionValue<N>& value2)
{
    CompressionValue<N> cprefix;

    const int max_trail_bit = std::max(value1.trail_bit_, value2.trail_bit_);
    const int maximum_prefix_length = N * CHAR_BIT - max_trail_bit;
    int current_bit_length = 0;

    std::array<uint8_t, N> prefix_bits;
    prefix_bits.fill(0);

    Prefix<N> xor_bits = value1.prefix_ ^ value2.prefix_;

    for (int byte_id = GetMSBByteStart(value1.prefix_); MSBContinueIteration(byte_id, value1.prefix_); MSBByteIncrement(byte_id))
    {
        if (xor_bits[byte_id] == 0)
        {
            if (current_bit_length + CHAR_BIT <= maximum_prefix_length)
            {
                prefix_bits[byte_id] = value1[byte_id];
                current_bit_length += CHAR_BIT;
            } else
            {
                const int bit_accessor = maximum_prefix_length - current_bit_length;
                prefix_bits[byte_id] = value1[byte_id] & ~(xor_bits[byte_id] | LowBitMask[bit_accessor]);

                return CompressionValue<N>(std::move(prefix_bits), max_trail_bit);
            }
        } else
        {
            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if ((xor_bits[byte_id] & (kHighBit >> bit_index)) || maximum_prefix_length <= current_bit_length)
                {
                    /* The first unequal bit has been found */
                    prefix_bits[byte_id] = value1[byte_id] & ~(xor_bits[byte_id] | LowBitMask[bit_index]);

                    return CompressionValue<N>(std::move(prefix_bits), N * CHAR_BIT - current_bit_length);
                } else
                {
                    ++current_bit_length;
                }
            }
        }
    }

    return (value1.trail_bit_ >= value2.trail_bit_ ? value1 : value2);
}


template<int N>
bool
CompressionValue<N>::IsEmpty() const
{
    if (front_bit_ + trail_bit_ >= N * CHAR_BIT)
    {
        return true;
    } else
    {
        return false;
    }
}


template<int N>
void
CompressionValue<N>::SetFrontBit(const int front_bit)
{
    cmc_assert(CHAR_BIT * N >= front_bit + trail_bit_);
    front_bit_ = front_bit;
}

template<int N>
int
CompressionValue<N>::GetTrailBit() const
{
    return trail_bit_;
}

template<int N>
int
CompressionValue<N>::GetFrontBit() const
{
    return front_bit_;
}

template<int N>
void
CompressionValue<N>::ToggleTailUntilNextUnsetBit()
{
    int byte_id = GetLSBByteStart(prefix_);
    /* Get to the current start byte */
    const int byte_increments = static_cast<int>(trail_bit_ / CHAR_BIT);
    if (byte_increments == N)
    {
        /* The value has already been completely worked out */
        return;
    }

    for (int increment = 0; increment < byte_increments; ++increment)
    {
        LSBByteIncrement(byte_id);
    }

    bool has_encountered_one = false;

    int bit_accessor = trail_bit_ - byte_increments * CHAR_BIT;

    for (int byte_index = byte_id; LSBContinueIteration(byte_index, prefix_); LSBByteIncrement(byte_index))
    {
        for (int bit_index = bit_accessor; bit_index < CHAR_BIT; ++bit_index)
        {
            if (prefix_[byte_index] & (kLowBit << bit_index))
            {
                /* If there is a one at this position */
                has_encountered_one = true;
                /* We toggle the bit (unset it) and continue iterating until we have found leading zero */
                prefix_[byte_index] &= ~(kLowBit << bit_index);
                continue;
            } else
            {
                /* If there is a zero at this position */
                if (has_encountered_one)
                {
                    /* If we have aready seen a one before, we will set the current bit to one and are finished */
                    prefix_[byte_index] |= (kLowBit << bit_index);
                    trail_bit_ += bit_index - bit_accessor;
                    return;
                } else
                {
                    /* If we have not seen a one before, we just keep iterating until we will find one */
                    continue;
                }
            }
        }

        /* Update the trainling end */
        trail_bit_ += CHAR_BIT - bit_accessor;

        /* Reset the bit_accesor, for the next byte */
        bit_accessor = 0;
    }
}

template<int N>
template<typename T>
CompressionValue<N>::CompressionValue(const T& value)
: prefix_(value), front_bit_{0}, trail_bit_{0} {};

template<int N>
CompressionValue<N>::CompressionValue(std::array<uint8_t, N>&& serialized_value)
: prefix_{std::move(serialized_value)} {};

template<int N>
CompressionValue<N>::CompressionValue(std::array<uint8_t, N>&& serialized_value, const int trail_bit)
: prefix_{std::move(serialized_value)}, trail_bit_{trail_bit}{};

template<int N>
inline
int
CompressionValue<N>::GetNumberTrailingZeros() const
{
    return prefix_.GetNumberTrailingZeros();
}

template<int N>
inline
int
CompressionValue<N>::GetNumberLeadingZeros() const
{
    return prefix_.GetNumberLeadingZeros();
}

template<int N>
inline
void
CompressionValue<N>::UpdateTrailBitCount()
{
    trail_bit_ = prefix_.GetNumberTrailingZeros();
}

template<int N>
inline
void
CompressionValue<N>::UpdateFrontBitCount()
{
    front_bit_ = prefix_.GetNumberLeadingZeros();
}

template<int N>
void
CompressionValue<N>::ClearNextSetBitFromTail()
{
    int byte_id = GetLSBByteStart(prefix_);
    /* Get to the current start byte */
    const int byte_increments = static_cast<int>(trail_bit_ / CHAR_BIT);
    if (byte_increments == N)
    {
        /* The value has already been completely worked out */
        return;
    }
    for (int increment = 0; increment < byte_increments; ++increment)
    {
        LSBByteIncrement(byte_id);
    }

    int bit_accessor = trail_bit_ - byte_increments * CHAR_BIT;

    for (int byte_index = byte_id; LSBContinueIteration(byte_index, prefix_); LSBByteIncrement(byte_index))
    {
        for (int bit_index = bit_accessor; bit_index < CHAR_BIT; ++bit_index)
        {
            if (prefix_[byte_index] & (kLowBit << bit_index))
            {
                /* If there is a one, clear this position and we are finished, else continue */
                prefix_[byte_index] &= ~(kLowBit << bit_index);
                trail_bit_ += bit_index - bit_accessor;   
                return;
            }
        }

        /* Update the trainling end */
        trail_bit_ += CHAR_BIT - bit_accessor;

        /* Reset the bit accesor */
        bit_accessor = 0;
    }
}

template<int N>
inline
int
CompressionValue<N>::GetCountOfSignificantBits() const
{
    //cmc_debug_msg("front bit: ", front_bit_, " und trail_bit: ", trail_bit_);
    cmc_assert(front_bit_ + trail_bit_ <= N * CHAR_BIT);
    return N * CHAR_BIT - front_bit_ - trail_bit_;
}


template<int N>
std::vector<uint8_t>
CompressionValue<N>::GetSignificantBitsInBigEndianOrdering(const int bit_start_offset) const
{
    cmc_assert(bit_start_offset >= 0 && bit_start_offset < CHAR_BIT);

    const int num_bits = GetCountOfSignificantBits();
    //if (num_bits == 0) {cmc_err_msg("Num significant bits is zero (in setsignificant bits inbigendian...)");}
    const int num_bytes = (num_bits + bit_start_offset) / CHAR_BIT + (((num_bits + bit_start_offset) % CHAR_BIT) != 0 ? 1 : 0);
    //cmc_debug_msg("NumyBytes: ", num_bytes);

    if (num_bits == 0)
    {
        /* In case there are no significant bits, we are returning an empty vector */
        return std::vector<uint8_t>();
    }
    
    std::vector<uint8_t> bytes(num_bytes, 0);

    const int front_bit_position = front_bit_ - static_cast<int>(front_bit_ / CHAR_BIT) * CHAR_BIT;

    int bits_written = 0;

    const int offset_shifting = front_bit_position - bit_start_offset;

    bool do_offset_shifiting = true;

    const int bit_shift_to_front = (offset_shifting < 0 ? std::abs(offset_shifting) : offset_shifting);
    const int bit_shift_to_back = (offset_shifting < 0 ? CHAR_BIT - std::abs(offset_shifting) : CHAR_BIT - offset_shifting);

    const int front_bits_written = (offset_shifting < 0 ? bit_shift_to_front : CHAR_BIT - bit_shift_to_front);
    const int tail_bits_written = (offset_shifting < 0 ? bit_shift_to_back : CHAR_BIT - bit_shift_to_back);

    int vec_offset = (offset_shifting < 0 ? 1 : 0);

    for (int byte_id = GetMSBByteStart(prefix_), iter = 0, vec_index = 0;
         MSBContinueIteration(byte_id, prefix_);
         MSBByteIncrement(byte_id), ++iter)
    {
        /* Iterate until we are in the byte holding the start bit */
        if (front_bit_ >= (iter + 1) * CHAR_BIT) {continue;}

        if (offset_shifting < 0 && do_offset_shifiting)
        {
            /* If not all significant bits have space within the first byte */
            bytes[vec_index] = (prefix_[byte_id] >> std::abs(offset_shifting));
            bits_written = CHAR_BIT - bit_start_offset;
            do_offset_shifiting = false;
            if (bits_written >= num_bits) {break;}
        }

        /* Fill the back of the previous byte, if the start_bit has an offset */
        if (bit_shift_to_front > 0 && ((offset_shifting >= 0 && bits_written > 0) ||
                                       (offset_shifting < 0 && bits_written > CHAR_BIT - bit_start_offset)))
        {
            bytes[vec_index + vec_offset - 1] |= (prefix_[byte_id] >> (bit_shift_to_back));
            bits_written += tail_bits_written;
            if (bits_written >= num_bits) {break;}
        }

        bytes[vec_index + vec_offset] |= (prefix_[byte_id] << (bit_shift_to_front));
        bits_written += front_bits_written;    
        if (bits_written >= num_bits) {break;}

        ++vec_index;
    }

    return bytes;
}

template<int N>
std::vector<uint8_t>
CompressionValue<N>::GetSignificantBitsInBigEndianOrdering_Offset_4_Stride_12() const
{
    const int offset = 4;
    const int max_contiguous_bit_count = 12;
    if (max_contiguous_bit_count >= GetCountOfSignificantBits())
    {
        return GetSignificantBitsInBigEndianOrdering(offset);
    }

    const uint8_t nullify_front_bits = 0x0F;

    const std::vector<uint8_t> significant_bits = GetSignificantBitsInBigEndianOrdering();

    const int num_significant_bits = GetCountOfSignificantBits();

    const int num_bytes = (num_significant_bits + offset * ((num_significant_bits / max_contiguous_bit_count) + 1)) / CHAR_BIT + ((((num_significant_bits + offset * ((num_significant_bits / max_contiguous_bit_count) + 1))) % CHAR_BIT) != 0 ? 1 : 0);
    
    std::vector<uint8_t> strided_vec;
    strided_vec.reserve(num_bytes);

    int bit_count = 0;

    const int processed_bits = 4;

    int shift_to_back = processed_bits;
    int shift_to_front = processed_bits;

    bool switch_position = true;
    strided_vec.push_back(0);

    for (int index = 0; index < significant_bits.size(); ++index)
    {
        if (bit_count + processed_bits > max_contiguous_bit_count)
        {
            /* Offset the current bits */
            strided_vec.push_back(significant_bits[index] >> shift_to_back);
            bit_count = processed_bits;
            switch_position = !switch_position;
        } else
        {
            if (switch_position)
            {
                strided_vec.back() |= significant_bits[index] >> shift_to_back;
                bit_count += processed_bits;
            } else
            {
                strided_vec.push_back(significant_bits[index] >> shift_to_back);
                bit_count += processed_bits;
            }
        }

        if (bit_count + processed_bits > max_contiguous_bit_count)
        {
            /* Offset the current bits */
            strided_vec.push_back((significant_bits[index] & nullify_front_bits));
            bit_count = processed_bits;
            switch_position = !switch_position;
        } else
        {
            if (switch_position)
            {
                strided_vec.push_back(significant_bits[index] << shift_to_front);
                bit_count += processed_bits;   
            } else
            {
                strided_vec.back() |= (significant_bits[index] & nullify_front_bits);
                bit_count += processed_bits;
            }
        }
    }

    return strided_vec;
}

}

#endif /* !CMC_PREFIX_HXX */
