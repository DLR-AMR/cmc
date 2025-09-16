#ifndef CMC_BYTE_COMPRESSION_VALUES_HXX
#define CMC_BYTE_COMPRESSION_VALUES_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_vector_view.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_byte_value.hxx"
#include "utilities/cmc_byte_iteration.hxx"

#include <array>
#include <vector>
#include <climits>
#include <cstring>
#include <bitset>
#include <cstddef>

namespace cmc
{

template<int N>
class SerializedCompressionValue;

template<typename T>
std::vector<T> ConvertCompressionValues(const VectorView<SerializedCompressionValue<sizeof(T)>>& compression_values);

template<int N>
struct SignificantBitsIndicator
{
    uint8_t front_bit_{0};
    uint8_t tail_bit_{N * CHAR_BIT};
};

template<int N>
class SerializedCompressionValue
{
public:
    SerializedCompressionValue(){
        bytes_.fill(0);
    }
    template<typename T> SerializedCompressionValue(const T& value);
    SerializedCompressionValue(std::array<uint8_t, N>&& serialized_value);
    SerializedCompressionValue(std::array<uint8_t, N>&& serialized_value, const uint8_t tail_bit);
    SerializedCompressionValue(const std::vector<uint8_t>& serialized_prefix, const int num_bits);

    ~SerializedCompressionValue() = default;

    SerializedCompressionValue(const SerializedCompressionValue& other) = default;
    SerializedCompressionValue& operator=(const SerializedCompressionValue& other) = default;
    SerializedCompressionValue(SerializedCompressionValue&& other) = default;
    SerializedCompressionValue& operator=(SerializedCompressionValue&& other) = default;

    void ToggleTailUntilNextUnsetBit();

    void ClearNextSetBitFromTail();

    int GetNumberTrailingZeros() const;
    int GetNumberLeadingZeros() const;
    int GetNumberTrailingOnes() const;
    void UpdateTailBitCount();
    void UpdateFrontBitCount();
    void SetFrontBit(const uint8_t front_bit);
    void SetTailBit(const uint8_t tail_bit);
    uint8_t GetTailBit() const;
    uint8_t GetFrontBit() const;
    int GetCountOfSignificantBits() const;
    bool IsEmpty() const;
    bool IsExactlyEmpty() const;
    std::vector<uint8_t> GetSignificantBitsInBigEndianOrdering() const;
    void ApplySuffix(const std::vector<uint8_t>& serialized_suffix, const int num_bits);

    //void AddIntegerResidual(const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits);
    //void AddXORResidualWithoutImplicitOneBit(const uint32_t lzc, const std::vector<uint8_t>& residual_bits);
    void _TestPCP4ReverseXORResidual(const SerializedCompressionValue<N>& residual);
    int GatherEqualTailBitCount() const;
    void TrimEqualTailSequence();
    void RevertEqualTailSequence();
    int GetSizeBits() const {return N * CHAR_BIT;}
    bool IsTailBitValueSet() const;

    template<typename T> 
    auto ReinterpretDataAs() const
        -> std::enable_if_t<sizeof(T) == N, T>
    {
        T value;
        std::memcpy(&value, bytes_.data(), sizeof(T));
        return value;
    }

    friend SerializedCompressionValue
    GetCommonPrefix <> (const SerializedCompressionValue<N>& value1, const SerializedCompressionValue<N>& value2);

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

    const std::array<uint8_t, N>& GetMemoryForReading() const
    {
        return bytes_;
    }

    bool IsBitSetAtPosition(const int bit_position) const;
    void ToggleBitAtPosition(const int bit_position);
    void ToggleBitAtPositionAndResetPreviousBits(const int bit_position);


    int GetTypeNumBits() const {return N * CHAR_BIT;}

    SerializedCompressionValue& operator^=(const SerializedCompressionValue& rhs)
    {
        for (int i = 0; i < N; ++i)
        {
            bytes_[i] ^= rhs.bytes_[i];
        }
        return *this;
    }

    void NullifyNonSignificantFrontBits();
    int GetLeadingZeroCountInSignificantBits() const;

    void PerformIntegerSubtraction(const SerializedCompressionValue& residual);
    void PerformIntegerAddition(const SerializedCompressionValue& residual);

    friend SerializedCompressionValue operator^(const SerializedCompressionValue& val1, const SerializedCompressionValue& val2)
    {
        SerializedCompressionValue xor_val = val1;
        for (int i = 0; i < N; ++i)
        {
            xor_val[i] ^= val2.bytes_[i];
        }
        return xor_val;
    };
private:
    std::array<uint8_t, N> bytes_{};
    SignificantBitsIndicator<N> indicators_;
};

template<int N>
inline
SerializedCompressionValue<N>
GetNullCompressionValue()
{
    SerializedCompressionValue<N> null_value;
    null_value.SetTailBit(uint8_t{0});
    return null_value;
}

template<int N>
SerializedCompressionValue<N>::SerializedCompressionValue(const std::vector<uint8_t>& serialized_prefix, const int num_bits)
{
    cmc_assert(N >= serialized_prefix.size());

    int bits_set = 0;
    for (int byte_id = GetMSBByteStart<N>(), index = 0; MSBContinueIteration<N>(byte_id); MSBByteIncrement(byte_id), ++index)
    {
        bytes_[byte_id] = serialized_prefix[index];
        bits_set += CHAR_BIT;
        if (bits_set >= num_bits) {break;}
    }

    indicators_.tail_bit_ = N * CHAR_BIT - num_bits;
};

template<int N>
int
SerializedCompressionValue<N>::GetLeadingZeroCountInSignificantBits() const
{
    //TODO: Only iterate until tail bit is reached

    if (this->IsEmpty()) {return 0;}
    cmc_assert(indicators_.tail_bit_ < N * CHAR_BIT);
    
    int num_leading_zeros = 0;
    int front_bit = indicators_.front_bit_;

    int counter = indicators_.front_bit_;

    for (int byte_id = GetMSBByteStart<N>(), iter = 0;
         MSBContinueIteration<N>(byte_id);
         MSBByteIncrement(byte_id), ++iter)
    {
        if (counter + indicators_.tail_bit_ >= N * CHAR_BIT)
        {
            return num_leading_zeros;
        }
        /* Iterate until we are in the byte holding the bit after the current front bit */
        if (front_bit >= (iter + 1) * CHAR_BIT) {continue;}

        for (int bit_index = front_bit % CHAR_BIT; bit_index < CHAR_BIT; ++bit_index, ++counter)
        {
            if (bytes_[byte_id] & (kHighBit >> bit_index))
            {
                /* If true, the prefix does hold a zero at this position */
                return num_leading_zeros;
            } else
            {
                /* In this case a zero is present at the given position */
                ++num_leading_zeros;
            }
        }

        /* We reset the front bit in order to start at the beginning of the next byte */
        front_bit = 0;
    }
    return num_leading_zeros;
}

template<int N>
void
SerializedCompressionValue<N>::NullifyNonSignificantFrontBits()
{
    /* Get the number of bits for this prefix */
    const int num_bits = GetCountOfSignificantBits();

    if (num_bits == N * CHAR_BIT)
    {
        /* In case there are no non-significant bits, we return immediately */
        return;
    }

    int front_bits_to_nullify = indicators_.front_bit_;

    /* Nullify the first part of the prefix */
    for (int byte_id = GetMSBByteStart<N>();
         MSBContinueIteration<N>(byte_id);
         MSBByteIncrement(byte_id))
    {
        if (front_bits_to_nullify <= 0) {break;}

        if (front_bits_to_nullify >= CHAR_BIT)
        {
            /* Nullify the whole byte */
            bytes_[byte_id]  = uint8_t{0};
            front_bits_to_nullify -= CHAR_BIT;
        } else 
        {
            const uint8_t bit_mask = LowBitMask[front_bits_to_nullify];
            bytes_[byte_id] &= bit_mask;
            break;
        }
    }

}


/* This functions needs the serialized suffix to be aligned at the high bits ( e.g. four bit prefix: 0b(p1 p2 p3 p4 0 0 0 0) */
template<int N>
void
SerializedCompressionValue<N>::ApplySuffix(const std::vector<uint8_t>& serialized_suffix, const int num_bits)
{
    if (num_bits <= 0) {return;}

    cmc_assert(indicators_.tail_bit_ > 0);
    cmc_assert(static_cast<int>(serialized_suffix.size() * CHAR_BIT) >= num_bits);

    int bits_written = 0;

    const int intern_trail_bit = indicators_.tail_bit_ - 1;
    const int offset_shifting_ = (indicators_.tail_bit_ - static_cast<int>(indicators_.tail_bit_ / CHAR_BIT) * CHAR_BIT);
    const int offset_shifting = (offset_shifting_ > 0 ? CHAR_BIT - offset_shifting_ : 0);

    for (int byte_id = GetMSBByteStart<N>(), iter = 0, vec_index = 0;
         MSBContinueIteration<N>(byte_id);
         MSBByteIncrement(byte_id), ++iter)
    {
        /* Iterate until we are in the byte holding the bit after the current tail bit */
        if (intern_trail_bit < (N - iter - 1) * CHAR_BIT) {continue;}

        if (offset_shifting > 0 && bits_written != 0)
        {
            bytes_[byte_id] |= (serialized_suffix[vec_index - 1] << (CHAR_BIT - offset_shifting));
            bits_written += offset_shifting;
            if (bits_written >= num_bits) {break;}
        }

        bytes_[byte_id] |= (serialized_suffix[vec_index] >> offset_shifting);
        bits_written += CHAR_BIT - offset_shifting;
        if (bits_written >= num_bits) {break;}
        ++vec_index;
    }

    indicators_.tail_bit_ -= num_bits;
}

template<int N>
SerializedCompressionValue<N>
GetCommonPrefix(const SerializedCompressionValue<N>& value1, const SerializedCompressionValue<N>& value2)
{
    const int max_tail_bit = std::max(value1.indicators_.tail_bit_, value2.indicators_.tail_bit_);
    const int maximum_bytes_length = N * CHAR_BIT - max_tail_bit;
    int current_bit_length = 0;

    std::array<uint8_t, N> bytes_bits;
    bytes_bits.fill(0);

    SerializedCompressionValue<N> cprefix = value1 ^ value2;

    for (int byte_id = GetMSBByteStart<N>(); MSBContinueIteration<N>(byte_id); MSBByteIncrement(byte_id))
    {
        if (cprefix.bytes_[byte_id] == 0)
        {
            if (current_bit_length + CHAR_BIT < maximum_bytes_length)
            {
                bytes_bits[byte_id] = value1[byte_id];
                current_bit_length += CHAR_BIT;
            } else if (current_bit_length + CHAR_BIT == maximum_bytes_length)
            {
                bytes_bits[byte_id] = value1[byte_id];
                return SerializedCompressionValue<N>(std::move(bytes_bits), max_tail_bit);
            }
            else
            {
                const int bit_accessor = maximum_bytes_length - current_bit_length;
                bytes_bits[byte_id] = value1[byte_id] & ~(LowBitMask[bit_accessor]);

                return SerializedCompressionValue<N>(std::move(bytes_bits), max_tail_bit);
            }
        } else
        {
            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if (maximum_bytes_length <= current_bit_length)
                {
                    /* The maximum length has been reached */
                    bytes_bits[byte_id] = value1[byte_id] & ~(LowBitMask[bit_index]);

                    return SerializedCompressionValue<N>(std::move(bytes_bits), N * CHAR_BIT - current_bit_length);
                }
                else if ((cprefix.bytes_[byte_id] & (kHighBit >> bit_index)))
                {
                    /* The first unequal bit has been found */
                    bytes_bits[byte_id] = value1[byte_id] & ~(LowBitMask[bit_index]);

                    return SerializedCompressionValue<N>(std::move(bytes_bits), N * CHAR_BIT - current_bit_length);
                } else
                {
                    ++current_bit_length;
                }
            }
        }
    }

    return (value1.indicators_.tail_bit_ >= value2.indicators_.tail_bit_ ? value1 : value2);
}


template<int N>
bool
SerializedCompressionValue<N>::IsEmpty() const
{
    if (indicators_.front_bit_ + indicators_.tail_bit_ >= static_cast<uint8_t>(N * CHAR_BIT))
    {
        return true;
    } else
    {
        return false;
    }
}

template<int N>
bool
SerializedCompressionValue<N>::IsExactlyEmpty() const
{
    if (indicators_.front_bit_ + indicators_.tail_bit_ == static_cast<uint8_t>(N * CHAR_BIT))
    {
        return true;
    } else
    {
        return false;
    }
}

template<int N>
bool
SerializedCompressionValue<N>::IsBitSetAtPosition(const int bit_position) const
{
    const int byte_pos = bit_position / CHAR_BIT;
    const int bit_pos = bit_position - byte_pos * CHAR_BIT;

    int byte_id = GetLSBByteStart<N>();
    for (int counter = 0; LSBContinueIteration<N>(byte_id) && counter < byte_pos; LSBByteIncrement(byte_id), ++counter)
    {
        //Iterate until the byte_id is correctly set
    }

    return (bytes_[byte_id] & (kLowBit << bit_pos) ? true : false);
}

template<int N>
void
SerializedCompressionValue<N>::ToggleBitAtPosition(const int bit_position)
{
    const int byte_pos = bit_position / CHAR_BIT;
    const int bit_pos = bit_position - byte_pos * CHAR_BIT;

    int byte_id = GetLSBByteStart<N>();
    for (int counter = 0; LSBContinueIteration<N>(byte_id) && counter < byte_pos; LSBByteIncrement(byte_id), ++counter)
    {
        //Iterate until the byte_id is correctly set
    }

    bytes_[byte_id] ^= (kLowBit << bit_pos);
}

template<int N>
void
SerializedCompressionValue<N>::ToggleBitAtPositionAndResetPreviousBits(const int bit_position)
{
    const bool bit_value = IsBitSetAtPosition(bit_position);

    const int byte_pos = bit_position / CHAR_BIT;
    const int bit_pos = bit_position - byte_pos * CHAR_BIT;

    int byte_id = GetLSBByteStart<N>();
    for (int counter = 0; LSBContinueIteration<N>(byte_id) && counter <= byte_pos; LSBByteIncrement(byte_id), ++counter)
    {
        //Iterate until the byte_id is correctly set
        if  (counter == byte_pos)
        {
            bytes_[byte_id] = (bytes_[byte_id] & ~(LowBitMask[CHAR_BIT - 1 - bit_pos])) ||
                              ((bit_value == true) ?  uint8_t{0xFF} >> (CHAR_BIT - bit_pos) : uint8_t{0x01} << bit_pos);
        } else
        {
            /* Fill the whole byte */
            bytes_[byte_id] = (bit_value ==  true ? uint8_t{0x00} : uint8_t{0xFF});
        }
    }
}

template<int N>
void
SerializedCompressionValue<N>::SetFrontBit(const uint8_t front_bit)
{
    //cmc_assert(CHAR_BIT * N >= front_bit + indicators_.tail_bit_);
    indicators_.front_bit_ = front_bit;
}

template<int N>
void
SerializedCompressionValue<N>::SetTailBit(const uint8_t tail_bit)
{
    //cmc_assert(CHAR_BIT * N >= front_bit + indicators_.tail_bit_);
    indicators_.tail_bit_ = tail_bit;
}

template<int N>
uint8_t
SerializedCompressionValue<N>::GetTailBit() const
{
    return indicators_.tail_bit_;
}

template<int N>
uint8_t
SerializedCompressionValue<N>::GetFrontBit() const
{
    return indicators_.front_bit_;
}

template<int N>
void
SerializedCompressionValue<N>::ToggleTailUntilNextUnsetBit()
{
    int byte_id = GetLSBByteStart<N>();
    /* Get to the current start byte */
    const int byte_increments = static_cast<int>(indicators_.tail_bit_ / CHAR_BIT);
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

    int bit_accessor = indicators_.tail_bit_ - byte_increments * CHAR_BIT;

    for (int byte_index = byte_id; LSBContinueIteration<N>(byte_index); LSBByteIncrement(byte_index))
    {
        for (int bit_index = bit_accessor; bit_index < CHAR_BIT; ++bit_index)
        {
            if (bytes_[byte_index] & (kLowBit << bit_index))
            {
                /* If there is a one at this position */
                has_encountered_one = true;
                /* We toggle the bit (unset it) and continue iterating until we have found leading zero */
                bytes_[byte_index] &= ~(kLowBit << bit_index);
                continue;
            } else
            {
                /* If there is a zero at this position */
                if (has_encountered_one)
                {
                    /* If we have aready seen a one before, we will set the current bit to one and are finished */
                    bytes_[byte_index] |= (kLowBit << bit_index);
                    indicators_.tail_bit_ += bit_index - bit_accessor;
                    return;
                } else
                {
                    /* If we have not seen a one before, we just keep iterating until we will find one */
                    continue;
                }
            }
        }

        /* Update the tail */
        indicators_.tail_bit_ += CHAR_BIT - bit_accessor;

        /* Reset the bit_accesor, for the next byte */
        bit_accessor = 0;
    }
}

template<int N>
template<typename T>
SerializedCompressionValue<N>::SerializedCompressionValue(const T& value)
{
    cmc_assert(N == static_cast<int>(sizeof(T)));
    std::memcpy(this->bytes_.data(), &value, sizeof(T));

    indicators_.front_bit_ = 0;
    indicators_.tail_bit_ = 0;
};

template<int N>
SerializedCompressionValue<N>::SerializedCompressionValue(std::array<uint8_t, N>&& serialized_value)
: bytes_{std::move(serialized_value)}
{
    indicators_.front_bit_ = 0;
    indicators_.tail_bit_ = 0;
};

template<int N>
SerializedCompressionValue<N>::SerializedCompressionValue(std::array<uint8_t, N>&& serialized_value, const uint8_t tail_bit)
: bytes_{std::move(serialized_value)}
{
    indicators_.tail_bit_ = 0;
    indicators_.tail_bit_ = tail_bit;
}

template<int N>
inline
int
SerializedCompressionValue<N>::GetNumberTrailingZeros() const
{
    int num_trailing_zeros = 0;
    for (int byte_id = GetLSBByteStart<N>(); LSBContinueIteration<N>(byte_id); LSBByteIncrement(byte_id))
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
                    /* If true, the value does not hold a zero-bit at this position */
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
inline
int
SerializedCompressionValue<N>::GetNumberTrailingOnes() const
{
    int num_trailing_ones = 0;
    for (int byte_id = GetLSBByteStart<N>(); LSBContinueIteration<N>(byte_id); LSBByteIncrement(byte_id))
    {
        /* Check whether the full byte is zero */
        if (bytes_[byte_id] == kOneByte)
        {
            num_trailing_ones += CHAR_BIT;
        } else
        {
            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if (bytes_[byte_id] & (kLowBit << bit_index))
                {
                    /* If true, the value does hold a one-bit at this position */
                    ++num_trailing_ones;
                } else
                {
                    return num_trailing_ones;
                }
            }
        }
    }
    return num_trailing_ones;
}

template<int N>
inline
int
SerializedCompressionValue<N>::GetNumberLeadingZeros() const
{
    int num_leading_zeros = 0;
    for (int byte_id = GetMSBByteStart<N>(); MSBContinueIteration<N>(byte_id); MSBByteIncrement(byte_id))
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
                    /* If true, the value does hold not a zero-bit at this position */
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
inline
bool
SerializedCompressionValue<N>::IsTailBitValueSet() const
{
    const int lsb_byte_pos = indicators_.tail_bit_ / CHAR_BIT;
    int byte_pos = GetLSBByteStart<N>();
    for (int byte_id = 0; byte_id < lsb_byte_pos; LSBByteIncrement(byte_pos)){}

    const int bit_pos = indicators_.tail_bit_ - lsb_byte_pos * CHAR_BIT;
    return ((bytes_[byte_pos] >> bit_pos) & uint8_t{1});
}

template<int N>
inline
int
SerializedCompressionValue<N>::GatherEqualTailBitCount() const
{
    const int num_trailing_ones = this->GetNumberTrailingOnes();
    const int num_trailing_zeros = this->GetNumberTrailingZeros();
    const int max_bit_seq = std::max(num_trailing_ones, num_trailing_zeros);
    return max_bit_seq;
}


template<int N>
inline
void
SerializedCompressionValue<N>::TrimEqualTailSequence()
{
    const int num_trailing_ones = this->GetNumberTrailingOnes();
    const int num_trailing_zeros = this->GetNumberTrailingZeros();
    const int max_bit_seq = std::max(num_trailing_ones, num_trailing_zeros);
    if (max_bit_seq < this->GetSizeBits())
    {
        indicators_.tail_bit_ = static_cast<uint8_t>(max_bit_seq);
    }
}

template<int N>
inline
void
SerializedCompressionValue<N>::RevertEqualTailSequence()
{
    if (indicators_.tail_bit_ == 0)
    {
        return;
    }

    /* Gather whether the tail bit is set or not */
    const bool is_tail_value_set = IsTailBitValueSet();

    /* Depending on whether the value is set or not, we create a contrary sequence which will be appended */
    const uint8_t sequence = is_tail_value_set ? 0x00 : 0xFF;
    const std::vector<uint8_t> max_tail_sequence(N, sequence);

    /* Append the previously trimmed tail sequence */
    this->ApplySuffix(max_tail_sequence,  indicators_.tail_bit_);

    cmc_assert(indicators_.tail_bit_ == 0);
}

template<int N>
inline
void
SerializedCompressionValue<N>::UpdateTailBitCount()
{
    indicators_.tail_bit_ = this->GetNumberTrailingZeros();
}

template<int N>
inline
void
SerializedCompressionValue<N>::UpdateFrontBitCount()
{
    indicators_.front_bit_ = this->GetNumberLeadingZeros();
}

template<int N>
void
SerializedCompressionValue<N>::ClearNextSetBitFromTail()
{
    int byte_id = GetLSBByteStart<N>();
    /* Get to the current start byte */
    const int byte_increments = static_cast<int>(indicators_.tail_bit_ / CHAR_BIT);
    if (byte_increments == N)
    {
        /* The value has already been completely worked out */
        return;
    }
    for (int increment = 0; increment < byte_increments; ++increment)
    {
        LSBByteIncrement(byte_id);
    }

    int bit_accessor = indicators_.tail_bit_ - byte_increments * CHAR_BIT;

    for (int byte_index = byte_id; LSBContinueIteration<N>(byte_index); LSBByteIncrement(byte_index))
    {
        for (int bit_index = bit_accessor; bit_index < CHAR_BIT; ++bit_index)
        {
            if (bytes_[byte_index] & (kLowBit << bit_index))
            {
                /* If there is a one, clear this position and we are finished, else continue */
                bytes_[byte_index] &= ~(kLowBit << bit_index);
                indicators_.tail_bit_ += bit_index - bit_accessor;   
                return;
            }
        }

        /* Update the tail */
        indicators_.tail_bit_ += CHAR_BIT - bit_accessor;

        /* Reset the bit accesor */
        bit_accessor = 0;
    }
}

template<int N>
inline
int
SerializedCompressionValue<N>::GetCountOfSignificantBits() const
{
    cmc_assert(indicators_.front_bit_ + indicators_.tail_bit_ <= N * CHAR_BIT);
    return N * CHAR_BIT - indicators_.front_bit_ - indicators_.tail_bit_;
}

template<int N>
std::vector<uint8_t>
SerializedCompressionValue<N>::GetSignificantBitsInBigEndianOrdering() const
{
    /* Get the number of bits for this prefix */
    const int num_bits = GetCountOfSignificantBits();

    if (num_bits == 0)
    {
        /* In case there are no significant bits, we are returning an empty vector */
        return std::vector<uint8_t>();
    }

    /* Get the number of bytes needed to accomodate the prefix */
    const int num_bytes = (num_bits / CHAR_BIT) + ((num_bits % CHAR_BIT) != 0 ? 1 : 0);
    
    /* Allocate the needed memory for the prefix */
    std::vector<uint8_t> bytes(num_bytes, 0);

    const int offset_shifting = indicators_.front_bit_ - static_cast<int>(indicators_.front_bit_ / CHAR_BIT) * CHAR_BIT;

    int bits_written = 0;

    const int bit_shift_to_front = offset_shifting;
    const int bit_shift_to_back = CHAR_BIT - offset_shifting;

    const int front_bits_written = CHAR_BIT - bit_shift_to_front;
    const int tail_bits_written = CHAR_BIT - bit_shift_to_back;

    for (int byte_id = GetMSBByteStart<N>(), iter = 0, vec_index = 0;
         MSBContinueIteration<N>(byte_id);
         MSBByteIncrement(byte_id), ++iter)
    {
        /* Iterate until we are in the byte holding the start bit */
        if (indicators_.front_bit_ >= (iter + 1) * CHAR_BIT) {continue;}

        /* Fill the back of the previous byte, if the start bit has an offset */
        if ((bits_written > 0 && offset_shifting != 0))
        {
            /* Assign the front bits to the the back of the previous byte */
            bytes[vec_index - 1] |= (bytes_[byte_id] >> (bit_shift_to_back));
            bits_written += tail_bits_written;
            if (bits_written >= num_bits) {break;}
        }

        /* Assign the the bits to the front part of the current  byte */
        bytes[vec_index] |= (bytes_[byte_id] << (bit_shift_to_front));
        bits_written += front_bits_written;
        if (bits_written >= num_bits) {break;}

        ++vec_index;
    }

    return bytes;
}

template<int N, typename T>
inline
auto
SerializeIntegerValueNatively(const T& value)
 -> std::enable_if_t<std::is_integral_v<T>, std::array<uint8_t, N>>
{
    cmc_assert(static_cast<size_t>(N) == sizeof(T));
    const uint8_t* mem_ptr = reinterpret_cast<const uint8_t*>(&value);
    std::array<uint8_t, N> serialized;
    for (int i = 0; i < N; ++i, ++mem_ptr)
    {
        serialized[i] = *mem_ptr;
    }
    return serialized;
}

template<>
inline void
SerializedCompressionValue<1>::PerformIntegerAddition(const SerializedCompressionValue<1>& residual)
{
    /* Get the value and the residual */
    uint8_t value = this->template ReinterpretDataAs<uint8_t>();
    const uint8_t res = residual.template ReinterpretDataAs<uint8_t>();
    /* Perform the integer addition */
    value += res;
    /* Store the value */
    bytes_ = std::array<uint8_t, 1>{value};
}

template<>
inline void
SerializedCompressionValue<2>::PerformIntegerAddition(const SerializedCompressionValue<2>& residual)
{
    /* Get the value and the residual */
    uint16_t value = this->template ReinterpretDataAs<uint16_t>();
    const uint16_t res = residual.template ReinterpretDataAs<uint16_t>();
    /* Perform the integer addition */
    value += res;
    /* Store the value */
    bytes_ = SerializeIntegerValueNatively<2, uint16_t>(value);
}

template<>
inline void
SerializedCompressionValue<4>::PerformIntegerAddition(const SerializedCompressionValue<4>& residual)
{
    /* Get the value and the residual */
    uint32_t value = this->template ReinterpretDataAs<uint32_t>();
    const uint32_t res = residual.template ReinterpretDataAs<uint32_t>();
    /* Perform the integer addition */
    value += res;
    /* Store the value */
    bytes_ = SerializeIntegerValueNatively<4, uint32_t>(value);
}

template<>
inline void
SerializedCompressionValue<8>::PerformIntegerAddition(const SerializedCompressionValue<8>& residual)
{
    /* Get the value and the residual */
    uint64_t value = this->template ReinterpretDataAs<uint64_t>();
    const uint64_t res = residual.template ReinterpretDataAs<uint64_t>();
    /* Perform the integer addition */
    value += res;
    /* Store the value */
    bytes_ = SerializeIntegerValueNatively<8, uint64_t>(value);
}

template<>
inline void
SerializedCompressionValue<1>::PerformIntegerSubtraction(const SerializedCompressionValue<1>& residual)
{
    /* Get the value and the residual */
    uint8_t value = this->template ReinterpretDataAs<uint8_t>();
    const uint8_t res = residual.template ReinterpretDataAs<uint8_t>();
    /* Perform the integer addition */
    value -= res;
    /* Store the value */
    bytes_ = std::array<uint8_t, 1>{value};
}

template<>
inline void
SerializedCompressionValue<2>::PerformIntegerSubtraction(const SerializedCompressionValue<2>& residual)
{
    /* Get the value and the residual */
    uint16_t value = this->template ReinterpretDataAs<uint16_t>();
    const uint16_t res = residual.template ReinterpretDataAs<uint16_t>();
    /* Perform the integer addition */
    value -= res;
    /* Store the value */
    bytes_ = SerializeIntegerValueNatively<2, uint16_t>(value);
}

template<>
inline void
SerializedCompressionValue<4>::PerformIntegerSubtraction(const SerializedCompressionValue<4>& residual)
{
    /* Get the value and the residual */
    uint32_t value = this->template ReinterpretDataAs<uint32_t>();
    const uint32_t res = residual.template ReinterpretDataAs<uint32_t>();
    /* Perform the integer addition */
    value -= res;
    /* Store the value */
    bytes_ = SerializeIntegerValueNatively<4, uint32_t>(value);
}

template<>
inline void
SerializedCompressionValue<8>::PerformIntegerSubtraction(const SerializedCompressionValue<8>& residual)
{
    /* Get the value and the residual */
    uint64_t value = this->template ReinterpretDataAs<uint64_t>();
    const uint64_t res = residual.template ReinterpretDataAs<uint64_t>();
    /* Perform the integer addition */
    value -= res;
    /* Store the value */
    bytes_ = SerializeIntegerValueNatively<8, uint64_t>(value);
}

template<int N>
void
SerializedCompressionValue<N>::_TestPCP4ReverseXORResidual(const SerializedCompressionValue<N>& residual)
{
    *this ^= residual;
}

/* Convert a view of ComrpessionValues to a vector of values of a certain data type */
template<typename T>
std::vector<T>
ConvertCompressionValues(const VectorView<SerializedCompressionValue<sizeof(T)>>& compression_values)
{
    std::vector<T> vals;
    vals.reserve(compression_values.size()); 

    for (auto val_iter = compression_values.begin(); val_iter != compression_values.end(); ++val_iter)
    {
        vals.emplace_back(val_iter->template ReinterpretDataAs<T>());
    }

    return vals;
}


}

#endif /* !CMC_BYTE_COMPRESSION_VALUES_HXX */
