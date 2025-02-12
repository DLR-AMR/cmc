#ifndef CMC_PREFIX_HXX
#define CMC_PREFIX_HXX

#include "utilities/cmc_serialized_compression_value_forward.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"
#include "utilities/cmc_arithmetic_encoder.hxx"
#include "utilities/cmc_bytes.hxx"

#include <array>
#include <vector>
#include <climits>
#include <cstring>
#include <bitset>

namespace cmc
{

template <typename T>
std::vector<T>
GetCompressionValuesAs(const std::vector<CompressionValue<sizeof(T)>>& serialized_values, const int start_index, const int num_elements);

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
    void SetTailBit(const int tail_bit);
    int GetTrailBit() const;
    int GetFrontBit() const;
    int GetCountOfSignificantBits() const;
    bool IsEmpty() const;
    std::vector<uint8_t> GetSignificantBitsInBigEndianOrdering() const;
    //std::vector<uint8_t> GetSignificantBitsInBigEndianOrdering_Offset_4_Stride_12() const;
    void ApplyPrefix(const std::vector<uint8_t>& serialized_prefix, const int num_bits);
    void AddIntegerResidual(const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits);
    void AddXORResidualWithoutImplicitOneBit(const uint32_t lzc, const std::vector<uint8_t>& residual_bits);

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

    CompressionValue& operator^=(const CompressionValue& rhs)
    {
        prefix_ ^= rhs.prefix_;
        return *this;
    }

    void NullifyNonSignificantFrontBits();
    int GetLeadingZeroCountInSignificantBits() const;

    void PerformIntegerSubtraction(const CompressionValue& residual);
    void PerformIntegerAddition(const CompressionValue& residual);
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

template<int N>
int
CompressionValue<N>::GetLeadingZeroCountInSignificantBits() const
{
    //TODO only iterate until tail bit is reached

    if (this->IsEmpty()) {return 0;}
    cmc_assert(trail_bit_ == 0);
    
    int num_leading_zeros = 0;
    int front_bit = front_bit_;

    int counter = front_bit_;

    for (int byte_id = GetMSBByteStart(prefix_), iter = 0, vec_index = 0;
         MSBContinueIteration(byte_id, prefix_);
         MSBByteIncrement(byte_id), ++iter)
    {
        if (counter + trail_bit_ >= N * CHAR_BIT)
        {
            return num_leading_zeros;
        }
        /* Iterate until we are in the byte holding the bit after the current front bit */
        if (front_bit >= (iter + 1) * CHAR_BIT) {continue;}

        for (int bit_index = front_bit % CHAR_BIT; bit_index < CHAR_BIT; ++bit_index, ++counter)
        {
            if (prefix_[byte_id] & (kHighBit >> bit_index))
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
CompressionValue<N>::NullifyNonSignificantFrontBits()
{
    /* Get the number of bits for this prefix */
    const int num_bits = GetCountOfSignificantBits();

    if (num_bits == N * CHAR_BIT)
    {
        /* In case there are no non-significant bits, we return immediately */
        return;
    }

    int front_bits_to_nullify = front_bit_;

    /* Nullify the first part of the prefix */
    for (int byte_id = GetMSBByteStart(prefix_);
         MSBContinueIteration(byte_id, prefix_);
         MSBByteIncrement(byte_id))
    {
        if (front_bits_to_nullify <= 0) {break;}

        if (front_bits_to_nullify >= CHAR_BIT)
        {
            /* Nullify the whole byte */
            prefix_[byte_id]  = uint8_t{0};
            front_bits_to_nullify -= CHAR_BIT;
        } else 
        {
            const uint8_t bit_mask = LowBitMask[front_bits_to_nullify];
            prefix_[byte_id] &= bit_mask;
            break;
        }
    }

}

/* This functions needs the serialized prefix to be aligned at the high bits ( e.g. four bit prefix: 0b(p1 p2 p3 p4 0 0 0 0) */
template<int N>
void
CompressionValue<N>::ApplyPrefix(const std::vector<uint8_t>& serialized_prefix, const int num_bits)
{
    if (num_bits <= 0) {return;}
    cmc_assert(trail_bit_ > 0);

    int bits_written = 0;

    const int intern_trail_bit = trail_bit_ - 1;
    const int offset_shifting_ = (trail_bit_ - static_cast<int>(trail_bit_ / CHAR_BIT) * CHAR_BIT);
    const int offset_shifting = (offset_shifting_ > 0 ? CHAR_BIT - offset_shifting_ : 0);

    for (int byte_id = GetMSBByteStart(prefix_), iter = 0, vec_index = 0;
         MSBContinueIteration(byte_id, prefix_);
         MSBByteIncrement(byte_id), ++iter)
    {
        /* Iterate until we are in the byte holding the bit after the current tail bit */
        if (intern_trail_bit < (N - iter - 1) * CHAR_BIT) {continue;}

        if (offset_shifting > 0 && bits_written != 0)
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
            if (current_bit_length + CHAR_BIT < maximum_prefix_length)
            {
                prefix_bits[byte_id] = value1[byte_id];
                current_bit_length += CHAR_BIT;
            } else if (current_bit_length + CHAR_BIT == maximum_prefix_length)
            {
                prefix_bits[byte_id] = value1[byte_id];
                return CompressionValue<N>(std::move(prefix_bits), max_trail_bit);
            }
            else
            {
                const int bit_accessor = maximum_prefix_length - current_bit_length;
                prefix_bits[byte_id] = value1[byte_id] & ~(LowBitMask[bit_accessor]);

                return CompressionValue<N>(std::move(prefix_bits), max_trail_bit);
            }
        } else
        {
            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if (maximum_prefix_length <= current_bit_length)
                {
                    /* The maximum length has been reached */
                    prefix_bits[byte_id] = value1[byte_id] & ~(LowBitMask[bit_index]);

                    return CompressionValue<N>(std::move(prefix_bits), N * CHAR_BIT - current_bit_length);
                }
                else if ((xor_bits[byte_id] & (kHighBit >> bit_index)))
                {
                    /* The first unequal bit has been found */
                    prefix_bits[byte_id] = value1[byte_id] & ~(LowBitMask[bit_index]);

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
    //cmc_assert(CHAR_BIT * N >= front_bit + trail_bit_);
    front_bit_ = front_bit;
}

template<int N>
void
CompressionValue<N>::SetTailBit(const int tail_bit)
{
    //cmc_assert(CHAR_BIT * N >= front_bit + trail_bit_);
    trail_bit_ = tail_bit;
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
CompressionValue<N>::GetSignificantBitsInBigEndianOrdering() const
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

    const int offset_shifting = front_bit_ - static_cast<int>(front_bit_ / CHAR_BIT) * CHAR_BIT;

    int bits_written = 0;

    const int bit_shift_to_front = offset_shifting;
    const int bit_shift_to_back = CHAR_BIT - offset_shifting;

    const int front_bits_written = CHAR_BIT - bit_shift_to_front;
    const int tail_bits_written = CHAR_BIT - bit_shift_to_back;

    for (int byte_id = GetMSBByteStart(prefix_), iter = 0, vec_index = 0;
         MSBContinueIteration(byte_id, prefix_);
         MSBByteIncrement(byte_id), ++iter)
    {
        /* Iterate until we are in the byte holding the start bit */
        if (front_bit_ >= (iter + 1) * CHAR_BIT) {continue;}

        /* Fill the back of the previous byte, if the start bit has an offset */
        if ((bits_written > 0 && offset_shifting != 0))
        {
            /* Assign the front bits to the the back of the previous byte */
            bytes[vec_index - 1] |= (prefix_[byte_id] >> (bit_shift_to_back));
            bits_written += tail_bits_written;
            if (bits_written >= num_bits) {break;}
        }

        /* Assign the the bits to the front part of the current  byte */
        bytes[vec_index] |= (prefix_[byte_id] << (bit_shift_to_front));
        bits_written += front_bits_written;
        if (bits_written >= num_bits) {break;}

        ++vec_index;
    }

    return bytes;
}

template<>
inline void
CompressionValue<1>::PerformIntegerAddition(const CompressionValue<1>& residual)
{
    /* Get the value and the residual */
    uint8_t value = this->template ReinterpretDataAs<uint8_t>();
    const uint8_t res = residual.template ReinterpretDataAs<uint8_t>();
    /* Perform the integer addition */
    value += res;
    /* Store the value */
    prefix_ = Prefix<1>(value);
}

template<>
inline void
CompressionValue<2>::PerformIntegerAddition(const CompressionValue<2>& residual)
{
    /* Get the value and the residual */
    uint16_t value = this->template ReinterpretDataAs<uint16_t>();
    const uint16_t res = residual.template ReinterpretDataAs<uint16_t>();
    /* Perform the integer addition */
    value += res;
    /* Store the value */
    prefix_ = Prefix<2>(value);
}

template<>
inline void
CompressionValue<4>::PerformIntegerAddition(const CompressionValue<4>& residual)
{
    /* Get the value and the residual */
    uint32_t value = this->template ReinterpretDataAs<uint32_t>();
    const uint32_t res = residual.template ReinterpretDataAs<uint32_t>();
    /* Perform the integer addition */
    value += res;
    /* Store the value */
    prefix_ = Prefix<4>(value);
}

template<>
inline void
CompressionValue<8>::PerformIntegerAddition(const CompressionValue<8>& residual)
{
    /* Get the value and the residual */
    uint64_t value = this->template ReinterpretDataAs<uint64_t>();
    const uint64_t res = residual.template ReinterpretDataAs<uint64_t>();
    /* Perform the integer addition */
    value += res;
    /* Store the value */
    prefix_ = Prefix<8>(value);
}

template<>
inline void
CompressionValue<1>::PerformIntegerSubtraction(const CompressionValue<1>& residual)
{
    /* Get the value and the residual */
    uint8_t value = this->template ReinterpretDataAs<uint8_t>();
    const uint8_t res = residual.template ReinterpretDataAs<uint8_t>();
    /* Perform the integer addition */
    value -= res;
    /* Store the value */
    prefix_ = Prefix<1>(value);
}

template<>
inline void
CompressionValue<2>::PerformIntegerSubtraction(const CompressionValue<2>& residual)
{
    /* Get the value and the residual */
    uint16_t value = this->template ReinterpretDataAs<uint16_t>();
    const uint16_t res = residual.template ReinterpretDataAs<uint16_t>();
    /* Perform the integer addition */
    value -= res;
    /* Store the value */
    prefix_ = Prefix<2>(value);
}

template<>
inline void
CompressionValue<4>::PerformIntegerSubtraction(const CompressionValue<4>& residual)
{
    /* Get the value and the residual */
    uint32_t value = this->template ReinterpretDataAs<uint32_t>();
    const uint32_t res = residual.template ReinterpretDataAs<uint32_t>();
    /* Perform the integer addition */
    value -= res;
    /* Store the value */
    prefix_ = Prefix<4>(value);
}

template<>
inline void
CompressionValue<8>::PerformIntegerSubtraction(const CompressionValue<8>& residual)
{
    /* Get the value and the residual */
    uint64_t value = this->template ReinterpretDataAs<uint64_t>();
    const uint64_t res = residual.template ReinterpretDataAs<uint64_t>();
    /* Perform the integer addition */
    value -= res;
    /* Store the value */
    prefix_ = Prefix<8>(value);
}

template<int N>
void
CompressionValue<N>::AddIntegerResidual(const uint32_t encoded_lzc, const std::vector<uint8_t>& residual_bits)
{
    //cmc_assert(not residual_bits.empty());
    cmc_assert(residual_bits.size() <= static_cast<size_t>(N));

    /* Decode the encoded LZC */
    const auto [signum, lzc] = arithmetic_encoding::DecodeLZC(encoded_lzc);

    //cmc_debug_msg("Signum: ", signum, ", lzc: ", lzc);
    if (lzc == N * CHAR_BIT)
    {
        /* The residual is empty, therefore nothing has to be added */
        return;
    }

    std::array<uint8_t, N> serialized_val;
    serialized_val.fill(uint8_t{0});

    /* Set up an compression value holding the residual in the front bits */
    CompressionValue<N> residual(serialized_val);
    residual.SetTailBit(N * CHAR_BIT - lzc);

    /* Afterwards, we add the the implicit one bit */
    residual.ApplyPrefix(std::vector<uint8_t>{0x80}, 1);
    const int remaining_bits = N * CHAR_BIT - lzc - 1;
    //const int remaining_bits = N * CHAR_BIT - lzc;

    cmc_assert((not residual_bits.empty()) || ((remaining_bits == 0) && residual_bits.empty()));

    if (residual.GetTrailBit() > 0)
    {
        /* And finally, we combine it with the actual remaining residual bits */
        residual.ApplyPrefix(residual_bits, remaining_bits);
    }

    /* Now, we should have a full CompressionValue resembling the residual */
    cmc_assert(residual.GetFrontBit() == 0 && residual.GetTrailBit() == 0);

    //if constexpr (N == 4)
    //{
    //    uint32_t eeee = residual.template ReinterpretDataAs<uint32_t>();
    //    cmc_debug_msg("The residual is: ", std::bitset<32>(eeee));
    //}
    /* Dependeing on the signum, either add or subtract the residual */
    if (signum == true)
    {
        this->PerformIntegerAddition(residual);
    } else
    {
        this->PerformIntegerSubtraction(residual);
    }
}

template<int N>
void
CompressionValue<N>::AddXORResidualWithoutImplicitOneBit(const uint32_t lzc, const std::vector<uint8_t>& residual_bits)
{
    //cmc_assert(not residual_bits.empty());
    cmc_assert(residual_bits.size() <= static_cast<size_t>(N));

    //cmc_debug_msg("Signum: ", signum, ", lzc: ", lzc);
    if (lzc >= N * CHAR_BIT)
    {
        /* The residual is empty, therefore nothing has to be added */
        return;
    }

    std::array<uint8_t, N> serialized_val;
    serialized_val.fill(uint8_t{0});

    /* Set up an compression value holding the residual in the front bits */
    CompressionValue<N> residual(serialized_val);
    residual.SetTailBit(N * CHAR_BIT - lzc);
    

    const int remaining_bits = N * CHAR_BIT - lzc;

    if (residual.GetTrailBit() > 0)
    {
        /* And finally, we combine it with the actual remaining residual bits */
        residual.ApplyPrefix(residual_bits, remaining_bits);
    }

    /* Now, we should have a full CompressionValue resembling the residual */
    cmc_assert(residual.GetFrontBit() == 0 && residual.GetTrailBit() == 0);

    /* Reverse the XOR operation */
    prefix_ ^= residual.prefix_;
}

template <typename T>
std::vector<T>
GetCompressionValuesAs(const std::vector<CompressionValue<sizeof(T)>>& serialized_values, const int start_index, const int num_elements)
{
    cmc_assert(start_index >= 0 && num_elements > 0);
    cmc_assert(static_cast<size_t>(start_index + num_elements) <= serialized_values.size());

    std::vector<T> vals;
    vals.reserve(num_elements);

    const int end_index = start_index + num_elements;

    for (int idx = start_index; idx < end_index; ++idx)
    {
        vals.push_back(serialized_values[idx].template ReinterpretDataAs<T>());
    }

    return vals;
}

}

#endif /* !CMC_PREFIX_HXX */
