#ifndef CMC_BITS_SPAN_HXX
#define CMC_BITS_SPAN_HXX

#include "cmc.hxx"
#include "utilities/cmc_bits_vector.hxx"

namespace cmc::bits
{

/**
 * A viewer class to extract bit sequences of a (padded) big-ednain stream serialized by cmc::bits::vector
 */
class vector_view
{
public:
    vector_view() = default;
    vector_view(const uint64_t* data)
    : data_{data}, pos_{0}, current_value_{ConvertBigEndianToNativeEndianness(*data)}, bit_position_{kBitIndexStart} {};

    void MoveToNextBit();
    void MoveToNextByteStart();
    bool IsCurrentBitSet() const;
    void SkipNumberOfBits(const size_t num_bits);

    void MoveToOffsetBitInStream(const size_t global_bit_position_bigendian_stream);
    template<UnsignedIntegerType T> T GetNextBitSequence(const int num_bits);
    bool GetNextBit();
    void SetStart(const uint64_t* data);
private:
    const uint64_t* data_{nullptr};
    int64_t pos_{0};
    uint64_t current_value_{0};
    int64_t bit_position_{kBitIndexStart};
};

inline void
vector_view::MoveToOffsetBitInStream(const size_t global_bit_position_bigendian_stream)
{
    /* Determine the value index in which the bit lies */
    pos_ = global_bit_position_bigendian_stream >> 6;
    const int be_bit_pos = global_bit_position_bigendian_stream - (pos_ << 6);
    current_value_ = ConvertBigEndianToNativeEndianness(*(data_ + pos_));
    bit_position_ = kBitIndexStart - (global_bit_position_bigendian_stream - pos_ * sizeof(uint64_t) * kCharBit); 
}

inline bool
vector_view::IsCurrentBitSet() const
{
    return ((current_value_ >> bit_position_) & uint64_t{1});
}

inline void
vector_view::MoveToNextBit()
{
    --bit_position_;
    if (bit_position_ < 0) [[unlikely]]
    {
        ++pos_;
        current_value_ = ConvertBigEndianToNativeEndianness(*(data_ + pos_));
        bit_position_ = kBitIndexStart;
    }
}
   
/* In case, the view already points to a start of a full byte, the pointer remains unchanged */
inline void
vector_view::MoveToNextByteStart()
{
    bit_position_ = ((bit_position_ + 1) / kCharBit) * kCharBit - 1;
    if (bit_position_ < 0) [[unlikely]]
    {
        ++pos_;
        current_value_ = ConvertBigEndianToNativeEndianness(*(data_ + pos_));
        bit_position_ = kBitIndexStart;
    }
}

inline bool
vector_view::GetNextBit()
{
    const bool bit = (current_value_ >> bit_position_) & uint64_t{1};
    --bit_position_;
    if (bit_position_ < 0) [[unlikely]]
    {
        ++pos_;
        current_value_ = ConvertBigEndianToNativeEndianness(*(data_ + pos_));
        bit_position_ = kBitIndexStart;
    }
    return bit;
}

template<UnsignedIntegerType T>
T
vector_view::GetNextBitSequence(const int num_bits)
{
    /* Check if there are bits to extarct */
    if (num_bits <= 0) [[unlikely]] 
    {
        return T();
    }

    /* Check if the type is capable of holding the bit-sequence */
    if (sizeof(T) * kCharBit < num_bits) [[unlikely]]
    {
        cmc_err_msg("The specified number of bits (", num_bits, ") does not fit in the requested type (maximum number of bits: ", sizeof(T) * kCharBit, ")!");
    }

    /* Check if the requested bit sequence is within the current value */
    const bool are_bits_entriely_in_current_val = ((bit_position_ + 1) >= num_bits);

    if (are_bits_entriely_in_current_val) [[likely]]
    {
        /* Compute the shift the sequence needs to have towards the least significant bit position */
        const int shift_to_lsb = (bit_position_ + 1) - num_bits;

        /* Nullify the remaining bits, that are not part of the sequence */
        const uint64_t bit_sequence = (current_value_ >> shift_to_lsb) & (~uint64_t{0} >> (64 - num_bits));
        
        /* Update the bit position */
        bit_position_ -= num_bits;
        if (bit_position_ < 0) [[unlikely]]
        {
            ++pos_;
            bit_position_ = kBitIndexStart;
            current_value_ = ConvertBigEndianToNativeEndianness(*(data_ + pos_));
        }

        /* The bit-sequence has been extracted in big endian, potentially, we need to swap the ordering*/
        return static_cast<T>(bit_sequence);
    } else
    {
        /* If the bit sequence does not lay entirely in the current value */
        //(bit_position_ + 1) bits are in the current value and (num_bits - (bit_position_ + 1)) are in the next value 
        /* Fill the sequence at the correct positions with the bits from the current value */
        uint64_t bit_sequence = ((*(data_ + pos_)) << (num_bits - (bit_position_ + 1))) & (~uint64_t{0} >> (64 - num_bits));
        
        /* Move to next value */
        ++pos_;
        current_value_ = ConvertBigEndianToNativeEndianness(*(data_ + pos_));

        /* Fill the sequence with the second part from the next value */
        bit_sequence |= (current_value_ >> (64 - (num_bits - (bit_position_ + 1))));

        /* Update the bit position */
        bit_position_ = kBitIndexStart - (num_bits - (bit_position_ + 1));

        //There is no need to check whether the bit_position_ becomes negative,
        //because due to the design (the bit_position always points to the start of the next bit-sequence,
        //and, therefore, at least one bit needs to reside in the preivous value, such that at max 63 bits
        //can be extarcted in the succeeding value)
        
        /* The bit-sequence has been extracted in big endian, potentially, we need to swap the ordering*/
        return static_cast<T>(bit_sequence);
    }
}

inline void
vector_view::SkipNumberOfBits(const size_t num_bits)
{
    cmc_assert(num_bits > 0);
    if (num_bits > 64) [[unlikely]]
    {
        cmc_err_msg("At maximum 64 bits can be skipped in the view!");
    }

    bit_position_ -= num_bits;

    if (bit_position_ < 0) [[unlikely]]
    {
        /* Rotate the bit position and load the next value from the stream */
        bit_position_ = 64 + bit_position_;
        ++pos_;
        current_value_ = ConvertBigEndianToNativeEndianness(*(data_ + pos_));
    }
}

inline void
vector_view::SetStart(const uint64_t* data)
{
    data_ = data;
    pos_ = 0;
    current_value_ = ConvertBigEndianToNativeEndianness(*(data_ + pos_));
    bit_position_ = kBitIndexStart;
}

}


#endif /* !CMC_BITS_SPAN_HXX */
