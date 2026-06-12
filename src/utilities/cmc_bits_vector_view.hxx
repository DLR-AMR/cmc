#ifndef CMC_BITS_SPAN_HXX
#define CMC_BITS_SPAN_HXX

#include "cmc.hxx"
#include "utilities/cmc_bits_vector.hxx"

namespace cmc::bits
{

using vector_view = vector_view_base<false>;
using vector_view_in_memory = vector_view_base<true>;

constexpr int kBitReadingStartIDX = 64;

/**
 * A viewer class to extract bit sequences of a (padded) big-endian stream serialized by cmc::bits::vector
 */
template<bool InMemory = false>
class vector_view_base
{
public:
    vector_view_base() = default;
    vector_view_base(const uint64_t* const data)
    : data_{data}, pos_{0}, current_value_{ConvertBigEndianToNativeEndianness(*data)}, bit_position_{kBitReadingStartIDX} {};
    vector_view_base(const cmc::bits::vector& vector)
    : data_{vector.vector_.data()}, pos_{-1}, current_value_{}, bit_position_{0} {
        if constexpr (not InMemory) {cmc_err_msg("A bits_vector_view on an active bits_vector can only be constructed InMemory!");}
    };


    void MoveToNextBit();
    void MoveToNextByteStart();
    void MoveToNextVectorValueStart();
    bool IsCurrentBitSet();
    void SkipNumberOfBits(const size_t num_bits);
    void SkipArbitraryNumberOfBits(const size_t num_bits);

    void MoveToOffsetBitInStream(const size_t global_bit_position_bigendian_stream);
    template<UnsignedIntegerType T> T GetNextBitSequence(const int num_bits);
    bool GetNextBit();

private:
    void GetValueAtPos();
    const uint64_t* data_{nullptr};
    int64_t pos_{-1};
    uint64_t current_value_{0};
    int64_t bit_position_{0};
};

template<bool InMemory>
inline void
vector_view_base<InMemory>::GetValueAtPos()
{
    if constexpr (not InMemory)
    {
        current_value_ = ConvertBigEndianToNativeEndianness(*(data_ + pos_));
    } else
    {
        current_value_ = *(data_ + pos_);
    }
}

template<bool InMemory>
inline void
vector_view_base<InMemory>::MoveToOffsetBitInStream(const size_t global_bit_position_bigendian_stream)
{
    if (global_bit_position_bigendian_stream % 64 == 0)
    {
        /* Set the pointer to the bit before */
        pos_ = (global_bit_position_bigendian_stream / 64) - 1;
        bit_position_ = 0;
    } else
    {
        /* Determine the value index in which the bit lies */
        pos_ = global_bit_position_bigendian_stream >> 6;
        this->GetValueAtPos();
        bit_position_ = kBitReadingStartIDX - (global_bit_position_bigendian_stream - pos_ * sizeof(uint64_t) * kCharBit); 
    }
}

template<bool InMemory>
inline bool
vector_view_base<InMemory>::IsCurrentBitSet()
{
    if (bit_position_ <= 0) [[unlikely]]
    {
        ++pos_;
        this->GetValueAtPos();
        bit_position_ = kBitReadingStartIDX;
    }

    return ((current_value_ >> (bit_position_ - 1)) & uint64_t{1});
}

template<bool InMemory>
inline void
vector_view_base<InMemory>::MoveToNextBit()
{
    if (bit_position_ <= 0) [[unlikely]]
    {
        ++pos_;
        this->GetValueAtPos();
        bit_position_ = kBitReadingStartIDX - 1;
    } else
    {
        --bit_position_;
    }
}
   
/* In case, the view already points to a start of a full byte, the pointer remains unchanged */
template<bool InMemory>
inline void
vector_view_base<InMemory>::MoveToNextByteStart()
{
    if (bit_position_ <= 0) [[unlikely]]
    {
        ++pos_;
        this->GetValueAtPos();
        bit_position_ = kBitReadingStartIDX;
    } else
    {
        bit_position_ = ((bit_position_) / kCharBit) * kCharBit;
    }
}

/* In case, the view already points tot he start of a value, the pointer remains unchanged */
template<bool InMemory>
inline void
vector_view_base<InMemory>::MoveToNextVectorValueStart()
{
    if (bit_position_ != kBitReadingStartIDX)
    {
        bit_position_ = kBitReadingStartIDX;
        ++pos_;
        this->GetValueAtPos();
    }
}

template<bool InMemory>
inline bool
vector_view_base<InMemory>::GetNextBit()
{
    if (bit_position_ <= 0) [[unlikely]]
    {
        ++pos_;
        this->GetValueAtPos();
        bit_position_ = kBitReadingStartIDX;
    }

    const bool bit = (current_value_ >> (bit_position_ - 1)) & uint64_t{1};
    --bit_position_;

    return bit;
}

/**
 * Get the specified bit sequence supplied in the given type T,
 * such that the sequence is aligned to the least significant bit position
 * Example: Get 5-Bit sequence in 32-Bit type: 0b00000000000000000000000000000XXXXX
 */
template<bool InMemory>
template<UnsignedIntegerType T>
T
vector_view_base<InMemory>::GetNextBitSequence(const int num_bits)
{
    /* Check if there are bits to extarct */
    if (num_bits <= 0) [[unlikely]] 
    {
        return T{};
    }

    /* Check if the type is capable of holding the bit-sequence */
    if (static_cast<int>(sizeof(T) * kCharBit) < num_bits) [[unlikely]]
    {
        cmc_err_msg("The specified number of bits (", num_bits, ") does not fit in the requested type (maximum number of bits: ", sizeof(T) * kCharBit, ")!");
    }

    if (bit_position_ <= 0) [[unlikely]]
    {
        ++pos_;
        this->GetValueAtPos();
        bit_position_ = kBitReadingStartIDX;
    }

    /* Check if the requested bit sequence is within the current value */
    const bool are_bits_entriely_in_current_val = (bit_position_ >= num_bits);

    if (are_bits_entriely_in_current_val) [[likely]]
    {
        /* Compute the shift the sequence needs to have towards the least significant bit position */
        const int shift_to_lsb = bit_position_ - num_bits;

        /* Nullify the remaining bits, that are not part of the sequence */
        const uint64_t bit_sequence = (current_value_ >> shift_to_lsb) & (~uint64_t{0} >> (64 - num_bits));
        
        /* Update the bit position */
        bit_position_ -= num_bits;

        return static_cast<T>(bit_sequence);
    } else
    {
        /* If the bit sequence does not lay entirely in the current value */
        //(bit_position_ + 1) bits are in the current value and (num_bits - (bit_position_ + 1)) are in the next value 
        /* Fill the sequence at the correct positions with the bits from the current value */
        uint64_t bit_sequence = (current_value_ << (num_bits - bit_position_)) & (~uint64_t{0} >> (64 - num_bits));
        
        /* Move to next value */
        ++pos_;
        this->GetValueAtPos();

        /* Fill the sequence with the second part from the next value */
        bit_sequence |= (current_value_ >> (64 - (num_bits - bit_position_)));

        /* Update the bit position */
        bit_position_ = kBitReadingStartIDX - (num_bits - bit_position_);

        //There is no need to check whether the bit_position_ becomes negative,
        //because due to the design (the bit_position always points to the start of the next bit-sequence,
        //and, therefore, at least one bit needs to reside in the preivous value, such that at max 63 bits
        //can be extracted in the succeeding value)
        
        return static_cast<T>(bit_sequence);
    }
}

template<bool InMemory>
inline void
vector_view_base<InMemory>::SkipNumberOfBits(const size_t num_bits)
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
        this->GetValueAtPos();
    }
}

template<bool InMemory>
inline void
vector_view_base<InMemory>::SkipArbitraryNumberOfBits(const size_t num_bits)
{
    cmc_assert(num_bits > 0);
    int num_bits_to_skip = num_bits;

    while (num_bits_to_skip > 64)
    {
        this->SkipNumberOfBits(64);
        num_bits_to_skip -= 64;
    }

    /* Skip the remainder as well */
    this->SkipNumberOfBits(num_bits_to_skip);
}

}

#endif /* !CMC_BITS_SPAN_HXX */
