#ifndef CMC_BIT_VECTOR_HXX
#define CMC_BIT_VECTOR_HXX

#include "utilities/cmc_log_functions.hxx"

#include <vector>
#include <climits>

namespace cmc
{

namespace bit_vector
{

constexpr size_t kCharBit = CHAR_BIT;
constexpr size_t kBitIndexStart = kCharBit - 1;
constexpr int kCharBitInteger = CHAR_BIT;

constexpr uint8_t kNullifyAllExceptTwoTailBits = (uint8_t) 0x03;
constexpr uint8_t kNullifyFourFrontBits = (uint8_t) 0x0F;
constexpr uint8_t kNullifyAllExceptTailBit = (uint8_t) 0x01;


/**
 * @brief The vector traverses the bits from high to low [7,...,0]b bits. It always starts at the "most significant bit" of a byte.
 * 
 * The assumption is, that the layout of the bytes/bits looks like:
 * 
 *     Bit Positions: 0 to 7      
 *  [ 7, 6, 5, 4, 3, 2, 1, 0 ]  [ 7, 6, 5, 4, 3, 2, 1, 0 ]  [ 7, 6, 5, 4, 3, 2, 1, 0 ]  ...
 *  |-------- Byte 0 --------|  |-------- Byte 1 --------|  |-------- Byte 2 --------|  ...
 *    ^ Traversal starts here an goes onwards ->
 *  
 */
class BitVector
{
public:
    using iterator = std::vector<uint8_t>::iterator;
    using const_iterator = std::vector<uint8_t>::const_iterator;

    BitVector() = default;
    BitVector(const size_t num_bytes)
    : vector_(num_bytes, 0){};
    BitVector(std::vector<uint8_t>&& bytes)
    : vector_{std::move(bytes)}{};

    BitVector(const BitVector& other) = default;
    BitVector& operator=(const BitVector& other) = default;
    BitVector(BitVector&& other) = default;
    BitVector& operator=(BitVector&& other) = default;

    ~BitVector() = default;
    void Reserve(const size_t num_bytes)
    {
        vector_.reserve(num_bytes);
    }

    void AppendBit(const uint8_t byte);
    void AppendBit(const bool is_bit_set);
    void AppendTwoBits(const uint8_t byte);
    void AppendFourBits(const uint8_t byte);
    void AppendBits(const uint8_t byte, const int num_bits);
    void AppendBits(const std::vector<uint8_t>& bytes, const int num_bits);

    void AppendBits(const BitVector& bit_vector);

    template<int N> void AppendBytes(const std::array<uint8_t, N>& byte_stream);
    void AppendBytes(const std::vector<uint8_t>& byte_stream);

    void AddPaddingToFullByte();

    void MoveDataInto(std::vector<uint8_t>& vector);

    #if 0
    void IncrementBitPosition(size_t& iterator, const int diff = 1);
    size_t GetFirstBitPositionOfByte() const;
    bool IsEndOfByteReached(size_t& iterator) const;
    #endif

    iterator begin() { return vector_.begin(); };
    iterator end() { return vector_.end(); };
    const_iterator begin() const { return vector_.begin(); };
    const_iterator end() const { return vector_.end(); };
    const_iterator cbegin() const { return vector_.cbegin(); };
    const_iterator cend() const { return vector_.cend(); };

    size_t size() const { return vector_.size();};
    size_t size_bits() const { return vector_.size() * kCharBit - bit_position_ - 1;};

    const uint8_t* data() const {return vector_.data();};
    std::pair<std::vector<uint8_t>, size_t> GetBits() const {return std::make_pair(vector_, this->size_bits());};
    const std::vector<uint8_t>& GetData() const {return vector_;};

    void TrimToContent() {if (bit_position_ == kBitIndexStart) {/* If the last byte is empty, we remove it */ vector_.pop_back();}};

private:
    size_t bit_position_{kBitIndexStart};
    size_t byte_position_{0};
    std::vector<uint8_t> vector_{0};
};

inline
uint8_t
GetNullifyBitMaskLSB(const int num_nonzero_bits)
{
    cmc_assert(num_nonzero_bits >= 0 && num_nonzero_bits <= CHAR_BIT);
    return (uint8_t{0xFF} << num_nonzero_bits);
}

class BitVectorView
{
public:
    BitVectorView()
    : data_{nullptr}, num_bytes_{0} {};
    BitVectorView(const uint8_t* data, const size_t num_bytes)
    : data_{data}, num_bytes_{num_bytes} {};

    const uint8_t* begin() { return data_; };
    const uint8_t* end() { return data_ + num_bytes_; };

    void MoveToStartBit(const size_t byte_position, const size_t bit_position)
    {
        byte_position_ = byte_position;
        bit_position_ = bit_position;
    }

    bool IsCurrentBitSet() const
    {
        return ((data_[byte_position_] >> bit_position_) & uint8_t{1});
    }

    void MoveToNextBit()
    {
        if (bit_position_ > 0)
        {
            --bit_position_;
        } else
        {
            ++byte_position_;
            bit_position_  = kBitIndexStart;
        }
    }
    
    void MoveToNextByte()
    {
        ++byte_position_;
        bit_position_ = kBitIndexStart;   
    }
    
    std::vector<uint8_t> GetNextBitSequence(const size_t num_bits)
    {
        if (num_bits == 0) {return std::vector<uint8_t>();}
        
        /* Number of bytes and bits needed to be extracted */
        const size_t num_bytes_ = num_bits / kCharBit;
        const size_t num_bits_ = num_bits % kCharBit;

        /* Save the start byte position */
        const int starting_byte_position = byte_position_;

        /* Number of bytes to allocate */
        const size_t num_bytes = num_bytes_ + (num_bits_ != 0 ? 1 : 0);
        std::vector<uint8_t> bit_sequence(num_bytes, 0);

        const int bit_shift_to_front = kBitIndexStart - bit_position_;
        const int bit_shift_to_back = kCharBitInteger - bit_shift_to_front;

        for (size_t num_bits_extracted = 0, byte_id = 0; num_bits_extracted < num_bits; ++byte_id, ++byte_position_)
        {
            if (bit_shift_to_front > 0 && num_bits_extracted != 0)
            {
                bit_sequence[byte_id - 1] |= (data_[byte_position_] >> bit_shift_to_back);
                num_bits_extracted += (kCharBitInteger - bit_shift_to_back);

                if (num_bits_extracted >= num_bits){break;}
            }

            bit_sequence[byte_id] |= (data_[byte_position_] << bit_shift_to_front);

            num_bits_extracted += (kCharBitInteger - bit_shift_to_front);
        }

        /* Update the byte and bit position correctly  */
        byte_position_ = starting_byte_position + (num_bits_ <= bit_position_ ? 0 : 1) + num_bytes_;
        //bit_position_ += (num_bits_ <= bit_position_ ? -num_bits_ : kCharBit - (num_bits_ - bit_position_));

        if (num_bits_ <= bit_position_)
        {
            bit_position_ = bit_position_ - num_bits_;
        } else
        {
            bit_position_ = kCharBit - (num_bits_ - bit_position_);
        }

        /* Potentially, we need to erase the lower bits from the last extracted byte */
        if (num_bits_ != 0)
        {
            const int bits_to_nullify = kCharBitInteger - num_bits_;
            bit_sequence.back() &= GetNullifyBitMaskLSB(bits_to_nullify);
        }

        return bit_sequence;
    }

private:
    const uint8_t* data_;
    std::size_t byte_position_{0};
    std::size_t bit_position_{kBitIndexStart};
    size_t num_bytes_;
};

/**
 * @brief Append a single bit to the BitVector from the given \a byte
 * 
 * @param byte The byte from which the single bit (0b0000000x) will be appended
 */
inline void
BitVector::AppendBit(const uint8_t byte)
{
    if (bit_position_ > 0)
    {
        /* If the bit fits into the current byte */
        vector_.back() |= ((kNullifyAllExceptTailBit & byte) << bit_position_);
        bit_position_ -= 1;
    } else
    {
        /* If the single bit fills the current byte */
        vector_.back() |= (kNullifyAllExceptTailBit & byte);
        vector_.emplace_back(0);
        bit_position_ = kBitIndexStart;
        ++byte_position_;
    }
}

inline void
BitVector::AppendBit(const bool is_bit_set)
{
    if (bit_position_ > 0)
    {
        /* If the bit fits into the current byte */
        if (is_bit_set)
        {
            vector_.back() |= (uint8_t{1} << bit_position_);
        }
        --bit_position_;
    } else
    {
        /* If the single bit fills the current byte */
        if (is_bit_set)
        {
            vector_.back() |= uint8_t{1};
        }
        vector_.emplace_back(0);
        bit_position_ = kBitIndexStart;
        ++byte_position_;
    }
}

/**
 * @brief Append two bits from the given \a byte to the vector 
 * 
 * @param byte The byte from which the two bits (0b000000xx) are appended
 */
inline
void
BitVector::AppendTwoBits(const uint8_t byte)
{
    if (bit_position_ > 1)
    {
        /* The two bits fit into the current byte */
        vector_.back() |= ((kNullifyAllExceptTwoTailBits & byte) << (bit_position_ - 1));
        bit_position_ -= 2;
    } else if (bit_position_ == 1)
    {
        /* The two bits fit into the current byte */
        vector_.back() |= (kNullifyAllExceptTwoTailBits & byte);
        vector_.push_back(0);
        bit_position_ = kBitIndexStart;
        ++byte_position_;
    }
    else if (bit_position_ == 0)
    {
        /* Only a single bit will fit into the current byte, the other bit has to set in the next byte */
        const uint8_t in_byte = (kNullifyAllExceptTwoTailBits & byte);
        vector_.back() |= (in_byte >> 1);
        vector_.push_back(0);
        vector_.back() |= (in_byte << kBitIndexStart);
        bit_position_ = kBitIndexStart - 1;
        ++byte_position_;
    } else
    {
        /* A new byte needs to be added and the bits will be assigned there */
        vector_.push_back(0);
        vector_.back() |= ((kNullifyAllExceptTwoTailBits & byte) << (kBitIndexStart - 1));
        bit_position_ = kBitIndexStart - 2;
        ++byte_position_;
    }
}

/**
 * @brief Append four bits from the given \a byte to the vector 
 * 
 * @param byte The byte from which the four bits (0b0000xxxx) are appended
 */
inline
void 
BitVector::AppendFourBits(const uint8_t byte)
{
    if (bit_position_ > 3)
    {
        /* The four bits fit into the current byte */
        vector_.back() |= ((kNullifyFourFrontBits & byte) << (bit_position_ - 3));
        bit_position_ -= 4;
    } else if (bit_position_ == 3)
    {
        /* The four bits fit into the current byte */
        vector_.back() |= ((kNullifyFourFrontBits & byte) << (bit_position_ - 3));
        vector_.push_back(0);
        bit_position_ = kBitIndexStart;
        ++byte_position_;
    }
    else if (bit_position_ < 3)
    {
        /* The bits of the four bits needs to be split across byte borders */
        const uint8_t in_byte = (kNullifyFourFrontBits & byte);
        vector_.back() |= (in_byte >> (3 - bit_position_));
        vector_.push_back(0);
        vector_.back() |= (in_byte << (4 + bit_position_ + 1));
        bit_position_ = kBitIndexStart - (3 - bit_position_);
        ++byte_position_;
    } else
    {
        /* A new byte needs to be added and the bits will be assigned there */
        vector_.push_back(0);
        vector_.back() |= ((kNullifyAllExceptTwoTailBits & byte) << (kBitIndexStart - 3));
        bit_position_ = kBitIndexStart - 4;
        ++byte_position_;
    }
}

inline void
BitVector::AppendBits(const BitVector& bit_vec)
{
    //cmc_assert(bit_vec.size_bits() > 0);
    if (bit_vec.size_bits() == 0) {return;}

    //cmc_debug_msg("Size of bit vector to be appended: ", bit_vec.size(), " in bits: ", bit_vec.size_bits());
    /* We copy all complete bytes */
    for (auto byte_iter = bit_vec.begin(); byte_iter != std::prev(bit_vec.end()); ++byte_iter)
    {
        const uint8_t byte = *byte_iter;
        this->AppendBits(byte, kCharBitInteger);
    }

    /* Additionally, we copy the remaining bits from the last byte */
    if (bit_vec.bit_position_ != kBitIndexStart)
    {
        const uint8_t last_byte = (*std::prev(bit_vec.end())) >> bit_vec.bit_position_;
        const int num_bits = kBitIndexStart - bit_vec.bit_position_;
        //cmc_debug_msg("IN AppendBits(BitVec): num_bits = ", num_bits, ", und bit position: ", bit_vec.bit_position_); 
        this->AppendBits(last_byte, num_bits);
    }
}

inline
uint8_t
GetNullifyBitMaskMSB(const int num_nonzero_bits)
{
    //cmc_debug_msg("In GetNullifyBitMASkMSB: num_nonzero_bits: ", num_nonzero_bits);
    cmc_assert(num_nonzero_bits >= 0 && num_nonzero_bits <= CHAR_BIT);
    return (uint8_t{0xFF} >> (kCharBit - num_nonzero_bits));
}

/**
 * @brief Append a given number of bits from the given \a byte to the vector 
 * 
 * @param byte The byte from which the four bits (0b00xxxxxx) are appended
 * @param num_bits The number of bits to append  (^ e.g. above 6 Bits)
 */
inline
void
BitVector::AppendBits(const uint8_t byte, const int num_bits)
{
    cmc_assert(num_bits > 0);
    if (num_bits <= 0) {return;}

    const size_t num_bits_decremented = static_cast<size_t>(num_bits - 1);

    if (bit_position_ > num_bits_decremented)
    {
        /* The bits fit into the current byte */
        vector_.back() |= ((GetNullifyBitMaskMSB(num_bits) & byte) << (bit_position_ - num_bits_decremented));
        bit_position_ -= num_bits;
    } else if (bit_position_ == num_bits_decremented)
    {
        vector_.back() |= ((GetNullifyBitMaskMSB(num_bits) & byte) << (bit_position_ - num_bits_decremented));
        vector_.push_back(0);
        bit_position_ = kBitIndexStart;
        ++byte_position_;
    } else if (bit_position_ < num_bits_decremented)
    {
        /* The bits need to be split across byte borders */
        const uint8_t in_byte = (GetNullifyBitMaskMSB(num_bits) & byte);
        vector_.back() |= (in_byte >> (num_bits_decremented - bit_position_));
        vector_.push_back(0);
        vector_.back() |= (in_byte << (kBitIndexStart - num_bits_decremented + bit_position_ + 1));
        bit_position_ = kBitIndexStart - (num_bits_decremented - bit_position_);
        ++byte_position_;
    }

    //cmc_debug_msg("In AppendedBits(byte, num_bits): num_bits was: ", num_bits, ", new bit position_ = ", bit_position_);
}

/**
 * @brief It is assumed that the vector hold bits that are aligned to the "most significant" bit.
 * The bits are appended to the BitVector as follows:
 * 
 *                  BitVector           \a bytes
 * Bit Position: [7,6,5,...,1,0]  |  [ 7,6,5,...,1,0] [7,6,5,...,1,0]
 *                                ^ append \a bytes here
 * 
 * @param bytes A vector of (front-shifted) bytes which shall be appended to the vector
 * @param num_bits The amount of bits to be appended
 */
inline
void
BitVector::AppendBits(const std::vector<uint8_t>& bytes, const int num_bits)
{
    cmc_assert(num_bits > 0);
    if (num_bits <= 0) {return;}

    int num_bits_appended = 0;

    for (auto byte_iter = bytes.begin(); byte_iter != bytes.end(); ++byte_iter)
    {
        /* Compute the amount of bits we are appending to the vector within this iteration */
        const int bits_to_append = (num_bits - num_bits_appended) >= CHAR_BIT ? CHAR_BIT : (num_bits - num_bits_appended);

        /* It may be possible that the vector is larger than the number of bits to append */
        if (bits_to_append > 0)
        {
            /* We need to tail-align the bits in order to use the given AppendBits(...)-function.
             * In case we are appending a full byte, the shift-operation is a no-op */
            const uint8_t byte_to_append = (*byte_iter) >> (CHAR_BIT - bits_to_append);

            /* Append the bits from this iteration */
            AppendBits(byte_to_append, bits_to_append);

            /* Add that we have assigned the bits from this iteration */
            num_bits_appended += bits_to_append;
        } else
        {
            break;
        }
    }
}

template<int N>
inline void
BitVector::AppendBytes(const std::array<uint8_t, N>& byte_stream)
{
    for (int idx = 0; idx < N; ++idx)
    {
        AppendBits(byte_stream[idx], kCharBit);
    }
}

inline void
BitVector::AppendBytes(const std::vector<uint8_t>& byte_stream)
{
    for (auto byte_iter = byte_stream.begin(); byte_iter != byte_stream.end(); ++byte_iter)
    {
        AppendBits(*byte_iter, kCharBit);
    }
}

inline void
BitVector::AddPaddingToFullByte()
{
    if (bit_position_ != kBitIndexStart)
    {
        vector_.emplace_back(0);
        bit_position_ = kBitIndexStart;
        ++byte_position_;
    }
}

inline void
BitVector::MoveDataInto(std::vector<uint8_t>& vector)
{
    vector  = std::move(vector_);
    *this = BitVector();
}


}

}

#endif /* !CMC_BIT_VECTOR_HXX */
