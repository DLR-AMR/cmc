#ifndef CMC_BIT_VECTOR_HXX
#define CMC_BIT_VECTOR_HXX

#include <vector>
#include <climits>

namespace cmc
{

namespace bit_vector
{

constexpr size_t kCharBit = CHAR_BIT;
constexpr size_t kBitIndexStart = kCharBit - 1;

constexpr uint8_t kOneBitPrefix = 0x00;
constexpr uint8_t kTwoBitPrefix = 0x01;
constexpr uint8_t kNullifyTwoTailBits = ~((uint8_t) 0x03);
constexpr uint8_t kNullifyThirdAndFourthTailBits = ~((uint8_t) 0x0C);
constexpr uint8_t kPrefixContinues = 0x0F;


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
    void AppendTwoBits(const uint8_t byte);
    void AppendFourBits(const uint8_t byte);
    void AppendBits(const uint8_t byte, const int num_bits);
    void AppendBits(const std::vector<uint8_t>& bytes, const int num_bits);

    void ToogleBit(const size_t& byte_position, const size_t& bit_position);
    void ClearBit(const size_t& byte_position, const size_t& bit_position);
    void SetBit(const size_t& byte_position, const size_t& bit_position);
    bool IsBitSet(const size_t& byte_position, const size_t& bit_position);

    void IncrementBitPosition(size_t& iterator, const int diff = 1);
    size_t GetFirstBitPositionOfByte() const;
    bool IsEndOfByteReached(size_t& iterator) const;

    //TODO: Add byte iterator; ...and maybe Bit iterator? ...and PrefixLengthByteIterator?

private:
    size_t bit_position_{kBitIndexStart};
    size_t byte_position_{0};
    std::vector<uint8_t> vector_{0};
};


/**
 * @brief Append a single bit to the BitVector from the given \a byte
 * 
 * @param byte The byte from which the single bit (0b0000000x) will be appended
 */
inline void
BitVector::AppendBit(const uint8_t byte)
{
    if (bit_position_ >= 0)
    {
        /* If the bit fits into the current byte */
        vector_.back() |= ((kNullifyAllExceptTailBit & byte) << bit_position_);
        bit_position_ -= 1;
    } else
    {
        /* If a new byte has to be added which will then hold the given bit */
        vector_.emplace_back(0);
        vector_.back() |= ((kNullifyAllExceptTailBit & byte) << kBitIndexStart);
        bit_position_ = kBitIndexStart - 1;
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
    if (bit_position_ >= 1)
    {
        /* The two bits fit into the current byte */
        vector_.back() |= ((kNullifyAllExceptTwoTailBits & byte) << (bit_position_ - 1));
        bit_position_ -= 2;
    } else if (bit_position_ == 0)
    {
        /* Only a single bit will fit into the current byte, the other bit has to set in the next byte */
        const uint8_t in_byte = (kNullifyAllExceptTwoTailBits & byte);
        vector_.back() |= (in_byte >> 1);
        vector_.emplace_back(0);
        vector_.back() |= (in_byte << kBitIndexStart);
        bit_position_ = kBitIndexStart - 1;
        ++byte_position_;
    } else
    {
        /* A new byte needs to be added and the bits will be assigned there */
        vector_.emplace_back(0);
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
    if (bit_position_ >= 3)
    {
        /* The four bits fit into the current byte */
        vector_.back() |= ((kNullifyFourFrontBits & byte) << (bit_position_ - 3));
        bit_position_ -= 4;
    } else if (bit_position_ < 3 && bit_position_ >= 0)
    {
        /* The bits of the four bits needs to be split across byte borders */
        const uint8_t in_byte = (kNullifyFourFrontBits & byte);
        vector_.back() |= (in_byte >> (3 - bit_position_));
        vector_.emplace_back(0);
        vector_.back() |= (in_byte << (4 + bit_position_ + 1));
        bit_position_ = kBitIndexStart - (bit_position_ + 1);
        ++byte_position_;
    } else
    {
        /* A new byte needs to be added and the bits will be assigned there */
        vector_.emplace_back(0);
        vector_.back() |= ((kNullifyAllExceptTwoTailBits & byte) << (kBitIndexStart - 3));
        bit_position_ = kBitIndexStart - 4;
        ++byte_position_;
    }
}

inline
uint8_t
GetNullifyBitMask(const int num_nonzero_bits)
{
    cmc_assert(num_nonzero_bits >= 0 && num_nonzero_bits <= CHAR_BIT);
    return ((uint8_t) 0xFF >> (kCharBit - num_nonzero_bits));
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
    const size_t num_bits_decremented = static_cast<size_t>(num_bits - 1);

    if (bit_position_ >= num_bits_decremented)
    {
        /* The four bits fit into the current byte */
        vector_.back() |= ((GetNullifyBitMask(num_bits) & byte) << (bit_position_ - num_bits_decremented));
        bit_position_ -= num_bits;
    } else if (bit_position_ < num_bits_decremented && bit_position_ >= 0)
    {
        /* The bits of the four bits needs to be split across byte borders */
        const uint8_t in_byte = (GetNullifyBitMask(num_bits) & byte);
        vector_.back() |= (in_byte >> (num_bits_decremented - bit_position_));
        vector_.emplace_back(0);
        vector_.back() |= (in_byte << (kBitIndexStart - num_bits_decremented + bit_position_ + 1));
        bit_position_ = kBitIndexStart - (bit_position_ + 1);
        ++byte_position_;
    } else
    {
        /* A new byte needs to be added and the bits will be assigned there */
        vector_.emplace_back(0);
        vector_.back() |= ((GetNullifyBitMask(num_bits) & byte) << (kBitIndexStart - num_bits_decremented));
        bit_position_ = kBitIndexStart - num_bits;
        ++byte_position_;
    }
}

/**
 * @brief It is assumed that the vector hold bits that are aligned to the "most significant" Bit.
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
    int num_bits_appended = 0;

    for (auto byte_iter = bytes.begin(); byte_iter != bytes.end(); ++byte_iter)
    {
        /* Compute the amount of bits we are appending to the vector within this iteration */
        const int bits_to_append = (num_bits - num_bits_appended) > CHAR_BIT ? CHAR_BIT : (num_bits - num_bits_appended);

        /* We need to tail-align the bits in order to use the given AppendBits(...)-function.
         * In case we are appending a full byte, the shift-operation is a no-op */
        const uint8_t byte_to_append = (*byte_iter) >> (CHAR_BIT - bits_to_append);

        /* Append the bits from this iteration */
        AppendBits(byte_to_append, bits_to_append);

        /* Add that we have assigned the bits from this iteration */
        num_bits_appended += CHAR_BIT;
    }
}

inline
void
BitVector::IncrementBitPosition(size_t& iterator, const int diff)
{
    iterator -= diff;
}

inline
size_t
BitVector::GetFirstBitPositionOfByte() const
{
    return kBitIndexStart;
}

inline
bool
BitVector::IsEndOfByteReached(size_t& iterator) const
{
    return (iterator >= 0 ? false : true);
}

}

}

#endif /* !CMC_BIT_VECTOR_HXX */
