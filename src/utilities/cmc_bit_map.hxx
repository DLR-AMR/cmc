#ifndef CMC_BIT_MAP_HXX
#define CMC_BIT_MAP_HXX


#include "utilities/cmc_utilities.hxx"

#include <vector>
#include <limits>
#include <climits>

namespace cmc
{

namespace bit_map
{

constexpr size_t kCharBit = CHAR_BIT;


class BitMap
{
public:
    BitMap() = default;
    BitMap(const size_t num_bits)
    : vector_(num_bits / kCharBit + 1, 0){};

    BitMap(const BitMap& other) = default;
    BitMap& operator=(const BitMap& other) = default;
    BitMap(BitMap&& other) = default;
    BitMap& operator=(BitMap&& other) = default;

    ~BitMap() = default;
    void Reserve(const size_t num_bits)
    {
        vector_.reserve(num_bits / kCharBit + 1);
    }

    void AppendBit(const bool bit);
    void AppendSetBit();
    void AppendUnsetBit();

    void ToogleBit(const size_t& byte_position, const size_t& bit_position);
    void ClearBit(const size_t& byte_position, const size_t& bit_position);
    void SetBit(const size_t& byte_position, const size_t& bit_position);
    bool IsBitSet(const size_t& byte_position, const size_t& bit_position);

    void ToogleBit(const size_t global_bit_position);
    void ClearBit(const size_t global_bit_position);
    void SetBit(const size_t global_bit_position);
    bool IsBitSet(const size_t global_bit_position);


    //TODO: Add byte iterator; ...and maybe Bit iterator? ...and PrefixLengthByteIterator?

private:
    size_t bit_position_{0};
    size_t byte_position_{0};
    std::vector<uint8_t> vector_{0};
};


/**
 * @brief Append a single bit to the BitMap from the given \a byte
 * 
 * @param bit The bit to be appended
 */
inline void
BitMap::AppendBit(const bool bit)
{
    if (bit == true)
    {
        AppendSetBit();
    } else
    {
        AppendUnsetBit();
    }
}

inline void
BitMap::AppendSetBit()
{
    if (bit_position_ < kCharBit)
    {
        /* If the bit fits into the current byte */
        vector_.back() |= (uint8_t{1} << bit_position_);
        ++bit_position_;
    } else
    {
        /* If a new byte has to be added which will then hold the given bit */
        vector_.emplace_back(1);
        bit_position_ = 1;
        ++byte_position_;
    }
}

inline void
BitMap::AppendUnsetBit()
{
    if (bit_position_ < kCharBit)
    {
        /* If the bit fits into the current byte */
        ++bit_position_;
    } else
    {
        /* If a new byte has to be added which will then hold the given bit */
        vector_.emplace_back(0);
        bit_position_ = 1;
        ++byte_position_;
    }
}

inline void
BitMap::ToogleBit(const size_t& byte_position, const size_t& bit_position)
{
    cmc_assert(byte_position < vector_.size());
    vector_[byte_position] ^= uint8_t{1} << bit_position;
}

inline void
BitMap::ClearBit(const size_t& byte_position, const size_t& bit_position)
{
    cmc_assert(byte_position < vector_.size());
    vector_[byte_position] &= ~(uint8_t{1} << bit_position);
}

inline void
BitMap::SetBit(const size_t& byte_position, const size_t& bit_position)
{
    cmc_assert(byte_position < vector_.size());
    vector_[byte_position] |= (uint8_t{1} << bit_position);
}

inline bool
BitMap::IsBitSet(const size_t& byte_position, const size_t& bit_position)
{
    cmc_assert(byte_position < vector_.size());
    return ((vector_[byte_position] >> bit_position) & uint8_t{1});
}

inline void
BitMap::ToogleBit(const size_t global_bit_position)
{
    ToogleBit(global_bit_position / kCharBit, global_bit_position % kCharBit);
}

inline void
BitMap::ClearBit(const size_t global_bit_position)
{
    ClearBit(global_bit_position / kCharBit, global_bit_position % kCharBit);
}

inline void
BitMap::SetBit(const size_t global_bit_position)
{
    SetBit(global_bit_position / kCharBit, global_bit_position % kCharBit);
}

inline bool
BitMap::IsBitSet(const size_t global_bit_position)
{
    return IsBitSet(global_bit_position / kCharBit, global_bit_position % kCharBit);
}

}

}

#endif /* !CMC_BIT_MAP_HXX */
