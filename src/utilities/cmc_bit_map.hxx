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

    /**
     * @brief A read-only forward-iterator for the BitMap in order to read it's bit field
     * 
     */
    struct Iterator 
    {
        using iterator_category = std::forward_iterator_tag;
        using difference_type   = std::ptrdiff_t;
        using value_type        = bool;
        using pointer           = const uint8_t*;
        using reference         = const uint8_t&;

        Iterator(pointer byte) : byte_(byte) {};
        Iterator(pointer byte, const int bit_position_constraint) : byte_{byte}, bit_position_{bit_position_constraint} {
            if (bit_position_ >= kCharBit)
            {
                ++(this->byte_);
                this->bit_position_ = 0;
            }
        };

        value_type operator*() const { return (*byte_ >> bit_position_) & uint8_t{1}; }
        Iterator& operator++() { ++bit_position_ ; if(bit_position_ >= kCharBit){bit_position_ = 0; ++byte_;}; return *this; }  
        Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
        friend bool operator== (const Iterator& a, const Iterator& b) { return a.byte_ == b.byte_ && a.bit_position_ == b.bit_position_; };
        friend bool operator!= (const Iterator& a, const Iterator& b) { return a.byte_ != b.byte_ || a.bit_position_ != b.bit_position_; };  

    private:
        pointer byte_;
        int bit_position_{0};
    };

    BitMap() = default;
    BitMap(const size_t num_bits)
    : vector_(num_bits / kCharBit + 1, 0), num_bits_{num_bits}{};

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

    void ToggleBit(const size_t& byte_position, const size_t& bit_position);
    void ClearBit(const size_t& byte_position, const size_t& bit_position);
    void SetBit(const size_t& byte_position, const size_t& bit_position);
    bool IsBitSet(const size_t& byte_position, const size_t& bit_position);

    void ToggleBit(const size_t global_bit_position);
    void ClearBit(const size_t global_bit_position);
    void SetBit(const size_t global_bit_position);
    bool IsBitSet(const size_t global_bit_position);

    Iterator begin() const {return Iterator(vector_.data());};
    const Iterator end() const {return Iterator(vector_.data() + vector_.size() - 1, num_bits_ % kCharBit);}

    size_t size() const {return num_bits_;};

private:
    size_t bit_position_{0};
    size_t byte_position_{0};
    std::vector<uint8_t> vector_{0};
    size_t num_bits_{0};
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
    ++num_bits_;
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
    ++num_bits_;
}

inline void
BitMap::ToggleBit(const size_t& byte_position, const size_t& bit_position)
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
BitMap::ToggleBit(const size_t global_bit_position)
{
    ToggleBit(global_bit_position / kCharBit, global_bit_position % kCharBit);
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
