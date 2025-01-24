#ifndef CMC_BIT_MAP_HXX
#define CMC_BIT_MAP_HXX


#include "utilities/cmc_utilities.hxx"

#include <vector>
#include <climits>

namespace cmc
{

namespace bit_map
{

/* The amount of bits within a byte */
constexpr int kCharBitAsInt = CHAR_BIT;
constexpr size_t kCharBit = static_cast<size_t>(kCharBitAsInt);

/**
 * @brief A class representing a contiguous bit field with some convenience functions
 */
class BitMap
{
public:

    using byte_iterator = std::vector<uint8_t>::iterator;
    using const_byte_iterator = std::vector<uint8_t>::const_iterator;
    /**
     * @brief A read-only forward-iterator for the BitMap in order to read it's bit field
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
            if (bit_position_ >= kCharBitAsInt)
            {
                ++(this->byte_);
                this->bit_position_ = 0;
            }
        };

        value_type operator*() const { return (*byte_ >> bit_position_) & uint8_t{1}; }
        Iterator& operator++() { ++bit_position_ ; if(bit_position_ >= kCharBitAsInt){bit_position_ = 0; ++byte_;}; return *this; }  
        Iterator operator++(int) { Iterator tmp = *this; ++(*this); return tmp; }
        friend bool operator== (const Iterator& a, const Iterator& b) { return a.byte_ == b.byte_ && a.bit_position_ == b.bit_position_; };
        friend bool operator!= (const Iterator& a, const Iterator& b) { return a.byte_ != b.byte_ || a.bit_position_ != b.bit_position_; };  

    private:
        pointer byte_;
        int bit_position_{0};
    };

    BitMap() = default;
    BitMap(const size_t num_bits)
    : vector_(num_bits / kCharBit + (num_bits % kCharBit != 0 ? 1 : 0), 0), num_bits_{num_bits}{};

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
    bool IsBitSet(const size_t& byte_position, const size_t& bit_position) const;

    void ToggleBit(const size_t global_bit_position);
    void ClearBit(const size_t global_bit_position);
    void SetBit(const size_t global_bit_position);
    bool IsBitSet(const size_t global_bit_position) const;

    Iterator begin() const {return Iterator(vector_.data());};
    const Iterator end() const {return Iterator(vector_.data() + vector_.size() - (num_bits_ % kCharBit != 0 ? 1 : 0), num_bits_ % kCharBit);}

    byte_iterator begin_bytes() { return vector_.begin(); };
    byte_iterator end_bytes() { return vector_.end(); };
    const_byte_iterator begin_bytes() const { return vector_.begin(); };
    const_byte_iterator end_bytes() const { return vector_.end(); };
    const_byte_iterator cbegin_bytes() const { return vector_.cbegin(); };
    const_byte_iterator cend_bytes() const { return vector_.cend(); };

    size_t size() const {return num_bits_;};
    size_t size_bytes() const {return vector_.size();};

    const std::vector<uint8_t>& GetByteData() const {return vector_;};

    friend class BitMapView;
private:
    size_t bit_position_{0};
    size_t byte_position_{0};
    std::vector<uint8_t> vector_{0};
    size_t num_bits_{0};
};

class BitMapView
{
public:
    BitMapView()
    : data_{nullptr}, size_{0} {};
    BitMapView(const uint8_t* data, std::size_t num_bits)
    : data_{data}, size_{num_bits}{};
    BitMapView(const BitMap& bitmap)
    : data_{bitmap.vector_.data()}, size_{bitmap.num_bits_}{};

    bool GetBit(const size_t global_bit_position);

    void MoveToStartBit(const size_t global_bit_position);

    bool GetNextBit();

private:
    const uint8_t* data_;
    std::size_t byte_position_{0};
    std::size_t bit_position_{0};
    std::size_t size_{0};
};

/** BitMapView Member Functions **/
inline bool
BitMapView::GetNextBit()
{
    if (bit_position_ < kCharBit)
    {
        /* Get the current bit of the current byte */
        const bool next_bit = ((data_[byte_position_] >> bit_position_) & uint8_t{1});
        ++bit_position_;
        return next_bit;
    } else
    {
        /* Get the first bit of the next byte */
        ++byte_position_;
        const bool next_bit = (data_[byte_position_] & uint8_t{1});
        bit_position_ = 1;
        return next_bit;
    }
}

inline void
BitMapView::MoveToStartBit(const size_t global_bit_position)
{
    byte_position_ = global_bit_position / kCharBit;
    bit_position_ = global_bit_position % kCharBit;
}

inline bool
BitMapView::GetBit(const size_t global_bit_position)
{
    return ((data_[global_bit_position / kCharBit] >> (global_bit_position % kCharBit)) & uint8_t{1});
}

/** BitMap Member Functions **/
/**
 * @brief Append a single bit to the BitMap
 * 
 * @param bit The bit to be appended (either 'true' == 1 or 'false' == 0)
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

/**
 * @brief Append a 'true'-bit to the BitMap
 */
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

/**
 * @brief Append a 'false'-bit to the BitMap
 */
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

/**
 * @brief Toggle a bit at the given position
 * 
 * @param byte_position The byte-number of the bit to be toggled
 * @param bit_position The intra-byte position of the bit
 */
inline void
BitMap::ToggleBit(const size_t& byte_position, const size_t& bit_position)
{
    cmc_assert(byte_position < vector_.size());
    vector_[byte_position] ^= uint8_t{1} << bit_position;
}

/**
 * @brief Clear a bit at the given position
 * 
 * @param byte_position The byte-number of the bit to be cleared
 * @param bit_position The intra-byte position of the bit
 */
inline void
BitMap::ClearBit(const size_t& byte_position, const size_t& bit_position)
{
    cmc_assert(byte_position < vector_.size());
    vector_[byte_position] &= ~(uint8_t{1} << bit_position);
}

/**
 * @brief Set a bit at the given position
 * 
 * @param byte_position The byte-number of the bit to be set
 * @param bit_position The intra-byte position of the bit
 */
inline void
BitMap::SetBit(const size_t& byte_position, const size_t& bit_position)
{
    cmc_assert(byte_position < vector_.size());
    vector_[byte_position] |= (uint8_t{1} << bit_position);
}

/**
 * @brief Check whether the bit at the given position is set
 * 
 * @param byte_position The byte-number of the bit to be checked
 * @param bit_position The intra-byte position of the bit
 * @return true The bit at this position is set
 * @return false The bit at this position is **not** set
 */
inline bool
BitMap::IsBitSet(const size_t& byte_position, const size_t& bit_position) const
{
    cmc_assert(byte_position < vector_.size());
    return ((vector_[byte_position] >> bit_position) & uint8_t{1});
}

/**
 * @brief Toggle a bit at the given global bit position (0,...,N)
 * 
 * @param bit_position The global bit position
 */
inline void
BitMap::ToggleBit(const size_t global_bit_position)
{
    ToggleBit(global_bit_position / kCharBit, global_bit_position % kCharBit);
}

/**
 * @brief CLear a bit at the given global bit position (0,...,N)
 * 
 * @param bit_position The global bit position
 */
inline void
BitMap::ClearBit(const size_t global_bit_position)
{
    ClearBit(global_bit_position / kCharBit, global_bit_position % kCharBit);
}

/**
 * @brief Set a bit at the given global bit position (0,...,N)
 * 
 * @param bit_position The global bit position
 */
inline void
BitMap::SetBit(const size_t global_bit_position)
{
    SetBit(global_bit_position / kCharBit, global_bit_position % kCharBit);
}

/**
 * @brief Check whether a bit at the given global bit position (0,...,N) is set
 * 
 * @param bit_position The global bit position
 * @return true The bit at this position is set
 * @return false The bit at this position is **not** set
 */
inline bool
BitMap::IsBitSet(const size_t global_bit_position) const
{
    return IsBitSet(global_bit_position / kCharBit, global_bit_position % kCharBit);
}

}

}

#endif /* !CMC_BIT_MAP_HXX */
