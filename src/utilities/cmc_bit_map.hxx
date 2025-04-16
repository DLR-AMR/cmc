#ifndef CMC_BIT_MAP_HXX
#define CMC_BIT_MAP_HXX


#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_log_functions.hxx"

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

    BitMap(const std::vector<uint8_t>& bitmap_conform_bitfield, const size_t num_bits)
    : vector_(bitmap_conform_bitfield), num_bits_{num_bits}{
        if (bitmap_conform_bitfield.size() != (num_bits / kCharBit + (num_bits % kCharBit != 0 ? 1 : 0)) || num_bits == 0)
        {
            throw std::invalid_argument("The supplied vector size does not match the amount of bits or the amount of bits is zero.");
        }
        byte_position_ = vector_.size() - 1;
        bit_position_ = (num_bits % kCharBit == 0 ? kCharBit : (num_bits % kCharBit));
    };
    BitMap(std::vector<uint8_t>&& bitmap_conform_bitfield, const size_t num_bits)
    : vector_(std::move(bitmap_conform_bitfield)), num_bits_{num_bits}{
        if (bitmap_conform_bitfield.size() != (num_bits / kCharBit + (num_bits % kCharBit != 0 ? 1 : 0)) || num_bits == 0)
        {
            throw std::invalid_argument("The supplied vector size does not match the amount of bits or the amount of bits is zero.");
        }
        byte_position_ = vector_.size() - 1;
        bit_position_ = (num_bits % kCharBit == 0 ? kCharBit : (num_bits % kCharBit));
    };

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
    bool IsEmpty() const {return num_bits_ == 0;}; 

    const uint8_t* data() const {return vector_.data();};

    uint8_t* data_overwrite() {return &vector_[0];}

    void AppendBits(const BitMap& bitmap);

    const std::vector<uint8_t>& GetByteData() const {return vector_;};
    void MoveDataInto(std::vector<uint8_t>& vector);
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
    
    std::vector<uint8_t> GetNextNumberOfBits(const size_t num_bits);

    void MoveToNextByte();

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

inline
uint8_t
GetNullifyBitMaskMSB(const int num_nonzero_bits)
{
    cmc_assert(num_nonzero_bits >= 0 && num_nonzero_bits <= kCharBitAsInt);
    return (uint8_t{0xFF} >> (kCharBit - num_nonzero_bits));
}

inline std::vector<uint8_t>
BitMapView::GetNextNumberOfBits(const size_t num_bits)
{
    if (num_bits == 0) {return std::vector<uint8_t>();}
        
    /* Potentially, we need to set the pointer to the next bit correctly */
    if (bit_position_ == kCharBit) {++byte_position_; bit_position_ = 0;}

    /* Number of bytes and bits needed to be extracted */
    const size_t num_bytes_ = num_bits / kCharBit;
    const size_t num_bits_ = num_bits % kCharBit;

    /* Save the start byte position */
    const int starting_byte_position = byte_position_;

    /* Number of bytes to allocate */
    const size_t num_bytes = num_bytes_ + (num_bits_ != 0 ? 1 : 0);
    std::vector<uint8_t> bit_sequence(num_bytes, 0);

    /* Define the shift operations for the extraction */
    const int bit_shift_to_front = kCharBitAsInt - static_cast<int>(bit_position_);
    const int bit_shift_to_back = bit_position_;

    /* Get the bits in the correct order and store them in the output bit sequence */
    for (size_t num_bits_extracted = 0, byte_id = 0; num_bits_extracted < num_bits; ++byte_id, ++byte_position_)
    {
        if (bit_shift_to_back > 0 && num_bits_extracted != 0)
        {
            bit_sequence[byte_id - 1] |= (data_[byte_position_] << bit_shift_to_front);
            num_bits_extracted += (kCharBitAsInt - bit_shift_to_front);

            if (num_bits_extracted >= num_bits){break;}
        }

        bit_sequence[byte_id] |= (data_[byte_position_] >> bit_shift_to_back);

        num_bits_extracted += (kCharBitAsInt - bit_shift_to_back);
    }

    /* Update the byte and bit position of the view correctly */
    byte_position_ = starting_byte_position + num_bytes_ + (num_bits_ + bit_position_ < kCharBit ? 0 : 1);
    bit_position_ = (bit_position_ + num_bits_) % kCharBit;

    /* Potentially, we need to erase the upper bits from the last extracted byte */
    if (num_bits_ != 0)
    {
        bit_sequence.back() &= GetNullifyBitMaskMSB(num_bits_);
    }

    return bit_sequence;
}

inline void
BitMapView::MoveToNextByte()
{
    ++byte_position_;
    bit_position_ = 0;
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
 * @brief Append the bits from another bitmap 
 * 
 * @param bitmap The bitmap whose bits will be appended to this bitmap 
 */
inline void
BitMap::AppendBits(const BitMap& bitmap)
{
    if (&bitmap == this)
    {
        throw std::invalid_argument("The BitMap cannot append it's content to itself.");
    }

    if (bitmap.size() == 0){return;}

    /* Potentially, prepare the setup to start the copying process */
    if (bit_position_ >= kCharBit)
    {
        vector_.emplace_back(0);
        ++byte_position_;
        bit_position_ = 0;
    } else if (bit_position_ == 0 && num_bits_ == 0)
    {
        /* In this case, we remove the preset byte in the vector, because it will be added later */
        vector_.pop_back();
    }

    /* Define the shift operations */
    const int shift_to_back = kCharBit - bit_position_;
    const int shift_to_front = bit_position_;

    /* Append all bytes from the bitmap to the vector */
    for (auto byte_iter = bitmap.begin_bytes(); byte_iter != bitmap.end_bytes(); ++byte_iter)
    {
        /* Get the current byte */
        uint8_t byte = *byte_iter;

        if (shift_to_front == 0)
        {
            vector_.emplace_back(0);
        }

        vector_.back() |= (byte << shift_to_front);

        if (shift_to_front > 0)
        {
            vector_.emplace_back(0);
            vector_.back() |= (byte >> shift_to_back);
        }
    }

    /* Potentially, we need to remove the last byte, since it may be empty */
    if (shift_to_front > 0 && (((bitmap.size() % kCharBit) + bit_position_ <=  kCharBit)))
    {
        vector_.pop_back();
    }

    /* Set the byte and bit indicators correctly */
    byte_position_ = vector_.size() - 1;
    bit_position_ = ((bitmap.size() % kCharBit) + bit_position_  == kCharBit ? kCharBit : (bit_position_ + bitmap.size()) % kCharBit);

    /* Update the bit count */
    num_bits_ += bitmap.size();
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

inline void
BitMap::MoveDataInto(std::vector<uint8_t>& vector)
{
    vector  = std::move(vector_);
    *this = BitMap();
}

}

}

#endif /* !CMC_BIT_MAP_HXX */
