#ifndef CMC_BITS_VECTOR_HXX
#define CMC_BITS_VECTOR_HXX

#include "cmc.hxx"
#include "utilities/cmc_bits.hxx"

#include <vector>

namespace cmc::bits
{

/* Forward declaration of the vector_view */
class vector_view;

/* Define some global constants for the BitsVector*/
constexpr int64_t kBitIndexStart = 63;
constexpr int64_t kNextByteMissing = -1;

class vector
{
public:
    template <UnsignedIntegerType T, IntegerType U> void AppendBits(const T value, const U start_pos, const U end_pos);
    
    void AppendBit(const bool bit);
    void AppendSetBit();
    void AppendUnsetBit();

    std::vector<uint8_t> GetSerializedByteStream() const;
    std::vector<uint8_t> GetSerializedByteStreamPadded() const;

    void Reserve(const size_t num_bits);
    size_t size() const;
    size_t size_bytes() const;

    friend class vector_view;
private:
    std::vector<uint64_t> vector_;
    int64_t bit_position_{-1};
};

template<typename T, IntegerType U>
inline int
GetNumSignificantBits(const U start_pos, const U end_pos)
{
    cmc_assert(start_pos + end_pos <= sizeof(T) * kCharBit);
    return sizeof(T) * kCharBit - start_pos - end_pos;
}

inline void
vector::Reserve(const size_t num_bits)
{
    vector_.reserve(num_bits / 64 + 1);
}

/* Append bits from \var value specified by the start and end bit */
template <UnsignedIntegerType T, IntegerType U>
inline void
vector::AppendBits(const T value, const U start_pos, const U end_pos)
{
    cmc_assert(start_pos >= 0 && end_pos >= 0); //Additionally, in case of singed types
    cmc_assert(start_pos <= 64 && end_pos <= 64 && start_pos + end_pos <= 64);

    /* If the value is empty, we do not need to add anything */
    if (start_pos + end_pos >= sizeof(T) * kCharBit) [[unlikely]]
    {
        return;
    }

    /* Check if there needs to be an additional value pushed to the vector before the bits are appended */
    if (bit_position_ < 0) [[unlikely]]
    {
        vector_.emplace_back();
        bit_position_ = kBitIndexStart;
    }

    /* Compute the shift of the start position due to zero-extension during integer promotion */
    const int start_pos_shift = start_pos + (sizeof(uint64_t) - sizeof(T)) * kCharBit;

    /* Get the value as an promoted uint64_t and nullified insignificant bits */
    const uint64_t value_u64 = static_cast<uint64_t>(value) & ((~uint64_t{0} >> start_pos_shift) & (~uint64_t{0} << end_pos));

    /* Check if the significant bits fit into the current value */
    const bool do_bits_fit_into_current_val = ((bit_position_ + 1) >= GetNumSignificantBits<uint64_t>(start_pos_shift, end_pos));

    /* If the bits fit into the current value, we can just append it */
    if (do_bits_fit_into_current_val) [[likely]]
    {
        /* Align the value correctly */
        const uint64_t aligned_value_u64 = (kBitIndexStart - bit_position_ >= start_pos_shift ? 
                                            value_u64 >> ((kBitIndexStart - bit_position_) - start_pos_shift)
                                            : value_u64 << (start_pos_shift - (kBitIndexStart - bit_position_)));

        /* Append the bits from this value */
        vector_.back() |= aligned_value_u64;

        /* Offset the bit position by the number of appended bits */
        bit_position_ -= GetNumSignificantBits<uint64_t>(start_pos_shift, end_pos);
    } else
    {
        /* In case the bits do not fit into the current value, we need to add the significant bits in two phases */
        /* In the first phase, we fit the bits that fit into the current value in place */
        const int right_shift = (kBitIndexStart - bit_position_) - start_pos_shift;
        vector_.back() |= (value_u64 >> right_shift);

        /* In the next phase, we add the new value with the remaining bits set */
        const int left_shift = 64 - (GetNumSignificantBits<uint64_t>(start_pos_shift, end_pos) - (bit_position_ + 1));
        vector_.push_back(value_u64 << left_shift);

        /* We need to offset the position counters */
        bit_position_ = left_shift - 1;
    }
}

inline void
vector::AppendBit(const bool bit)
{
    if (bit_position_ < 0) [[unlikely]]
    {
        vector_.emplace_back();
        bit_position_ = kBitIndexStart;
    }

    vector_.back() |= (static_cast<uint64_t>(bit) << bit_position_);
    --bit_position_;
}

inline void
vector::AppendSetBit()
{
    if (bit_position_ < 0) [[unlikely]]
    {
        vector_.emplace_back();
        bit_position_ = kBitIndexStart;
    }

    vector_.back() |= (uint64_t{1} << bit_position_);
    --bit_position_;
}

inline void
vector::AppendUnsetBit()
{
    if (bit_position_ < 0) [[unlikely]]
    {
        vector_.emplace_back();
        bit_position_ = kBitIndexStart;
    }

    --bit_position_;
}

inline std::vector<uint8_t>
vector::GetSerializedByteStream() const
{
    constexpr size_t type_size = sizeof(uint64_t);

    /* Allocate the serialized byte stream */
    std::vector<uint8_t> serialized_encoding(vector_.size() * type_size);

    /* Iterate through the encoded stream and append the bytes in the correct endianness (big endian) */
    for (size_t idx{0}; idx < vector_.size(); ++idx)
    {
        /* Account for the endianness */
        const auto serialized = SerializeValueBE(vector_[idx]);
        /* Copy the bytes */
        std::copy_n(serialized.data(), type_size, serialized_encoding.data() + idx * type_size);
    }

    /* Potentially, pop back the bytes that are empty at the end */
    const size_t num_bits_to_pop = std::abs(bit_position_ + 1) / kCharBit;
    for (size_t idx{0}; idx < num_bits_to_pop; ++idx)
    {
        serialized_encoding.pop_back();
    }

    /* Return the serialized value */
    return serialized_encoding;
}

inline std::vector<uint8_t>
vector::GetSerializedByteStreamPadded() const
{
    constexpr size_t type_size = sizeof(uint64_t);

    /* Allocate the serialized byte stream */
    std::vector<uint8_t> serialized_encoding(vector_.size() * type_size);

    /* Iterate through the encoded stream and append the bytes in the correct endianness (big endian) */
    for (size_t idx{0}; idx < vector_.size(); ++idx)
    {
        /* Account for the endianness */
        const auto serialized = SerializeValueBE(vector_[idx]);
        /* Copy the bytes */
        std::copy_n(serialized.data(), type_size, serialized_encoding.data() + idx * type_size);
    }

    /* Return the serialized value */
    return serialized_encoding;
}

inline size_t
vector::size() const
{
    return (vector_.size() >= 1 ? (vector_.size() - 1) * sizeof(uint64_t) * kCharBit + (bit_position_ > 0 ? (64 - (bit_position_ + 1)) : 64) : 0);
}

inline size_t
vector::size_bytes() const
{
    return (vector_.size() * sizeof(uint64_t) - (std::abs(bit_position_ + 1) / kCharBit));
}

}

#endif /* !CMC_BITS_VECTOR_HXX */
