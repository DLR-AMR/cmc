#ifndef CMC_BITS_VECTOR_HXX
#define CMC_BITS_VECTOR_HXX

#include "cmc.hxx"
#include "utilities/cmc_endian.hxx"
#include "utilities/cmc_bits.hxx"

#include <vector>
#include <execution>

namespace cmc::bits
{

/* Forward declaration of the vector_view_base */
template<bool InMemory>
class vector_view_base;

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

    std::vector<uint8_t> GetSerializedByteStreamBETrimmed() const;
    std::vector<uint64_t> GetSerializedOffsetByteStreamBE(const int lsb_bit_offset) const;
    std::vector<uint64_t> GetSerializedByteStreamBE() const;

    void Reserve(const size_t num_bits);
    size_t size() const;
    size_t size_bytes() const;

    template<bool InMemory> friend class vector_view_base;
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
vector::GetSerializedByteStreamBETrimmed() const
{
    static_assert(std::endian::native == std::endian::big || std::endian::native == std::endian::little,
                  "Only little-endian and big-endian systems are suppoprted!");
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
    const size_t num_bytes_to_pop = std::abs(bit_position_ + 1) / kCharBit;
    for (size_t idx{0}; idx < num_bytes_to_pop; ++idx)
    {
        serialized_encoding.pop_back();
    }

    /* Return the serialized value */
    return serialized_encoding;
}

inline std::vector<uint64_t>
vector::GetSerializedOffsetByteStreamBE(const int lsb_bit_offset) const
{
    static_assert(std::endian::native == std::endian::big || std::endian::native == std::endian::little,
                  "Only little-endian and big-endian systems are suppoprted!");

    cmc_assert(lsb_bit_offset < sizeof(uint64_t) * cmc::bits::kCharBit && lsb_bit_offset >= 0);
    
    if (this->vector_.empty()) [[unlikely]]
    {
        return this->vector_;
    }

    if (lsb_bit_offset == 0) [[unlikely]]
    {
        return this->GetSerializedByteStreamBE();
    }

    /* Shift parameters */
    const int shift_right = lsb_bit_offset;
    const int shift_left = (cmc::bits::kCharBit * sizeof(uint64_t)) - lsb_bit_offset;

    /* Check if the vector needs to be extended */
    const bool has_vec_to_be_extended = not (bit_position_ + 1 - lsb_bit_offset >= 0);

    /** In this case the offset fits into the current vector as well **/
    /* Allocate the serialized byte stream */
    std::vector<uint64_t> offset_stream(this->vector_.size() + (has_vec_to_be_extended ? 1 : 0));

    for (size_t idx{0}; idx < vector_.size(); ++idx)
    {
        if (idx > 0)
        {
            /* Apply the left shift to the previous value */
            offset_stream[idx - 1] |= ConvertToBigEndian<uint64_t>(this->vector_[idx - 1] << shift_left);
        }
    
        /* Apply the right shfit to the current value */
        offset_stream[idx] = ConvertToBigEndian<uint64_t>(this->vector_[idx] >> shift_right);
    }

    /* Potentially set the new value that has been introduced due to the offset */
    if (has_vec_to_be_extended)
    {
        /* Apply the lkast left-shift to the newly added value*/
        offset_stream.back() = ConvertToBigEndian<uint64_t>(this->vector_.back() << shift_left);
    }

    return offset_stream;

    #if 0
    constexpr size_t type_size = sizeof(uint64_t);

    /* Allocate the serialized byte stream */
    std::vector<uint64_t> offset_stream((this->vector_.size() + 1) * type_size);

    /* Shift parameters */
    const int shift_right = lsb_bit_offset;
    const int shift_left = (cmc::bits::kCharBit * sizeof(uint64_t)) - lsb_bit_offset;

    const size_t num_vals = this->vector_.size();
    
    /** Offset the bit stream **/
    for (size_t idx{0}; idx < num_vals; ++idx)
    {
        /* Add the front part */
        offset_stream[idx] |= (this->vector_[idx] >> shift_right);

        /* Add the back part */
        offset_stream[idx + 1] |= (this->vector_[idx] << shift_left);
    }

    /* Compute the novel bit position with the offset */
    const int offset_bit_idx = this->bit_position_ - lsb_bit_offset;
    if (offset_bit_idx < 0)
    {
        /* In this case, the last value has not hed any data */
        offset_stream.pop_back();
    }

    /** Copy the bytes in big endian **/
    /* Allocate the serialized byte stream */
    std::vector<uint8_t> serialized_encoding(offset_stream.size() * type_size);

    /* Iterate through the encoded stream and append the bytes in the correct endianness (big endian) */
    for (size_t idx{0}; idx < offset_stream.size(); ++idx)
    {
        /* Account for the endianness */
        const auto serialized = SerializeValueBE(offset_stream[idx]);
        /* Copy the bytes */
        std::copy_n(serialized.data(), type_size, serialized_encoding.data() + idx * type_size);
    }

    /** Remove insignificant bytes **/
    const size_t num_bytes_to_pop = (offset_bit_idx + 1 < 0 ? sizeof(uint64_t) - 1 : (offset_bit_idx + 1)/ kCharBit);
    for (size_t idx{0}; idx < num_bytes_to_pop; ++idx)
    {
        serialized_encoding.pop_back();
    }

    /* Return the serialized value */
    return serialized_encoding;
    #endif
}

inline std::vector<uint64_t>
vector::GetSerializedByteStreamBE() const
{
    static_assert(std::endian::native == std::endian::big || std::endian::native == std::endian::little,
                  "Only little-endian and big-endian systems are suppoprted!");

    if constexpr (std::endian::native == std::endian::big)
    {
        /* On a big endian system, we can juste return the encoded vector */
        return vector_;
    } else 
    {
        /* On a little endian machine, we need to swap the bytes */
        std::vector<uint64_t> serialized_encoding(vector_.size());
        std::transform(std::execution::par_unseq, vector_.cbegin(), vector_.cend(), serialized_encoding.begin(), std::byteswap<uint64_t>);
        return serialized_encoding;
    }

    #if 0
    constexpr size_t type_size = sizeof(uint64_t);
    /* Allocate the serialized byte stream */
    std::vector<uint64_t> serialized_encoding(vector_.size() * type_size);

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
    #endif
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
