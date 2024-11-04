#ifndef CMC_PREFIX_ENCODING_HXX
#define CMC_PREFIX_ENCODING_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_prefix.hxx"
#include "utilities/cmc_span.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_span.hxx"

#include "utilities/cmc_log_functions.hxx"

#include <vector>
#include <utility>
#include <type_traits>
#include <cstring>


//Delete later
#include <bitset>

namespace cmc
{

constexpr int kTwoBitPrefixEncodingOffset = 0b01;
constexpr int kFourBitPrefixEncodingStart = 0b1000;
constexpr int kFourBitPrefixEncodingOffset = 0b0011;
constexpr uint8_t kIndicatePrefixLongerThanByte = 0b1110;
constexpr int kMoreThanFourBitPrefixEncodingOffset = 0b1000;
constexpr uint8_t kSkipToNextByteInEncodingStream = 0b1111;

constexpr int kShiftTwoBitPrefixLength = bit_vector::kCharBit - 2;
constexpr int kShiftFourBitPrefixLength = bit_vector::kCharBit - 4;
constexpr size_t kTwoBitCodeSize = 2;
constexpr size_t kFourBitCodeSize = 4;
constexpr uint8_t kMaxOneBytePrefixEncoding = 0b1101;

template<typename T>
class LevelwisePrefixData
{
public:
    LevelwisePrefixData() = delete;

    LevelwisePrefixData(const int size_hint_num_elements)
    {
        /* Allocate memory given the supplied size hint */
        refinement_indicator.Reserve(size_hint_num_elements);
        prefix_indicator.Reserve(size_hint_num_elements);
        prefixes.reserve(size_hint_num_elements);
    };

    /* Set the next bit for the prefix indication */
    void SetPrefixIndicatorBit(const bool has_prefix_been_found)
    {
        prefix_indicator.AppendBit(has_prefix_been_found);
    }

    /* Set the next bit for the refinement indication*/
    void SetRefinementIndicatorBit(const bool is_refined)
    {
        refinement_indicator.AppendBit(is_refined);
    }

    /* Append a common prefix */
    void SetPrefix(const CompressionValue<sizeof(T)>& prefix)
    {
        prefixes.push_back(prefix);
    }

    /* Append a common prefix */
    void SetPrefix(CompressionValue<sizeof(T)>&& prefix)
    {
        prefixes.push_back(std::move(prefix));
    }

    std::vector<uint8_t> EncodeLevelData() const;

    bit_map::BitMap refinement_indicator;
    bit_map::BitMap prefix_indicator;
    std::vector<CompressionValue<sizeof(T)>> prefixes;
};

inline void
EncodeAndAppendPrefixLength(bit_vector::BitVector& encoded_prefix_lengths, const int prefix_length)
{
    cmc_assert(prefix_length >= 1);

    //cmc_debug_msg("encode prefix length call with: ", prefix_length);

    if (prefix_length <= 2)
    {
        /* If the prefix is of length one or two, we need two bits to encode the prefix */
        encoded_prefix_lengths.AppendTwoBits(static_cast<uint8_t>(prefix_length - kTwoBitPrefixEncodingOffset));
    } else if (prefix_length <= kMoreThanFourBitPrefixEncodingOffset)
    {
        /* We need four bits to encode the prefix */
        encoded_prefix_lengths.AppendFourBits(static_cast<uint8_t>(kFourBitPrefixEncodingStart + prefix_length - kFourBitPrefixEncodingOffset));
    } else
    {
        //cmc_debug_msg("recursive encode");
        /* In case the prefix is longer than a byte, we indicate that by a four bit sequence and encode the actual prefix length afterwards */
        encoded_prefix_lengths.AppendFourBits(kIndicatePrefixLongerThanByte);
        /* We apply this encoding recursive until, the length has been encoded */
        EncodeAndAppendPrefixLength(encoded_prefix_lengths, prefix_length - kMoreThanFourBitPrefixEncodingOffset);
    }
}

inline void
EncodeAndAppendPrefixLengthByteBoundary(bit_vector::BitVector& encoded_prefix_lengths)
{
    /* We append the four bit code to skip to the next byte */
    encoded_prefix_lengths.AppendFourBits(kSkipToNextByteInEncodingStream);
}

inline void
EncodeAndAppendSuffixLengthByteBoundary(bit_vector::BitVector& encoded_suffix_lengths)
{
    /* We append the four bit code to skip to the next byte */
    encoded_suffix_lengths.AppendFourBits(kSkipToNextByteInEncodingStream);
}

inline void
EncodeAndAppendPrefix(bit_vector::BitVector& encoded_prefixes, const std::vector<uint8_t>& prefix_bits, const int prefix_length)
{
    //cmc_debug_msg("In EncodeAndAppendPrefix: size of vec: ", prefix_bits.size());
    //cmc_debug_msg("Size of encoded_prefix befor appending: ", encoded_prefixes.size());
    encoded_prefixes.AppendBits(prefix_bits, prefix_length);
    //cmc_debug_msg("Size of encoded_prefix after appending: ", encoded_prefixes.size());
}

template <typename T>
void
PushBackValueToByteStream(std::vector<uint8_t>& byte_stream, const T& value)
{
    /* Serialize the value to a collection of bytes (in big endian order) */
    const std::array<uint8_t, sizeof(T)> serialized_value = SerializeValue(value, Endian::Big);

    /* Append the the bytes to the byte stream */
    std::copy_n(serialized_value.begin(), sizeof(T), std::back_insert_iterator(byte_stream));
}

/* A value is extracted from the byte stream starting at position \a pos to which the iterator points to.
 * The serialized value is stored in big endian order.  */
template <typename T, typename Iter>
auto GetValueFromByteStream(Iter pos)
    -> std::enable_if_t<std::is_fundamental_v<T>, T>
{
    /* Get the number of bytes for this data type */
    const size_t type_length = sizeof(T);

    /* Deserialize the given */
    const std::array<uint8_t, sizeof(T)> deserialized_value = DeserializeValue<T>(pos, Endian::Big);

    T value;

    /* Copy the bytes in the correct endianness to the 'value' */
    std::memcpy(static_cast<void*>(&value), static_cast<const void*>(deserialized_value.data()), type_length);

    return value;
}

/* Encode and serialize the data of a level */
template<typename T>
std::vector<uint8_t>
LevelwisePrefixData<T>::EncodeLevelData() const
{
    /* The amount of bits for the refinement and prefix indication should be equal, since it equals the number of elements in the mesh */
    cmc_assert(refinement_indicator.size() == prefix_indicator.size() && refinement_indicator.size_bytes() == prefix_indicator.size_bytes());

    /* Get an estimate for the memory allocation of the serialized data */
    const size_t memory_estimate = prefixes.size();

    /* Allocate a new BitMap for the prefix indicators based on the indications from the coarsening process */
    bit_map::BitMap prefix_indications(prefix_indicator);
    
    /* Allocate a BitVector for the prefix lengths */
    bit_vector::BitVector prefix_lengths;
    prefix_lengths.Reserve(memory_estimate / 2 + 1);

    /* Allocate a BitVector for the prefixes */
    bit_vector::BitVector prefix_encodings;
    prefix_encodings.Reserve(memory_estimate + 1);

    int prefix_idx = 0;

    /* We are iterating over the prefix indications and check whether there is an actual prefix left corresponding to this ID or not.
     * Afterwards, we are encoding the length and the prefix itself and store them in a serialized way. */
    for (auto pfi_iter = prefix_indicator.begin(); pfi_iter != prefix_indicator.end(); ++pfi_iter, ++prefix_idx)
    {
        //cmc_debug_msg("Prefix index is: ", prefix_idx);
        /* Check if a prefix is indicated, if so it will be encoded if it is not empty; if no prefix is indicated, we just continue */
        if (*pfi_iter)
        {
            if (prefixes[prefix_idx].IsEmpty())
            {
                //cmc_debug_msg("Prefix is empty and will be cleared.");
                /* If the prefix is empty, it has been extracted to a coarser level.
                 * Therefore, we need to unset the indicator bit at this position. */
                prefix_indications.ClearBit(prefix_idx);
            } else
            {
                //cmc_debug_msg("Prefix Bit is set: Prefix hat significant bits: ", prefixes[prefix_idx].GetCountOfSignificantBits(), "\nsize of bytes in big endian: ", prefixes[prefix_idx].GetSignificantBitsInBigEndianOrdering().size());
                /* If the prefix is not empty, it's length has to be encoded and stored. */
                const int prefix_length = prefixes[prefix_idx].GetCountOfSignificantBits();
                EncodeAndAppendPrefixLength(prefix_lengths, prefix_length);


            //uint32_t bitssss = suff_iter->template ReinterpretDataAs<uint32_t>();
            //cmc_debug_msg("Suffix_Length: ", suffix_length, "\t Bitset: ", std::bitset<32>(bitssss));

            //std::vector<uint8_t> s = prefixes[prefix_idx].GetSignificantBitsInBigEndianOrdering();
            //uint32_t valsss{0};
            //memcpy(&valsss, s.data(), s.size());
            //cmc_debug_msg("Prefix_Length: ", prefix_length, "\t Bitset: ", std::bitset<32>(valsss));

                /* Afterwards, the actual prefix itself has to encoded and stored as well */
                EncodeAndAppendPrefix(prefix_encodings, prefixes[prefix_idx].GetSignificantBitsInBigEndianOrdering(), prefix_length);
            }
        }
        
    }
    cmc_debug_msg("prefix_idx = ", prefix_idx, "     sollte gleich ANzahl Elemente sein");
    /* We want to add a flag in the prefix lengths, that here (at the end of the process local encodings) is a boundary where we skip to the next byte in the prefix lengths as well as in the prefix encodings */
    EncodeAndAppendPrefixLengthByteBoundary(prefix_lengths);
    cmc_debug_msg("In encode level data");
    /* Pack all the data to a single stream */
    const size_t num_bits_indicator = prefix_indications.size();
    cmc_debug_msg("num_bits_indicator: ", num_bits_indicator);
    const size_t num_bytes_indicator = prefix_indications.size_bytes();
    cmc_debug_msg("num_bytes_indicator: ", num_bytes_indicator);
    const size_t num_bytes_prefix_length_encodings = prefix_lengths.size();
    cmc_debug_msg("num_bytes_prefix_length_encodings: ", num_bytes_prefix_length_encodings);
    const size_t num_bytes_prefix_encodings = prefix_encodings.size();
    cmc_debug_msg("num_bytes_prefix_encodings: ", num_bytes_prefix_encodings);

    /* Overall memory estimate (including the just computed number of bytes for this whole level) */
    const size_t num_bytes_level = 4 * sizeof(size_t) + 2 * num_bytes_indicator + num_bytes_prefix_length_encodings + num_bytes_prefix_encodings;
    cmc_debug_msg("num_bytes_level: ", num_bytes_level);

    /* Now, we fill a buffer for the whole level */
    std::vector<uint8_t> serialized_level_data;
    serialized_level_data.reserve(num_bytes_level);

    //TODO: This has to be adapted for parallel gathering of the data 
    /* At first, we store the number of bytes for the whole level */
    PushBackValueToByteStream(serialized_level_data, num_bytes_level);

    /* Afterwards, we store the number of bits for each of the both indicators, the number of bytes for the prefix lengths
     * and then the number of bytes for the encoding of the actual prefixes */
    PushBackValueToByteStream(serialized_level_data, num_bits_indicator);
    PushBackValueToByteStream(serialized_level_data, num_bytes_prefix_length_encodings);
    PushBackValueToByteStream(serialized_level_data, num_bytes_prefix_encodings);

    /* Afterwards, we will copy the indicators as well as the encodings */
    std::copy_n(refinement_indicator.begin_bytes(), refinement_indicator.size_bytes(), std::back_insert_iterator(serialized_level_data));
    std::copy_n(prefix_indications.begin_bytes(), prefix_indications.size_bytes(), std::back_insert_iterator(serialized_level_data));
    std::copy_n(prefix_lengths.begin(), prefix_lengths.size(), std::back_insert_iterator(serialized_level_data));
    std::copy_n(prefix_encodings.begin(), prefix_encodings.size(), std::back_insert_iterator(serialized_level_data));

    cmc_debug_msg("serialized_level_data size: ", serialized_level_data.size());
    return serialized_level_data;
}


template <int N>
inline
std::vector<uint8_t>
EncodeSuffixes(const std::vector<CompressionValue<N>>& suffixes)
{
    /* A BitMap to collect suffix indication bits */
    bit_map::BitMap suffix_indications;
    suffix_indications.Reserve(suffixes.size());

    /* A BitVector to store the encoded suffix lengths */
    bit_vector::BitVector suffix_lengths;
    suffix_lengths.Reserve(suffixes.size() / 2 + 1);
    
    /* A BitVector to store the actual encoded suffix */
    bit_vector::BitVector suffix_encodings;
    suffix_encodings.Reserve(suffixes.size() + 1);

    /* Iterate over all suffixes and check whether they exist (and have not been fully extracted), 
     * if so, they will be encoded */
    int all_suffix_bits = 0;
    for (auto suff_iter = suffixes.begin(); suff_iter != suffixes.end(); ++suff_iter)
    {
        /* If the suffix is empty, we indicate this and continue with the next suffix */
        if (suff_iter->IsEmpty())
        {
            suffix_indications.AppendBit(false);
            continue;
        }

        /* In case the suffix is not empty, we will encode it's length as well as the suffix itself  */
        suffix_indications.AppendBit(true);

        /* If the prefix is not empty, it's length has to be encoded and stored. */
        const int suffix_length = suff_iter->GetCountOfSignificantBits();
        EncodeAndAppendPrefixLength(suffix_lengths, suffix_length);

        //cmc_debug_msg("Suffix TailBitPos: ", suff_iter->GetTrailBit());
        //if constexpr (N == 4)
        //{
        ////uint32_t bitssss = suff_iter->template ReinterpretDataAs<uint32_t>();
        ////cmc_debug_msg("Suffix_Length: ", suffix_length, "\t Bitset: ", std::bitset<32>(bitssss));
        ////
        //std::vector<uint8_t> s = suff_iter->GetSignificantBitsInBigEndianOrdering();
        //std::vector<uint8_t> sssss = s;
        //sssss[0] = s[3]; sssss[1] = s[2]; sssss[2] = s[1]; sssss[3] = s[0];
        ////
        //uint32_t valsss{0};
        //memcpy(&valsss, s.data(), s.size());
        //cmc_debug_msg("Bitset: ", std::bitset<32>(valsss), "    Suffix_Length: ", suffix_length);
        //}
        //cmc_debug_msg("Suffix length: ", suffix_length);
        all_suffix_bits += suffix_length;

        /* Afterwards, the actual prefix itself has to encoded and stored as well */
        EncodeAndAppendPrefix(suffix_encodings, suff_iter->GetSignificantBitsInBigEndianOrdering(), suffix_length);
    }

    cmc_debug_msg("\nAll leftover Suffix Bits: ", all_suffix_bits, "\n");

    /* We want to add a flag in the prefix lengths, that here (at the end of the process local encodings) is a boundary where we skip to the next byte in the prefix lengths as well as in the prefix encodings */
    EncodeAndAppendSuffixLengthByteBoundary(suffix_lengths);

    /* Pack all the data to a single stream */
    const size_t num_bits_indicator = suffix_indications.size();
    cmc_debug_msg("num_bits_indicator: ", num_bits_indicator);
    const size_t num_bytes_indicator = suffix_indications.size_bytes();
    cmc_debug_msg("num_bytes_indicator: ", num_bytes_indicator);
    const size_t num_bytes_suffix_length_encodings = suffix_lengths.size();
    cmc_debug_msg("num_bytes_suffix_length_encodings: ", num_bytes_suffix_length_encodings);
    const size_t num_bytes_suffix_encodings = suffix_encodings.size();
    cmc_debug_msg("num_bytes_suffix_encodings: ", num_bytes_suffix_encodings);

    /* Calculate the amount of bytes needed to encode the suffixes */
    const size_t num_bytes_suffixes = 4 * sizeof(size_t) + num_bytes_indicator + num_bytes_suffix_length_encodings + num_bytes_suffix_encodings;
    cmc_debug_msg("num_bytes_suffixes: ", num_bytes_suffixes);
    /* Now, we fill a buffer for the suffixes */
    std::vector<uint8_t> serialized_suffix_data;
    serialized_suffix_data.reserve(num_bytes_suffixes);

    //TODO: This has to be adapted for parallel gathering of the data 
    /* At first, we store the number of bytes for the suffix-level */
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffixes);

    /* Afterwards, we store the number of bits for the suffix indicators, the number of bytes for the suffix lengths
     * and then the number of bytes for the encoding of the actual suffixes */
    PushBackValueToByteStream(serialized_suffix_data, num_bits_indicator);
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffix_length_encodings);
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffix_encodings);

    /* Copy the serialized data to the buffer */
    std::copy_n(suffix_indications.begin_bytes(), num_bytes_indicator, std::back_insert_iterator(serialized_suffix_data));
    std::copy_n(suffix_lengths.begin(), num_bytes_suffix_length_encodings, std::back_insert_iterator(serialized_suffix_data));
    std::copy_n(suffix_encodings.begin(), num_bytes_suffix_encodings, std::back_insert_iterator(serialized_suffix_data));

    return serialized_suffix_data;
}



template <int N>
inline
std::vector<uint8_t>
EncodePlainSuffixes(const std::vector<CompressionValue<N>>& suffixes)
{
    /* A BitVector to store the actual encoded suffix */
    bit_vector::BitVector suffix_encodings;
    suffix_encodings.Reserve(3 * suffixes.size() + 1);

    /* Iterate over all suffixes and check whether they exist (and have not been fully extracted), 
     * if so, they will be encoded */
    int all_suffix_bits = 0;
    for (auto suff_iter = suffixes.begin(); suff_iter != suffixes.end(); ++suff_iter)
    {
        /* If the suffix is empty, we indicate this and continue with the next suffix */
        if (suff_iter->IsEmpty())
        {
            continue;
        }

        all_suffix_bits += suff_iter->GetCountOfSignificantBits();

        /* Afterwards, the actual prefix itself has to encoded and stored as well */
        EncodeAndAppendPrefix(suffix_encodings, suff_iter->GetSignificantBitsInBigEndianOrdering(), suff_iter->GetCountOfSignificantBits());
    }

    cmc_debug_msg("\nAll leftover Suffix Bits: ", all_suffix_bits, "\n");

    const size_t num_bytes_suffix_encodings = suffix_encodings.size();
    cmc_debug_msg("num_bytes_suffix_encodings: ", num_bytes_suffix_encodings);

    /* Calculate the amount of bytes needed to encode the suffixes */
    const size_t num_bytes_suffixes = sizeof(size_t) + num_bytes_suffix_encodings;
    cmc_debug_msg("num_bytes_suffixes: ", num_bytes_suffixes);

    /* Now, we fill a buffer for the suffixes */
    std::vector<uint8_t> serialized_suffix_data;
    serialized_suffix_data.reserve(num_bytes_suffixes);

    //TODO: This has to be adapted for parallel gathering of the data 
    /* At first, we store the number of bytes for the suffix-level */
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffixes);

    /* Copy the serialized data to the buffer */
    std::copy_n(suffix_encodings.begin(), num_bytes_suffix_encodings, std::back_insert_iterator(serialized_suffix_data));

    return serialized_suffix_data;
}

class PrefixDecoder
{
public:
    PrefixDecoder() = delete;

    PrefixDecoder(const std::vector<uint8_t>& serialized_data)
    : serialized_data_{serialized_data} {
        /* Offset for the first encoded values in the stream */
        const size_t offset = sizeof(size_t);

        /* Get the following amount of bytes for the current level */
        current_level_bytes_ = GetValueFromByteStream<size_t>(serialized_data_.data());
        cmc_debug_msg("current_level_bytes_: ", current_level_bytes_);
        /* Get the length of all arrays */
        const size_t num_bits_indicator = GetValueFromByteStream<size_t>(serialized_data_.data() + offset);
        cmc_debug_msg("num_bits_indicator: ", num_bits_indicator);
        const size_t num_bytes_prefix_length = GetValueFromByteStream<size_t>(serialized_data_.data() + 2 * offset);
        cmc_debug_msg("num_bytes_prefix_length: ", num_bytes_prefix_length);
        const size_t num_bytes_prefixes = GetValueFromByteStream<size_t>(serialized_data_.data() + 3 * offset);
        cmc_debug_msg("num_bytes_prefixes: ", num_bytes_prefixes);

        /* Get the number of bytes that are encoded for each indicator */
        const size_t num_bytes_indicator = num_bits_indicator / bit_map::kCharBit + (num_bits_indicator % bit_map::kCharBit != 0 ? 1 : 0);
        cmc_debug_msg("num_bytes_indicator: ", num_bytes_indicator);
        /* Create views on the first level */
        size_t start_offset = 4 * offset;
        refinement_indications_ = bit_map::BitMapView(serialized_data_.data() + start_offset, num_bytes_indicator);

        start_offset += num_bytes_indicator;
        prefix_indications_ =  bit_map::BitMapView(serialized_data_.data() + start_offset, num_bytes_indicator);

        start_offset += num_bytes_indicator;
        prefix_lengths_ = bit_vector::BitVectorView(serialized_data_.data() + start_offset, num_bytes_prefix_length);

        start_offset += num_bytes_prefix_length;
        prefixes_ = bit_vector::BitVectorView(serialized_data_.data() + start_offset, num_bytes_prefixes);
    };

    void MoveToNextLevel();
    void MoveToSuffixLevel();
    bool GetNextPrefixIndicatorBit();

    bool GetNextRefinementIndicatorBit();

    std::pair<std::vector<uint8_t>, int> GetNextPrefix();

    void MoveToGlobalPositionWithinLevel(const size_t global_start_position = 0);

    void MoveToPlainEncodedSuffixLevel();
    std::vector<uint8_t> GetNextPlainSuffix(const size_t num_bits);

private:
    int DecodeTwoBitCode(const uint8_t two_bit_encoding) const;
    int DecodeFourBitCode(const uint8_t four_bit_encoding) const;
    bool IsPrefixLengthEncodedInFourBitCode(const uint8_t four_bit_encoding) const;
    bool IsMoveToNextByteIndicator(const uint8_t four_bit_encoding) const;
    int GetNextPrefixLength();

    const std::vector<uint8_t>& serialized_data_;
    size_t current_level_bytes_{0};
    size_t processed_bytes_{0};

    bit_map::BitMapView refinement_indications_;
    bit_map::BitMapView prefix_indications_;

    bit_vector::BitVectorView prefix_lengths_;
    bit_vector::BitVectorView prefixes_;
};

inline
void
PrefixDecoder::MoveToGlobalPositionWithinLevel(const size_t global_start_position)
{
    //TODO: In parallel, we need to move to the first local element. Unfortunately, this has to be done serial from the start of the byte stream
}

inline
bool
PrefixDecoder::GetNextPrefixIndicatorBit()
{
    return prefix_indications_.GetNextBit();
}

inline
bool
PrefixDecoder::GetNextRefinementIndicatorBit()
{
    return refinement_indications_.GetNextBit();
}

inline
int
PrefixDecoder::DecodeTwoBitCode(const uint8_t two_bit_encoding) const
{
    return static_cast<int>((two_bit_encoding >> kShiftTwoBitPrefixLength) + kTwoBitPrefixEncodingOffset);
}

inline
int
PrefixDecoder::DecodeFourBitCode(const uint8_t four_bit_encoding) const
{
    return static_cast<int>((four_bit_encoding >> kShiftFourBitPrefixLength) - (kFourBitPrefixEncodingStart - kFourBitPrefixEncodingOffset));
}

inline
bool
PrefixDecoder::IsPrefixLengthEncodedInFourBitCode(const uint8_t four_bit_encoding) const
{
    return ((four_bit_encoding >> kShiftFourBitPrefixLength) <= kMaxOneBytePrefixEncoding);
}

inline
bool
PrefixDecoder::IsMoveToNextByteIndicator(const uint8_t four_bit_encoding) const
{
    return ((four_bit_encoding >> kShiftFourBitPrefixLength) == kSkipToNextByteInEncodingStream);
}

inline
void
PrefixDecoder::MoveToPlainEncodedSuffixLevel()
{
    const size_t offset = sizeof(size_t);

    /* Add the amount of bytes of the current level to the counter of all processed bytes */
    processed_bytes_ += current_level_bytes_;

    /* Get the following amount of bytes for the current level */
    current_level_bytes_ = GetValueFromByteStream<size_t>(serialized_data_.data() + processed_bytes_);
    cmc_debug_msg("Current lvl bbytes (in plain suff) is: ", current_level_bytes_);
    prefixes_ = bit_vector::BitVectorView(serialized_data_.data() + processed_bytes_ + offset, current_level_bytes_ - offset);
}

inline
void
PrefixDecoder::MoveToSuffixLevel()
{
    const size_t offset = sizeof(size_t);

    /* Add the amount of bytes of the current level to the counter of all processed bytes */
    processed_bytes_ += current_level_bytes_;

    /* Get the following amount of bytes for the current level */
    current_level_bytes_ = GetValueFromByteStream<size_t>(serialized_data_.data() + processed_bytes_);

    /* Get the length of all arrays */
    const size_t num_bits_indicator = GetValueFromByteStream<size_t>(serialized_data_.data() + processed_bytes_ + offset);
    cmc_debug_msg("num bits indicator: ", num_bits_indicator);
    const size_t num_bytes_prefix_length = GetValueFromByteStream<size_t>(serialized_data_.data() + processed_bytes_ + 2 * offset);
    cmc_debug_msg("num_bytes_prefix_length: ", num_bytes_prefix_length);
    const size_t num_bytes_prefixes = GetValueFromByteStream<size_t>(serialized_data_.data() + processed_bytes_ + 3 * offset);
    cmc_debug_msg("num_bytes_prefixes: ", num_bytes_prefixes);

    /* Get the number of bytes that are encoded for the indicator */
    const size_t num_bytes_indicator = num_bits_indicator / bit_map::kCharBit + (num_bits_indicator % bit_map::kCharBit != 0 ? 1 : 0);
    cmc_debug_msg("num_bytes_indicator: ", num_bytes_indicator);

    /* Create views on the first level */
    size_t start_offset = processed_bytes_ + 4 * offset;
    refinement_indications_ = bit_map::BitMapView(serialized_data_.data() + start_offset, 0);
    cmc_debug_msg("refinement_indications_ start offset: ", start_offset);

    prefix_indications_ = bit_map::BitMapView(serialized_data_.data() + start_offset, num_bytes_indicator);
    cmc_debug_msg("prefix_indications_ start offset: ", start_offset);

    start_offset += num_bytes_indicator;
    prefix_lengths_ = bit_vector::BitVectorView(serialized_data_.data() + start_offset, num_bytes_prefix_length);
    cmc_debug_msg("prefix_lengths_ start offset: ", start_offset);

    start_offset += num_bytes_prefix_length;
    prefixes_ = bit_vector::BitVectorView(serialized_data_.data() + start_offset, num_bytes_prefixes);
    cmc_debug_msg("prefixes_ start offset: ", start_offset);

}

inline
void
PrefixDecoder::MoveToNextLevel()
{
    const size_t offset = sizeof(size_t);

    /* Add the amount of bytes of the current level to the counter of all processed bytes */
    processed_bytes_ += current_level_bytes_;
    cmc_debug_msg("processed_bytes_: ", processed_bytes_);
    /* Get the following amount of bytes for the current level */
    current_level_bytes_ = GetValueFromByteStream<size_t>(serialized_data_.data() + processed_bytes_);

    /* Get the length of all arrays */
    const size_t num_bits_indicator = GetValueFromByteStream<size_t>(serialized_data_.data() + processed_bytes_ + offset);
    cmc_debug_msg("num bits indicator: ", num_bits_indicator);
    const size_t num_bytes_prefix_length = GetValueFromByteStream<size_t>(serialized_data_.data() + processed_bytes_ + 2 * offset);
    cmc_debug_msg("num_bytes_prefix_length: ", num_bytes_prefix_length);
    const size_t num_bytes_prefixes = GetValueFromByteStream<size_t>(serialized_data_.data() + processed_bytes_ + 3 * offset);
    cmc_debug_msg("num_bytes_prefixes: ", num_bytes_prefixes);

    /* Get the number of bytes that are encoded for each indicator */
    const size_t num_bytes_indicator = num_bits_indicator / bit_map::kCharBit + (num_bits_indicator % bit_map::kCharBit != 0 ? 1 : 0);
    cmc_debug_msg("num_bytes_indicator: ", num_bytes_indicator);

    /* Create views on the first level */
    size_t start_offset = processed_bytes_ + 4 * offset;
    refinement_indications_ = bit_map::BitMapView(serialized_data_.data() + start_offset, num_bytes_indicator);
    cmc_debug_msg("refinement_indications_ start offset: ", start_offset);

    start_offset += num_bytes_indicator;
    prefix_indications_ = bit_map::BitMapView(serialized_data_.data() + start_offset, num_bytes_indicator);
    cmc_debug_msg("prefix_indications_ start offset: ", start_offset);

    start_offset += num_bytes_indicator;
    prefix_lengths_ = bit_vector::BitVectorView(serialized_data_.data() + start_offset, num_bytes_prefix_length);
    cmc_debug_msg("prefix_lengths_ start offset: ", start_offset);

    start_offset += num_bytes_prefix_length;
    prefixes_ = bit_vector::BitVectorView(serialized_data_.data() + start_offset, num_bytes_prefixes);
    cmc_debug_msg("prefixes_ start offset: ", start_offset);
}

inline
int
PrefixDecoder::GetNextPrefixLength()
{
    /* Check if a two or four bit code is encoded */
    bool is_four_bit_code = prefix_lengths_.IsCurrentBitSet();

    int prefix_length{0};
    
    if (is_four_bit_code)
    {
        /* In this case, we encoded a four bit sequence */
        std::vector<uint8_t> pref_length = prefix_lengths_.GetNextBitSequence(kFourBitCodeSize);

        if (IsPrefixLengthEncodedInFourBitCode(pref_length.front()))
        {
            /* If a prefix length is encoded, we decode it */
            prefix_length += DecodeFourBitCode(pref_length.front());
        } else
        {
            /* There are two special codes which need to be handled seperately */
            /* Check if we need to move to next full byte */
            if (IsMoveToNextByteIndicator(pref_length.front()))
            {
                /* If so, we need to move wihtin the arrays of the prefix_length as well as in the prefix encoding 
                 * to the start of the next byte */
                prefix_lengths_.MoveToNextByte();
                prefixes_.MoveToNextByte();
            } else
            {
                /* If the code is not a "Move"-Indicator, it means that the prefix is longer than a single byte
                 * and we call the prefix length decoding recursively */
                prefix_length += kMoreThanFourBitPrefixEncodingOffset + GetNextPrefixLength();
            }
        }
    } else
    {
        /* If it is not a four bit code, we encoded a two bit sequence */
        std::vector<uint8_t> pref_length = prefix_lengths_.GetNextBitSequence(kTwoBitCodeSize);

        /* Decode the prefix length */
        prefix_length += DecodeTwoBitCode(pref_length.front());
    }

    return prefix_length;
}

inline
std::pair<std::vector<uint8_t>, int>
PrefixDecoder::GetNextPrefix()
{
    /* Get the length of the next prefix */
    const int prefix_length = GetNextPrefixLength();


    //uint32_t bitssss = suff_iter->template ReinterpretDataAs<uint32_t>();
    //cmc_debug_msg("Suffix_Length: ", suffix_length, "\t Bitset: ", std::bitset<32>(bitssss));

    //std::vector<uint8_t> s = prefixes_.GetNextBitSequence(static_cast<size_t>(prefix_length));
    //uint32_t valsss{0};
    //memcpy(&valsss, s.data(), s.size());
    //cmc_debug_msg("Decoded Prefix_Length: ", prefix_length, "\t Bitset: ", std::bitset<32>(valsss));

    /* Get the actual prefix of the decoded length */
    return std::make_pair(prefixes_.GetNextBitSequence(static_cast<size_t>(prefix_length)), prefix_length);
}

inline
std::vector<uint8_t>
PrefixDecoder::GetNextPlainSuffix(const size_t num_bits)
{
    return prefixes_.GetNextBitSequence(num_bits);
}

}

#endif /* !CMC_PREFIX_ENCODING_HXX */
