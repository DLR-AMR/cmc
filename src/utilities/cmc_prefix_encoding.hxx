#ifndef CMC_PREFIX_ENCODING_HXX
#define CMC_PREFIX_ENCODING_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_prefix.hxx"
#include "utilities/cmc_vector_view.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_vector_view.hxx"
#include "utilities/cmc_huffman.hxx"
#include "utilities/cmc_ac_model.hxx"
#include "utilities/cmc_arithmetic_encoder.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <vector>
#include <utility>
#include <type_traits>
#include <cstring>


//Delete later
#include <bitset>

namespace cmc
{

enum SuffixEncoding {Plain, LengthEncoding, ArithmeticLengthEncoding, EncodeFirstOne};

template <int N> using SuffixEncodingFunc = std::vector<uint8_t> (const std::vector<CompressionValue<N>>&);

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
    std::vector<uint8_t> EncodeLevelData(const huffman::HuffmanTree<int>& huffman_encoder) const;

    std::vector<uint8_t> EncodeLevelDataAsLeadingZeroCount() const;

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

    cmc_debug_msg("Reconstructed Value From Stream: ", value);
    
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

/* Encode and serialize the data of a level */
template<typename T>
std::vector<uint8_t>
LevelwisePrefixData<T>::EncodeLevelData(const huffman::HuffmanTree<int>& huffman_encoder) const
{
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
        /* Check if a prefix is indicated, if so it will be encoded if it is not empty; if no prefix is indicated, we just continue */
        if (*pfi_iter)
        {
            if (prefixes[prefix_idx].IsEmpty())
            {
                /* If the prefix is empty, it has been extracted to a coarser level.
                 * Therefore, we need to unset the indicator bit at this position. */
                prefix_indications.ClearBit(prefix_idx);
            } else
            {
                /* If the prefix is not empty, it's length has to be encoded and stored. */
                const uint8_t num_bits_prefix = static_cast<uint8_t>(prefixes[prefix_idx].GetCountOfSignificantBits());

                /* Encode the prefix length with the huffman tree */
                //cmc_debug_msg("Symbol to huffman encode: ", static_cast<int>(num_bits_prefix));
                const auto [code, code_length] = huffman_encoder.EncodeSymbol(num_bits_prefix);

                if (code_length > 0)
                {
                    /* The code length is only encoded when there is an actual code. A code of length zero occurs, when the whole block is 
                    * is coarsened to only a single value (the optimal code in that case is of length zero) */
                    prefix_lengths.AppendBits(code, static_cast<int>(code_length));
                }

                if (num_bits_prefix != 0)
                {
                    /* If there is a prefix to encode, we get the prefix and set it in the encoded stream */
                    prefix_encodings.AppendBits(prefixes[prefix_idx].GetSignificantBitsInBigEndianOrdering(), num_bits_prefix);
                }

                // const int prefix_length = prefixes[prefix_idx].GetCountOfSignificantBits();
                //EncodeAndAppendPrefixLength(prefix_lengths, prefix_length);
                /* Afterwards, the actual prefix itself has to encoded and stored as well */
                //EncodeAndAppendPrefix(prefix_encodings, prefixes[prefix_idx].GetSignificantBitsInBigEndianOrdering(), static_cast<int>(num_bits_prefix));
            }
        }
        
    }
    cmc_debug_msg("prefix_idx = ", prefix_idx, " sollte gleich Anzahl Elemente sein");

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
    const size_t num_bytes_level = 4 * sizeof(size_t) + num_bytes_indicator + num_bytes_prefix_length_encodings + num_bytes_prefix_encodings;
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
EncodeSuffixLengths(const std::vector<CompressionValue<N>>& suffixes)
{
    /* A BitVector to store the encoded suffix lengths */
    bit_vector::BitVector suffix_lengths;
    suffix_lengths.Reserve(suffixes.size() / 2 + 1);
    
    /* A BitVector to store the actual encoded suffix */
    bit_vector::BitVector suffix_encodings;
    suffix_encodings.Reserve(suffixes.size() + 1);

    std::vector<int> frequency(N * bit_vector::kCharBit + 1, int{0});

    /* Iterate over all suffixes on this level */
    for (auto val_iter = suffixes.begin(); val_iter != suffixes.end(); ++val_iter)
    {
        /* Get the count of significant bits, i.e. the length of the prefix */
        const int suff_length = static_cast<uint16_t>(val_iter->GetCountOfSignificantBits());

        /* Increment the counter for the given prefix length */
        ++frequency[suff_length];
    }

    /* Create a huffman tree out of it */
    huffman::HuffmanTree<int> huffman_encoder(frequency);
    cmc_debug_msg("Huffmann Coder has been instantiated");

    /* Iterate over all suffixes and check whether they exist (and have not been fully extracted), 
     * if so, they will be encoded */
    int all_suffix_bits = 0;

    for (auto suff_iter = suffixes.begin(); suff_iter != suffixes.end(); ++suff_iter)
    {
        /* If the prefix is not empty, it's length has to be encoded and stored. */
        const int suffix_length = suff_iter->GetCountOfSignificantBits();

        /* Encode the prefix length with the huffman tree */
        const auto [code, code_length] = huffman_encoder.EncodeSymbol(static_cast<uint8_t>(suffix_length));

        if (code_length > 0)
        {
            /* The code length is only encoded when there is an actual code. A code of length zero occurs, when the whole block is 
             * is coarsened to only a single value (the optimal code in that case is of length zero) */
            suffix_lengths.AppendBits(code, static_cast<int>(code_length));
        }

        if (suffix_length != 0)
        {
            /* If there is a prefix to encode, we get the prefix and set it in the encoded stream */
            const std::vector<uint8_t> suffix = suff_iter->GetSignificantBitsInBigEndianOrdering();
            suffix_encodings.AppendBits(suffix, suffix_length);
        }


        all_suffix_bits += suffix_length;
    }

    cmc_debug_msg("\nAll leftover Suffix Bits: ", all_suffix_bits, "\n");

    /* We want to add a flag in the prefix lengths, that here (at the end of the process local encodings) is a boundary where we skip to the next byte in the prefix lengths as well as in the prefix encodings */
    //EncodeAndAppendSuffixLengthByteBoundary(suffix_lengths);

    /* Pack all the data to a single stream */
    const size_t num_bytes_suffix_length_encodings = suffix_lengths.size();
    cmc_debug_msg("num_bytes_suffix_length_encodings: ", num_bytes_suffix_length_encodings);
    const size_t num_bytes_suffix_encodings = suffix_encodings.size();
    cmc_debug_msg("num_bytes_suffix_encodings: ", num_bytes_suffix_encodings);

    /* Calculate the amount of bytes needed to encode the suffixes */
    const size_t num_bytes_suffixes = 3 * sizeof(size_t) + sizeof(uint16_t) + num_bytes_suffix_length_encodings + num_bytes_suffix_encodings;
    cmc_debug_msg("num_bytes_suffixes: ", num_bytes_suffixes);
    /* Now, we fill a buffer for the suffixes */
    std::vector<uint8_t> serialized_suffix_data;
    serialized_suffix_data.reserve(num_bytes_suffixes);

    //TODO: This has to be adapted for parallel gathering of the data 
    /* At first, we store the number of bytes for the suffix-level */
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffixes);

    /* Afterwards, we store the number of bits for the suffix indicators, the number of bytes for the suffix lengths
     * and then the number of bytes for the encoding of the actual suffixes */
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffix_length_encodings);
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffix_encodings);

    /* Serialioze the huffman tree and store it */
    const size_t symbol_size = 6;
    bit_vector::BitVector serialized_huffman_tree = huffman_encoder.SerializeHuffmanTree(symbol_size);
    cmc_debug_msg("Huffmann tree has been serialized");
    
    const uint16_t huffman_tree_size = serialized_huffman_tree.size();
    PushBackValueToByteStream(serialized_suffix_data, huffman_tree_size);
    std::copy_n(serialized_huffman_tree.begin(), huffman_tree_size, std::back_insert_iterator(serialized_suffix_data));

    /* Copy the serialized data to the buffer */
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

    int i = 0;
    /* Iterate over all suffixes and check whether they exist (and have not been fully extracted), 
     * if so, they will be encoded */
    int all_suffix_bits = 0;
    std::vector<int> null_bits(33,0);
    int null_bit_counter = 0;
    for (auto suff_iter = suffixes.begin(); suff_iter != suffixes.end(); ++suff_iter)
    {
        /* If the suffix is empty, we indicate this and continue with the next suffix */
        if (suff_iter->IsEmpty())
        {
            continue;
        }

        CompressionValue<N> val = *suff_iter;
        val.NullifyNonSignificantFrontBits();
        int num_null_bits = val.GetNumberLeadingZeros() - val.GetFrontBit();
   
        null_bit_counter += num_null_bits;
        ++null_bits[num_null_bits];



        #if 0
        if constexpr (N == 4)
        {
        if (i < 1000)
        {
            cmc_debug_msg("\n");
            CompressionValue<N> vall = *suff_iter;
            cmc_debug_msg("Count significant bits: ", suff_iter->GetCountOfSignificantBits());
            cmc_debug_msg("Start Bit: ", vall.GetFrontBit(), ", Tail Bit: ", vall.GetTrailBit());
            cmc_debug_msg("Count leading zeros: ", suff_iter->GetNumberLeadingZeros());
            cmc_debug_msg("Count trailing zeros: ", suff_iter->GetNumberTrailingZeros());
            CompressionValue<N> prev_val = *std::prev(suff_iter);

            uint32_t val = vall.template ReinterpretDataAs<uint32_t>();
            cmc_debug_msg(std::bitset<32>(val));
            vall.NullifyNonSignificantFrontBits();
            prev_val.NullifyNonSignificantFrontBits();
            val = vall.template ReinterpretDataAs<uint32_t>();

            cmc_debug_msg(std::bitset<32>(val));
            uint32_t diff = val - prev_val.template ReinterpretDataAs<uint32_t>();
            vall ^= prev_val;
            uint32_t xor_val = vall.template ReinterpretDataAs<uint32_t>();
            cmc_debug_msg(std::bitset<32>(xor_val), " = xor", " mit significant leading zeros: ", vall.GetLeadingZeroCountInSignificantBits());
            cmc_debug_msg(std::bitset<32>(diff), " = diff");
            uint32_t adjusted_diff = (val >= prev_val.template ReinterpretDataAs<uint32_t>() ? val - prev_val.template ReinterpretDataAs<uint32_t>() : prev_val.template ReinterpretDataAs<uint32_t>() -val);
            cmc_debug_msg(std::bitset<32>(adjusted_diff), " = diff_adjusted (needs explicit storage of sign)");
            ++i;
        }
        }
        #endif
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

    for (auto iter = null_bits.begin(); iter != null_bits.end(); ++iter)
    {
        cmc_debug_msg("Null bits of length: ", std::distance(null_bits.begin(), iter), " frequency: ", *iter);
    } 
    cmc_debug_msg("Num null bits not extracted: ", null_bit_counter);
    return serialized_suffix_data;
}


template <int N>
inline
std::vector<uint8_t>
EncodeFirstOneInSuffixes(const std::vector<CompressionValue<N>>& suffixes)
{
    /* A BitVector to store the actual encoded suffix */
    bit_vector::BitVector suffix_encodings;
    suffix_encodings.Reserve(3 * suffixes.size() + 1);

    int i = 0;
    /* Iterate over all suffixes and check whether they exist (and have not been fully extracted), 
     * if so, they will be encoded */
    int all_suffix_bits = 0;

    std::vector<arithmetic_encoding::Letter> alphabet;
    const uint32_t alphabet_size = static_cast<uint32_t>(N * bit_vector::kCharBit + 1);
    alphabet.reserve(alphabet_size);

    /* Construct the alphabet */
    for (uint32_t letter = 0; letter < alphabet_size; ++letter)
    {
        alphabet.emplace_back(arithmetic_encoding::Letter{letter, 1});
    }

    /* Iterate over all suffixes on this level */
    for (auto val_iter = suffixes.begin(); val_iter != suffixes.end(); ++val_iter)
    {
        /* Get the leading zero count in the suffix */
        const int lzc_in_suffix = val_iter->GetLeadingZeroCountInSignificantBits();

        /* Increment the counter for the given prefix length */
        ++(alphabet[lzc_in_suffix].frequency);
    }

    /* Construct the frequency model from the length information */
    arithmetic_encoding::StaticFrequencyModel frequency_model(alphabet);

    /* Construct the arithmetic encoder */
    arithmetic_encoding::Encoder arm_encoder(std::make_unique<cmc::arithmetic_encoding::StaticFrequencyModel>(frequency_model));

    for (auto suff_iter = suffixes.begin(); suff_iter != suffixes.end(); ++suff_iter)
    {
        CompressionValue<N> suffix = *suff_iter;

        /* Get the leading zero count in the suffix */
        const int lzc_in_suffix = suffix.GetLeadingZeroCountInSignificantBits();

        /* The leading zero count is processed by the arithmetic encoder */
        arm_encoder.EncodeSymbol(lzc_in_suffix);

        /* The remaining suffix is encoded */
        if (lzc_in_suffix == 0 && not suffix.IsEmpty())
        {
            /* Afterwards, the actual prefix itself has to encoded and stored as well */
            EncodeAndAppendPrefix(suffix_encodings, suffix.GetSignificantBitsInBigEndianOrdering(), suffix.GetCountOfSignificantBits());
        } else
        {
            /* When its greater than zero, we can trim the zeros */
            suffix.SetFrontBit(suffix.GetFrontBit() + lzc_in_suffix);

            if (not suffix.IsEmpty())
            {
                EncodeAndAppendPrefix(suffix_encodings, suffix.GetSignificantBitsInBigEndianOrdering(), suffix.GetCountOfSignificantBits());
            }
        }

        all_suffix_bits += suffix.GetCountOfSignificantBits();
    }

    bit_map::BitMap encoded_lzc_stream = arm_encoder.GetEncodedBitStream();

    const size_t num_bytes_suffix_lzc = encoded_lzc_stream.size_bytes();
    cmc_debug_msg("num_bytes_suffix_lzc: ", encoded_lzc_stream.size_bytes());

    cmc_debug_msg("\nAll leftover Suffix Bits: ", all_suffix_bits, "\n");

    const size_t num_bytes_suffix_encodings = suffix_encodings.size();
    cmc_debug_msg("num_bytes_suffix_encodings: ", num_bytes_suffix_encodings);

    /* Calculate the amount of bytes needed to encode the suffixes */
    const size_t num_bytes_suffixes = sizeof(size_t) + num_bytes_suffix_encodings + num_bytes_suffix_lzc;
    cmc_debug_msg("num_bytes_suffixes: ", num_bytes_suffixes);

    /* Now, we fill a buffer for the suffixes */
    std::vector<uint8_t> serialized_suffix_data;
    serialized_suffix_data.reserve(num_bytes_suffixes);

    //TODO: This has to be adapted for parallel gathering of the data 
    /* At first, we store the number of bytes for the suffix-level */
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffixes);

    /* Copy the serialized data to the buffer */
    std::copy_n(encoded_lzc_stream.begin_bytes(), num_bytes_suffix_lzc, std::back_insert_iterator(serialized_suffix_data));
    std::copy_n(suffix_encodings.begin(), num_bytes_suffix_encodings, std::back_insert_iterator(serialized_suffix_data));

    //for (auto iter = null_bits.begin(); iter != null_bits.end(); ++iter)
    //{
    //    cmc_debug_msg("Null bits of length: ", std::distance(null_bits.begin(), iter), " frequency: ", *iter);
    //} 
    //cmc_debug_msg("Num null bits not extracted: ", null_bit_counter);
    return serialized_suffix_data;
}


template <int N>
inline
std::vector<uint8_t>
EncodeSuffixLengthsArithmeticEncoder(const std::vector<CompressionValue<N>>& suffixes)
{
    /* A BitVector to store the encoded suffix lengths */
    //bit_vector::BitVector suffix_lengths;
    //suffix_lengths.Reserve(suffixes.size() / 2 + 1);
    
    /* A BitVector to store the actual encoded suffix */
    bit_vector::BitVector suffix_encodings;
    suffix_encodings.Reserve(suffixes.size() + 1);

    std::vector<arithmetic_encoding::Letter> alphabet;
    const uint32_t alphabet_size = static_cast<uint32_t>(N * bit_vector::kCharBit + 1);
    alphabet.reserve(alphabet_size);

    /* Construct the alphabet */
    for (uint32_t letter = 0; letter < alphabet_size; ++letter)
    {
        alphabet.emplace_back(arithmetic_encoding::Letter{letter, 1});
    }

    /* Iterate over all suffixes on this level */
    for (auto val_iter = suffixes.begin(); val_iter != suffixes.end(); ++val_iter)
    {
        /* Get the count of significant bits, i.e. the length of the prefix */
        const int suff_length = val_iter->GetCountOfSignificantBits();

        /* Increment the counter for the given prefix length */
        ++(alphabet[suff_length].frequency);
    }

    /* Construct the frequency model from the length information */
    arithmetic_encoding::StaticFrequencyModel frequency_model(alphabet);

    /* Construct the arithmetic encoder */
    arithmetic_encoding::Encoder arm_encoder(std::make_unique<cmc::arithmetic_encoding::StaticFrequencyModel>(frequency_model));

    /* Iterate over all suffixes and check whether they exist (and have not been fully extracted), 
     * if so, they will be encoded */
    int all_suffix_bits = 0;

    for (auto suff_iter = suffixes.begin(); suff_iter != suffixes.end(); ++suff_iter)
    {
        /* If the prefix is not empty, it's length has to be encoded and stored. */
        uint32_t suffix_length = static_cast<uint32_t>(suff_iter->GetCountOfSignificantBits());
        
        /* Encode the suffix length */
        arm_encoder.EncodeSymbol(suffix_length);

        if (suffix_length > 0)
        {
           /* If the suffix is not empty, we can leave the out the last bit of the suffix, because it will be a "one".
            * Since all zeros at the tail have been truncated by "PerformTailTruncation". */
           suffix_length -= 1;
        }

        if (suffix_length != 0)
        {
            /* If there is a prefix to encode, we get the prefix and set it in the encoded stream */
            const std::vector<uint8_t> suffix = suff_iter->GetSignificantBitsInBigEndianOrdering();
            suffix_encodings.AppendBits(suffix, suffix_length);
        }


        all_suffix_bits += suffix_length;
    }

    cmc_debug_msg("\nAll leftover Suffix Bits: ", all_suffix_bits, "\n");

    /* We want to add a flag in the prefix lengths, that here (at the end of the process local encodings) is a boundary where we skip to the next byte in the prefix lengths as well as in the prefix encodings */
    //EncodeAndAppendSuffixLengthByteBoundary(suffix_lengths);

    /* Get the encoded bit stream from the arithmetic encoder */
    bit_map::BitMap suffix_lengths = arm_encoder.GetEncodedBitStream();

    /* Pack all the data to a single stream */
    const size_t num_bytes_suffix_length_encodings = suffix_lengths.size_bytes();
    cmc_debug_msg("num_bytes_suffix_length_encodings: ", num_bytes_suffix_length_encodings);
    const size_t num_bytes_suffix_encodings = suffix_encodings.size();
    cmc_debug_msg("num_bytes_suffix_encodings: ", num_bytes_suffix_encodings);

    /* Calculate the amount of bytes needed to encode the suffixes */
    const size_t num_bytes_suffixes = 3 * sizeof(size_t) + sizeof(uint16_t) + num_bytes_suffix_length_encodings + num_bytes_suffix_encodings;
    cmc_debug_msg("num_bytes_suffixes: ", num_bytes_suffixes);
    /* Now, we fill a buffer for the suffixes */
    std::vector<uint8_t> serialized_suffix_data;
    serialized_suffix_data.reserve(num_bytes_suffixes);

    //TODO: This has to be adapted for parallel gathering of the data 
    /* At first, we store the number of bytes for the suffix-level */
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffixes);

    /* Afterwards, we store the number of bits for the suffix indicators, the number of bytes for the suffix lengths
     * and then the number of bytes for the encoding of the actual suffixes */
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffix_length_encodings);
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffix_encodings);

    /* Serialioze the huffman tree and store it */
    //const size_t symbol_size = 6;
    //bit_vector::BitVector serialized_huffman_tree = huffman_encoder.SerializeHuffmanTree(symbol_size);
    //cmc_debug_msg("Huffmann tree has been serialized");
    cmc_debug_msg("Encode alphabet of arithmetic coder and append (however this is only ~100 bytes)");

    //const uint16_t huffman_tree_size = serialized_huffman_tree.size();
    //PushBackValueToByteStream(serialized_suffix_data, huffman_tree_size);
    //std::copy_n(serialized_huffman_tree.begin(), huffman_tree_size, std::back_insert_iterator(serialized_suffix_data));

    /* Copy the serialized data to the buffer */
    std::copy_n(suffix_lengths.begin_bytes(), num_bytes_suffix_length_encodings, std::back_insert_iterator(serialized_suffix_data));
    std::copy_n(suffix_encodings.begin(), num_bytes_suffix_encodings, std::back_insert_iterator(serialized_suffix_data));

    return serialized_suffix_data;
}

template <int N>
inline
std::vector<int>
GetLeadingZeroFrequency(const std::vector<CompressionValue<N>>& values)
{
    std::vector<int> frequency(N * bit_vector::kCharBit + 1, int{0});

    /* Iterate over all prefixes on this level */
    for (auto val_iter = values.begin(); val_iter != values.end(); ++val_iter)
    {
        /* Get the count of significant bits, i.e. the length of the prefix */
        const int leading_zero_count = val_iter->GetLeadingZeroCountInSignificantBits();

        /* Increment the counter for the given prefix length */
        ++frequency[leading_zero_count];
    }

    int i = 0;
    for (auto fiter = frequency.begin(); fiter != frequency.end(); ++fiter, ++i)
    {
        cmc_debug_msg("Frequency of LeadingZero Count: ", i, " is ", *fiter);
    }

    return frequency;
}

#if 1

template <int N>
inline
std::vector<uint8_t>
EncodeSuffixesXORScheme(const std::vector<CompressionValue<N>>& suffixes)
{
    /* A BitVector to store the actual encoded suffix */
    bit_vector::BitVector suffix_encodings;
    suffix_encodings.Reserve(3 * suffixes.size() + 1);

    //int i = 5000000;
    std::vector<CompressionValue<N>> xor_suffixes;
    xor_suffixes.reserve(suffixes.size());

    /* Just copy the first value, since it cannot be xor'ed */
    xor_suffixes.push_back(suffixes.front());
    xor_suffixes.back().NullifyNonSignificantFrontBits();

    for (auto suff_iter = std::next(suffixes.begin()); suff_iter != suffixes.end(); ++suff_iter)
    {
        const CompressionValue<N>& previous_value = *std::prev(suff_iter);
        xor_suffixes.push_back(*suff_iter);
        xor_suffixes.back() ^= previous_value;
        xor_suffixes.back().NullifyNonSignificantFrontBits();
    }

    /* Get the frequency of the leading zero counts */
    std::vector<int> leading_zero_frequency = GetLeadingZeroFrequency<N>(xor_suffixes);

    /* Create a huffman tree out of it */
    huffman::HuffmanTree<int> huffman_encoder(leading_zero_frequency);
    cmc_debug_msg("Huffmann Coder has been instantiated");
    /* Serialioze the huffman tree and store it */
    const size_t symbol_size = 6;
    bit_vector::BitVector serialized_huffman_tree = huffman_encoder.SerializeHuffmanTree(symbol_size);
    cmc_debug_msg("Huffmann tree has been serialized");
    const uint16_t num_bits_huffman_tree = serialized_huffman_tree.size_bits();

    cmc_debug_msg("Size of leading zero count huffman tree in suffixes: ", num_bits_huffman_tree);

    //TODO: We just push to placeholders bacl 
    suffix_encodings.AppendBits(uint8_t{0}, 8);
    suffix_encodings.AppendBits(uint8_t{0}, 8);

    /* Store the serialized huffman tree */
    suffix_encodings.AppendBits(serialized_huffman_tree);
    cmc_debug_msg("Suffix Encoding begins");
    /* Iterate over all suffixes and encode them */
    for (auto suff_iter = xor_suffixes.begin(); suff_iter != xor_suffixes.end(); ++suff_iter)
    {
        /* Get the leading zero count in the significant part of the value */
        const int leading_zero_count = suff_iter->GetLeadingZeroCountInSignificantBits();
        //cmc_debug_msg("Old front bit: ", suff_iter->GetFrontBit(), "LZC: ", leading_zero_count, ", new front bit: ", suff_iter->GetFrontBit() + leading_zero_count, ", tail_bit: ", suff_iter->GetTrailBit());
        /* Only indicate the non zero part as significant */
        suff_iter->SetFrontBit(suff_iter->GetFrontBit() + leading_zero_count);
        /* Encode the leading zero count and append it to the suffix encodings */
        const auto [code, code_length] = huffman_encoder.EncodeSymbol(leading_zero_count);
        /* Set the leading zero count length */
        if (code_length > 0)
        {
            /* The code length is only encoded when there is an actual code. A code of length zero occurs, when the whole data
            * has only a sngle code (the optimal code in that case is of length zero) */
            suffix_encodings.AppendBits(code, static_cast<int>(code_length));
        }
        const int significant_bits = suff_iter->GetCountOfSignificantBits();
        if (significant_bits > 0)
        {
            /* Encode the remaining suffix */
            EncodeAndAppendPrefix(suffix_encodings, suff_iter->GetSignificantBitsInBigEndianOrdering(), significant_bits); 
        }
    }

    /* Calculate the amount of bytes needed to encode the suffixes */
    const size_t num_bytes_suffixes = sizeof(size_t) + suffix_encodings.size();
    cmc_debug_msg("num_bytes_suffixes: ", num_bytes_suffixes);

    /* Now, we fill a buffer for the suffixes */
    std::vector<uint8_t> serialized_suffix_data;
    serialized_suffix_data.reserve(num_bytes_suffixes);

    //TODO: This has to be adapted for parallel gathering of the data 
    /* At first, we store the number of bytes for the suffix-level */
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffixes);

    /* Copy the serialized data to the buffer */
    std::copy_n(suffix_encodings.begin(), suffix_encodings.size(), std::back_insert_iterator(serialized_suffix_data));

    return serialized_suffix_data;
}

#endif


template <int N>
inline
std::vector<uint8_t>
EncodeSuffixesXORSchemeNew(const std::vector<CompressionValue<N>>& suffixes)
{
    /* A BitVector to store the actual encoded suffix */
    bit_vector::BitVector suffix_encodings;
    suffix_encodings.Reserve(3 * suffixes.size() + 1);

    //int i = 5000000;
    std::vector<CompressionValue<N>> xor_suffixes;
    xor_suffixes.reserve(suffixes.size());

    /* Just copy the first value, since it cannot be xor'ed */
    xor_suffixes.push_back(suffixes.front());
    xor_suffixes.back().NullifyNonSignificantFrontBits();

    for (auto suff_iter = std::next(suffixes.begin()); suff_iter != suffixes.end(); ++suff_iter)
    {
        const CompressionValue<N>& previous_value = *std::prev(suff_iter);
        xor_suffixes.push_back(*suff_iter);
        xor_suffixes.back() ^= previous_value;
        xor_suffixes.back().NullifyNonSignificantFrontBits();
    }

    /* Get the frequency of the leading zero counts */
    std::vector<int> leading_zero_frequency = GetLeadingZeroFrequency<N>(xor_suffixes);

    /* We nullify the first encodings */
    leading_zero_frequency[0] = 0;
    leading_zero_frequency[1] = 0;
    //leading_zero_frequency[2] = 0;

    /* Create a huffman tree out of it */
    huffman::HuffmanTree<int> huffman_encoder(leading_zero_frequency);
    cmc_debug_msg("Huffmann Coder has been instantiated");
    /* Serialioze the huffman tree and store it */
    const size_t symbol_size = 6;
    bit_vector::BitVector serialized_huffman_tree = huffman_encoder.SerializeHuffmanTree(symbol_size);
    cmc_debug_msg("Huffmann tree has been serialized");
    const uint16_t num_bits_huffman_tree = serialized_huffman_tree.size_bits();

    cmc_debug_msg("Size of leading zero count huffman tree in suffixes: ", num_bits_huffman_tree);

    //TODO: We just push to placeholders bacl 
    suffix_encodings.AppendBits(uint8_t{0}, 8);
    suffix_encodings.AppendBits(uint8_t{0}, 8);

    /* Store the serialized huffman tree */
    suffix_encodings.AppendBits(serialized_huffman_tree);
    cmc_debug_msg("Suffix Encoding begins");
    /* Iterate over all suffixes and encode them */
    for (auto suff_iter = xor_suffixes.begin(); suff_iter != xor_suffixes.end(); ++suff_iter)
    {
        /* Get the leading zero count in the significant part of the value */
        const int leading_zero_count = suff_iter->GetLeadingZeroCountInSignificantBits();
        if (leading_zero_count >= 2)
        {
            /* Push back a flag that we encode a LZC */
            suffix_encodings.AppendBit(true);
            /* Only indicate the non zero part as significant */
            suff_iter->SetFrontBit(suff_iter->GetFrontBit() + leading_zero_count);
            /* Encode the leading zero count and append it to the suffix encodings */
            const auto [code, code_length] = huffman_encoder.EncodeSymbol(leading_zero_count);
            /* Set the leading zero count length */
            if (code_length > 0)
            {
                /* The code length is only encoded when there is an actual code. A code of length zero occurs, when the whole data
                * has only a sngle code (the optimal code in that case is of length zero) */
                suffix_encodings.AppendBits(code, static_cast<int>(code_length));
            }
        } else
        {
            /* Push back a flag that we do not encode a LZC */
            suffix_encodings.AppendBit(false);
        }
    
        const int significant_bits = suff_iter->GetCountOfSignificantBits();
        if (significant_bits > 0)
        {
            /* Encode the remaining suffix */
            EncodeAndAppendPrefix(suffix_encodings, suff_iter->GetSignificantBitsInBigEndianOrdering(), significant_bits); 
        }
    }

    /* Calculate the amount of bytes needed to encode the suffixes */
    const size_t num_bytes_suffixes = sizeof(size_t) + suffix_encodings.size();
    cmc_debug_msg("num_bytes_suffixes: ", num_bytes_suffixes);

    /* Now, we fill a buffer for the suffixes */
    std::vector<uint8_t> serialized_suffix_data;
    serialized_suffix_data.reserve(num_bytes_suffixes);

    //TODO: This has to be adapted for parallel gathering of the data 
    /* At first, we store the number of bytes for the suffix-level */
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffixes);

    /* Copy the serialized data to the buffer */
    std::copy_n(suffix_encodings.begin(), suffix_encodings.size(), std::back_insert_iterator(serialized_suffix_data));

    return serialized_suffix_data;
}


template <int N>
inline
std::vector<uint8_t>
EncodeSuffixesXORSchemeNew2(const std::vector<CompressionValue<N>>& suffixes)
{
    /* A BitVector to store the actual encoded suffix */
    bit_vector::BitVector suffix_encodings;
    suffix_encodings.Reserve(3 * suffixes.size() + 1);

    //int i = 5000000;
    std::vector<CompressionValue<N>> xor_suffixes;
    xor_suffixes.reserve(suffixes.size());

    /* Just copy the first value, since it cannot be xor'ed */
    xor_suffixes.push_back(suffixes.front());
    xor_suffixes.back().NullifyNonSignificantFrontBits();

    for (auto suff_iter = std::next(suffixes.begin()); suff_iter != suffixes.end(); ++suff_iter)
    {
        const CompressionValue<N>& previous_value = *std::prev(suff_iter);
        xor_suffixes.push_back(*suff_iter);
        xor_suffixes.back() ^= previous_value;
        xor_suffixes.back().NullifyNonSignificantFrontBits();
    }

    if constexpr (N == 4)
    {
        for (int i =12000000; i < 12001000; ++i)
        {
            CompressionValue<N> val = xor_suffixes[i];
            val.NullifyNonSignificantFrontBits();
            uint32_t int_val = val.template ReinterpretDataAs<uint32_t>();

            cmc_debug_msg("index ", i, " hat xor_suff int: ", std::bitset<32>(int_val), " = ", int_val, " mit num significant bits: ", val.GetCountOfSignificantBits());
        }
    }
    #if 0
    if constexpr (N == 4)
    {
        std::vector<uint32_t> int_suffixes;
        int_suffixes.reserve(xor_suffixes.size());

        for (int i =12000000; i < 12001000; ++i)
        {
            CompressionValue<N> val = suffixes[i];
            val.NullifyNonSignificantFrontBits();
            uint32_t int_val = val.template ReinterpretDataAs<uint32_t>();

            cmc_debug_msg("index ", i, " hat suff int: ", std::bitset<32>(int_val), " = ", int_val, " mit num significant bits: ", val.GetCountOfSignificantBits());
        }
    }
    #endif

    bit_vector::BitVector lzc_lengths;
    lzc_lengths.Reserve(xor_suffixes.size() / 2 + 1);

    cmc_debug_msg("Suffix Encoding begins");
    /* Iterate over all suffixes and encode them */
    for (auto suff_iter = xor_suffixes.begin(); suff_iter != xor_suffixes.end(); ++suff_iter)
    {
        /* Get the leading zero count in the significant part of the value */
        uint8_t leading_zero_count = static_cast<uint8_t>(suff_iter->GetLeadingZeroCountInSignificantBits());
        if (leading_zero_count > 15)
        {
            //cmc_debug_msg("LZC of over 15: length: ", leading_zero_count);
            leading_zero_count = 15;
        }
        lzc_lengths.AppendFourBits(leading_zero_count);
        lzc_lengths.AppendBits(leading_zero_count, 8);
        suff_iter->SetFrontBit(suff_iter->GetFrontBit() + leading_zero_count);
        
        const int significant_bits = suff_iter->GetCountOfSignificantBits();
        if (significant_bits > 0)
        {
            /* Encode the remaining suffix */
            EncodeAndAppendPrefix(suffix_encodings, suff_iter->GetSignificantBitsInBigEndianOrdering(), significant_bits); 
        }
    }
    cmc_debug_msg("Num suffixes: ", xor_suffixes.size(), " und size lzc vector: ", lzc_lengths.size());
    FILE* file = fopen("lzc_suffices.bin", "wb");
    fwrite(lzc_lengths.data(), sizeof(uint8_t), lzc_lengths.size(), file);
    fclose(file);

    /* Calculate the amount of bytes needed to encode the suffixes */
    const size_t num_bytes_suffixes = sizeof(size_t) + suffix_encodings.size();
    cmc_debug_msg("num_bytes_suffixes: ", num_bytes_suffixes);
    cmc_debug_msg("Num bytes to encode lzc: ", lzc_lengths.size());

    /* Now, we fill a buffer for the suffixes */
    std::vector<uint8_t> serialized_suffix_data;
    serialized_suffix_data.reserve(num_bytes_suffixes);

    //TODO: This has to be adapted for parallel gathering of the data 
    /* At first, we store the number of bytes for the suffix-level */
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffixes);

    /* Copy the serialized data to the buffer */
    std::copy_n(suffix_encodings.begin(), suffix_encodings.size(), std::back_insert_iterator(serialized_suffix_data));

    return serialized_suffix_data;
}


template <int N>
inline
std::vector<uint8_t>
EncodeSuffixesBitPlaneRLE(const std::vector<CompressionValue<N>>& suffixes)
{
    /* A BitVector to store the actual encoded suffix */
    bit_vector::BitVector suffix_encodings;
    suffix_encodings.Reserve(3 * suffixes.size() + 1);

    //int i = 5000000;
    std::vector<CompressionValue<N>> xor_suffixes;
    xor_suffixes.reserve(suffixes.size());

    /* Just copy the first value, since it cannot be xor'ed */
    xor_suffixes.push_back(suffixes.front());
    xor_suffixes.back().NullifyNonSignificantFrontBits();

    for (auto suff_iter = std::next(suffixes.begin()); suff_iter != suffixes.end(); ++suff_iter)
    {
        const CompressionValue<N>& previous_value = *std::prev(suff_iter);
        xor_suffixes.push_back(*suff_iter);
        xor_suffixes.back() ^= previous_value;
        xor_suffixes.back().NullifyNonSignificantFrontBits();
    }

    const size_t stride = 8;
    const size_t num_full_iterations = xor_suffixes.size() / stride;
    const size_t end_iteration_idx = num_full_iterations * stride;
    const size_t remaining_vals = xor_suffixes.size() % stride;

    for (auto iter = 0; iter < end_iteration_idx;)
    {
        /* Check if an encoding is possible */
        bool encode_bit_plane_rle_possible = true;
        const int num_significant_bits = xor_suffixes[iter].GetCountOfSignificantBits();
        for (auto idx = iter + 1; idx < iter + stride; ++idx)
        {
            if (num_significant_bits != xor_suffixes[idx].GetCountOfSignificantBits())
            {
                encode_bit_plane_rle_possible = false;
                break;
            }
        }

        /* Check if it possible to encode */
        if (encode_bit_plane_rle_possible)
        {
            suffix_encodings.AppendBit(true);

            if (num_significant_bits > 0)
            {
                /* Build bit plane */
                std::vector<std::vector<uint8_t>> significant_bits;
                significant_bits.reserve(stride);

                for (auto idx = 0; idx < stride; ++idx)
                {
                    significant_bits.emplace_back(xor_suffixes[iter + idx].GetSignificantBitsInBigEndianOrdering());
                }

                bit_vector::BitVector bit_plane;
                bit_plane.Reserve((num_significant_bits * stride) / 8 + 1);

                /* Iterate over all significant bit positions */
                for (auto bit_pos = 0; bit_pos < num_significant_bits; ++bit_pos)
                {
                    const int byte_id = bit_pos / bit_vector::kCharBit;
                    const int bit_id = bit_vector::kBitIndexStart - (bit_pos % bit_vector::kCharBit);

                    /* Iterate over all values in the stride */
                    for (auto val_iter = significant_bits.begin(); val_iter != significant_bits.end(); ++val_iter)
                    {
                        /* Append the current bit from this value */
                        bit_plane.AppendBit(static_cast<uint8_t>(((*val_iter)[byte_id] >> bit_id) & uint8_t{1}));
                    }
                }

                /* After the bit plane has been constructed, we try an RLE */
                bit_vector::BitVectorView view(bit_plane.data(), bit_plane.size());

                const size_t num_bits = stride * num_significant_bits;

                bool was_previous_bit_set = view.IsCurrentBitSet();
                view.MoveToNextBit();

                uint8_t run_length_counter = 0;

                for (auto bit_index = 1; bit_index < num_bits; ++bit_index)
                {
                    const bool is_current_bit_set = view.IsCurrentBitSet();
                    view.MoveToNextBit();

                    /* We can create a maxiumum run length of 16 which is resembled by the counter == 15 */
                    if (was_previous_bit_set == is_current_bit_set && run_length_counter < 15 && bit_index < num_bits - 1)
                    {
                        /* If there is a run length */
                        ++run_length_counter;
                    } else
                    {
                        /* If the bits differ */
                        /* Write out the previous run length */
                        suffix_encodings.AppendBit(was_previous_bit_set);
                        suffix_encodings.AppendFourBits(run_length_counter);
                        was_previous_bit_set = is_current_bit_set;
                        run_length_counter = 0;
                    }
                }
            }


            iter += stride;
        } else
        {
            suffix_encodings.AppendBit(false);
            if (num_significant_bits > 0)
            {
            EncodeAndAppendPrefix(suffix_encodings, xor_suffixes[iter].GetSignificantBitsInBigEndianOrdering(), num_significant_bits); 
            }
            ++iter;
        }
    }

    /* Just append the remaining suffixes */
    for (auto iter = num_full_iterations * stride; iter < xor_suffixes.size(); ++iter)
    {
        suffix_encodings.AppendBit(false);
        const int num_significant_bits = xor_suffixes[iter].GetCountOfSignificantBits();
        if (num_significant_bits > 0)
        {
            EncodeAndAppendPrefix(suffix_encodings, xor_suffixes[iter].GetSignificantBitsInBigEndianOrdering(), num_significant_bits); 
        }
    }

    /* Calculate the amount of bytes needed to encode the suffixes */
    const size_t num_bytes_suffixes = sizeof(size_t) + suffix_encodings.size();
    cmc_debug_msg("num_bytes_suffixes: ", num_bytes_suffixes);

    /* Now, we fill a buffer for the suffixes */
    std::vector<uint8_t> serialized_suffix_data;
    serialized_suffix_data.reserve(num_bytes_suffixes);

    //TODO: This has to be adapted for parallel gathering of the data 
    /* At first, we store the number of bytes for the suffix-level */
    PushBackValueToByteStream(serialized_suffix_data, num_bytes_suffixes);

    /* Copy the serialized data to the buffer */
    std::copy_n(suffix_encodings.begin(), suffix_encodings.size(), std::back_insert_iterator(serialized_suffix_data));

    return serialized_suffix_data;
}


template<typename T>
std::vector<uint8_t>
LevelwisePrefixData<T>::EncodeLevelDataAsLeadingZeroCount() const
{
    return std::vector<uint8_t>();
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
