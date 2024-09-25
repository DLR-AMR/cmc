#ifndef CMC_T8_PREFIX_ENCODING_HXX
#define CMC_T8_PREFIX_ENCODING_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_prefix.hxx"
#include "utilities/cmc_span.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"

#include <vector>

namespace cmc
{

constexpr int kTwoBitPrefixEncodingOffset = 0b01;
constexpr int kFourBitPrefixEncodingStart = 0b1000;
constexpr int kFourBitPrefixEncodingOffset = 0b0011;
constexpr uint8_t kIndicatePrefixLongerThanByte = 0b1110;
constexpr int kMoreThanFourBitPrefixEncodingOffset = 0b1000;
constexpr uint8_t kSkipToNextByteInEncodingStream = 0b1111;

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
    encoded_prefixes.AppendBits(prefix_bits, prefix_length);
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
                const int prefix_length = prefixes[prefix_idx].GetCountOfSignificantBits();
                EncodeAndAppendPrefixLength(prefix_lengths, prefix_length);

                /* Afterwards, the actual prefix itself has to encoded and stored as well */
                EncodeAndAppendPrefix(prefix_encodings, prefixes[prefix_idx].GetSignificantBitsInBigEndianOrdering(), prefix_length);
            }
        }
        
    }
    
    /* We want to add a flag in the prefix lengths, that here (at the end of the process local encodings) is a boundary where we skip to the next byte in the prefix lengths as well as in the prefix encodings */
    EncodeAndAppendPrefixLengthByteBoundary(prefix_lengths);

    /* Pack all the data to a single stream */
    const size_t num_bits_indicator = prefix_indicator.size();
    const size_t num_bytes_indicator = prefix_indicator.size_bytes();
    const size_t num_bytes_prefix_length_encodings = prefix_lengths.size();
    const size_t num_bytes_prefix_encodings = prefix_encodings.size();

    /* Overall memory estimate (including the just computed number of bytes for this whole level) */
    const size_t num_bytes_level = 4 * sizeof(size_t) + 2 * num_bytes_indicator + num_bytes_prefix_length_encodings + num_bytes_prefix_encodings;

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
    std::copy_n(prefix_indicator.begin_bytes(), prefix_indicator.size_bytes(), std::back_insert_iterator(serialized_level_data));
    std::copy_n(prefix_lengths.begin(), prefix_lengths.size(), std::back_insert_iterator(serialized_level_data));
    std::copy_n(prefix_encodings.begin(), prefix_encodings.size(), std::back_insert_iterator(serialized_level_data));

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

        /* Afterwards, the actual prefix itself has to encoded and stored as well */
        EncodeAndAppendPrefix(suffix_encodings, suff_iter->GetSignificantBitsInBigEndianOrdering(), suffix_length);
    }

    /* We want to add a flag in the prefix lengths, that here (at the end of the process local encodings) is a boundary where we skip to the next byte in the prefix lengths as well as in the prefix encodings */
    EncodeAndAppendSuffixLengthByteBoundary(suffix_lengths);

    /* Pack all the data to a single stream */
    const size_t num_bits_indicator = suffix_indications.size();
    const size_t num_bytes_indicator = suffix_indications.size_bytes();
    const size_t num_bytes_suffix_length_encodings = suffix_lengths.size();
    const size_t num_bytes_suffix_encodings = suffix_encodings.size();

    /* Calculate the amount of bytes needed to encode the suffixes */
    const size_t num_bytes_suffixes = 4 * sizeof(size_t) + num_bytes_indicator + num_bytes_suffix_length_encodings + num_bytes_suffix_encodings;

    /* Now, we fill a buffer for the suffixes */
    std::vector<uint8_t> serialized_suffix_data;
    serialized_suffix_data.reserve(num_bytes_suffixes);

    //TODO: This has to be adapted for parallel gathering of the data 
    /* At first, we store the number of bytes for the suffix-level*/
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

}

#endif /* !CMC_T8_PREFIX_ENCODING_HXX */
