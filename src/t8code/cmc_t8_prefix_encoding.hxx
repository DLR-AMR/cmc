#ifndef CMC_T8_PREFIX_ENCODING_HXX
#define CMC_T8_PREFIX_ENCODING_HXX

#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_prefix.hxx"
#include "utilities/cmc_span.hxx"

#include <vector>

namespace cmc
{

//constexpr int kEncodingLength = 4;


template<typename T>
class LevelwisePrefixes
{
public:
    LevelwisePrefixes() = delete;

    LevelwisePrefixes(const int size_hint)
    {
        prefixes_.reserve(size_hint);
        prefix_indicator_.reserve(size_hint / CHAR_BIT + 1);
        prefix_indicator_.push_back(uint8_t{0});
    };

    void SetPrefixIndicatorBit(const bool has_prefix_found)
    {
        if (current_bit_position_ == CHAR_BIT)
        {
            current_bit_position_ = 0;
            prefix_indicator_.push_back(uint8_t{0});
        }
        if (has_prefix_found)
        {
            /* Set a one for this element indicating that this element was corresponds to a prefix */
            prefix_indicator_.back() |= (uint8_t{1} << current_bit_position_);
        } else
        {
            /* Set a zero for this element indicating that there is not any prefix associated */
            //TODO: Can we just move the bit position, since per bdefault all bits are zero
            //prefix_indicator_.back() &= ~(uint8_t{1} << current_bit_position_);
        }
        ++current_bit_position_;
    }

    void SetPrefix(const CompressionValue<sizeof(T)>& prefix)
    {
        prefixes_.push_back(prefix);
    }
    void SetPrefix(CompressionValue<sizeof(T)>&& prefix)
    {
        prefixes_.push_back(std::move(prefix));
    }


    std::vector<CompressionValue<sizeof(T)>> prefixes_;
    std::vector<uint8_t> prefix_indicator_;
    int current_bit_position_{0};
    //int num_elements_per_level{0};
};

    // We will choose this encoding:
    //Write function get significant bytes in big endian ordern with offset by one nibble
    //4-Bit Code -> 12 Möglichkeiten für Bits: [1 - 12] für ANzahl Bits
    //4 Leftover-Möglichkeiten:
    // 0000 -> Continues with the succeeding byte 1001 1011
    // 3 Leftover Möglichkeiten
    // => z.B. Prefix Length 14 : [0000 yyyy][yyyy yyyy][0010 yy--]
    // -> Prefix_length | Num bytes needed
    //      [1-4]              1 Byte
    //      [5-12]             2 Bytes
    //      [13-16]            3 Bytes
    //      [17-24]            4 Bytes
    //      [25-28]            5 Bytes
    //      [29-36]            6 Bytes

static inline
bool IsEven(const int i)
{
    return (i % 2 == 0 ? true : false);
}

static inline
void
SetPrefixLengthInEncoding_Offset_4_Stride_12(const int num_significant_bits, std::vector<uint8_t>& strided_variable)
{
    int num_remaining_bits = num_significant_bits;
    uint8_t encoded_bit_length;
    for (int i = 0; i < strided_variable.size(); ++i)
    {
        /* Set the prefix length in every other byte */
        if (IsEven(i))
        {
            /* Set the prefix length */
            if (num_remaining_bits > 12)
            {
                /* Set the front four bits to zero (they already should be zero from the offset-space integration) */
                strided_variable[i] &= 0x0F;
                num_remaining_bits -= 12;
            } else
            {
                /* Set the number of remaining bits in the front four bits */
                const uint8_t remaining_bits = static_cast<uint8_t>(num_remaining_bits);
                strided_variable[i] |= (remaining_bits << 4);
            }
        }
    }
}


#if 0
//Previous Version; a new test is below

/* serialized_variable in-out variable*/
template<typename T>
int
EncodeAndAppendLevelwisePrefixes(std::vector<uint8_t>& serialized_variable, const std::vector<LevelwisePrefixes<T>>& prefixes)
{
    int num_bytes = 0;
    int num_set_bits = 0;
    int num_prefixes = 0;
    //TODO: Revert the reverse iterator to normal ones in the level-wise for-loop (eventually not)

    int avg_pref_length = 0;
    int avg_pref_length_count = 0;
    int count_zero_prefix = 0;
    int count_smaller_two = 0;
    int count_smaller_four = 0;
    int count_smaller_eight = 0;
    int count_smaller_12 = 0;
    int count_smaller_16 = 0;
    int count_greater_16 = 0;
    /* Encode the prefixes of the data */
    for (auto lw_pfi_iter = prefixes.rbegin(); lw_pfi_iter != prefixes.rend(); ++lw_pfi_iter)
    {
        int index = 0;
        const std::vector<CompressionValue<sizeof(T)>>& prefixes = lw_pfi_iter->prefixes_;
        num_prefixes += prefixes.size();
        //cmc_debug_msg("\n\n Num prefixes on levl: ", num_prefixes);
        for (auto pfi_iter = lw_pfi_iter->prefix_indicator_.begin(); pfi_iter != lw_pfi_iter->prefix_indicator_.end(); ++pfi_iter)
        {
            /* Check if a prefix is given and if so encode the prefix and store it within the serialization */
            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if (CheckIfBitIsSet(*pfi_iter, bit_index))
                {
                    ++num_set_bits;
                    /* In this case, there is a prefix present, we need to encode it */
                    //cmc_assert(prefixes[index].GetCountOfSignificantBits() >= 0 && prefixes[index].GetCountOfSignificantBits() <= kMaxEncodingLength);
                    //if (prefixes[index].GetCountOfSignificantBits() > 32 || prefixes[index].GetCountOfSignificantBits() < 0)
                    //{
                    //    cmc_debug_msg(" Hier ist ein komischer Fall: ", index,", bit index: ", bit_index);
                    //}
                    const int num_significant_bits = prefixes[index].GetCountOfSignificantBits();


                    if (num_significant_bits == 0) {++count_zero_prefix;}
                    if (num_significant_bits <= 2) {++count_smaller_two;}
                    else if (num_significant_bits <= 4) {++count_smaller_four;}
                    else if (num_significant_bits <= 8) {++count_smaller_eight;}
                    else if (num_significant_bits <= 12) {++count_smaller_12;}
                    else if (num_significant_bits <= 16) {++count_smaller_16;}
                    else {++count_greater_16;}
                    //cmc_debug_msg("PrefixLength: ", num_significant_bits);
                    avg_pref_length += num_significant_bits;
                    ++avg_pref_length_count;


                    std::vector<uint8_t> strided_significant_bits = prefixes[index].GetSignificantBitsInBigEndianOrdering_Offset_4_Stride_12();

                    SetPrefixLengthInEncoding_Offset_4_Stride_12(num_significant_bits, strided_significant_bits);

                    num_bytes += strided_significant_bits.size();
                    /* Append the encoded prefix to the serialized variable */
                    std::copy_n(strided_significant_bits.begin(), strided_significant_bits.size(), std::back_inserter(serialized_variable));
                }

                ++index;
            }
        }
    }
    cmc_debug_msg("\n\n\nnum set bits: ", num_set_bits, "und num_prefixes: ", num_prefixes);
    cmc_debug_msg("The encoding of the levelwise prefixes took: ", num_bytes, " Bytes.");
    cmc_debug_msg("Average pref length: ", static_cast<float>(avg_pref_length) / static_cast<float>(avg_pref_length_count));
    cmc_debug_msg("Count zero prefixes: ", count_zero_prefix);
    cmc_debug_msg("Count smaller two prefixes: ", count_smaller_two);
    cmc_debug_msg("Count smaller four prefixes: ", count_smaller_four);
    cmc_debug_msg("Count smaller eight prefixes: ", count_smaller_eight);
    cmc_debug_msg("Count smaller 12 prefixes: ", count_smaller_12);
    cmc_debug_msg("Count smaller 16 prefixes: ", count_smaller_16);
    cmc_debug_msg("Count greater 16 prefixes: ", count_greater_16);

    return num_bytes;
}


#else

constexpr uint8_t kOneBitPrefix = 0x00;
constexpr uint8_t kTwoBitPrefix = 0x01;
constexpr uint8_t kNullifyTwoTailBits = ~((uint8_t) 0x03);
constexpr uint8_t kNullifyThirdAndFourthTailBits = ~((uint8_t) 0x0C);
constexpr uint8_t kPrefixContinues = 0x0F;



/* serialized_variable in-out variable*/
template<typename T>
std::pair<int, int>
EncodeAndAppendLevelwisePrefixes(std::vector<uint8_t>& serialized_variable, const std::vector<LevelwisePrefixes<T>>& prefixes)
{
    int num_bytes = 0;
    int num_prefixes = 0;

    int avg_pref_length = 0;
    int avg_pref_length_count = 0;
    int count_zero_prefix = 0;
    int count_smaller_two = 0;
    int count_smaller_four = 0;
    int count_smaller_eight = 0;
    int count_smaller_12 = 0;
    int count_smaller_16 = 0;
    int count_greater_16 = 0;

    int memory_estimate = 0;
    for (auto lw_pfi_iter = prefixes.rbegin(); lw_pfi_iter != prefixes.rend(); ++lw_pfi_iter)
    {
        memory_estimate += lw_pfi_iter->prefixes_.size();
    }

    std::vector<uint8_t> prefix_lengths;
    prefix_lengths.reserve(memory_estimate / 2 + 1);

    std::vector<uint8_t> prefix_encodings;
    prefix_encodings.reserve(memory_estimate + 1);

    //int start_offset_lengths = CHAR_BIT;
    int start_offset_encodings = CHAR_BIT;

    int tail_offset_lengths = 0;
    //int tail_offset_encodings = 0;

    /* Encode the prefixes of the data */
    for (auto lw_pfi_iter = prefixes.rbegin(); lw_pfi_iter != prefixes.rend(); ++lw_pfi_iter)
    {
        int index = 0;
        const std::vector<CompressionValue<sizeof(T)>>& prefixes = lw_pfi_iter->prefixes_;
        num_prefixes += prefixes.size();
        cmc_debug_msg("Next level\n\n mit num prefixes: ", prefixes.size());
        for (auto pfi_iter = lw_pfi_iter->prefix_indicator_.begin(); pfi_iter != lw_pfi_iter->prefix_indicator_.end(); ++pfi_iter)
        {
            /* Check if a prefix is given and if so encode the prefix and store it within the serialization */
            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if (index >= prefixes.size())
                {
                    break;
                }
                if (CheckIfBitIsSet(*pfi_iter, bit_index))
                {
                    const int num_significant_bits = prefixes[index].GetCountOfSignificantBits();
                    if (num_significant_bits == 0){cmc_debug_msg("Zero significant bits in pref encoding");}
                    std::vector<uint8_t> significant_bits = prefixes[index].GetSignificantBitsInBigEndianOrdering();

                    if (num_significant_bits == 0) {++count_zero_prefix;}
                    if (num_significant_bits <= 2) {++count_smaller_two;}
                    else if (num_significant_bits <= 4) {++count_smaller_four;}
                    else if (num_significant_bits <= 8) {++count_smaller_eight;}
                    else if (num_significant_bits <= 12) {++count_smaller_12;}
                    else if (num_significant_bits <= 16) {++count_smaller_16;}
                    else {++count_greater_16;}
                    //cmc_debug_msg("PrefixLength: ", num_significant_bits);
                    avg_pref_length += num_significant_bits;
                    ++avg_pref_length_count;

                    if (num_significant_bits <= 2)
                    {
                        //The prefix is at most two bits long
                        /* Encode the length of the prefix */
                        if (tail_offset_lengths <= 0)
                        {
                            prefix_lengths.push_back(uint8_t{0});
                            tail_offset_lengths = CHAR_BIT;
                        }

                        if (num_significant_bits == 1)
                        {
                            prefix_lengths.back() |= (kOneBitPrefix << (tail_offset_lengths - 2));
                        } else
                        {
                            prefix_lengths.back() |= (kTwoBitPrefix << (tail_offset_lengths - 2));
                        }

                        tail_offset_lengths -= 2;

                        /* Encode the prefix itself */
                        if (start_offset_encodings >= CHAR_BIT)
                        {
                            prefix_encodings.push_back((significant_bits[0] >> (CHAR_BIT - 2 - start_offset_encodings)));
                            start_offset_encodings = 2;
                        } else if (start_offset_encodings >= CHAR_BIT - 1)
                        {
                            const uint8_t inter_pref_encoding = significant_bits[0] >> (CHAR_BIT - 1);
                            prefix_encodings.back() |= (inter_pref_encoding << (CHAR_BIT - 1));
                            const uint8_t inter_pref_encoding1 = significant_bits[0] << 1;
                            prefix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - 1));

                            //prefix_encodings.back() |= ((significant_bits[0] >> (CHAR_BIT - 1)) << (CHAR_BIT - 1));
                            //prefix_encodings.push_back(((significant_bits[0] << 1) >> (CHAR_BIT - 1)));
                            start_offset_encodings = 1;
                        } else
                        {
                            prefix_encodings.back() |= (significant_bits[0] >> (CHAR_BIT - 2 - start_offset_encodings));
                            start_offset_encodings += 2;
                        }

                    } else if (num_significant_bits <= CHAR_BIT)
                    {
                        //The prefix is not longer than one byte 
                        const uint8_t pref_length = static_cast<uint8_t>(0x08) + (num_significant_bits - 3);

                        /* Encode the length of the prefix */
                        if (tail_offset_lengths <= 0)
                        {
                            prefix_lengths.push_back(pref_length);
                            tail_offset_lengths = 4;
                        } else if (tail_offset_lengths <= 2)
                        {
                            /* The prefix length is split via two bytes */
                            prefix_lengths.back() |= (pref_length >> 2);
                            prefix_lengths.push_back((pref_length & kNullifyThirdAndFourthTailBits) << 6);
                            tail_offset_lengths = CHAR_BIT - 2;
                        } else
                        {
                            prefix_lengths.back() |= (pref_length << (tail_offset_lengths - 4));
                            tail_offset_lengths -= 4;
                        }

                        /* Encode the prefix itself */
                        if (num_significant_bits <= CHAR_BIT - start_offset_encodings)
                        {
                            prefix_encodings.back() |= (significant_bits[0] >> (CHAR_BIT - num_significant_bits - start_offset_encodings)); //invalid write
                        } else
                        {
                            /* Needs to be split */
                            const uint8_t inter_pref_encoding = significant_bits[0] >> start_offset_encodings;
                            prefix_encodings.back() |= inter_pref_encoding << start_offset_encodings;
                            const uint8_t inter_pref_encoding1 = significant_bits[0] << (num_significant_bits - (CHAR_BIT - start_offset_encodings));
                            prefix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - (num_significant_bits - (CHAR_BIT - start_offset_encodings))));


                            //prefix_encodings.back() |= ((significant_bits[0] >> start_offset_encodings) << start_offset_encodings);
                            //prefix_encodings.push_back((significant_bits[0] << (num_significant_bits - (CHAR_BIT - start_offset_encodings))) >> (CHAR_BIT - (num_significant_bits - (CHAR_BIT - start_offset_encodings))));
                            start_offset_encodings = num_significant_bits - (CHAR_BIT - start_offset_encodings);
                        }
                    } else
                    {
                        //The prefix is longer than a single byte
                        int num_bits_encoded = 0;
                        int remaining_bits = num_significant_bits;

                        /* Encode the length (indicating that the prefix coiontinues) */
                        while (remaining_bits > CHAR_BIT)
                        {
                            const uint8_t pref_length = kPrefixContinues;
                            /* Encode the length of the prefix */
                            if (tail_offset_lengths <= 0)
                            {
                                prefix_lengths.push_back(pref_length);
                                tail_offset_lengths = 4;
                            } else if (tail_offset_lengths <= 2)
                            {
                                /* The prefix length is split via two bytes */
                                prefix_lengths.back() |= (pref_length >> 2);
                                prefix_lengths.push_back((pref_length & kNullifyThirdAndFourthTailBits) << 6);
                                tail_offset_lengths = CHAR_BIT - 2;
                            } else
                            {
                                prefix_lengths.back() |= (pref_length << (tail_offset_lengths - 4));
                                tail_offset_lengths -= 4;
                            }
                            remaining_bits -= CHAR_BIT;
                            
                            //Maybe todo:
                            /* We store a whole nibble for the length after a prefix indication (if possible, otherwise another "prefix continues" sequence) */
                            //if(remaining_bits <= 16) {break;}
                        }

                        /* Store the remaing bit length as a nibble */
                        //const uint8_t pref_length = remaining_bits - 1;
                        const uint8_t pref_length = remaining_bits; //We need to store the remaining bits as they are, (otherwise it may be problematic for example 9,16,17,24,25,... bit prefixes)
                        if (tail_offset_lengths <= 0)
                        {
                            prefix_lengths.push_back(pref_length);
                            tail_offset_lengths = 4;
                        } else if (tail_offset_lengths <= 2)
                        {
                            /* The prefix length is split via two bytes */
                            prefix_lengths.back() |= (pref_length >> 2);
                            prefix_lengths.push_back((pref_length & kNullifyThirdAndFourthTailBits) << 6);
                            tail_offset_lengths = CHAR_BIT - 2;
                        } else
                        {
                            prefix_lengths.back() |= (pref_length << (tail_offset_lengths - 4));
                            tail_offset_lengths -= 4;
                        }
                        
                        /* Encode the prefix itself */
                        int pref_index = 0;
                        while (num_bits_encoded < num_significant_bits)
                        {
                            const int num_considered_bits = (num_significant_bits - num_bits_encoded > CHAR_BIT ? CHAR_BIT : num_significant_bits - num_bits_encoded);
                            /* The bit hsift is performed in order to ensure that there are zeros at the end (however, the not considered bits in the byte should already been zero) */
                            const uint8_t significant_byte = (num_considered_bits == CHAR_BIT ? significant_bits[pref_index] : ((significant_bits[pref_index] >> (CHAR_BIT - num_considered_bits)) << (CHAR_BIT - num_considered_bits)));
                            
                            if (num_considered_bits <= CHAR_BIT - start_offset_encodings)
                            {
                                prefix_encodings.back() |= (significant_bits[0] >> (CHAR_BIT - num_considered_bits - start_offset_encodings)); 
                            } else
                            {
                                /* Needs to be split */
                                const uint8_t inter_pref_encoding = significant_bits[0] >> start_offset_encodings;
                                prefix_encodings.back() |= inter_pref_encoding << start_offset_encodings;
                                const uint8_t inter_pref_encoding1 = significant_bits[0] << (num_considered_bits - (CHAR_BIT - start_offset_encodings));
                                prefix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings))));


                                //prefix_encodings.back() |= ((significant_bits[0] >> start_offset_encodings) << start_offset_encodings);
                                //prefix_encodings.push_back((significant_bits[0] << (num_considered_bits - (CHAR_BIT - start_offset_encodings))) >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings))));
                                start_offset_encodings = num_considered_bits - (CHAR_BIT - start_offset_encodings);
                            }
                            num_bits_encoded += num_considered_bits;
                            ++pref_index;
                        }
                    }
                }

                ++index;
            }
        }
    }

    //Copy the prefix length as well as the prefix encodings */
    /* Append the encoded prefix to the serialized variable */
    std::copy_n(prefix_lengths.begin(), prefix_lengths.size(), std::back_inserter(serialized_variable));
    std::copy_n(prefix_encodings.begin(), prefix_encodings.size(), std::back_inserter(serialized_variable));

    cmc_debug_msg("\n\n\nnum_prefixes: ", num_prefixes);
    cmc_debug_msg("Average pref length: ", static_cast<float>(avg_pref_length) / static_cast<float>(avg_pref_length_count));
    cmc_debug_msg("Count zero prefixes: ", count_zero_prefix);
    cmc_debug_msg("Count smaller two prefixes: ", count_smaller_two);
    cmc_debug_msg("Count smaller four prefixes: ", count_smaller_four);
    cmc_debug_msg("Count smaller eight prefixes: ", count_smaller_eight);
    cmc_debug_msg("Count smaller 12 prefixes: ", count_smaller_12);
    cmc_debug_msg("Count smaller 16 prefixes: ", count_smaller_16);
    cmc_debug_msg("Count greater 16 prefixes: ", count_greater_16);

    return std::make_pair(prefix_lengths.size(), prefix_encodings.size());
}





/* serialized_variable in-out variable*/
template<typename T>
std::pair<int, int>
EncodeAndAppendLevelwisePrefixesEGU(std::vector<uint8_t>& serialized_variable, const std::vector<LevelwisePrefixes<T>>& lvl_ws_prefixes)
{
    int num_bytes = 0;
    int num_prefixes = 0;

    int memory_estimate = 0;
    for (auto lw_pfi_iter = lvl_ws_prefixes.rbegin(); lw_pfi_iter != lvl_ws_prefixes.rend(); ++lw_pfi_iter)
    {
        memory_estimate += lw_pfi_iter->prefixes_.size();
    }

    std::vector<uint8_t> prefix_lengths;
    prefix_lengths.reserve(memory_estimate / 2 + 1);

    std::vector<uint8_t> prefix_encodings;
    prefix_encodings.reserve(memory_estimate + 1);

    /* Encode the prefixes of the data */
    for (auto lw_pfi_iter = lvl_ws_prefixes.rbegin(); lw_pfi_iter != lvl_ws_prefixes.rend(); ++lw_pfi_iter)
    {
        int index = 0;
        const std::vector<CompressionValue<sizeof(T)>>& prefixes = lw_pfi_iter->prefixes_;
        num_prefixes += prefixes.size();
        /* We start each level with a whole new byte */
        int start_offset_encodings = CHAR_BIT;
        int tail_offset_lengths = 0;
        cmc_debug_msg("Next level\n\n mit num prefixes: ", prefixes.size());
        for (auto pfi_iter = lw_pfi_iter->prefix_indicator_.begin(); pfi_iter != lw_pfi_iter->prefix_indicator_.end(); ++pfi_iter)
        {
            /* Check if a prefix is given and if so encode the prefix and store it within the serialization */
            for (int bit_index = 0; bit_index < CHAR_BIT; ++bit_index)
            {
                if (index >= prefixes.size())
                {
                    break;
                }
                if (CheckIfBitIsSet(*pfi_iter, bit_index))
                {
                    const int num_significant_bits = prefixes[index].GetCountOfSignificantBits();
                    if (num_significant_bits == 0){cmc_debug_msg("Zero significant bits in pref encoding");}
                    if (num_significant_bits == 0){continue;}
                    std::vector<uint8_t> significant_bits = prefixes[index].GetSignificantBitsInBigEndianOrdering();

                    if (num_significant_bits <= 2)
                    {
                        /* Encode the length of the prefix */
                        if (tail_offset_lengths <= 0)
                        {
                            prefix_lengths.push_back(uint8_t{0});
                            tail_offset_lengths = CHAR_BIT;
                        }

                        if (num_significant_bits == 1)
                        {
                            prefix_lengths.back() |= (kOneBitPrefix << (tail_offset_lengths - 2));
                        } else
                        {
                            prefix_lengths.back() |= (kTwoBitPrefix << (tail_offset_lengths - 2));
                        }

                        tail_offset_lengths -= 2;


                        /* Encode the prefix itself */
                        if (start_offset_encodings >= CHAR_BIT)
                        {
                            //prefix_encodings.push_back((significant_bits[0] >> (CHAR_BIT - 2 - start_offset_encodings)));
                            prefix_encodings.push_back((significant_bits[0] >> (CHAR_BIT - num_significant_bits)));
                            start_offset_encodings = num_significant_bits;
                        } else if (start_offset_encodings >= CHAR_BIT - 1 && num_significant_bits == 2)
                        {
                            const uint8_t inter_pref_encoding = significant_bits[0] >> (CHAR_BIT - 1);
                            prefix_encodings.back() |= (inter_pref_encoding << (CHAR_BIT - 1));
                            const uint8_t inter_pref_encoding1 = significant_bits[0] << 1;
                            prefix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - 1));

                            start_offset_encodings = 1;
                        } else
                        {
                            prefix_encodings.back() |= (significant_bits[0] >> (CHAR_BIT - num_significant_bits - start_offset_encodings));
                            start_offset_encodings += num_significant_bits;
                        }

                    } else if (num_significant_bits <= CHAR_BIT)
                    {
                        //The prefix is not longer than one byte 
                        const uint8_t pref_length = static_cast<uint8_t>(0x08) + (num_significant_bits - 3);

                        /* Encode the length of the prefix */
                        if (tail_offset_lengths <= 0)
                        {
                            prefix_lengths.push_back(pref_length << 4);
                            tail_offset_lengths = 4;
                        } else if (tail_offset_lengths <= 2)
                        {
                            /* The prefix length is split via two bytes */
                            prefix_lengths.back() |= (pref_length >> 2);
                            prefix_lengths.push_back((pref_length & kNullifyThirdAndFourthTailBits) << 6);
                            //prefix_lengths.push_back(pref_length << 6);
                            tail_offset_lengths = CHAR_BIT - 2;
                        } else
                        {
                            prefix_lengths.back() |= (pref_length << (tail_offset_lengths - 4));
                            tail_offset_lengths -= 4;
                        }

                        /* Encode the prefix itself */
                        if (start_offset_encodings >= CHAR_BIT)
                        {
                            prefix_encodings.push_back(significant_bits[0] >> (CHAR_BIT - num_significant_bits));
                            start_offset_encodings = num_significant_bits;
                        } else if (num_significant_bits <= CHAR_BIT - start_offset_encodings)
                        {
                            prefix_encodings.back() |= (significant_bits[0] >> (CHAR_BIT - num_significant_bits - start_offset_encodings));
                            start_offset_encodings += num_significant_bits;
                        } else
                        {
                            /* Needs to be split */
                            //Here is an update needed, if the prefix is longer, we need the middle bytes need to be aligned front first (and not the back)
                            const uint8_t inter_pref_encoding = significant_bits[0] >> start_offset_encodings;
                            prefix_encodings.back() |= inter_pref_encoding << start_offset_encodings;
                            //const uint8_t inter_pref_encoding1 = significant_bits[0] << (num_significant_bits - (CHAR_BIT - start_offset_encodings));
                            //prefix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - (num_significant_bits - (CHAR_BIT - start_offset_encodings))));
                            const uint8_t inter_pref_encoding1 = significant_bits[0] << (CHAR_BIT - start_offset_encodings);
                            prefix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - (num_significant_bits - (CHAR_BIT - start_offset_encodings))));

                            //prefix_encodings.back() |= ((significant_bits[0] >> start_offset_encodings) << start_offset_encodings);
                            //prefix_encodings.push_back((significant_bits[0] << (num_significant_bits - (CHAR_BIT - start_offset_encodings))) >> (CHAR_BIT - (num_significant_bits - (CHAR_BIT - start_offset_encodings))));
                            start_offset_encodings = num_significant_bits - (CHAR_BIT - start_offset_encodings);
                        }
                    } else
                    {
                        //The prefix is longer than a single byte
                        int remaining_bits = num_significant_bits;

                        /* Encode the length (indicating that the prefix continues) */
                        while (remaining_bits > CHAR_BIT)
                        {
                            const uint8_t pref_length = kPrefixContinues;
                            /* Encode the length of the prefix */
                            if (tail_offset_lengths <= 0)
                            {
                                prefix_lengths.push_back(pref_length << 4);
                                tail_offset_lengths = 4;
                            } else if (tail_offset_lengths <= 2)
                            {
                                /* The prefix length is split via two bytes */
                                prefix_lengths.back() |= (pref_length >> 2);
                                prefix_lengths.push_back((pref_length & kNullifyThirdAndFourthTailBits) << 6);
                                //prefix_lengths.push_back(pref_length << static_cast<uint8_t>(6));
                                tail_offset_lengths = CHAR_BIT - 2;
                            } else
                            {
                                prefix_lengths.back() |= (pref_length << (tail_offset_lengths - 4));
                                tail_offset_lengths -= 4;
                            }
                            remaining_bits -= CHAR_BIT;
                        }

                        /* Store the remaing bit length as a nibble */
                        const uint8_t pref_length = remaining_bits; //We need to store the remaining bits as they are, (otherwise it may be problematic for example 9,16,17,24,25,... bit prefixes)
                        if (tail_offset_lengths <= 0)
                        {
                            prefix_lengths.push_back(pref_length << 4);
                            tail_offset_lengths = 4;
                        } else if (tail_offset_lengths <= 2)
                        {
                            /* The prefix length is split via two bytes */
                            prefix_lengths.back() |= (pref_length >> 2);
                            prefix_lengths.push_back((pref_length & kNullifyThirdAndFourthTailBits) << 6);
                            //prefix_lengths.push_back(pref_length << 6);
                            tail_offset_lengths = CHAR_BIT - 2;
                        } else
                        {
                            prefix_lengths.back() |= (pref_length << (tail_offset_lengths - 4));
                            tail_offset_lengths -= 4;
                        }
                        
                        /* Encode the prefix itself */
                        int num_bits_encoded = 0;
                        int pref_index = 0;
                        while (num_bits_encoded < num_significant_bits)
                        {
                            const int num_considered_bits = (num_significant_bits - num_bits_encoded > CHAR_BIT ? CHAR_BIT : num_significant_bits - num_bits_encoded);
                            /* The bit hsift is performed in order to ensure that there are zeros at the end (however, the not considered bits in the byte should already been zero) */
                            //const uint8_t significant_byte = (num_considered_bits == CHAR_BIT ? significant_bits[index] : ((significant_bits[index] >> (CHAR_BIT - num_considered_bits)) << (CHAR_BIT - num_considered_bits)));
                            //if (pref_index >= significant_bits.size()){cmc_debug_msg("Pref index too large: ", pref_index, " und size of significant bits = ", significant_bits.size(), ", current prefixes.size() = ", prefixes.size());}
                            const uint8_t significant_byte = significant_bits[pref_index];

                            if (start_offset_encodings >= CHAR_BIT)
                            {
                                prefix_encodings.push_back(uint8_t{0});
                                start_offset_encodings = 0;
                            }

                            if (num_considered_bits <= CHAR_BIT - start_offset_encodings)
                            {
                                prefix_encodings.back() |= (significant_byte >> ((CHAR_BIT - start_offset_encodings) - num_considered_bits));
                                start_offset_encodings += num_considered_bits;
                            } else
                            {
                                /* Needs to be split */
                                const uint8_t inter_pref_encoding = significant_byte >> start_offset_encodings;
                                prefix_encodings.back() |= inter_pref_encoding << start_offset_encodings;
                                const uint8_t inter_pref_encoding1 = significant_byte << (CHAR_BIT - start_offset_encodings);
                                prefix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings))));


                                //prefix_encodings.back() |= ((significant_bits[0] >> start_offset_encodings) << start_offset_encodings);
                                //prefix_encodings.push_back((significant_bits[0] << (num_considered_bits - (CHAR_BIT - start_offset_encodings))) >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings))));
                                start_offset_encodings = num_considered_bits - (CHAR_BIT - start_offset_encodings);
                            }
                            num_bits_encoded += num_considered_bits;
                            ++pref_index;
                        }
                    }
                }

                ++index;
            }
        }
        
    }

    //Copy the prefix length as well as the prefix encodings */
    /* Append the encoded prefix to the serialized variable */
    std::copy_n(prefix_lengths.begin(), prefix_lengths.size(), std::back_inserter(serialized_variable));
    std::copy_n(prefix_encodings.begin(), prefix_encodings.size(), std::back_inserter(serialized_variable));

    return std::make_pair(prefix_lengths.size(), prefix_encodings.size());
}

#endif



//DO the above for the suffix
//...




/* serialized_variable in-out variable*/
template<int N>
std::tuple<int, int, int>
EncodeAndAppendSuffixesEGU(std::vector<uint8_t>& serialized_variable, const std::vector<CompressionValue<N>>& suffixes)
{
    std::vector<uint8_t> suffix_lengths;
    suffix_lengths.reserve(suffixes.size() / 2 + 1);

    std::vector<uint8_t> suffix_encodings;
    suffix_encodings.reserve(suffixes.size() * 2);

    std::vector<uint8_t> rle;
    rle.reserve(suffixes.size());

    uint8_t current_rl = 0;
    int last_state = 0;
    int avg_length = 0;
    int rl_count = 0;
    bool front_indicator = false;

    cmc_debug_msg("Size of sffixes before encoding: ", suffixes.size());
    /* Generate run length encoding */
    for (auto iter = suffixes.begin(); iter != suffixes.end(); ++iter)
    {
        const int num_significant_bits = iter->GetCountOfSignificantBits();

        if (last_state == 1)
        {
            if (current_rl == 8 && num_significant_bits > 0)
            {
                if (front_indicator)
                {
                    rle.back() |= ((uint8_t)0x08 & (current_rl - 1)) << 4;
                    front_indicator = !front_indicator;
                } else
                {
                    rle.push_back((uint8_t)0x08 & (current_rl - 1));
                    front_indicator = !front_indicator;
                }

                avg_length += current_rl;
                ++rl_count;
                current_rl = 1;
                last_state = 1;
            }
            else if (num_significant_bits > 0)
            {
                ++current_rl;
            } else
            {
                if (front_indicator)
                {
                    rle.back() |= ((uint8_t)0x08 & (current_rl - 1)) << 4;
                    front_indicator = !front_indicator;
                } else
                {
                    rle.push_back((uint8_t)0x08 & (current_rl - 1));
                    front_indicator = !front_indicator;
                }

                avg_length += current_rl;
                ++rl_count;
                current_rl = 1;
                last_state = 0;
            }
        } else
        {
            if (current_rl == 8 && num_significant_bits == 0)
            {
                if (front_indicator)
                {
                    rle.back() |= ((current_rl - 1)) << 4;
                    front_indicator = !front_indicator;
                } else
                {
                    rle.push_back((current_rl - 1));
                    front_indicator = !front_indicator;
                }

                avg_length += current_rl;
                ++rl_count;
                current_rl = 1;
                last_state = 0;
            }
            else if (num_significant_bits == 0)
            {
                ++current_rl;
            } else
            {
                if (front_indicator)
                {
                    rle.back() |= ((current_rl - 1)) << 4;
                    front_indicator = !front_indicator;
                } else
                {
                    rle.push_back((current_rl - 1));
                    front_indicator = !front_indicator;
                }

                avg_length += current_rl;
                ++rl_count;
                current_rl = 1;
                last_state = 1;
            }
        }
    }

    cmc_debug_msg("Size of rle: ", rle.size());
    std::copy_n(rle.begin(), rle.size(), std::back_inserter(serialized_variable));


    int start_offset_encodings = CHAR_BIT;
    int tail_offset_lengths = 0;

    /* Encode the suffix and it's length */
    for (auto iter = suffixes.begin(); iter != suffixes.end(); ++iter)
    {
        const int num_significant_bits = iter->GetCountOfSignificantBits();
        //cmc_debug_msg("Num significant bits in suffix: ", num_significant_bits);
        if (num_significant_bits > 0)
        {
            /* There is a suffix left */
            std::vector<uint8_t> significant_bits = iter->GetSignificantBitsInBigEndianOrdering();

            if (num_significant_bits <= 2)
            {
                //The prefix is at most two bits long
                /* Encode the length of the prefix */
                if (tail_offset_lengths <= 0)
                {
                    suffix_lengths.push_back(uint8_t{0});
                    tail_offset_lengths = CHAR_BIT;
                }

                if (num_significant_bits == 1)
                {
                    suffix_lengths.back() |= (kOneBitPrefix << (tail_offset_lengths - 2));
                } else
                {
                    suffix_lengths.back() |= (kTwoBitPrefix << (tail_offset_lengths - 2));
                }

                tail_offset_lengths -= 2;

                /* Encode the prefix itself */
                if (start_offset_encodings >= CHAR_BIT)
                {
                    //prefix_encodings.push_back((significant_bits[0] >> (CHAR_BIT - 2 - start_offset_encodings)));
                    suffix_encodings.push_back((significant_bits[0] >> (CHAR_BIT - num_significant_bits)));
                    start_offset_encodings = num_significant_bits;
                } else if (start_offset_encodings >= CHAR_BIT - 1 && num_significant_bits == 2)
                {
                    const uint8_t inter_pref_encoding = significant_bits[0] >> (CHAR_BIT - 1);
                    suffix_encodings.back() |= (inter_pref_encoding << (CHAR_BIT - 1));
                    const uint8_t inter_pref_encoding1 = significant_bits[0] << 1;
                    suffix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - 1));

                    start_offset_encodings = 1;
                } else
                {
                    suffix_encodings.back() |= (significant_bits[0] >> (CHAR_BIT - num_significant_bits - start_offset_encodings));
                    start_offset_encodings += num_significant_bits;
                }

            } else if (num_significant_bits <= CHAR_BIT)
            {
                //The prefix is not longer than one byte 
                const uint8_t pref_length = static_cast<uint8_t>(0x08) + (num_significant_bits - 3);

                /* Encode the length of the prefix */
                if (tail_offset_lengths <= 0)
                {
                    suffix_lengths.push_back(pref_length << 4);
                    tail_offset_lengths = 4;
                } else if (tail_offset_lengths <= 2)
                {
                    /* The prefix length is split via two bytes */
                    suffix_lengths.back() |= (pref_length >> 2);
                    suffix_lengths.push_back((pref_length & kNullifyThirdAndFourthTailBits) << 6);
                    //prefix_lengths.push_back(pref_length << 6);
                    tail_offset_lengths = CHAR_BIT - 2;
                } else
                {
                    suffix_lengths.back() |= (pref_length << (tail_offset_lengths - 4));
                    tail_offset_lengths -= 4;
                }

                /* Encode the prefix itself */
                if (start_offset_encodings >= CHAR_BIT)
                {
                    suffix_encodings.push_back(significant_bits[0] >> (CHAR_BIT - num_significant_bits));
                    start_offset_encodings = num_significant_bits;
                } else if (num_significant_bits <= CHAR_BIT - start_offset_encodings)
                {
                    suffix_encodings.back() |= (significant_bits[0] >> (CHAR_BIT - num_significant_bits - start_offset_encodings));
                    start_offset_encodings += num_significant_bits;
                } else
                {
                    /* Needs to be split */
                    //Here is the error, we cannot do "<< start_offset_encodings", because the front gets truncated !!
                    const uint8_t inter_pref_encoding = significant_bits[0] >> start_offset_encodings;
                    suffix_encodings.back() |= (inter_pref_encoding << start_offset_encodings);
                    //const uint8_t inter_pref_encoding1 = significant_bits[0] << (num_significant_bits - (CHAR_BIT - start_offset_encodings));
                    //prefix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - (num_significant_bits - (CHAR_BIT - start_offset_encodings))));
                    const uint8_t inter_pref_encoding1 = significant_bits[0] << (CHAR_BIT - start_offset_encodings);
                    suffix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - (num_significant_bits - (CHAR_BIT - start_offset_encodings))));

                    //prefix_encodings.back() |= ((significant_bits[0] >> start_offset_encodings) << start_offset_encodings);
                    //prefix_encodings.push_back((significant_bits[0] << (num_significant_bits - (CHAR_BIT - start_offset_encodings))) >> (CHAR_BIT - (num_significant_bits - (CHAR_BIT - start_offset_encodings))));
                    start_offset_encodings = num_significant_bits - (CHAR_BIT - start_offset_encodings);
                }
            } else
            {
                //[cmc] DEBUG: Next encoding bytes: 10110101111100010111010101110001
                //[cmc] DEBUG: Num bytes for prefix allocated: 2 und Pref length is: 15 und byte pos: 0 und bit position: 0
                //[cmc] DEBUG:  Encodign byte pos: 0 und bit position: 0
                //[cmc] DEBUG: Decoded prefix hat num_bits: 15
                //[cmc] DEBUG: Before Application of prefix single: trail bit pos: 17 und num significant bits: 15
                //[cmc] DEBUG: Previous Value: num_sigificant bits: 15, in big endian order: 010000111000000 00000000000000000
                //[cmc] DEBUG: SUffix: num_sigificant bits: 15, in big endian order: 00000000000000001110101001110001
                //[cmc] DEBUG: After Application of prefix single: trail bit pos: 2 und num significant bits: 30
                //[cmc] DEBUG: Prefixed Value: num_sigificant bits: 30, in big endian order: 01000011100000001110001111010100
                //[cmc] DEBUG: Next encoding bytes: 00100001101101011111000101110101
                //[cmc] DEBUG: Num bytes for prefix allocated: 2 und Pref length is: 15 und byte pos: 1 und bit position: 0
                //[cmc] DEBUG:  Encodign byte pos: 1 und bit position: 7
                //[cmc] DEBUG: Decoded prefix hat num_bits: 15
                //[cmc] DEBUG: Before Application of prefix single: trail bit pos: 17 und num significant bits: 15
                //[cmc] DEBUG: Previous Value: num_sigificant bits: 15, in big endian order: 01000011100000000000000000000000
                //[cmc] DEBUG: SUffix: num_sigificant bits: 15, in big endian order: 00000000000000001110101001111000
                //[cmc] DEBUG: After Application of prefix single: trail bit pos: 2 und num significant bits: 30
                //[cmc] DEBUG: Prefixed Value: num_sigificant bits: 30, in big endian order: 01000011100000001111000111010100

                //[cmc] DEBUG: Decompressed Data: 257.78, initial data: 257.78 und Error: 0
                //[cmc] DEBUG: Initial Value: 010000111000000 01110001111010100
                //[cmc] DEBUG: Decompr Value: 010000111000000 01110001111010100
                //[cmc] DEBUG: Decompressed Data: 257.889, initial data: 257.78 und Error: 0.109375
                //[cmc] DEBUG: Initial Value: 010000111000000 01110001111010100
                //[cmc] DEBUG: Decompr Value: 010000111000000 01111000111010100

                //[cmc] DEBUG: PrefixENcodings: 0, bitset: 00000000 01110001
                //[cmc] DEBUG: The actual suffix in big endian: 11101010 01110001
                //[cmc] DEBUG: PrefixENcodings: 1, bitset: 00000000 01110101
                //[cmc] DEBUG: The actual suffix in big endian: 11101010 01110001
                //[cmc] DEBUG: PrefixENcodings: 2, bitset: 00000000 11110001 -> This should be: 1110001x
                //[cmc] DEBUG: The actual suffix in big endian: 00110111 10100001
                //[cmc] DEBUG: PrefixENcodings: 3, bitset: 000000001 0110101
                //[cmc] DEBUG: The actual suffix in big endian: 00110111 10100001


                //[cmc] DEBUG: PrefixENcodings: 0, bitset: 00000000 01110001
                //[cmc] DEBUG: PrefixENcodings: 1, bitset: 00000000 0 1110101
                //[cmc] DEBUG: PrefixENcodings: 2, bitset: 00000000 11110001
                //[cmc] DEBUG: PrefixENcodings: 3, bitset: 00000000 10 110101

                //The prefix is longer than a single byte
                int remaining_bits = num_significant_bits;

                /* Encode the length (indicating that the prefix continues) */
                while (remaining_bits > CHAR_BIT)
                {
                    const uint8_t pref_length = kPrefixContinues;
                    /* Encode the length of the prefix */
                    if (tail_offset_lengths <= 0)
                    {
                        suffix_lengths.push_back(pref_length << 4);
                        tail_offset_lengths = 4;
                    } else if (tail_offset_lengths <= 2)
                    {
                        /* The prefix length is split via two bytes */
                        suffix_lengths.back() |= (pref_length >> 2);
                        suffix_lengths.push_back((pref_length & kNullifyThirdAndFourthTailBits) << 6);
                        //prefix_lengths.push_back(pref_length << static_cast<uint8_t>(6));
                        tail_offset_lengths = CHAR_BIT - 2;
                    } else
                    {
                        suffix_lengths.back() |= (pref_length << (tail_offset_lengths - 4));
                        tail_offset_lengths -= 4;
                    }
                    remaining_bits -= CHAR_BIT;
                    
                    //Maybe todo:
                    /* We store a whole nibble for the length after a prefix indication (if possible, otherwise another "prefix continues" sequence) */
                    //if(remaining_bits <= 16) {break;}
                }

                /* Store the remaing bit length as a nibble */
                //const uint8_t pref_length = remaining_bits - 1;
                const uint8_t pref_length = remaining_bits; //We need to store the remaining bits as they are, (otherwise it may be problematic for example 9,16,17,24,25,... bit prefixes)
                if (tail_offset_lengths <= 0)
                {
                    suffix_lengths.push_back(pref_length << 4);
                    tail_offset_lengths = 4;
                } else if (tail_offset_lengths <= 2)
                {
                    /* The prefix length is split via two bytes */
                    suffix_lengths.back() |= (pref_length >> 2);
                    suffix_lengths.push_back((pref_length & kNullifyThirdAndFourthTailBits) << 6);
                    //prefix_lengths.push_back(pref_length << 6);
                    tail_offset_lengths = CHAR_BIT - 2;
                } else
                {
                    suffix_lengths.back() |= (pref_length << (tail_offset_lengths - 4));
                    tail_offset_lengths -= 4;
                }
                
                /* Encode the suffix itself */
                int num_bits_encoded = 0;
                int pref_index = 0;
                // 01110001 1110101 |0 Start: byte 1, bit 7, Significant Bits: 15 (0xxx xxxx) (111000)
                ////// Old try for tests
                while (num_bits_encoded < num_significant_bits)
                {
                    const int num_considered_bits = (num_significant_bits - num_bits_encoded > CHAR_BIT ? CHAR_BIT : num_significant_bits - num_bits_encoded);
                    /* The bit hsift is performed in order to ensure that there are zeros at the end (however, the not considered bits in the byte should already been zero) */
                    //const uint8_t significant_byte = (num_considered_bits == CHAR_BIT ? significant_bits[index] : ((significant_bits[index] >> (CHAR_BIT - num_considered_bits)) << (CHAR_BIT - num_considered_bits)));
                    //if (pref_index >= significant_bits.size()){cmc_debug_msg("Pref index too large: ", pref_index, " und size of significant bits = ", significant_bits.size(), ", current prefixes.size() = ", prefixes.size());}
                    const uint8_t significant_byte = significant_bits[pref_index];

                    if (start_offset_encodings >= CHAR_BIT)
                    {
                        //suffix_encodings.push_back(uint8_t{0});
                        //start_offset_encodings = 0;
                        suffix_encodings.push_back(significant_byte >> (CHAR_BIT - num_considered_bits));
                        start_offset_encodings = num_considered_bits;
                    } else if (num_considered_bits <= CHAR_BIT - start_offset_encodings)
                    {
                        suffix_encodings.back() |= (significant_byte >> ((CHAR_BIT - start_offset_encodings) - num_considered_bits));
                        start_offset_encodings += num_considered_bits;
                    } else
                    {
                        //The problem: if the the prefix/suffx strechtes over three or more bytes, al "middle" needs to be reversly filled
                        /* Needs to be split */
                        
                            const uint8_t inter_pref_encoding = significant_byte >> start_offset_encodings;
                            suffix_encodings.back() |= inter_pref_encoding << start_offset_encodings;

                        //if (iter == suffixes.begin() + 1)
                        //{
                        //    const uint8_t ival = inter_pref_encoding << start_offset_encodings;
                        //    uint16_t interimval{0};
                        //    std::memcpy(&interimval, &ival, 1);
                        //    cmc_debug_msg("Iterim val1 is: ", std::bitset<16>(interimval));
                        //}
                        const uint8_t inter_pref_encoding1 = significant_byte << (CHAR_BIT - start_offset_encodings);
                        suffix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings))));
                        //suffix_encodings.push_back(inter_pref_encoding1);
                        
                        //if (iter == suffixes.begin() + 1)
                        //{
                        //    const uint8_t ival = inter_pref_encoding1 >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings)));
                        //    uint16_t interimval{0};
                        //    std::memcpy(&interimval, &ival, 1);
                        //    cmc_debug_msg("Iterim val2 is: ", std::bitset<16>(interimval));
                        //}

                        //prefix_encodings.back() |= ((significant_bits[0] >> start_offset_encodings) << start_offset_encodings);
                        //prefix_encodings.push_back((significant_bits[0] << (num_considered_bits - (CHAR_BIT - start_offset_encodings))) >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings))));
                        start_offset_encodings = num_considered_bits - (CHAR_BIT - start_offset_encodings);
                        
                    }
                    num_bits_encoded += num_considered_bits;
                    ++pref_index;
                }
                /***** New TRY *********************/
                #if 0
                //Looks better but not completely correct
                if (start_offset_encodings >= CHAR_BIT)
                {
                    suffix_encodings.push_back(uint8_t{0});
                    start_offset_encodings = 0;
                }
                /* Part of previous encoding byte filled */
                const uint8_t significant_byte = significant_bits[pref_index];
                /* Add the first part right aligned */
                const uint8_t inter_pref_start_encoding = significant_byte >> start_offset_encodings;
                suffix_encodings.back() |= inter_pref_start_encoding << start_offset_encodings;
                num_bits_encoded += (CHAR_BIT - start_offset_encodings);

                int leftover_bits = start_offset_encodings;

                while (num_significant_bits - num_bits_encoded >= CHAR_BIT)
                {
                    if (leftover_bits > 0)
                    {
                        const uint8_t inter_pref_enc = significant_bits[pref_index] << (CHAR_BIT - leftover_bits);
                        suffix_encodings.push_back(inter_pref_enc);
                        num_bits_encoded += leftover_bits;
                    } else
                    {
                        suffix_encodings.push_back(uint8_t{0});
                    }

                    ++pref_index;

                    suffix_encodings.back() |= significant_bits[pref_index] >> start_offset_encodings;
                    num_bits_encoded += (CHAR_BIT - start_offset_encodings);
                }

                int enc_remaining_bits = num_significant_bits - num_bits_encoded;

                if (enc_remaining_bits == 0)
                {
                    //Nothing to do
                    start_offset_encodings = CHAR_BIT;
                } else
                {
                    if (start_offset_encodings == 0)
                    {
                        ++pref_index;
                    }
                    /* Push the remainder tail aligned */
                    const uint8_t rem_suffix = significant_bits[pref_index] << (CHAR_BIT - start_offset_encodings);
                    suffix_encodings.push_back(rem_suffix >> (CHAR_BIT - enc_remaining_bits));
                    start_offset_encodings = enc_remaining_bits;
                }

                #endif
                #if 0
                /***********************************/
                //left adjust (if not first part) and das aktuelle beschrieben Byte noch komplett ausgefüllt wird
                //if (is_first_part_written && num_significant_bits - (num_bits_encoded + num_considered_bits) >= CHAR_BIT) -> left_adjust
                bool is_first_part_written = false;
                bool adjust_to_front = false;

                while (num_bits_encoded < num_significant_bits)
                {

                    const int num_considered_bits = (num_significant_bits - num_bits_encoded > CHAR_BIT ? CHAR_BIT : num_significant_bits - num_bits_encoded);
                    /* The bit hsift is performed in order to ensure that there are zeros at the end (however, the not considered bits in the byte should already been zero) */
                    //const uint8_t significant_byte = (num_considered_bits == CHAR_BIT ? significant_bits[index] : ((significant_bits[index] >> (CHAR_BIT - num_considered_bits)) << (CHAR_BIT - num_considered_bits)));
                    //if (pref_index >= significant_bits.size()){cmc_debug_msg("Pref index too large: ", pref_index, " und size of significant bits = ", significant_bits.size(), ", current prefixes.size() = ", prefixes.size());}
                    const uint8_t significant_byte = significant_bits[pref_index];

                    if (num_considered_bits < CHAR_BIT)
                    {
                        /* We are filling the last part */
                        if (start_offset_encodings >= CHAR_BIT)
                        {
                            suffix_encodings.push_back(significant_byte >> (CHAR_BIT - num_considered_bits));
                            start_offset_encodings = num_considered_bits;
                            num_bits_encoded += num_considered_bits;
                        } else if (num_considered_bits <== CHAR_BIT - start_offset_encodings)
                        {
                            //Has to be Left aligned
                            suffix_encodings.back() |= (significant_byte >> (CHAR_BIT - num_considered_bits));
                            num_bits_encoded += num_considered_bits;
                            start_offset_encodings += num_considered_bits;
                        }  else
                        {
                            /* Needs to be split */

                            //First part needs to be left aligned
                            const uint8_t inter_pref_encoding = significant_byte >> start_offset_encodings;
                            suffix_encodings.back() |= inter_pref_encoding;

                            //Second part needs to right aligned in the next byte
                            const uint8_t inter_pref_encoding1 = significant_byte << (CHAR_BIT - start_offset_encodings);
                            suffix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings))));
                            
                            start_offset_encodings = num_considered_bits - (CHAR_BIT - start_offset_encodings);
                            num_bits_encoded += num_considered_bits;
                        }
                    } else
                    {
                        //num_considered_bits == CHAR_BIT
                        if (start_offset_encodings >= CHAR_BIT)
                        {
                            suffix_encodings.push_back(significant_byte);
                            start_offset_encodings == CHAR_BIT;
                            num_bits_encoded += CHAR_BIT;
                        } else if (start_offset_encodings == 0)
                        {
                            suffix_encodings.back() |= significant_byte;
                            start_offset_encodings == CHAR_BIT;
                            num_bits_encoded += CHAR_BIT;
                        } else
                        {
                            /* The bits are split over two bytes */
                            if (!is_first_part_written)
                            {
                                //Right align first part
                                const uint8_t inter_pref_encoding = significant_byte >> start_offset_encodings;
                                suffix_encodings.back() |= inter_pref_encoding << start_offset_encodings;
                                is_first_part_written = true;
                                num_bits_encoded += (CHAR_BIT - start_offset_encodings);
                            } else
                            {
                                //Left align other parts

                            }

                            /* Second part of significant bits */
                            const uint8_t lal_inter_pref_encoding = significant_byte << (CHAR_BIT - start_offset_encodings);
                            const uint8_t lal_pref_length = (num_considered_bits - (CHAR_BIT - start_offset_encodings));
                            //const uint8_t shift_to_fit = ((CHAR_BIT - lal_pref_length) > (num_significant_bits - num_bits_encoded - lal_pref_length) ? ((CHAR_BIT - lal_pref_length) - (num_significant_bits - num_bits_encoded - lal_pref_length)) : 0);
                            suffix_encodings.push_back(lal_inter_pref_encoding);
                            num_bits_encoded += lal_pref_length;
                            start_offset_encodings = lal_pref_length;
                        }

                    }

                    if (start_offset_encodings >= CHAR_BIT)
                    {
                        //suffix_encodings.push_back(uint8_t{0});
                        //start_offset_encodings = 0;
                        if (num_considered_bits < CHAR_BIT)
                        {
                            suffix_encodings.push_back(significant_byte >> (CHAR_BIT - num_considered_bits));
                            start_offset_encodings = num_considered_bits;
                        } else
                        {
                            suffix_encodings.push_back(significant_byte);
                            start_offset_encodings = CHAR_BIT;
                        }
                    } else if (num_considered_bits < CHAR_BIT - start_offset_encodings)
                    {
                        //Always needs to be right adjusted
                        suffix_encodings.back() |= (significant_byte >> ((CHAR_BIT - start_offset_encodings) - num_considered_bits));
                        start_offset_encodings += num_considered_bits;
                    } else if (num_considered_bits < CHAR_BIT && num_considered_bits == CHAR_BIT - start_offset_encodings)
                    {
                        /* Here we need to tail */
                    } else
                    {
                        //The problem: if the the prefix/suffx strechtes over three or more bytes, al "middle" needs to be reversly filled
                        /* Needs to be split */
                        
                        const uint8_t inter_pref_encoding = significant_byte >> start_offset_encodings;
                        suffix_encodings.back() |= inter_pref_encoding << start_offset_encodings;

                        if (iter == suffixes.begin() + 1)
                        {
                            const uint8_t ival = inter_pref_encoding << start_offset_encodings;
                            uint16_t interimval{0};
                            std::memcpy(&interimval, &ival, 1);
                            cmc_debug_msg("Iterim val1 is: ", std::bitset<16>(interimval));
                        }

                        const uint8_t inter_pref_encoding1 = significant_byte << (CHAR_BIT - start_offset_encodings);
                        suffix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings))));
                        
                        if (iter == suffixes.begin() + 1)
                        {
                            const uint8_t ival = inter_pref_encoding1 >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings)));
                            uint16_t interimval{0};
                            std::memcpy(&interimval, &ival, 1);
                            cmc_debug_msg("Iterim val2 is: ", std::bitset<16>(interimval));
                        }

                        start_offset_encodings = num_considered_bits - (CHAR_BIT - start_offset_encodings);
                    }
                    num_bits_encoded += num_considered_bits;
                    ++pref_index;
                }


                #if 0
                //Hier gilt sowieso num_significant_bits > CHAR_BIT
                bool adjust_left = false;
                bool first_time_tail_adjusted = true;
                while (num_bits_encoded < num_significant_bits)
                {
                    const int num_considered_bits = (num_significant_bits - num_bits_encoded > CHAR_BIT ? CHAR_BIT : num_significant_bits - num_bits_encoded);
                    const uint8_t significant_byte = significant_bits[pref_index];

                    if (start_offset_encodings >= CHAR_BIT)
                    {
                        suffix_encodings.push_back(uint8_t{0});
                        start_offset_encodings = 0;
                    }

                    if (num_considered_bits == CHAR_BIT)
                    {
                        if (start_offset_encodings != 0)
                        {
                            /* The significant byte is split over two bytes */
                            if (first_time_tail_adjusted)
                            {
                                //Tail adjust first part of the byte  
                                const uint8_t inter_pref_encoding = significant_byte >> start_offset_encodings;
                                suffix_encodings.back() |= inter_pref_encoding << start_offset_encodings;

                                start_offset_encodings += num_considered_bits;  
                                first_time_tail_adjusted = false;
                            } else
                            {
                                //Front adjust significant byte

                            }
                        }
                    }



                    if (num_considered_bits <= CHAR_BIT - start_offset_encodings)
                    {
                        if (first_time_tail_adjusted || num_significant_bits - num_bits_encoded - num_considered_bits < CHAR_BIT)
                        {
                            suffix_encodings.back() |= (significant_byte >> ((CHAR_BIT - start_offset_encodings) - num_considered_bits));
                            start_offset_encodings += num_considered_bits;
                            first_time_tail_adjusted = false;
                        }
                    } else
                    {}

                    //First Time Right adjust
                    if (first_time_tail_adjusted)
                    {
                        const uint8_t inter_pref_encoding = significant_byte >> start_offset_encodings;
                        suffix_encodings.back() |= inter_pref_encoding << start_offset_encodings;
                        num_bits_encoded += CHAR_BIT - start_offset_encodings;
                        first_time_tail_adjusted = false;
                    } else
                    {

                    }
                    //If remaining bits larger equal than 8, left assign
                    if (num_significant_bits - num_bits_encoded >= 8)
                    {
                        //Left assign

                    }
                    //Danach left assign

                    if (num_considered_bits <= CHAR_BIT - start_offset_encodings)
                    {
                        if (num_bits_encoded + num_considered_bits >= num_significant_bits)
                        {
                            //Adjust to tail, since we are finished

                        } else
                        {
                            //Adjust to left
                        }
                        suffix_encodings.back() |= (significant_byte >> ((CHAR_BIT - start_offset_encodings) - num_considered_bits));
                        start_offset_encodings += num_considered_bits;
                    }
                }
                #endif
                #if 0

                while (num_bits_encoded < num_significant_bits)
                {
                    const int num_considered_bits = (num_significant_bits - num_bits_encoded > CHAR_BIT ? CHAR_BIT : num_significant_bits - num_bits_encoded);
                    /* The bit hsift is performed in order to ensure that there are zeros at the end (however, the not considered bits in the byte should already been zero) */
                    //const uint8_t significant_byte = (num_considered_bits == CHAR_BIT ? significant_bits[index] : ((significant_bits[index] >> (CHAR_BIT - num_considered_bits)) << (CHAR_BIT - num_considered_bits)));
                    //if (pref_index >= significant_bits.size()){cmc_debug_msg("Pref index too large: ", pref_index, " und size of significant bits = ", significant_bits.size(), ", current prefixes.size() = ", prefixes.size());}
                    const uint8_t significant_byte = significant_bits[pref_index];

                    if (start_offset_encodings >= CHAR_BIT)
                    {
                        //suffix_encodings.push_back(uint8_t{0});
                        //start_offset_encodings = 0;
                        suffix_encodings.push_back(significant_byte >> (CHAR_BIT - num_considered_bits));
                        start_offset_encodings = num_considered_bits;
                    } else if (num_considered_bits <= CHAR_BIT - start_offset_encodings)
                    {
                        suffix_encodings.back() |= (significant_byte >> ((CHAR_BIT - start_offset_encodings) - num_considered_bits));
                        start_offset_encodings += num_considered_bits;
                    } else
                    {
                        //The problem: if the the prefix/suffx strechtes over three or more bytes, al "middle" needs to be reversly filled
                        /* Needs to be split */
                        
                        const uint8_t inter_pref_encoding = significant_byte >> start_offset_encodings;
                        suffix_encodings.back() |= inter_pref_encoding << start_offset_encodings;
                        if (iter == suffixes.begin() + 1)
                        {
                            const uint8_t ival = inter_pref_encoding << start_offset_encodings;
                            uint16_t interimval{0};
                            std::memcpy(&interimval, &ival, 1);
                            cmc_debug_msg("Iterim val1 is: ", std::bitset<16>(interimval));
                        }
                        const uint8_t inter_pref_encoding1 = significant_byte << (CHAR_BIT - start_offset_encodings);
                        suffix_encodings.push_back(inter_pref_encoding1 >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings))));
                        if (iter == suffixes.begin() + 1)
                        {
                            const uint8_t ival = inter_pref_encoding1 >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings)));
                            uint16_t interimval{0};
                            std::memcpy(&interimval, &ival, 1);
                            cmc_debug_msg("Iterim val2 is: ", std::bitset<16>(interimval));
                        }

                        //prefix_encodings.back() |= ((significant_bits[0] >> start_offset_encodings) << start_offset_encodings);
                        //prefix_encodings.push_back((significant_bits[0] << (num_considered_bits - (CHAR_BIT - start_offset_encodings))) >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings))));
                        start_offset_encodings = num_considered_bits - (CHAR_BIT - start_offset_encodings);
                        
                        #if 0
                        const uint8_t inter_pref_encoding = significant_byte >> start_offset_encodings;
                        if (was_left_adjusted)
                        {
                            suffix_encodings.back() |= inter_pref_encoding;
                        } else
                        {
                            suffix_encodings.back() |= inter_pref_encoding << start_offset_encodings;
                        }
                        num_bits_encoded += CHAR_BIT - start_offset_encodings;

                        const uint8_t inter_pref_encoding1 = significant_byte << (CHAR_BIT - start_offset_encodings);
                        const uint8_t tail_adjusted_pref_encoding = inter_pref_encoding1 >> (CHAR_BIT - (num_considered_bits - (CHAR_BIT - start_offset_encodings)));
                        if (num_bits_encoded + (num_considered_bits - (CHAR_BIT - start_offset_encodings)) < num_considered_bits)
                        {
                            /* Needs to added left adjusted */
                            suffix_encodings.push_back(tail_adjusted_pref_encoding << (CHAR_BIT - start_offset_encodings));
                            was_left_adjusted = true;
                        } else
                        {
                            /* It can be added right adjusted */
                            suffix_encodings.push_back(tail_adjusted_pref_encoding);
                            was_left_adjusted = false;
                        }

                        start_offset_encodings = num_considered_bits - (CHAR_BIT - start_offset_encodings);
                        #endif
                    }
                    num_bits_encoded += num_considered_bits;
                    ++pref_index;
                }
                #endif

                #endif
            }

        }
    }
    
    std::copy_n(suffix_lengths.begin(), suffix_lengths.size(), std::back_inserter(serialized_variable));
    std::copy_n(suffix_encodings.begin(), suffix_encodings.size(), std::back_inserter(serialized_variable));

    return std::make_tuple(rle.size(), suffix_lengths.size(), suffix_encodings.size());
}



constexpr uint8_t kUnsetFrontNibble = 0x0F;
constexpr uint8_t kUnsetTailNibble = 0xF0;
constexpr int kNibbleSize = 4;


constexpr uint8_t kPrefixContinuationIndicator = 0x0F;
constexpr uint8_t kFourBitLengthEncodingOffset = 0x03;
constexpr uint8_t kTwoBitLengthEncodingOffset = 0x01;
constexpr uint8_t kNullifyFourthBit = ~(0x08);


//Currently, new version
inline
std::tuple<std::vector<uint8_t>, int>
GetNextPrefix(const VectorView<uint8_t>& prefix_lengths, int& length_byte_position, int& length_bit_position,
              const VectorView<uint8_t>& prefix_encodings, int& encoding_byte_position, int& encoding_bit_position)
{
    if (length_bit_position <= 0)
    {
        ++length_byte_position;
        length_bit_position = CHAR_BIT;
    }

    const uint8_t current_length_byte = prefix_lengths[length_byte_position];

    uint8_t pref_length{0};
    bool was_prefix_continuation_before = false;

    //Bit shoft of eight is applied which obviously leads zero; CHar_BIt -1 had to be cheked
    if (CheckIfBitIsSet(current_length_byte,  length_bit_position - 1))
    {
        /* A four bit length code is present */
        bool pref_length_is_not_found = true;

        while (pref_length_is_not_found)
        {
            if (length_bit_position <= 0)
            {
                ++length_byte_position;
                length_bit_position = CHAR_BIT;
            }

            if (length_bit_position - 4 >= 0)
            {
                /* The whole length is encoded within this byte */
                const uint8_t inter_pref_encoding = (prefix_lengths[length_byte_position] << (CHAR_BIT - length_bit_position));

                const uint8_t pref_length_encoding = inter_pref_encoding >> 4;

                /* Check if there is a prefix length continuation indicator */
                if (pref_length_encoding != kPrefixContinuationIndicator)
                {
                    /* The encoded prefix needs to be transformed for retrieving the actual length.
                     * The fourth bit needs (indicating the bit size) needs to be cleared and the offset needs to be added */
                    if (!was_prefix_continuation_before)
                    {
                        pref_length += (pref_length_encoding & kNullifyFourthBit) + kFourBitLengthEncodingOffset;
                    } else
                    {
                        /* If a continuation flag has been given, the following prefix_length is not offseted and directly corresponds to the remaining length */
                        pref_length += pref_length_encoding;
                    }
                    pref_length_is_not_found = false;
                    was_prefix_continuation_before = false;
                } else
                {
                    was_prefix_continuation_before = true;
                    pref_length += CHAR_BIT;
                }
                
                length_bit_position -= 4;
            } else
            {
                /* The length code is split between bytes */

                const uint8_t inter_pref_encoding1 = (prefix_lengths[length_byte_position] << (CHAR_BIT - length_bit_position));
                const uint8_t inter_pref_encoding2 = inter_pref_encoding1 >> 4;
                const uint8_t inter_pref_encoding3 = prefix_lengths[length_byte_position + 1] >> (CHAR_BIT - (4 - length_bit_position));
                const uint8_t pref_length_encoding = inter_pref_encoding2 | inter_pref_encoding3;

                /* Check if there is a prefix length continuation indicator */
                if (pref_length_encoding != kPrefixContinuationIndicator)
                {
                    /* The encoded prefix needs to be transformed for retrieving the actual length.
                     * The fourth bit needs (indicating the bit size) needs to be cleared and the offset needs to be added */
                    if (!was_prefix_continuation_before)
                    {
                        pref_length += (pref_length_encoding & kNullifyFourthBit) + kFourBitLengthEncodingOffset;
                    } else
                    {
                        /* If a continuation flag has bee given, the following prefix_length is not offseted and directly corresponds to the remaining length */
                        pref_length += pref_length_encoding;
                    }
                    pref_length_is_not_found = false;
                    was_prefix_continuation_before = false;
                } else
                {
                    was_prefix_continuation_before = true;
                    pref_length += CHAR_BIT;
                }

                ++length_byte_position;
                length_bit_position = CHAR_BIT - (4 - length_bit_position);
            }
        }
    }  else
    {
        /* A two bit length code is present */
        const uint8_t inter_pref_encoding = current_length_byte << (CHAR_BIT - length_bit_position);
        const uint8_t pref_length_encoding = inter_pref_encoding >> 6;
        pref_length = pref_length_encoding + kTwoBitLengthEncodingOffset;
        length_bit_position -= 2;
    }

    /* The length of the next prefix */
    const int prefix_length = static_cast<int>(pref_length);

    uint32_t following_encoding_bytes{0};
    std::memcpy(&following_encoding_bytes, prefix_encodings.begin() + encoding_byte_position, 4);

    const int num_bytes = static_cast<int>(prefix_length / CHAR_BIT) + (prefix_length % CHAR_BIT != 0 ? 1 : 0);

    std::vector<uint8_t> prefix(num_bytes, uint8_t{0});
    int prefix_tail_position = CHAR_BIT;
    int prefix_byte_posititon = 0;

    int num_bits_extracted = 0;

    /* Extract all prefix bits */
    while (num_bits_extracted < prefix_length)
    {
        if (encoding_bit_position >= CHAR_BIT)
        {
            ++encoding_byte_position;
            encoding_bit_position = 0;
        }

        /* Check if the remaining part of the prefix lays is the current byte */
        const int current_pref_length = std::min(prefix_length - num_bits_extracted, CHAR_BIT - encoding_bit_position);
        
        const uint8_t inter_extracted_bit_sequence = (prefix_encodings[encoding_byte_position] >> encoding_bit_position);
        const uint8_t extracted_bit_sequence = inter_extracted_bit_sequence << (CHAR_BIT - current_pref_length);

        if (prefix_tail_position <= 0)
        {
            prefix_tail_position = CHAR_BIT;
            ++prefix_byte_posititon;
        }

        /* Check whether the prefix length fits in the current byte of the extracted prefix */
        if (current_pref_length <= prefix_tail_position)
        {
            prefix[prefix_byte_posititon] |= (extracted_bit_sequence >> (CHAR_BIT - prefix_tail_position));
            prefix_tail_position -= current_pref_length;
        } else
        {
            prefix[prefix_byte_posititon] |= (extracted_bit_sequence >> (CHAR_BIT - prefix_tail_position));
            ++prefix_byte_posititon;
            prefix[prefix_byte_posititon] |= (extracted_bit_sequence << prefix_tail_position);

            prefix_tail_position = CHAR_BIT - (current_pref_length - prefix_tail_position);
        }

        num_bits_extracted += current_pref_length;
        encoding_bit_position += current_pref_length;
    }

    return std::make_pair(prefix, prefix_length);
}

const uint8_t kNullifyAllExceptFourthTailBit = 0x08;
const uint8_t kNullifyExceptHighBit = 0x80;
const uint8_t kNullifyAllExceptThreeTailBits = 0x07;
const uint8_t kNullifyAllExceptFiveToSevenTailBits = 0x70;
const uint8_t kAllOnes = 0xFF;
const uint8_t kAllZeros = 0x00;
const uint8_t kRLEZero = 0x00;
const uint8_t kRLEOne = 0x08;

//TODO put this somewhere generally
//const uint8_t CharBit = CHAR_BIT;

/* Currently from low to high bit (This should change soon) */
inline static
void
InsertBitMultipleTimes(const bool bit_is_set, const uint8_t number_of_copies, std::vector<uint8_t>& vec_to_insert, uint8_t& current_bit_position)
{
    if (current_bit_position >= CHAR_BIT)
    {
        vec_to_insert.push_back(uint8_t{0});
        current_bit_position = 0;
    }

    if (bit_is_set)
    {
        /* Insert a one multiple times */
        if (number_of_copies <= CHAR_BIT - current_bit_position)
        {
            /* Nullify the tail bits up to the bit starting position */
            const uint8_t bit_insertion = kAllOnes << current_bit_position;
            
            /* Nullify the front bits */
            const uint8_t bit_nullifier = bit_insertion >> ((CHAR_BIT - current_bit_position) - number_of_copies);

            /* Set the ones in the correct place */
            vec_to_insert.back() |= (bit_insertion & bit_nullifier);

            current_bit_position += number_of_copies;
        } else
        {
            /* Nullify the tail bits up to the bit starting position */
            const uint8_t bit_insertion = kAllOnes << current_bit_position;
            vec_to_insert.back() |= bit_insertion;
            vec_to_insert.push_back(kAllOnes >> (CHAR_BIT - (number_of_copies - (CHAR_BIT - current_bit_position))));
            current_bit_position = (number_of_copies - (CHAR_BIT - current_bit_position));
        }
    } else
    {
        /* Insert a zero multiple times */
        if (number_of_copies <= CHAR_BIT - current_bit_position)
        {
            current_bit_position += number_of_copies;
        } else
        {
            vec_to_insert.push_back(uint8_t{0});
            current_bit_position = (number_of_copies - (CHAR_BIT - current_bit_position));
        }
    }
}

inline
std::pair<std::vector<uint8_t>, int>
DecodeRunLengthEncoding(const VectorView<uint8_t>& rle)
{
    cmc_debug_msg("Start decode rle");
    uint8_t bit_position = CHAR_BIT;
    std::vector<uint8_t> decoded_rle;
    decoded_rle.reserve(rle.size() * 2);

    int num_bits = 0;

    cmc_debug_msg("Rle size at the beginning: ", rle.size());

    for (size_t index = 0; index < rle.size(); ++index)
    {
        /* Decode the back of the byte */

        /* Extract the number of copies from the first part of the byte */
        //Encoding: ((uint8_t)0x08 & (current_run_length - 1))
        const uint8_t num_copies_tail_encoding = (rle[index] & kNullifyAllExceptThreeTailBits) + uint8_t{1};

        //cmc_debug_msg("Tail encoding unum copies: ", (unsigned) num_copies_tail_encoding);
        if ((rle[index] & kNullifyAllExceptFourthTailBit) == kRLEZero)
        {
            /* There is a sequence of zeros given */
            InsertBitMultipleTimes(false, num_copies_tail_encoding, decoded_rle, bit_position);
        } else
        {
            /* There is a sequence of ones given */
            InsertBitMultipleTimes(true, num_copies_tail_encoding, decoded_rle, bit_position);
        }
        num_bits += num_copies_tail_encoding;



        /* Decode the front of the byte */

        /* Extract the number of copies from the first part of the byte */
        const uint8_t num_copies_front_encoding = ((rle[index] >> uint8_t{4}) & kNullifyAllExceptThreeTailBits) + uint8_t{1};

        if (rle[index] & kNullifyExceptHighBit == kRLEZero)
        {
            /* There is a sequence of zeros given */
            InsertBitMultipleTimes(false, num_copies_front_encoding, decoded_rle, bit_position);
        } else
        {
            /* There is a sequence of ones given */
            InsertBitMultipleTimes(true, num_copies_front_encoding, decoded_rle, bit_position);
        }
        num_bits += num_copies_front_encoding;
    }
    cmc_debug_msg("num bits directly after decoding: ", num_bits);
    return std::make_pair(decoded_rle, num_bits);
}

inline int
GetNumberBytesNeededForEncoding(const int num_significant_bits)
{
    cmc_assert(num_significant_bits >= 0);
    if (num_significant_bits == 0) {return 0;}

    const int num_significant_bits_ = num_significant_bits - 1;
    return (2 * (num_significant_bits_ / 12) + ((num_significant_bits_ % 12) >= 4 ? 2 : 1));
}

}

#endif /* !CMC_T8_PREFIX_ENCODING_HXX */
