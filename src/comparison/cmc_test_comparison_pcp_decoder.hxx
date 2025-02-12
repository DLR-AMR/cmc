#ifndef CMC_TEST_COMPARISON_PCP_DECODER_HXX
#define CMC_TEST_COMPARISON_PCP_DECODER_HXX


#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_prefix_encoding.hxx"

#include <vector>
#include <memory>

namespace cmc
{

namespace test_comparison
{

namespace light_amr_pcp
{

template <typename Iter>
CmcUniversalType
GetSerializedValueAtPosition(Iter pos, const CmcType type);

class _TestPCPDecoder
{
public:
    _TestPCPDecoder() = delete;

    _TestPCPDecoder(const std::vector<uint8_t>& serialized_data, const CmcType type)
    : serialized_data_{serialized_data}, type_{type} {
        cmc_debug_msg("Start decoding: Offset for decoding: ", byte_offset_);
        /* Initialize the the decoder by correctly setting up the correct view to the different encoded streams */
        num_bytes_ = GetValueFromByteStream<size_t>(serialized_data_.data());
        byte_offset_ += sizeof(size_t);
        cmc_debug_msg("after total byte count: Offset for decoding: ", byte_offset_);
 
        /* Get the root element's value that has been encoded in the encoded stream */
        root_elem_value_ = GetSerializedValueAtPosition(serialized_data_.data() + byte_offset_, type_);
        byte_offset_ += CmcTypeToBytes(type_);
        cmc_debug_msg("After root elem value: Offset for decoding: ", byte_offset_);
    };

    CmcUniversalType GetRootElementValue() const {return root_elem_value_;}

    void MoveToNextRefinementLevel()
    {
        cmc_debug_msg("move to next ref lvlv start");
        /* Get the number of bytes for the leading zero counts */
        const size_t bytes_lzc = GetValueFromByteStream<size_t>(serialized_data_.data() + byte_offset_);
        byte_offset_ += sizeof(size_t);
        cmc_debug_msg("bytes lzc: ", bytes_lzc);
        cmc_debug_msg("After LZC count on level: Offset for decoding: ", byte_offset_);

        /* Get the number of bytes for the actual encoded residuals */
        const size_t bytes_residual = GetValueFromByteStream<size_t>(serialized_data_.data() + byte_offset_);
        byte_offset_ += sizeof(size_t);
        cmc_debug_msg("bytes residual: ", bytes_residual);
        cmc_debug_msg("After residual encoding count: Offset for decoding: ", byte_offset_);


        /* Create view on the encoded LZC */
        encoded_lzc_ = bit_vector::BitVectorView(serialized_data_.data() + byte_offset_, bytes_lzc);
        byte_offset_ += bytes_lzc;
        cmc_debug_msg("encoded lzc view has been created");
        cmc_debug_msg("After LZC view: Offset for decoding: ", byte_offset_);

        /* Create a view on the encoded actual remaining residual bits */
        encoded_residuals_ = bit_vector::BitVectorView(serialized_data_.data() + byte_offset_, bytes_residual);
        byte_offset_ += bytes_residual;
        cmc_debug_msg("encoded residual view has been created");
        cmc_debug_msg("After esidual view: Offset for decoding: ", byte_offset_);
    }

    inline size_t GetNextEncodedLZC()
    {
        const std::vector<uint8_t> encoded_lzc = encoded_lzc_.GetNextBitSequence(4);
        const size_t decoded_lzc = static_cast<size_t>(encoded_lzc.front() >> 4);
        return decoded_lzc;
    }

    inline std::vector<uint8_t> GetNextResidualBitSequence(const size_t num_bits)
    {
        return encoded_residuals_.GetNextBitSequence(num_bits);
    }

private:
    const std::vector<uint8_t>& serialized_data_;
    CmcUniversalType root_elem_value_;
    CmcType type_;
    size_t num_bytes_{0};
    size_t byte_offset_{0};
    bit_vector::BitVectorView encoded_lzc_;
    bit_vector::BitVectorView encoded_residuals_;
};


template <typename Iter>
CmcUniversalType
GetSerializedValueAtPosition(Iter pos, const CmcType type)
{
    switch (type)
    {
        case CmcType::Int8_t:
        {
            const int8_t value = GetValueFromByteStream<int8_t>(pos);
            return CmcUniversalType{value};
        }
        break;
        case CmcType::Char:
        {
            const char value = GetValueFromByteStream<char>(pos);
            return CmcUniversalType{value};
        }
        break;
        case CmcType::Int16_t:
        {
            const int16_t value = GetValueFromByteStream<int16_t>(pos);
            return CmcUniversalType{value};
        }
        break;
        case CmcType::Int32_t:
        {
            const int32_t value = GetValueFromByteStream<int32_t>(pos);
            return CmcUniversalType{value};
        }
        break;
        case CmcType::Float:
        {
            const float value = GetValueFromByteStream<float>(pos);
            return CmcUniversalType{value};
        }
        break;
        case CmcType::Double:
        {
            const double value = GetValueFromByteStream<double>(pos);
            return CmcUniversalType{value};
        }
        break;
        case CmcType::Uint8_t:
        {
            const uint8_t value = GetValueFromByteStream<uint8_t>(pos);
            return CmcUniversalType{value};
        }
        break;
        case CmcType::Uint16_t:
        {
            const uint16_t value = GetValueFromByteStream<uint16_t>(pos);
            return CmcUniversalType{value};
        }
        break;
        case CmcType::Uint32_t:
        {
            const uint32_t value = GetValueFromByteStream<uint32_t>(pos);
            return CmcUniversalType{value};
        }
        break;
        case CmcType::Int64_t:
        {
            const int64_t value = GetValueFromByteStream<int64_t>(pos);
            return CmcUniversalType{value};
        }
        break;
        case CmcType::Uint64_t:
        {
            const uint64_t value = GetValueFromByteStream<uint64_t>(pos);
            return CmcUniversalType{value};
        }
        break;
        default:
            cmc_err_msg("An unknown data type has been supplied.");
    }

}


}

}

}

#endif /* !CMC_TEST_COMPARISON_PCP_DECODER_HXX */
