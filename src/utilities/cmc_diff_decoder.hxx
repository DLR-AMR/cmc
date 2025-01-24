#ifndef CMC_DIFF_DECODER_HXX
#define CMC_DIFF_DECODER_HXX

#include "utilities/cmc_ac_model.hxx"
#include "utilities/cmc_arithmetic_encoder.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_utilities.hxx"
#include "utilities/cmc_prefix_encoding.hxx"

#include <vector>
#include <memory>

namespace cmc
{

namespace diff
{

template <typename Iter>
CmcUniversalType
GetSerializedValueAtPosition(Iter pos, const CmcType type);

class DiffDecoder
{
public:
    DiffDecoder() = delete;

    DiffDecoder(const std::vector<uint8_t>& serialized_data, const CmcType type)
    : serialized_data_{serialized_data}, type_{type} {
        cmc_debug_msg("Start decoding: Offset for decoding: ", byte_offset_);
        /* Initialize the the decoder by correctly setting up the correct view to the different encoded streams */
        num_bytes_ = GetValueFromByteStream<size_t>(serialized_data_.data());
        byte_offset_ += sizeof(size_t);
        cmc_debug_msg("after total byte count: Offset for decoding: ", byte_offset_);
        /* Get the frequency model from the encoded stream */
        auto [frequency_model, num_bytes_encoding] = arithmetic_encoding::DecodeStaticFrequencyAlphabet(serialized_data_.data() + byte_offset_); 
        //frequency_model_ = std::make_unique<cmc::arithmetic_encoding::StaticFrequencyModel>(frequency_model);        
        byte_offset_ += num_bytes_encoding;
        cmc_debug_msg("After decoding alphabet: Offset for decoding: ", byte_offset_);

        cmc_debug_msg("Before arithmetic decoder is initialized");
        /* Create the arithmetic decoder */
        arm_decoder_ = std::make_unique<arithmetic_encoding::Decoder>(std::make_unique<cmc::arithmetic_encoding::StaticFrequencyModel>(frequency_model));
        cmc_debug_msg("after arithmetic decoder init");

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
        arm_decoder_->Reset(bit_map::BitMapView(serialized_data_.data() + byte_offset_, bytes_lzc * bit_map::kCharBit));
        //encoded_lzc_ = bit_map::BitMapView(serialized_data_.data() + byte_offset_, bytes_lzc * bit_map::kCharBit);
        byte_offset_ += bytes_lzc;
        cmc_debug_msg("encoded lzc view has been created");
        cmc_debug_msg("After LZC view: Offset for decoding: ", byte_offset_);
        /* Create a view on the encoded actual remaining residual bits */
        encoded_residuals_ = bit_vector::BitVectorView(serialized_data_.data() + byte_offset_, bytes_residual);
        byte_offset_ += bytes_residual;
        cmc_debug_msg("encoded residual view has been created");
        cmc_debug_msg("After esidual view: Offset for decoding: ", byte_offset_);

        /* Update the decoder for the new encoded LZC stream */
        //arm_decoder_->Reset(encoded_lzc_);
        cmc_debug_msg("arm decoder has been resetted");
    }

    inline uint32_t GetNextEncodedResidualLength()
    {
        return arm_decoder_->DecodeNextSymbol();
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
    //bit_map::BitMapView encoded_lzc_;
    bit_vector::BitVectorView encoded_residuals_;
    //std::unique_ptr<arithmetic_encoding::IACModel> frequency_model_;
    std::unique_ptr<arithmetic_encoding::Decoder> arm_decoder_;
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


#endif
