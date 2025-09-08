#ifndef CMC_RESIDUAL_QUANTIZATION_ENTROPY_CODER_HXX
#define CMC_RESIDUAL_QUANTIZATION_ENTROPY_CODER_HXX

#include "utilities/cmc_residual_quantization_entropy_alphabet.hxx"
#include "utilities/cmc_byte_compression_arithmetic_encoding_frequency_model.hxx"
#include "utilities/cmc_byte_compression_arithmetic_encoding.hxx"
#include "utilities/cmc_bit_map.hxx"

#include <memory>

namespace cmc::entropy_coding::arithmetic_coding
{

template <typename T>
class ResidualQuantizationEncoder : public Encoder<T>
{
public:
ResidualQuantizationEncoder()
    : Encoder<T>(std::make_unique<ResidualQuantizationAlphabet<T>>()) {};
    ResidualQuantizationEncoder(std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet)
    : Encoder<T>(std::move(alphabet)) {};
};

template <typename T>
class ResidualQuantizationDecoder : public Decoder
{
public:
    template <typename Iter> ResidualQuantizationDecoder(Iter alphabet_begin)
    : Decoder(DecodeByteCompressionStaticFrequencyAlphabet<T>(alphabet_begin, std::make_unique<ResidualQuantizationAlphabet<T>>())) {};

    template <typename Iter> ResidualQuantizationDecoder(Iter alphabet_begin, bit_map::BitMapView encoded_stream)
    : Decoder(DecodeByteCompressionStaticFrequencyAlphabet<T>(alphabet_begin, std::make_unique<ResidualQuantizationAlphabet<T>>()), encoded_stream) {};

    template <typename Iter> ResidualQuantizationDecoder(std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet, Iter alphabet_begin, bit_map::BitMapView encoded_stream)
    : Decoder(DecodeByteCompressionStaticFrequencyAlphabet<T>(alphabet_begin, std::move(alphabet)), encoded_stream) {};
};

}


#endif /* !CMC_RESIDUAL_QUANTIZATION_ENTROPY_CODER_HXX */
