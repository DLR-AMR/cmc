#ifndef CMC_PREFIX_ENTROPY_CODER_HXX
#define CMC_PREFIX_ENTROPY_CODER_HXX

#include "utilities/cmc_prefix_compression_entropy_alphabet.hxx"
#include "utilities/cmc_byte_compression_arithmetic_encoding_frequency_model.hxx"
#include "utilities/cmc_byte_compression_arithmetic_encoding.hxx"
#include "utilities/cmc_bit_map.hxx"

#include <memory>

namespace cmc::entropy_coding::arithmetic_coding
{

template <typename T>
class PrefixEncoder : public Encoder<T>
{
public:
    PrefixEncoder()
    : Encoder<T>(std::make_unique<PrefixCompressionAlphabet<T>>()) {};
};

template <typename T>
class PrefixDecoder : public Decoder
{
public:
    template <typename Iter> PrefixDecoder(Iter alphabet_begin)
    : Decoder(DecodeByteCompressionStaticFrequencyAlphabet<T>(alphabet_begin, std::make_unique<PrefixCompressionAlphabet<T>>())) {};

    template <typename Iter> PrefixDecoder(Iter alphabet_begin, bit_map::BitMapView encoded_stream)
    : Decoder(DecodeByteCompressionStaticFrequencyAlphabet<T>(alphabet_begin, std::make_unique<PrefixCompressionAlphabet<T>>()), encoded_stream) {};
};

}


#endif /* !CMC_PREFIX_ENTROPY_CODER_HXX */
