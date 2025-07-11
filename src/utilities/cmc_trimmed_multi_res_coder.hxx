#ifndef CMC_TRIMMED_MULTI_RES_ENTROPY_CODER_HXX
#define CMC_TRIMMED_MULTI_RES_ENTROPY_CODER_HXX

#include "utilities/cmc_trimmed_multi_res_compression_entropy_alphabet.hxx"
#include "utilities/cmc_byte_compression_arithmetic_encoding_frequency_model.hxx"
#include "utilities/cmc_byte_compression_arithmetic_encoding.hxx"
#include "utilities/cmc_bit_map.hxx"

#include <memory>

namespace cmc::entropy_coding::arithmetic_coding
{

template <typename T>
class TrimmedMultiResEncoder : public Encoder<T>
{
public:
    MultiResEncoder()
    : Encoder<T>(std::make_unique<MultiResCompressionAlphabet<T>>()) {};
    MultiResEncoder(std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet)
    : Encoder<T>(std::move(alphabet)) {};
};

template <typename T>
class TrimmedMultiResDecoder : public Decoder
{
public:
    template <typename Iter> MultiResDecoder(Iter alphabet_begin)
    : Decoder(DecodeByteCompressionStaticFrequencyAlphabet<T>(alphabet_begin, std::make_unique<TrimmedMultiResCompressionAlphabet<T>>())) {};

    template <typename Iter> MultiResDecoder(Iter alphabet_begin, bit_map::BitMapView encoded_stream)
    : Decoder(DecodeByteCompressionStaticFrequencyAlphabet<T>(alphabet_begin, std::make_unique<TrimmedMultiResCompressionAlphabet<T>>()), encoded_stream) {};

    template <typename Iter> MultiResDecoder(std::unique_ptr<IByteCompressionEntropyAlphabet> alphabet, Iter alphabet_begin, bit_map::BitMapView encoded_stream)
    : Decoder(DecodeByteCompressionStaticFrequencyAlphabet<T>(alphabet_begin, std::move(alphabet)), encoded_stream) {};
};

}


#endif /* !CMC_TRIMMED_MULTI_RES_ENTROPY_CODER_HXX */
