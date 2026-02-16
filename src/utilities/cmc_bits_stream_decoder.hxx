#ifndef CMC_BITS_STREAM_DECODER_HXX
#define CMC_BITS_STREAM_DECODER_HXX

#include "cmc.hxx"
#include "utilities/cmc_bits_vector_view.hxx"
#include "utilities/cmc_huffman_coder.hxx"

namespace cmc::bits
{

template<typename T>
class StreamDecoder 
{
public:
    StreamDecoder() = default;
    void StartHuffmanCodesDecoding(const uint8_t* start_huffman_codes_serialization);
    void StartDecoding(cmc::bits::vector_view encoding);
    T DecodeNextEntropySymbol();
    void SkipNextBits(const int num_bits_to_skip);
    template<UnsignedIntegerType U> U GetNextBitSequence(const int num_bits);
    bool GetNextBit();

    uint64_t GetNumberOfProcessedBytesForSymbolCodewordTable() const;

private:
    
    uint64_t num_processed_bytes_symbol_codes_{0};
    cmc::bits::vector_view encoded_stream_view_;
    cmc::entropy_coding::huffman::HuffmanDecodeMap<T> codes_;
};

template <typename T>
inline void
StreamDecoder<T>::StartDecoding(cmc::bits::vector_view encoding)
{
    encoded_stream_view_ = encoding;
}

template <typename T>
inline void
StreamDecoder<T>::StartHuffmanCodesDecoding(const uint8_t* start_huffman_codes_serialization)
{
    size_t offset{0};

    const cmc::entropy_coding::huffman::HuffmanCodeInfoType num_symbols = DeserializeValueBE<cmc::entropy_coding::huffman::HuffmanCodeInfoType>(start_huffman_codes_serialization);
    offset += sizeof(cmc::entropy_coding::huffman::HuffmanCodeInfoType);

    const cmc::entropy_coding::huffman::HuffmanCodeInfoType data_type = DeserializeValueBE<cmc::entropy_coding::huffman::HuffmanCodeInfoType>(start_huffman_codes_serialization + offset);
    offset += sizeof(cmc::entropy_coding::huffman::HuffmanCodeInfoType);

    if (static_cast<cmc::entropy_coding::huffman::HuffmanCodeInfoType>(ConvertToCmcType<T>()) != data_type) [[unlikely]]
    {
        cmc_err_msg("The template parameter does not coincide with the symbol type of the Huffman symbol frequency table.");
    }

    codes_.reserve(num_symbols);

    /* Iterate until the symbol frequency table has been re-created */
    for (cmc::entropy_coding::huffman::HuffmanCodeInfoType iter{0}; iter < num_symbols; ++iter)
    {
        /* De-Serialize the code word */
        const cmc::entropy_coding::huffman::HuffmanCodeWord deserialized_code_word = DeserializeValueBE<cmc::entropy_coding::huffman::HuffmanCodeWord>(start_huffman_codes_serialization + offset);
        offset += sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord);
    
        /* De-Serialize the symbol */
        const T deserialized_symbol = DeserializeValueBE<T>(start_huffman_codes_serialization + offset);
        offset += sizeof(T);
        
        /* Store the deserialized symbol with the code word */
        codes_[deserialized_code_word] = deserialized_symbol;
    }

    /* Store the processed bytes */
    num_processed_bytes_symbol_codes_ = offset;
}

template <typename T>
inline void
StreamDecoder<T>::SkipNextBits(const int num_bits_to_skip)
{
    encoded_stream_view_.SkipNumberOfBits(num_bits_to_skip);
}

template<typename T>
template<UnsignedIntegerType U>
inline U
StreamDecoder<T>::GetNextBitSequence(const int num_bits)
{
    return encoded_stream_view_.GetNextBitSequence<U>(num_bits);
}

template <typename T>
inline bool
StreamDecoder<T>::GetNextBit()
{
    return encoded_stream_view_.GetNextBit();
}

inline void 
AppendBit(uint64_t& code, const bool bit)
{
    code <<= 1;
    code |= uint64_t(bit);
}

template <typename T>
inline T
StreamDecoder<T>::DecodeNextEntropySymbol()
{
    /* Define the start codeword */
    uint64_t code = cmc::entropy_coding::huffman::CreateStartEncodedHuffmanCode(encoded_stream_view_.GetNextBit());

    /* Iterate until we found a corresponding leaf node from the codeword */
    while (codes_.find(code) == codes_.end())
    {
        /* Get the next bit from the stream */
        cmc::entropy_coding::huffman::CreateNextEncodedHuffmanCode(code, encoded_stream_view_.GetNextBit());
    }

    return codes_[code];
}

template <typename T>
inline uint64_t
StreamDecoder<T>::GetNumberOfProcessedBytesForSymbolCodewordTable() const
{
    return num_processed_bytes_symbol_codes_;
}

}

#endif /*! CMC_BITS_STREAM_DECODER_HXX */