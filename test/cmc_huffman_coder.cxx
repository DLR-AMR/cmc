#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_bits_vector.hxx"
#include "utilities/cmc_bits_vector_view.hxx"
#include "utilities/cmc_huffman_coder.hxx"
#include "utilities/cmc_bits_stream_decoder.hxx"

#include <vector>
#include <cstdint>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    /** EXAMPLE DATA GENERATION **/
    /* Create a symbol frequency table */
    std::vector<cmc::entropy_coding::huffman::HuffmanSymbol<int32_t>> int_sym_freq_table;

    std::vector<int32_t> int_symbols{0,2,0,0,3,5,6,3,2,4,5,6,1,1,0,0,0,0,0,0,5,4,5,3,2,1,3,4,5,6,6,3,2,3,6,0,0,0,0,3,2,2,1,5,4,3,1,1,1,4};
    int zero_freq{0}, one_freq{0}, two_freq{0}, three_freq{0}, four_freq{0}, five_freq{0}, six_freq{0}, seven_freq{0};
    for (auto sym_iter = int_symbols.begin(); sym_iter != int_symbols.end(); ++sym_iter)
    {
        switch (*sym_iter)
        {
            case 0:
            ++zero_freq;
            break;
            case 1:
            ++one_freq;
            break;
            case 2:
            ++two_freq;
            break;
            case 3:
            ++three_freq;
            break;
            case 4:
            ++four_freq;
            break;
            case 5:
            ++five_freq;
            break;
            case 6:
            ++six_freq;
            break;
        }
    }

    int_sym_freq_table.emplace_back(0, zero_freq);
    int_sym_freq_table.emplace_back(1, one_freq);
    int_sym_freq_table.emplace_back(2, two_freq);
    int_sym_freq_table.emplace_back(3, three_freq);
    int_sym_freq_table.emplace_back(4, four_freq);
    int_sym_freq_table.emplace_back(5, five_freq);
    int_sym_freq_table.emplace_back(6, six_freq);
    int_sym_freq_table.emplace_back(7, seven_freq);

    /** ENCODING **/
    /* Craete a Huffman encoder */
    cmc::entropy_coding::huffman::HuffmanCoder<int32_t> int_coder(int_sym_freq_table);

    cmc::bits::vector encoded_symbols;
    
    for (auto iter = int_symbols.begin(); iter != int_symbols.end(); ++iter)
    {
        /* Encode the current symbol */
        const cmc::entropy_coding::huffman::HuffmanCode code =  int_coder.EncodeSymbol(*iter);

        /* Store the bit sequence of the code */
        encoded_symbols.AppendBits(code.code_word, static_cast<int>(sizeof(cmc::entropy_coding::huffman::HuffmanCodeWord) * cmc::bits::kCharBit - code.code_length), 0);
    }

    /* Serialiaze the Huffman codes */
    const std::vector<uint8_t> serialized_huffman_codes = int_coder.SerializeHuffmanCodes();

    /* Serialize the stream encoded symbols */
    const std::vector<uint8_t> encoded_stream = encoded_symbols.GetSerializedByteStreamPadded();


    /** DECODING **/
    /* Construct a stream decoder for the given serialized Huffman codes */
    cmc::bits::StreamDecoder<int32_t> stream_decoder;
    stream_decoder.StartHuffmanCodesDecoding(serialized_huffman_codes.data());
    
    const uint64_t num_bytes_serialized_huff_codes = stream_decoder.GetNumberOfProcessedBytesForSymbolCodewordTable();
    cmc::ExpectTrue(num_bytes_serialized_huff_codes == 2 * sizeof(cmc::entropy_coding::huffman::HuffmanCodeInfoType)
                                                       + 7 * (sizeof(int32_t) + sizeof(uint64_t)));

    /* Define a view on the encoding of the symbols */
    cmc::bits::vector_view encoded_stream_view(reinterpret_cast<const uint64_t*>(encoded_stream.data()));

    /* Set the view as the start of the decoder */
    stream_decoder.StartDecoding(encoded_stream_view);

    /* Decode the symbols and compare for equality */
    for (size_t sym_idx{0}; sym_idx < int_symbols.size(); ++sym_idx)
    {
        /* Decode the next entropy symbol */
        const int32_t symbol = stream_decoder.DecodeNextEntropySymbol();
        
        /* Check for equality */
        cmc::ExpectTrue(symbol == int_symbols[sym_idx]);
    }

    /* Define a second stream without entropy codes */
    std::vector<uint8_t> stream2{0b01101101, 0b00010011, 0b10010010, 0b10111010,
                                 0b00011100, 0b11101011, 0b00000010, 0b11011001,
                                 0b00100111, 0b01010101, 0b11110000, 0b01001101,
                                 0b11001100, 0b11000010, 0b01010101, 0b11000000};
    /* Define a view on the stream */
    cmc::bits::vector_view encoded_stream_view2(reinterpret_cast<const uint64_t*>(stream2.data()));

    /* Set the view to the encoded data (without entropy codes) */
    cmc::bits::StreamDecoder<uint32_t> stream_decoder2;
    stream_decoder2.StartDecoding(encoded_stream_view2);

    const bool is_bit_at_pos_0_set = stream_decoder2.GetNextBit();
    cmc::ExpectTrue(is_bit_at_pos_0_set == false);
    
    stream_decoder2.SkipNextBits(5);

    const bool is_bit_at_pos_6_set = stream_decoder2.GetNextBit();
    cmc::ExpectTrue(is_bit_at_pos_6_set == false);

    const uint32_t byte_seq = stream_decoder2.GetNextBitSequence<uint32_t>(16);
    cmc::ExpectTrue(byte_seq == uint32_t{35273});

    stream_decoder2.SkipNextBits(48);

    const bool is_bit_at_pos_71_set = stream_decoder2.GetNextBit();
    cmc::ExpectTrue(is_bit_at_pos_71_set == true);

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
