#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_huffman_coder.hxx"

#include <vector>
#include <cstdint>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
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

    /* Craete a Huffman encoder */
    cmc::entropy_coding::huffman::HuffmanCoder<int32_t> int_coder(int_sym_freq_table);

    std::vector<uint8_t> int_sym_freq_table_serialized = int_coder.GetSerializedSymbolFrequencyTable();

    cmc::bit_vector::BitVector encoded_symbols;
    
    for (auto iter = int_symbols.begin(); iter != int_symbols.end(); ++iter)
    {
        const auto [code, num_bits] = int_coder.EncodeSymbol(*iter);

        encoded_symbols.AppendBits(code, num_bits);
    }

    encoded_symbols.TrimToContent();

    /* Create a HuffmanDecoder */
    cmc::entropy_coding::huffman::HuffmanDecoder<int32_t> int_decoder(int_sym_freq_table_serialized.data());
    
    const size_t kExpecetdNumBytesSerializedSymbolFrequencyTable = 64;
    cmc::ExpectTrue(int_sym_freq_table_serialized.size() == kExpecetdNumBytesSerializedSymbolFrequencyTable);

    /* Set a view on the decoding stream */
    cmc::bit_vector::BitVectorView encoded_view(encoded_symbols.data(), encoded_symbols.size());
    int_decoder.StartDecoding(encoded_view);

    /* Decode the symbols and compare for equality */
    for (size_t iter{0}; iter < int_symbols.size(); ++iter)
    {
        /* Decode the next symbol */
        const int32_t current_symbol = int_decoder.DecodeNextSymbol();

        /* Check for equality */
        cmc::ExpectTrue(current_symbol == int_symbols[iter]);
    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
