#include "cmc.hxx"
#include "test/cmc_test.hxx"

#include "amr/lossless/cmc_par_multi_res_extraction_util.hxx"
#include "utilities/cmc_bits_stream_decoder.hxx"
#include "utilities/cmc_huffman_coder.hxx"

#include <array>
#include <vector>
#include <span>
#include <numeric>
#include <algorithm>
#include <memory>
#include <limits>

void
Test2DFloatData()
{
    constexpr int32_t N = 16;
    constexpr int DIM = 2;

    /* We define some example data on one element */
    std::array<float, N> init_data{2.34,2.79,3.10,2.56,
                                   3.21,3.07,2.99,2.76,
                                   3.01,2.88,2.34,2.19,
                                   2.27,2.45,2.76,2.56};

    /* Compute the element encoding */
    const cmc::par::lossless::multi_res::IntraElementCoding<float, DIM, N> elem_coding = cmc::par::lossless::multi_res::ComputeElementEncoding<float, DIM, N>(std::span<float>(init_data));
    
    /* Get the coarse level predictor */
    const float coarse_predictor = elem_coding.GetCoarseValuePredictor();

    /* Setup of the Huffman Coder */
    ///////////////////////////////
    constexpr int num_entropy_symbols = cmc::par::lossless::multi_res::GetNumEntropySymbols<float>();
    std::vector<uint64_t> entropy_symbol_frequencies(num_entropy_symbols, 0);
    
    /* Iterate over the entropy codes and collect the frequencies */
    for (const auto& symbol : elem_coding.entropy_codes)
    {
        /* Update the frequency */
        ++entropy_symbol_frequencies[cmc::par::lossless::multi_res::MapEntropySymbolToArrayIndex<float>(symbol)];
    }

    std::vector<cmc::entropy_coding::huffman::EntropySymbol<cmc::par::lossless::multi_res::SymbolType>> global_symbol_frequencies;
    global_symbol_frequencies.reserve(num_entropy_symbols);

    for (int idx{0}; idx < num_entropy_symbols; ++idx)
    {
        /* Convert the index back to the entropy symbol */
        const cmc::par::lossless::multi_res::SymbolType entropy_symbol = cmc::par::lossless::multi_res::MapArrayIndexToEntropySymbol<float>(idx);

        /* Store the symbol with the global frequency */
        global_symbol_frequencies.emplace_back(entropy_symbol, entropy_symbol_frequencies[idx]);
    }

    /* We initialize the Huffman coder */
    cmc::entropy_coding::huffman::HuffmanCoder<cmc::par::lossless::multi_res::SymbolType> huffman_coder(global_symbol_frequencies);
    
    /* We serialize the Huffman coder */
    const std::vector<uint8_t> serialized_huffman_codes = huffman_coder.SerializeHuffmanCodes();
    ////////////////////////////////////
    /* End of the Huffman Coder setup */

    /* We encode the data */
    cmc::bits::vector elem_encoding;
    cmc::par::lossless::multi_res::PerformElementEncoding<float, DIM, N>(elem_encoding, huffman_coder, elem_coding);

    /* Serialize the stream of the encoded element data */
    const std::vector<uint64_t> encoded_elem_stream = elem_encoding.GetSerializedByteStreamBE();

    /* We setup a stream decoder for the data */
    /* Construct a stream decoder for the given serialized Huffman codes */
    cmc::bits::StreamDecoder<cmc::par::lossless::multi_res::SymbolType> stream_decoder;
    stream_decoder.StartHuffmanCodesDecoding(serialized_huffman_codes.data());

    /* Set the start of the decoder to the encoded stream */
    cmc::bits::vector_view encoded_stream_view(encoded_elem_stream.data());
    stream_decoder.StartDecoding(encoded_stream_view);

    /* We decompress the encoded data */
    std::array<float, N> decompressed_data = cmc::par::lossless::multi_res::PerformElementDecoding<float, DIM, N>(stream_decoder, coarse_predictor);

    /* We compare the initial and decompressed data for equality */
    for (int idx{0}; idx < N; ++idx)
    {
        cmc::ExpectTrue(init_data[idx] == decompressed_data[idx]); 
    }
}

void
Test3DDoubleData()
{
    constexpr int32_t N = 27;
    constexpr int DIM = 3;

    /* We define some example data on one element */
    std::array<double, N> init_data{2.54,2.49,3.03,2.26,
                                   3.41,3.12,2.87,2.66,
                                   3.04,2.89,2.44,2.23,
                                   2.17,2.25,2.66,2.45,
                                   3.15,3.27,3.87,3.99,
                                   4.05,4.44,4.53,4.78,
                                   4.99,4.89,4.34};

    /* Compute the element encoding */
    const cmc::par::lossless::multi_res::IntraElementCoding<double, DIM, N> elem_coding = cmc::par::lossless::multi_res::ComputeElementEncoding<double, DIM, N>(std::span<double>(init_data));
    
    /* Get the coarse level predictor */
    const double coarse_predictor = elem_coding.GetCoarseValuePredictor();

    /* Setup of the Huffman Coder */
    ///////////////////////////////
    constexpr int num_entropy_symbols = cmc::par::lossless::multi_res::GetNumEntropySymbols<double>();
    std::vector<uint64_t> entropy_symbol_frequencies(num_entropy_symbols, 0);
    
    /* Iterate over the entropy codes and collect the frequencies */
    for (const auto& symbol : elem_coding.entropy_codes)
    {
        /* Update the frequency */
        ++entropy_symbol_frequencies[cmc::par::lossless::multi_res::MapEntropySymbolToArrayIndex<double>(symbol)];
    }

    std::vector<cmc::entropy_coding::huffman::EntropySymbol<cmc::par::lossless::multi_res::SymbolType>> global_symbol_frequencies;
    global_symbol_frequencies.reserve(num_entropy_symbols);

    for (int idx{0}; idx < num_entropy_symbols; ++idx)
    {
        /* Convert the index back to the entropy symbol */
        const cmc::par::lossless::multi_res::SymbolType entropy_symbol = cmc::par::lossless::multi_res::MapArrayIndexToEntropySymbol<double>(idx);

        /* Store the symbol with the global frequency */
        global_symbol_frequencies.emplace_back(entropy_symbol, entropy_symbol_frequencies[idx]);
    }

    /* We initialize the Huffman coder */
    cmc::entropy_coding::huffman::HuffmanCoder<cmc::par::lossless::multi_res::SymbolType> huffman_coder(global_symbol_frequencies);
    
    /* We serialize the Huffman coder */
    const std::vector<uint8_t> serialized_huffman_codes = huffman_coder.SerializeHuffmanCodes();
    ////////////////////////////////////
    /* End of the Huffman Coder setup */

    /* We encode the data */
    cmc::bits::vector elem_encoding;
    cmc::par::lossless::multi_res::PerformElementEncoding<double, DIM, N>(elem_encoding, huffman_coder, elem_coding);

    /* Serialize the stream of the encoded element data */
    const std::vector<uint64_t> encoded_elem_stream = elem_encoding.GetSerializedByteStreamBE();

    /* We setup a stream decoder for the data */
    /* Construct a stream decoder for the given serialized Huffman codes */
    cmc::bits::StreamDecoder<cmc::par::lossless::multi_res::SymbolType> stream_decoder;
    stream_decoder.StartHuffmanCodesDecoding(serialized_huffman_codes.data());

    /* Set the start of the decoder to the encoded stream */
    cmc::bits::vector_view encoded_stream_view(encoded_elem_stream.data());
    stream_decoder.StartDecoding(encoded_stream_view);

    /* We decompress the encoded data */
    std::array<double, N> decompressed_data = cmc::par::lossless::multi_res::PerformElementDecoding<double, DIM, N>(stream_decoder, coarse_predictor);

    /* We compare the initial and decompressed data for equality */
    for (int idx{0}; idx < N; ++idx)
    {
        cmc::ExpectTrue(init_data[idx] == decompressed_data[idx]); 
    }
}

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {

    /* Transform some values and check their expecetd bit pattern */
    static_assert(sizeof(char) == sizeof(uint8_t));
    const uint8_t val_u8_trans = cmc::par::lossless::multi_res::TransformToUInteger<char>(static_cast<char>('a'));
    cmc::ExpectTrue(val_u8_trans == static_cast<uint8_t>(97));

    const uint16_t val_u16_trans = cmc::par::lossless::multi_res::TransformToUInteger<int16_t>(static_cast<int16_t>(12003));
    cmc::ExpectTrue(val_u16_trans == static_cast<uint16_t>(12003));

    static_assert(sizeof(float) == sizeof(uint32_t));
    const uint32_t val_u32_trans = cmc::par::lossless::multi_res::TransformToUInteger<float>(static_cast<float>(32.0));
    cmc::ExpectTrue(val_u32_trans == static_cast<uint32_t>(1107296256));

    static_assert(sizeof(double) == sizeof(uint64_t));
    const uint64_t val_u64_trans = cmc::par::lossless::multi_res::TransformToUInteger<double>(static_cast<double>(67.25));
    cmc::ExpectTrue(val_u64_trans == static_cast<uint64_t>(4634432714982817792));

    /* Create some entropy symbols and check their expeceted value */
    cmc::ExpectTrue(cmc::par::lossless::multi_res::GetNumEntropySymbols<float>() == 2 * (sizeof(float) * cmc::bits::kCharBit + 1));
    cmc::ExpectTrue(cmc::par::lossless::multi_res::GetNumEntropySymbols<double>() == 2 * (sizeof(double) * cmc::bits::kCharBit + 1));

    //The residual has 18 leading zeros and it is assumed to be greater than the approximation (expected: 128 + 18)
    cmc::par::lossless::multi_res::SymbolType s1 = cmc::par::lossless::multi_res::CreateEntropySymbol(true, static_cast<uint32_t>(14300));
    cmc::ExpectTrue(s1 == static_cast<cmc::par::lossless::multi_res::SymbolType>(146));
    cmc::ExpectTrue(cmc::par::lossless::multi_res::GetLZCFromEntropySymbol(s1) == 18);
    cmc::ExpectTrue(s1 == cmc::par::lossless::multi_res::MapArrayIndexToEntropySymbol<float>(cmc::par::lossless::multi_res::MapEntropySymbolToArrayIndex<float>(s1)));

    //The residual has 47 leading zeros and it is assumed to be smaller than the approximation (expected: 0 + 47)
    cmc::par::lossless::multi_res::SymbolType s2 = cmc::par::lossless::multi_res::CreateEntropySymbol(false, static_cast<uint64_t>(67210));
    cmc::ExpectTrue(s2 == static_cast<cmc::par::lossless::multi_res::SymbolType>(47));
    cmc::ExpectTrue(cmc::par::lossless::multi_res::GetLZCFromEntropySymbol(s2) == 47);
    cmc::ExpectTrue(s2 == cmc::par::lossless::multi_res::MapArrayIndexToEntropySymbol<double>(cmc::par::lossless::multi_res::MapEntropySymbolToArrayIndex<double>(s2)));

    /* Perform the calculation of a predictor*/
    const std::array<float, 4> exdata1{1.0,2.0,3.0,4.0};
    const float fmean1 = cmc::par::lossless::multi_res::ComputeArithmeticMean<float, 4>(exdata1);
    cmc::ExpectTrue(fmean1 >= 2.499969482421875 && fmean1 <= 2.5000152587890625);
    const float fmid1 = cmc::par::lossless::multi_res::ComputeMidRange<float, 4>(exdata1);
    cmc::ExpectTrue(fmid1 >= 2.499969482421875 && fmid1 <= 2.5000152587890625);

    /* Test an example intra element coding in 2D with float data */
    Test2DFloatData();
    
    /* Test an example intra element coding in 3D with double data */
    Test3DDoubleData();
    
    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
