#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_ac_model.hxx"
#include "utilities/cmc_entropy_alphabet.hxx"
#include "utilities/cmc_arithmetic_encoder.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <cstddef>
#include <vector>
#include <memory>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();

    /* An uncompressed message */
    std::vector<uint32_t> uncompressed_message{0,5,4,0,2,1,6,5,2,0,3,5,3,4,6,2,4,5,5,5,2,3,7,3,1,2,5,5,6,6,2,5,2,6,7,7};

    /* Total count of letters */
    const size_t num_symbols_in_message = uncompressed_message.size();

    /* The alphabet and the letter probabilities of the message */
    std::vector<cmc::entropy_coding::Letter> alphabet;
    alphabet.emplace_back(cmc::entropy_coding::Letter{0, 3});
    alphabet.emplace_back(cmc::entropy_coding::Letter{1, 2});
    alphabet.emplace_back(cmc::entropy_coding::Letter{2, 7});
    alphabet.emplace_back(cmc::entropy_coding::Letter{3, 4});
    alphabet.emplace_back(cmc::entropy_coding::Letter{4, 3});
    alphabet.emplace_back(cmc::entropy_coding::Letter{5, 9});
    alphabet.emplace_back(cmc::entropy_coding::Letter{6, 5});
    alphabet.emplace_back(cmc::entropy_coding::Letter{7, 3});

    /* Construct a frequency model based on the alphabet and the probabilities of the symbols in the message */
    cmc::arithmetic_encoding::StaticFrequencyModel frequency_model(alphabet);

    /* Construct an arithmetic encoder based on the alphabet with the given probabilities */
    cmc::arithmetic_encoding::Encoder encoder(std::make_unique<cmc::arithmetic_encoding::StaticFrequencyModel>(frequency_model));
    
    /* Iterate over all symbols in the message */
    for (size_t symbol_idx = 0; symbol_idx < num_symbols_in_message; ++symbol_idx)
    {
        /* Encode the current symbol */
        encoder.EncodeSymbol(uncompressed_message[symbol_idx]);
    }

    /* Complete the encoding process after all symbols have been encoded */
    encoder.FinishEncoding();

    /* Get the encoded bit stream */
    cmc::bit_map::BitMap encoded_stream = encoder.GetEncodedBitStream();

    /* Create a view on the encoded stream */
    cmc::bit_map::BitMapView encoded_stream_view(encoded_stream);

    /* Start the decoding based on the same frequency model */
    cmc::arithmetic_encoding::Decoder decoder(std::make_unique<cmc::arithmetic_encoding::StaticFrequencyModel>(frequency_model), encoded_stream_view);

    /* Set up a vector for the decompressed symbols */
    std::vector<uint32_t> decoded_symbols;
    decoded_symbols.reserve(num_symbols_in_message);

    /* Decode all symbols that have been encoded */
    for (size_t symbol_idx = 0; symbol_idx < num_symbols_in_message; ++symbol_idx)
    {
        /* Decode the next symbol */
        const uint32_t symbol = decoder.DecodeNextSymbol();
        decoded_symbols.push_back(symbol);
    }

    /* Compare the uncompressed and decompressed message for equality */
    for (size_t symbol_idx = 0; symbol_idx < num_symbols_in_message; ++symbol_idx)
    {
        //cmc::cmc_debug_msg("Symbol Index: ", symbol_idx, ", Uncompressed Symbol: ", uncompressed_message[symbol_idx], ", Decompressed Symbol: ", decoded_symbols[symbol_idx]);
        cmc::ExpectTrue(uncompressed_message[symbol_idx] == decoded_symbols[symbol_idx]);
    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
