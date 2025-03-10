#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_arithmetic_encoding.hxx"
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


    /* Create the arithmetic encoder with a static dictionary */
    cmc::entropy_coding::arithmetic_coding::Encoder encoder;

    /* Initilaize the alphabet */
    encoder.InitializeAlphabet();

    /* Create the alphabet with the given frequencies */
    for (auto sym_iter = uncompressed_message.begin(); sym_iter != uncompressed_message.end(); ++sym_iter)
    {
        encoder.UpdateSymbolFrequency(*sym_iter);
    }

    /* Setup the interior structures for encoding once the frequencies for the static dictionary have been collected */
    encoder.SetupEncoding();

    /* Encode each symbol of the message */
    for (size_t symbol_idx = 0; symbol_idx < num_symbols_in_message; ++symbol_idx)
    {
        /* Encode the current symbol */
        encoder.EncodeSymbol(uncompressed_message[symbol_idx]);
    }
    
    /* Complete the encoding process after all symbols have been encoded */
    encoder.FinishEncoding();

    /* Get the encoded bit stream */
    cmc::bit_map::BitMap encoded_stream = encoder.GetEncodedBitStream();

    cmc::cmc_debug_msg("Encoded stream length: ", encoded_stream.size());

    #if 0
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

    #endif

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
