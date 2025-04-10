#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_multi_res_entropy_coder.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_log_functions.hxx"

#include "mpi/cmc_mpi.hxx"

#include <cstddef>
#include <vector>
#include <memory>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();

    using variable_data_type_t = float; //!< The actual data type of the underlying variable (which may vary)
    using entropy_symbol_type_t = uint32_t; //!< The data type which is used for the symbols (which is fixed)

    const uint32_t sign = cmc::entropy_coding::arithmetic_coding::kByteCompressionSignumBit;
    /* An uncompressed message */
    //std::vector<entropy_symbol_type_t> uncompressed_message{8,5,4,0,2,1,6,sign+5,2,sign+0,3,5,sign+3,4,6,2,4,5,5,5,2,3,sign+7,3,1,2,5,5,6,sign+6,2,5,sign+2,6,7,7,cmc::entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte};
    std::vector<entropy_symbol_type_t> uncompressed_message{8,5,4,0,2,1,6,sign+5,2,sign+0,3,5,sign+3,4,6,2,4,5,5,5,2,3,sign+7,3,1,2,5,5,6,sign+6,2,5,sign+2,6,7,7,cmc::entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte, 8,5,4,0,2,1,6,sign+5,2,sign+0,3,5,sign+3,4,6,2,4,5,5,5,2,3,sign+7,3,1,2,5,5,5,5,5,2,3,sign+7,3,1,2,5,5,6,sign+6,2,5,sign+2,6,7,7,cmc::entropy_coding::arithmetic_coding::kByteCompressionSymbolJumpToNextByte, 8,5,4,0,2,1,6,sign+5,2,sign+0,3,5,sign+3,4,6,2,4,5,5,5,2,3,sign+7,3,1,2,5};

    /* Total count of letters */
    const size_t num_symbols_in_message = uncompressed_message.size();

    /* Setup the encoder */
    cmc::entropy_coding::arithmetic_coding::MultiResEncoder<variable_data_type_t> ac;

    /* Initialize the alphabet */
    ac.InitializeAlphabet();

    /* Update the frequencies accordingly */
    for (auto sym_iter = uncompressed_message.begin(); sym_iter != uncompressed_message.end(); ++sym_iter)
    {
        ac.UpdateSymbolFrequency(*sym_iter);
    }

    /* After the frequencies have been collected, we setup the encoding process */
    ac.SetupEncoding(MPI_COMM_SELF);

    /* Now, we start encoding the symbols */
    for (auto sym_iter = uncompressed_message.begin(); sym_iter != uncompressed_message.end(); ++sym_iter)
    {
        ac.EncodeSymbol(*sym_iter);
    }

    /* Afterwards, we flush the encoding stream, once all symbols have been encoded */
    ac.FinishEncoding();

    /* Encode the alphabet */
    cmc::bit_vector::BitVector encoded_alphabet = ac.EncodeAlphabet();

    /* Get the encoded bit stream */
    cmc::bit_map::BitMap encoded_stream = ac.GetEncodedBitStream();

    /* Create a view on the encoded stream */
    cmc::bit_map::BitMapView encoded_stream_view(encoded_stream);

    /* Decode the alphabet */
    //[[maybe_unused]] auto [frequency_model, num_alphabet_bytes] = cmc::entropy_coding::arithmetic_coding::DecodeByteCompressionStaticFrequencyAlphabet<variable_data_type_t>(encoded_alphabet.begin());

    /* Setup the entropy decoder */
    cmc::entropy_coding::arithmetic_coding::MultiResDecoder<variable_data_type_t> decoder(encoded_alphabet.begin(), encoded_stream_view);

    /* Set up the decoder */
    decoder.SetupDecoding();

    /* Set up a vector for the decompressed symbols */
    std::vector<entropy_symbol_type_t> decoded_symbols;
    decoded_symbols.reserve(num_symbols_in_message);

    /* Decode all symbols in the encoded message */
    for (size_t idx = 0; idx < num_symbols_in_message; ++idx)
    {
        /* Decode the next symbol */
        const entropy_symbol_type_t symbol = decoder.DecodeNextSymbol();
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
