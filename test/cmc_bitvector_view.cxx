#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_bit_vector.hxx"
#include "utilities/cmc_log_functions.hxx"

#include <vector>

int
main(void)
{

    /* Initialize cmc */
    cmc::CmcInitialize();

    const size_t num_bytes = 4;

    std::vector<uint8_t> bit_sequence;
    bit_sequence.reserve(num_bytes);

    /* Create a bit sequence */
    bit_sequence.push_back(uint8_t{0b11010111});
    bit_sequence.push_back(uint8_t{0b11000010});
    bit_sequence.push_back(uint8_t{0b00000111});///bei 4) : 00010000 001
    bit_sequence.push_back(uint8_t{0b01010101});
    
    /* Create a BitVectorView on the bit sequence */
    cmc::bit_vector::BitVectorView bv_view(bit_sequence.data(), num_bytes);

    std::vector<uint8_t> extracted_sequence;

    /* Retrieve the bit sequences MSB ('most significant bit'-first) fromt the sequence */
    extracted_sequence = bv_view.GetNextBitSequence(2);
    cmc::cmc_debug_msg("1 Value of bit sequence: ", static_cast<int>(extracted_sequence.front()));
    cmc::ExpectTrue(extracted_sequence.front() == uint8_t{0b11000000});

    extracted_sequence = bv_view.GetNextBitSequence(4);
    cmc::cmc_debug_msg("2 Value of bit sequence: ", static_cast<int>(extracted_sequence.front()));
    cmc::ExpectTrue(extracted_sequence.front() == uint8_t{0b01010000});

    extracted_sequence = bv_view.GetNextBitSequence(5);
    cmc::cmc_debug_msg("3 Value of bit sequence: ", static_cast<int>(extracted_sequence.front()));
    cmc::ExpectTrue(extracted_sequence.front() == uint8_t{0b11110000});

    extracted_sequence = bv_view.GetNextBitSequence(11);
    cmc::cmc_debug_msg("4 Values of bit sequence: ", static_cast<int>(extracted_sequence[0]), " und ", static_cast<int>(extracted_sequence[1]));
    cmc::ExpectTrue(extracted_sequence[0] == uint8_t{0b00010000});
    cmc::ExpectTrue(extracted_sequence[1] == uint8_t{0b00100000});

    extracted_sequence = bv_view.GetNextBitSequence(3);
    cmc::cmc_debug_msg("5 Value of bit sequence: ", static_cast<int>(extracted_sequence.front()));
    cmc::ExpectTrue(extracted_sequence.front() == uint8_t{0b11000000});

    extracted_sequence = bv_view.GetNextBitSequence(7);
    cmc::cmc_debug_msg("6 Value of bit sequence: ", static_cast<int>(extracted_sequence.front()));
    cmc::ExpectTrue(extracted_sequence.front() == uint8_t{0b10101010});

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
