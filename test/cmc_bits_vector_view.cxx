#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_bits_vector_view.hxx"

#include <cstddef>
#include <cstring>
#include <array>

int main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {

    /* Declare an example stream as serialized */
    const std::vector<uint8_t> serialized_byte_stream_padded{0b01000000, 0b00101101, 0b10001111, 0b10000000,
                                                             0b00011010, 0b11010000, 0b00000000, 0b00001110,
                                                             0b10001101, 0b01001010, 0b01101000, 0b01101110,
                                                             0b00000101, 0b01011101, 0b00000000, 0b00000000};

    /* Define a pointer to the data */
    const uint64_t* padded_data_ptr = reinterpret_cast<const uint64_t*>(serialized_byte_stream_padded.data());

    /* Define a view to the data */
    cmc::bits::vector_view view(padded_data_ptr);

    /* Extract the values from the stream and check them for correctness */
    const uint64_t extracted_val1 = view.GetNextBitSequence<uint64_t>(21);
    cmc::ExpectTrue(extracted_val1 == 525745);

    const uint32_t extracted_val2 = view.GetNextBitSequence<uint32_t>(19);
    cmc::ExpectTrue(extracted_val2 == 491546);

    const uint16_t extracted_val3 = view.GetNextBitSequence<uint16_t>(7);
    cmc::ExpectTrue(extracted_val3 == 104);

    view.SkipNumberOfBits(3);
    
    const uint16_t extracted_val4 = view.GetNextBitSequence<uint16_t>(10);
    cmc::ExpectTrue(extracted_val4 == 0);

    const uint32_t extracted_val5 = view.GetNextBitSequence<uint32_t>(0);
    cmc::ExpectTrue(extracted_val5 == 0);

    const uint8_t extracted_val6 = view.GetNextBitSequence<uint8_t>(4);
    cmc::ExpectTrue(extracted_val6 == 14);

    const uint64_t extracted_val7 = view.GetNextBitSequence<uint64_t>(48);
    cmc::ExpectTrue(extracted_val7 == 155350719137117);
}

/* Finalize cmc */
cmc::CmcFinalize();

return cmc::CMC_TEST_SUCCESS;
}
