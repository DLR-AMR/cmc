#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_bits_vector.hxx"

#include <cstddef>
#include <cstring>
#include <array>

int main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {

    /* Declare a new bits vector */
    cmc::bits::vector bits_vec;
    bits_vec.Reserve(68);

    const double fval1 = 14.765625;
    const uint64_t val1 = std::bit_cast<uint64_t>(fval1);
    //IEEE-754: 01000000 00101101 10001000 00000000 00000000 00000000 00000000 00000000 
    //          ^start_pos            ^end_pos
    const int start_pos1 = 0, end_pos1 = 43;

    const float fval2 = -1.9375;
    const uint32_t val2 = std::bit_cast<uint32_t>(fval2);
    ///IEEE-754: 10111111 11111000 00000000 00000000
    //                     ^start_pos ^end_pos
    const int start_pos2 = 9, end_pos2 = 12;

    const uint32_t val3 = 429;
    //BE: 00000000 00000000 00000001 10101101
    //                          ^start_pos  ^end_pos
    const int start_pos3 = 20, end_pos3 = 0; 

    const uint64_t val4 = 256000024567893;
    //BE: 00000000 00000000 11101000 11010100 10100110 10000110 11100000 01010101
    //    ^start_pos                                                            ^end_pos
    const int start_pos4 = 0, end_pos4 = 0;

    const uint8_t val5 = 104;
    //        01101000
    //start_pos^  ^end_pos
    const int start_pos5 = 1, end_pos5 = 3; 

    const uint16_t val6 = 0;
    const int start_pos6 = 1, end_pos6 = 15; //Empty

    /* Append the specified bit-sequences from the values */
    bits_vec.AppendBits(val1, start_pos1, end_pos1);
    bits_vec.AppendBits(val2, start_pos2, end_pos2);
    bits_vec.AppendBits(val3, start_pos3, end_pos3);
    bits_vec.AppendBits(val4, start_pos4, end_pos4);
    bits_vec.AppendBits(val5, start_pos5, end_pos5);
    bits_vec.AppendBits(val6, start_pos6, end_pos6);

    //Bits vector layout in BE should be:
    // 01000000 00101101 10001111 10000000 00011010 11010000 00000000 00001110 10001101 01001010 01101000 01101110 00000101 01011101
    const std::vector<uint8_t> expected_bytes{0b01000000, 0b00101101, 0b10001111, 0b10000000, 0b00011010, 0b11010000,
                                              0b00000000, 0b00001110, 0b10001101, 0b01001010, 0b01101000, 0b01101110,
                                              0b00000101, 0b01011101};
                                            
    const std::vector<uint8_t> serialized_byte_stream = bits_vec.GetSerializedByteStream();

    cmc::ExpectTrue(serialized_byte_stream.size() == expected_bytes.size());
    cmc::ExpectTrue(bits_vec.size_bytes() == expected_bytes.size());
    cmc::ExpectTrue(bits_vec.size() == 112);

    /* Compare each byte in the stream and check if the encoded value is correct */
    for (size_t i{0}; i < expected_bytes.size(); ++i)
    {
        cmc::ExpectTrue(expected_bytes[i] == serialized_byte_stream[i]);
    }
    
    /* Get the serialized stream padded to contain full uint64_t's */
    const std::vector<uint8_t> serialized_byte_stream_padded = bits_vec.GetSerializedByteStreamPadded();
    cmc::ExpectTrue(serialized_byte_stream_padded.size() == 16);

    /* Declare a new bits vector */
    cmc::bits::vector bits_vec2;

    /* Append bits to the vector such that the encoded stream in (BE should be): 10111000 00001000 11 */
    bits_vec2.AppendBit(true);
    bits_vec2.AppendBit(false);
    bits_vec2.AppendSetBit();
    bits_vec2.AppendSetBit();
    bits_vec2.AppendSetBit();
    bits_vec2.AppendUnsetBit();
    bits_vec2.AppendBits(uint8_t{2}, 0, 0);
    bits_vec2.AppendUnsetBit();
    bits_vec2.AppendUnsetBit();
    bits_vec2.AppendSetBit();
    bits_vec2.AppendSetBit();

    cmc::ExpectTrue(bits_vec2.size() == 18);
    cmc::ExpectTrue(bits_vec2.size_bytes() == 3);

    /* Define the bytes that are expected */
    std::vector<uint8_t> expected_bytes2{0b10111000, 0b00001000, 0b11000000};

    /* Get the serialized bit stream in big endian */
    const std::vector<uint8_t> serialized_byte_stream2 = bits_vec2.GetSerializedByteStream();

    cmc::ExpectTrue(serialized_byte_stream2.size() == 3);

    /* Check the expected bytes for equality with the serialization */
    for (size_t i{0}; i < expected_bytes2.size(); ++i)
    {
        cmc::ExpectTrue(expected_bytes2[i] == serialized_byte_stream2[i]);
    }
}

/* Finalize cmc */
cmc::CmcFinalize();

return cmc::CMC_TEST_SUCCESS;
}
