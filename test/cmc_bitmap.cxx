#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_bit_map.hxx"

int
main(void)
{

    /* Initialize cmc */
    cmc::CmcInitialize();

    /* Create a BitMap object */
    cmc::bit_map::BitMap bitmap;

    /* Allocate some memory for the bitmap */
    const size_t num_bits = 11;
    bitmap.Reserve(num_bits);

    /* Set the bits */
    bitmap.AppendBit(true); //0b00000001
    bitmap.AppendUnsetBit(); //0b00000001
    bitmap.AppendSetBit(); //0b00000101

    /* Check if the size of the bitmap coincides with the appended bits */
    cmc::ExpectTrue(bitmap.size() == 3);

    /* Set some bits */
    bitmap.AppendBit(true); //0b00001101
    bitmap.AppendBit(true); //0b00011101
    bitmap.AppendBit(true); //0b00111101

    bitmap.ToggleBit(3); //0b00110101
    bitmap.ClearBit(5); //0b00010101

    bitmap.AppendUnsetBit(); //0b00010101
    bitmap.SetBit(6); //0b01010101

    bitmap.AppendBit(false); //0b01010101 (first byte is full)

    /* Toggling a bit an checking whther its true */
    bitmap.ToggleBit(0, 5); //0b01110101
    cmc::ExpectTrue(bitmap.IsBitSet(5));
    cmc::ExpectTrue(bitmap.IsBitSet(0, 5));
    bitmap.ToggleBit(5); //0b01010101

    /* Fill a few bits of the second byte */
    bitmap.AppendBit(true); //0b01010101 00000001
    bitmap.AppendBit(true); //0b01010101 00000011
    bitmap.AppendBit(true); //0b01010101 00000111

    bitmap.ClearBit(1, 1); //0b01010101 00000101
    bitmap.SetBit(1, 2); //0b01010101 00000101 (stays unchanged)

    /* Check the final size */
    cmc::ExpectTrue(bitmap.size() == 11);

    int idx = 0;
    /* Check the alternating bit pattern that should have been employed */
    for (auto iter = bitmap.begin(); iter != bitmap.end(); ++iter)
    {
        if (idx % 2 == 0)
        {
            cmc::ExpectTrue(*iter);
        } else
        {
            cmc::ExpectFalse(*iter);
        }
        ++idx;
    }

    /* Fill a BitMap */
    cmc::bit_map::BitMap bitmap2;
    bitmap2.AppendBit(true); //0b00000001
    bitmap2.AppendBit(false); //0b00000001
    bitmap2.AppendBit(true); //0b00000101
    bitmap2.AppendBit(true); //0b00001101

    /* Copy constructor */
    cmc::bit_map::BitMap bitmap3 = bitmap2; //0b00001101

    /* Append a whole BitMap */
    bitmap3.AppendBits(bitmap2); //0b11011101

    cmc::ExpectTrue(bitmap3.size() == 8);

    /* Check whether the bitmap has been appended correctly */
    const uint8_t encoded_byte = *(bitmap3.data());
    cmc::ExpectTrue(encoded_byte == 221);

    /* Append some more */
    bitmap3.AppendBit(true);  //0b11011101 00000001
    bitmap3.AppendBit(false); //0b11011101 00000001
    bitmap3.AppendBit(false); //0b11011101 00000001
    bitmap3.AppendBit(false); //0b11011101 00000001
    bitmap3.AppendBit(false); //0b11011101 00000001
    bitmap3.AppendBit(true);  //0b11011101 00100001


    /* Append one to bitmap2 */
    bitmap2.AppendBit(true); //0b00011101

    /* And again append bitap2 to bitmap3 */
    bitmap3.AppendBits(bitmap2); //0b11011101 01100001 00000111

    cmc::ExpectTrue(bitmap3.size() == 19);

    const uint8_t first_encoded_byte = *(bitmap3.data());
    const uint8_t second_encoded_byte = *(bitmap3.data() + 1);
    const uint8_t third_encoded_byte = *(bitmap3.data() + 2);

    cmc::ExpectTrue(first_encoded_byte == 221);
    cmc::ExpectTrue(second_encoded_byte == 97);
    cmc::ExpectTrue(third_encoded_byte == 7);

    /* Append the BitMap itself */
    cmc::bit_map::BitMap bitmap4 = bitmap3; //0b00001101
    bitmap3.AppendBits(bitmap4); //0b11011101 01100001 11101111 00001110 00111011

    cmc::ExpectTrue(bitmap3.size() == 38);

    const uint8_t again_first_encoded_byte = *(bitmap3.data());
    const uint8_t again_second_encoded_byte = *(bitmap3.data() + 1);
    const uint8_t again_third_encoded_byte = *(bitmap3.data() + 2);
    const uint8_t again_fourth_encoded_byte = *(bitmap3.data() + 3);
    const uint8_t again_fifth_encoded_byte = *(bitmap3.data() + 4);

    cmc::ExpectTrue(again_first_encoded_byte == 221);
    cmc::ExpectTrue(again_second_encoded_byte == 97);
    cmc::ExpectTrue(again_third_encoded_byte == 239);
    cmc::ExpectTrue(again_fourth_encoded_byte == 14);
    cmc::ExpectTrue(again_fifth_encoded_byte == 59);

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
