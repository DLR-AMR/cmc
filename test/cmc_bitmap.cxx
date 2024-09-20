#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_bit_map.hxx"
#include "utilities/cmc_log_functions.hxx"

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
        cmc::cmc_debug_msg("Bit Value at pos: ", idx, " is: ", *iter);
        ++idx;
    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
