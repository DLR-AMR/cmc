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

    /* Append some more */
    bitmap.AppendBit(true);  //0b01111101
    bitmap.AppendBit(false); //0b01111101
    bitmap.AppendBit(false); //0b01111101 00000000
    bitmap.AppendBit(false); //0b01111101 00000000
    bitmap.AppendBit(false); //0b01111101 00000000
    bitmap.AppendBit(false); //0b01111101 00000000
    bitmap.AppendBit(true);  //0b01111101 00010000
    bitmap.AppendBit(false);  //0b01111101 00010000
    bitmap.AppendBit(false);  //0b01111101 00010000
    bitmap.AppendBit(true);  //0b01111101 10010000
    bitmap.AppendBit(false);  //0b01111101 00010000
    bitmap.AppendBit(true);  //0b01111101 10010000 00000010

    /* Create a BitMapView object */
    cmc::bit_map::BitMapView bitmap_viewer(bitmap);

    cmc::ExpectTrue(bitmap_viewer.GetNextBit() == true);
    cmc::ExpectTrue(bitmap_viewer.GetNextBit() == false);
    cmc::ExpectTrue(bitmap_viewer.GetNextBit() == true);
    cmc::ExpectTrue(bitmap_viewer.GetNextBit() == true);
    cmc::ExpectTrue(bitmap_viewer.GetNextBit() == true);

    std::vector<uint8_t> next_five_bits = bitmap_viewer.GetNextNumberOfBits(5);

    cmc::ExpectTrue(next_five_bits.size() == 1);
    cmc::ExpectTrue(next_five_bits.front() == 3);

    cmc::ExpectTrue(bitmap_viewer.GetNextBit() == false);
    cmc::ExpectTrue(bitmap_viewer.GetNextBit() == false);
    cmc::ExpectTrue(bitmap_viewer.GetNextBit() == true);
    cmc::ExpectTrue(bitmap_viewer.GetNextBit() == false);

    std::vector<uint8_t> next_four_bits = bitmap_viewer.GetNextNumberOfBits(4);

    cmc::ExpectTrue(next_four_bits.size() == 1);
    cmc::ExpectTrue(next_four_bits.front() == 10);

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
