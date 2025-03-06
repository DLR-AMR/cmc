#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_byte_compression_values.hxx"


#include <cstddef>
#include <vector>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {

        uint32_t initial_value = 1213;
        constexpr int N = static_cast<int>(sizeof(uint32_t));

        if constexpr (N != 4)
        {
            return cmc::CMC_TEST_SKIP;
        }

        /* Construct a compression value from the initial value */
        cmc::SerializedCompressionValue<N> scv(initial_value);

        /* Clear the last set bit */
        scv.ClearNextSetBitFromTail();

        /* Check if the result is expected */
        const uint32_t reinterp_val1 = scv.ReinterpretDataAs<uint32_t>();
        cmc::ExpectTrue(reinterp_val1 == uint32_t{1212});

        /* We update the tail bit count by skipping the tail bits until the next set bit */
        scv.UpdateTailBitCount();

        /* Check if the value consists of the expected significant bits */
        const int num_signinficant_bits1 = scv.GetCountOfSignificantBits();
        cmc::ExpectTrue(num_signinficant_bits1 == 30);

        /* Toggle the next set bits */
        scv.ToggleTailUntilNextUnsetBit();

        /* Check if the result is expected */
        const uint32_t reinterp_val2 = scv.ReinterpretDataAs<uint32_t>();
        cmc::ExpectTrue(reinterp_val2 == uint32_t{1216});

        /* Update to the front bit count to the first set bit */
        scv.UpdateFrontBitCount();

        /* Check if the value consists of the expected significant bits */
        const int num_signinficant_bits2 = scv.GetCountOfSignificantBits();
        cmc::ExpectTrue(num_signinficant_bits2 == 5);

        cmc::ExpectTrue(not scv.IsEmpty());

        /* Get the significant bits in big endian order */
        const std::vector<uint8_t> bits = scv.GetSignificantBitsInBigEndianOrdering();

        /* Check if the bits are as expected */
        cmc::ExpectTrue(bits.size() == 1);
        cmc::ExpectTrue(bits.front() == uint8_t{0b10011000});
        cmc::ExpectTrue(bits.front() == uint8_t{152});

        /* Define some sufffix bits that will be appended to the serialized compression value */
        const std::vector<uint8_t> suffix{uint8_t{0b11001000}};
        scv.ApplySuffix(suffix, 5);

        /* Check if the value consists of the expected significant bits */
        const int num_signinficant_bits3 = scv.GetCountOfSignificantBits();
        cmc::ExpectTrue(num_signinficant_bits3 == 10);

        /* Check if the result is expected */
        const uint32_t reinterp_val3 = scv.ReinterpretDataAs<uint32_t>();
        cmc::ExpectTrue(reinterp_val3 == uint32_t{1266});

        /* Get the significant bits in big endian order */
        const std::vector<uint8_t> bits2 = scv.GetSignificantBitsInBigEndianOrdering();

        /* Check if the bits are as expected */
        cmc::ExpectTrue(bits2.size() == 2);
        cmc::ExpectTrue(bits2.front() == uint8_t{0b10011110});
        cmc::ExpectTrue(bits2.back() == uint8_t{0b01000000});

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
