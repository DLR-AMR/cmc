#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_bytes.hxx"

#include <cstddef>
#include <cstring>
#include <array>

int main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    
    {
        /* Define an integer value */
        const uint32_t value1 = 234;

        /* Get the value in little as well as big endian representation */
        std::array<uint8_t, sizeof(uint32_t)> bvalue1_little_endian = cmc::SerializeValue(value1, cmc::Endian::Little);
        std::array<uint8_t, sizeof(uint32_t)> bvalue1_big_endian = cmc::SerializeValue(value1, cmc::Endian::Big);

        /* Deserialize the values and compare them in the native endianness */
        std::array<uint8_t, sizeof(uint32_t)> deserialized_value1_le = cmc::DeserializeValue<uint32_t>(bvalue1_little_endian.data(), cmc::Endian::Little);
        std::array<uint8_t, sizeof(uint32_t)> deserialized_value1_be = cmc::DeserializeValue<uint32_t>(bvalue1_big_endian.data(), cmc::Endian::Big);

        /* Copy the ordeed bytes to an actual type of the given data */
        uint32_t value1_from_le;
        std::memcpy(&value1_from_le, deserialized_value1_le.data(), sizeof(uint32_t));

        uint32_t value1_from_be;
        std::memcpy(&value1_from_be, deserialized_value1_be.data(), sizeof(uint32_t));

        /* Comapre them for equality with the initial value */
        cmc::ExpectTrue(value1 == value1_from_le);
        cmc::ExpectTrue(value1 == value1_from_be);

        /* Directly reconstruct the values from a pointer and the serialized endianness */
        const uint32_t reconstructed_value1_le = cmc::ReconstructSerializedValue<uint32_t>(bvalue1_little_endian.data(), cmc::Endian::Little);
        const uint32_t reconstructed_value1_be = cmc::ReconstructSerializedValue<uint32_t>(bvalue1_big_endian.data(), cmc::Endian::Big);

        /* Comapre them for equality with the initial value */
        cmc::ExpectTrue(value1 == reconstructed_value1_le);
        cmc::ExpectTrue(value1 == reconstructed_value1_be);


        /* Define a double  value */
        const double value2 = 234.0;

        /* Get the value in little as well as big endian representation */
        std::array<uint8_t, sizeof(double)> bvalue2_little_endian = cmc::SerializeValue(value2, cmc::Endian::Little);
        std::array<uint8_t, sizeof(double)> bvalue2_big_endian = cmc::SerializeValue(value2, cmc::Endian::Big);

        /* Deserialize the values and compare them in the native endianness */
        std::array<uint8_t, sizeof(double)> deserialized_value2_le = cmc::DeserializeValue<double>(bvalue2_little_endian.data(), cmc::Endian::Little);
        std::array<uint8_t, sizeof(double)> deserialized_value2_be = cmc::DeserializeValue<double>(bvalue2_big_endian.data(), cmc::Endian::Big);

        /* Copy the ordeed bytes to an actual type of the given data */
        double value2_from_le;
        std::memcpy(&value2_from_le, deserialized_value2_le.data(), sizeof(double));

        double value2_from_be;
        std::memcpy(&value2_from_be, deserialized_value2_be.data(), sizeof(double));

        /* Comapre them for equality with the initial value */
        cmc::ExpectTrue(value2 == value2_from_le);
        cmc::ExpectTrue(value2 == value2_from_be);

        /* Directly reconstruct the values from a pointer and the serialized endianness */
        const double reconstructed_value2_le = cmc::ReconstructSerializedValue<double>(bvalue2_little_endian.data(), cmc::Endian::Little);
        const double reconstructed_value2_be = cmc::ReconstructSerializedValue<double>(bvalue2_big_endian.data(), cmc::Endian::Big);

        /* Comapre them for equality with the initial value */
        cmc::ExpectTrue(value2 == reconstructed_value2_le);
        cmc::ExpectTrue(value2 == reconstructed_value2_be);

    }

    /* Finalize cmc */
    cmc::CmcFinalize();

    return cmc::CMC_TEST_SUCCESS;
}
