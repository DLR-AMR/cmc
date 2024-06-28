#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "utilities/cmc_prefix.hxx"
#include "t8code/cmc_t8_prefix_encoding.hxx"
#include "utilities/cmc_span.hxx"

#include <iostream>
#include <bitset>

int main(void)
{
    /* Initialize cmc */
    cmc_initialize();
    
    {

    return cmc::CMC_TEST_SKIP;

    float test_float1 = 5.877471754111437539843683E-39;
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&test_float1)) <<" fuer " << test_float1 << std::endl;
    cmc::Prefix<sizeof(float)> bin_float1(test_float1);

    std::cout << "Num leading zeros: " << bin_float1.GetNumberLeadingZeros() << std::endl;
    std::cout << "Num trailing zeros: " << bin_float1.GetNumberTrailingZeros() << std::endl;

    //float test_float2 = 10.625;
    //float test_float3 = 12.625;
    //cmc::Prefix<sizeof(float)> bin_float2(test_float2);
    //cmc::Prefix<sizeof(float)> bin_float3(test_float3);
    //cmc::CommonPrefix<4> cp = cmc::GetCommonPrefix(bin_float2, bin_float3);
//
    //std::cout << "Common prefix length: " << cp.bit_length << "It reads: " << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&cp.prefix)) << std::endl;
    
    float test_float4 = 0.036186754703521728515625;
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&test_float4)) <<" fuer " << test_float4 << std::endl;

    cmc::CompressionValue<4> cv(test_float4);

    cv.ToggleTailUntilNextUnsetBit();
    float ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;
    
    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;

    cv.ToggleTailUntilNextUnsetBit();
    ret_val = cv.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val)) <<" fuer " << ret_val << std::endl;


    std::cout << "\n\n Jetzt CLear Bits\n" << std::endl;


    float test_float5 = 0.036186754703521728515625;
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&test_float5)) <<" fuer " << test_float5 << std::endl;

    cmc::CompressionValue<4> cv2(test_float5);
    float ret_val2;
    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;

    cv2.ClearNextSetBitFromTail();
    ret_val2 = cv2.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val2)) <<" fuer " << ret_val2 << std::endl;


    float test_float6 = 0.150634765625;
    cmc::CompressionValue<4> cv3(test_float6);
    float ret_val3 = cv3.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val3)) <<" fuer " << ret_val3 << std::endl;    

    cv3.UpdateTrailBitCount();
    cv3.UpdateFrontBitCount();
    std::vector<uint8_t> vec = cv3.GetSignificantBitsInBigEndianOrdering();
    
    for (auto iter = vec.begin(); iter != vec.end(); ++iter)
    {
        unsigned val = *iter;
        std::cout << std::bitset<8*sizeof(val)>(val) << std::endl;
    }
#if 0
    std::cout << "Fuer cv2:"<< std::endl;
    std::vector<uint8_t> vec2 = cv2.GetSignificantBitsInBigEndianOrdering();
    
    for (auto iter = vec2.begin(); iter != vec2.end(); ++iter)
    {
        unsigned val = *iter;
        std::cout << std::bitset<8*sizeof(val)>(val) << std::endl;
    }

    float test_float7 = 8.644016075870922788090416E-39;// 9.330486409651907094501846E-38;
    cmc::CompressionValue<4> cv4(test_float7);

    float ret_val4 = cv4.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val4)) <<" fuer " << ret_val4 << std::endl; 

    cv4.UpdateTrailBitCount();
    cv4.UpdateFrontBitCount();

    std::vector<uint8_t> vec3 = cv4.GetSignificantBitsInBigEndianOrdering();
    
    for (auto iter = vec3.begin(); iter != vec3.end(); ++iter)
    {
        unsigned val = *iter;
        std::cout << std::bitset<8*sizeof(val)>(val) << std::endl;
    }
    

    float test_float19 = -267.0;
    cmc::CompressionValue<4> cv19(test_float19);

    float ret_val19 = cv19.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val19)) <<" fuer " << ret_val19 << std::endl;

    for (int h = 0; h < 4; ++h)
    {
        std::cout << "index: " << h << unsigned(cv19[h])<< std::endl;
    }

    std::cout << "After modificatopn: " << std::endl;

    cv19[0] = (uint8_t) 0xFF;
    cv19[1] = (uint8_t) 0xF0;
    for (int h = 0; h < 4; ++h)
    {
        std::cout << "index: " << h << unsigned(cv19[h])<< std::endl;
    }

    cv19.ToggleTailUntilNextUnsetBit();

    ret_val19 = cv19.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val19)) <<" fuer " << ret_val19 << std::endl;
    std::cout << "New test" << std::endl;

    std::cout << "OFFSET: 0" << std::endl;
    std::vector<uint8_t> sfbo = cv19.GetSignificantBitsInBigEndianOrdering();
    uint64_t iiii{0};
    std::memcpy(&iiii, sfbo.data(), sfbo.size());
    std::cout << std::bitset<8*sizeof(uint64_t)>(iiii) <<" fuer " << ret_val19 << " size of vec: " << sfbo.size() << std::endl;
    sfbo.clear();

    std::cout << "OFFSET: 2" << std::endl;
     sfbo = cv19.GetSignificantBitsInBigEndianOrdering(2);
    iiii = 0;
    std::memcpy(&iiii, sfbo.data(), sfbo.size());
    std::cout << std::bitset<8*sizeof(uint64_t)>(iiii) <<" fuer " << ret_val19 << " size of vec: " << sfbo.size() << std::endl;
    sfbo.clear();

    std::cout << "OFFSET: 4" << std::endl;
     sfbo = cv19.GetSignificantBitsInBigEndianOrdering(4);
    iiii = 0;
    std::memcpy(&iiii, sfbo.data(), sfbo.size());
    std::cout << std::bitset<8*sizeof(uint64_t)>(iiii) <<" fuer " << ret_val19 << " size of vec: " << sfbo.size() << std::endl;
    sfbo.clear();

    std::cout << "OFFSET: 6" << std::endl;
    sfbo = cv19.GetSignificantBitsInBigEndianOrdering(6);
    iiii = 0;
    std::memcpy(&iiii, sfbo.data(), sfbo.size());
    std::cout << std::bitset<8*sizeof(uint64_t)>(iiii) <<" fuer " << ret_val19 << " size of vec: " << sfbo.size() << std::endl;
    sfbo.clear();

    std::cout << "OFFSET: 7" << std::endl;
    sfbo = cv19.GetSignificantBitsInBigEndianOrdering(7);
    iiii = 0;
    std::memcpy(&iiii, sfbo.data(), sfbo.size());
    std::cout << std::bitset<8*sizeof(uint64_t)>(iiii) <<" fuer " << ret_val19 << " size of vec: " << sfbo.size() << std::endl;

    std::cout << "OFFSETed and strided " << std::endl;
    sfbo = cv19.GetSignificantBitsInBigEndianOrdering_Offset_4_Stride_12();
    iiii = 0;
    std::memcpy(&iiii, sfbo.data(), sfbo.size());
    std::cout << std::bitset<8*sizeof(uint64_t)>(iiii) <<" fuer " << ret_val19 << " size of vec: " << sfbo.size() << std::endl;

    cv19[0] = (uint8_t) 0x00;
    cv19[1] = (uint8_t) 0xF0;
    cv19[2] = (uint8_t) 0xFF;
    cv19[3] = (uint8_t) 0x00;

    cv19.UpdateTrailBitCount();
    cv19.UpdateFrontBitCount();

    ret_val19 = cv19.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val19)) <<" fuer " << ret_val19 << std::endl;
    std::cout << "New test" << std::endl;

    std::cout << "OFFSET: 0" << std::endl;
    sfbo = cv19.GetSignificantBitsInBigEndianOrdering();
    iiii = 0;
    std::memcpy(&iiii, sfbo.data(), sfbo.size());
    std::cout << std::bitset<8*sizeof(uint64_t)>(iiii) <<" fuer " << ret_val19 << " size of vec: " << sfbo.size() << std::endl;
    sfbo.clear();

    std::cout << "OFFSETed and strided22 " << std::endl;
    sfbo = cv19.GetSignificantBitsInBigEndianOrdering_Offset_4_Stride_12();
    iiii = 0;
    std::memcpy(&iiii, sfbo.data(), sfbo.size());
    std::cout << std::bitset<8*sizeof(uint64_t)>(iiii) <<" fuer " << ret_val19 << " size of vec: " << sfbo.size() << std::endl;

#endif
//00000000000000000000000000000000 01011111 00111000 00111100 00001100
//-------------------------------- 11000011 10000101 11110000 11111111


    cmc::CompressionValue<4> cv_test_33;
    std::vector<uint8_t> preffff{0xF0};
    std::vector<uint8_t> preffff2{0xF1, 0x00, 0x00, 0x00};
    uint32_t jjj{0};
    cv_test_33.ApplyPrefix(preffff, 6);
    float ret_val33 = cv_test_33.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val33)) << std::endl;

    cv_test_33.ApplyPrefix(preffff2, 12);
    ret_val33 = cv_test_33.ReinterpretDataAs<float>();
    std::cout << std::bitset<8*sizeof(float)>(*reinterpret_cast<uint32_t*>(&ret_val33)) << std::endl;


    std::cout << "Test for cmc_t8_prefix_encoding" << std::endl;

    //Length: 11100101  00000111  (Resulting from forward: 01 1001 1111 0001) 01100111 11000100   01 00001001 1100
    //Encoded Prefixes: 00010111 00011111 (Resulting from forawrd: 11 0101 000011111) 
    std::vector<uint8_t> pref_length{0x67, 0xC4};
    cmc::VectorView<uint8_t> pl_vv(pref_length.data(), pref_length.size());
    std::vector<uint8_t> pref_encodings{0x17, 0x1F};
    cmc::VectorView<uint8_t> pe_vv(pref_encodings.data(), pref_encodings.size());

    std::cout << "Size of pref_length: " << pref_length.size() << " und size of pref_encodings: " << pref_encodings.size() << std::endl;

    int length_bit_position = CHAR_BIT;
    int length_byte_position = 0;
    int encoding_bit_position = 0;
    int encoding_byte_position = 0;

    auto [preff1, preff1_length] = cmc::GetNextPrefix(pl_vv, length_byte_position, length_bit_position, pe_vv, encoding_byte_position, encoding_bit_position);
    uint16_t fffff = 0;
    std::memcpy(&fffff, preff1.data(), preff1.size());
    std::cout << std::bitset<8*sizeof(uint16_t)>(fffff) << " size of vec: " << preff1.size() << std::endl;
    std::cout << "Size of prefix: " << preff1_length << std::endl;
    std::cout << "Nach GetPrefix: length_byte_position: " << length_byte_position << ", length_bit_position: " << length_bit_position << ", encoding_byte_position: " << encoding_byte_position << ", encoding_bit_position: " << encoding_bit_position << std::endl; 
    auto [preff2, preff2_length] = cmc::GetNextPrefix(pl_vv, length_byte_position, length_bit_position, pe_vv, encoding_byte_position, encoding_bit_position);
    fffff = 0;
    std::memcpy(&fffff, preff2.data(), preff2.size());
    std::cout << std::bitset<8*sizeof(uint16_t)>(fffff) << " size of vec: " << preff2.size() << std::endl;
    std::cout << "Size of prefix: " << preff2_length << std::endl;
    std::cout << "Nach GetPrefix: length_byte_position: " << length_byte_position << ", length_bit_position: " << length_bit_position << ", encoding_byte_position: " << encoding_byte_position << ", encoding_bit_position: " << encoding_bit_position << std::endl;

    auto [preff3, preff3_length] = cmc::GetNextPrefix(pl_vv, length_byte_position, length_bit_position, pe_vv, encoding_byte_position, encoding_bit_position);
    fffff = 0;
    std::memcpy(&fffff, preff3.data(), preff3.size());
    std::cout << std::bitset<8*sizeof(uint16_t)>(fffff) << " size of vec: " << preff3.size() << std::endl;
    std::cout << "Size of prefix: " << preff3_length << std::endl;
    std::cout << "Nach GetPrefix: length_byte_position: " << length_byte_position << ", length_bit_position: " << length_bit_position << ", encoding_byte_position: " << encoding_byte_position << ", encoding_bit_position: " << encoding_bit_position << std::endl;


    //0111 1001   0010 1100 
    //Erwartung: 1 0000000|0 1111 00 |0
    std::vector<uint8_t> rle_encoding{0x79, 0x2C};
    cmc::VectorView<uint8_t> vvrle(rle_encoding.data(), rle_encoding.size());
    auto[decoded_rle, num_bits] = cmc::DecodeRunLengthEncoding(vvrle);

    cmc_debug_msg("Size of decoded rle: ", decoded_rle.size(), " and num Bits: ", (unsigned)num_bits);
    uint16_t dddddd = 0;
    std::memcpy(&dddddd, decoded_rle.data(), decoded_rle.size());
    std::cout << std::bitset<8*sizeof(uint16_t)>(dddddd) << std::endl;

    }

    /* Finalize cmc */
    cmc_finalize();

    return cmc::CMC_TEST_SUCCESS;
}
