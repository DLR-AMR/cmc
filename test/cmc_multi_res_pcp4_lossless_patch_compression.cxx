#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "test/cmc_test_patch_compression_util.hxx"
#include "patch/lossless/cmc_multi_res_pcp4_extraction.hxx"
#include "patch/lossless/cmc_multi_res_pcp4_decompression.hxx"

#include <vector>
#include <cmath>
#include <numeric>
#include <cstdint>
#include <array>
#include <span>
#include <string>
#include <bitset>

int
main(void)
{
    /* Initialize cmc */
    cmc::CmcInitialize();
    {

    //2D Test case
    {
    cmc::cmc_global_msg("2D Lossless Patch-based Compression Test");
    constexpr int32_t DIM = 2;
    constexpr std::array<int32_t, DIM> kDimLength{41, 57};
    const std::vector<float> init_data = cmc::test::GenerateExampleData_2D<kDimLength[0], kDimLength[1]>();
    const std::span<const float> initial_data(init_data);

    const std::string file_name("cmc_out_test_multi_res_patch_lossless_pcp4_serial.cmc");

    cmc::patch::lossless::multi_res::pcp4::CompressionVariable<float, DIM> compression_variable(initial_data, kDimLength);
 
    compression_variable.Compress();

    compression_variable.WriteCompressedData(file_name.c_str());

    cmc::patch::lossless::multi_res::pcp4::ReadCompressionInfo(file_name);
  
    cmc::patch::lossless::multi_res::pcp4::DecompressionVariable<float, DIM> decompression_variable(file_name);
 
    decompression_variable.Decompress();

    std::vector<float> decompressed_data = decompression_variable.GetDecompressedData();
    
    const int32_t num_elems = decompressed_data.size();

    for (int32_t elem_idx{0}; elem_idx < num_elems; ++elem_idx)
    {
        if (init_data[elem_idx] != decompressed_data[elem_idx])
        {
            cmc::cmc_debug_msg("init_data[", elem_idx, "] = ", init_data[elem_idx], ", decompressed[", elem_idx, "] = ", decompressed_data[elem_idx]);
            cmc::cmc_debug_msg("Init Bitset: ", std::bitset<32>(std::bit_cast<uint32_t>(init_data[elem_idx])), ", Decompressed Bitset: ", std::bitset<32>(std::bit_cast<uint32_t>(decompressed_data[elem_idx])));
        }

        cmc::ExpectTrue(init_data[elem_idx] == decompressed_data[elem_idx]);
    }

    }

    //3D Test case
    {
    cmc::cmc_global_msg("3D Lossless Patch-based Compression Test");
    constexpr int32_t DIM = 3;
    constexpr std::array<int32_t, DIM> kDimLength{128,128,128};
    const std::vector<float> init_data = cmc::test::GenerateExampleData_3D<kDimLength[0], kDimLength[1], kDimLength[2]>();
    const std::span<const float> initial_data(init_data);

    const std::string file_name("cmc_out_test_multi_res_patch_lossless_pcp4_serial.cmc");

    cmc::patch::lossless::multi_res::pcp4::CompressionVariable<float, DIM> compression_variable(initial_data, kDimLength);
 
    compression_variable.Compress();

    compression_variable.WriteCompressedData(file_name.c_str());

    cmc::patch::lossless::multi_res::pcp4::ReadCompressionInfo(file_name);

    cmc::patch::lossless::multi_res::pcp4::DecompressionVariable<float, DIM> decompression_variable(file_name);
 
    decompression_variable.Decompress();

    std::vector<float> decompressed_data = decompression_variable.GetDecompressedData();
    
    const int32_t num_elems = decompressed_data.size();

    for (int32_t elem_idx{0}; elem_idx < num_elems; ++elem_idx)
    {
        if (init_data[elem_idx] != decompressed_data[elem_idx])
        {
            cmc::cmc_debug_msg("init_data[", elem_idx, "] = ", init_data[elem_idx], ", decompressed[", elem_idx, "] = ", decompressed_data[elem_idx]);
            cmc::cmc_debug_msg("Init Bitset: ", std::bitset<32>(std::bit_cast<uint32_t>(init_data[elem_idx])), ", Decompressed Bitset: ", std::bitset<32>(std::bit_cast<uint32_t>(decompressed_data[elem_idx])));
        }

        cmc::ExpectTrue(init_data[elem_idx] == decompressed_data[elem_idx]);
    }

    }

    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
