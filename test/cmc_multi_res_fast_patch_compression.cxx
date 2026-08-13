#include "cmc.hxx"
#include "test/cmc_test.hxx"
#include "test/cmc_test_patch_compression_util.hxx"
#include "patch/lossless/cmc_fast_multi_res_extraction.hxx"
#include "patch/lossless/cmc_fast_multi_res_decompression.hxx"

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

    /* Check the corrected mean computation */
    std::vector<float> values{1.3f, 2.42f, 1.714f, 2.12f};
    std::vector<float> partial_values{values[0], values[1], values[2]};
    const float N = static_cast<float>(values.size());
    const float partial_sum = cmc::patch::lossless::multi_res::NeumaierKahanBabuskaSum<float>(partial_values);
    const float value_to_match = values[3];
    const float mmean = cmc::patch::lossless::multi_res::FindMatchingMeanValue(N, partial_sum, value_to_match);
    const float eval_mmean = cmc::patch::lossless::multi_res::ComputeFMAForImplicitValue(N, mmean, partial_sum);
    cmc::ExpectTrue(std::bit_cast<uint32_t>(value_to_match) == std::bit_cast<uint32_t>(eval_mmean));

    //2D Test case
    {
    cmc::cmc_global_msg("2D Fast Patch-based Compression Test");
    constexpr int32_t DIM = 2;
    constexpr std::array<int32_t, DIM> kDimLength{32, 32};
    const std::vector<float> init_data = cmc::test::GenerateExampleData_2D<kDimLength[0], kDimLength[1]>();

    const std::span<const float> initial_data(init_data);

    const std::string file_name("cmc_out_test_fast_multi_res_patch_lossless_serial.cmc");

    cmc::patch::lossless::multi_res::fast::CompressionVariable<float, DIM> compression_variable(initial_data, kDimLength);
 
    compression_variable.Compress();

    compression_variable.WriteCompressedData(file_name.c_str());

    cmc::patch::lossless::multi_res::fast::ReadCompressionInfo(file_name);
  
    cmc::patch::lossless::multi_res::fast::DecompressionVariable<float, DIM> decompression_variable(file_name);
 
    decompression_variable.Decompress();

    std::vector<float> decompressed_data = decompression_variable.GetDecompressedData();
    
    const int32_t num_elems = decompressed_data.size();

    constexpr float max_error_tol = 5E-07f;

    for (int32_t elem_idx{0}; elem_idx < num_elems; ++elem_idx)
    {
        if (init_data[elem_idx] != decompressed_data[elem_idx])
        {
            cmc::cmc_debug_msg("init_data[", elem_idx, "] = ", init_data[elem_idx], ", decompressed[", elem_idx, "] = ", decompressed_data[elem_idx]);
            cmc::cmc_debug_msg("Init Bitset: ", std::bitset<32>(std::bit_cast<uint32_t>(init_data[elem_idx])), ", Decompressed Bitset: ", std::bitset<32>(std::bit_cast<uint32_t>(decompressed_data[elem_idx])));
            cmc::cmc_debug_msg("Abs Residual: ", std::abs(init_data[elem_idx] - decompressed_data[elem_idx]));
        }
        cmc::ExpectTrue(std::abs(init_data[elem_idx] - decompressed_data[elem_idx]) <= max_error_tol);
    }
    }

    //3D Test case
    {
    cmc::cmc_global_msg("3D Fast Patch-based Compression Test");
    constexpr int32_t DIM = 3;
    constexpr std::array<int32_t, DIM> kDimLength{128,128,128};
    const std::vector<float> init_data = cmc::test::GenerateExampleData_3D<kDimLength[0], kDimLength[1], kDimLength[2]>();
    const std::span<const float> initial_data(init_data);

    const std::string file_name("cmc_out_test_fast_multi_res_patch_lossless_serial.cmc");

    cmc::patch::lossless::multi_res::fast::CompressionVariable<float, DIM> compression_variable(initial_data, kDimLength);
 
    compression_variable.Compress();

    compression_variable.WriteCompressedData(file_name.c_str());

    cmc::patch::lossless::multi_res::fast::ReadCompressionInfo(file_name);
  
    cmc::patch::lossless::multi_res::fast::DecompressionVariable<float, DIM> decompression_variable(file_name);
 
    decompression_variable.Decompress();

    std::vector<float> decompressed_data = decompression_variable.GetDecompressedData();
    
    const int32_t num_elems = decompressed_data.size();

    constexpr float max_error_tol = 1E-05f;

    for (int32_t elem_idx{0}; elem_idx < num_elems; ++elem_idx)
    {
        if (init_data[elem_idx] != decompressed_data[elem_idx])
        {
            cmc::cmc_debug_msg("init_data[", elem_idx, "] = ", init_data[elem_idx], ", decompressed[", elem_idx, "] = ", decompressed_data[elem_idx]);
            cmc::cmc_debug_msg("Init Bitset: ", std::bitset<32>(std::bit_cast<uint32_t>(init_data[elem_idx])), ", Decompressed Bitset: ", std::bitset<32>(std::bit_cast<uint32_t>(decompressed_data[elem_idx])));
            cmc::cmc_debug_msg("Abs Residual: ", std::abs(init_data[elem_idx] - decompressed_data[elem_idx]));
        }
        cmc::ExpectTrue(std::abs(init_data[elem_idx] - decompressed_data[elem_idx]) <= max_error_tol);
    }
    }


    }
    /* Finalize cmc */
    cmc::CmcFinalize();

    return 0;
}
