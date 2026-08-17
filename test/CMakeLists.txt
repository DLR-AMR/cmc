if (CMC_WITH_T8CODE)
if (CMC_ENABLE_MPI)

# Add the test for the intra element lossless multi res compression with multiple data
set(CMC_TEST_MULTI_RES_LOSSLESS_MULTI_DATA_COMPRESSION cmc_test_multi_res_lossless_multi_data_compression)
add_executable(${CMC_TEST_MULTI_RES_LOSSLESS_MULTI_DATA_COMPRESSION} ${CMAKE_CURRENT_LIST_DIR}/cmc_multi_res_lossless_compression_multi_data.cxx)
target_link_libraries(${CMC_TEST_MULTI_RES_LOSSLESS_MULTI_DATA_COMPRESSION} PUBLIC cmc)
target_include_directories(${CMC_TEST_MULTI_RES_LOSSLESS_MULTI_DATA_COMPRESSION} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
if(CMC_ENABLE_ALL_WARNINGS)
    target_compile_options(${CMC_TEST_MULTI_RES_LOSSLESS_MULTI_DATA_COMPRESSION} PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME ${CMC_TEST_MULTI_RES_LOSSLESS_MULTI_DATA_COMPRESSION} COMMAND ${CMC_TEST_MULTI_RES_LOSSLESS_MULTI_DATA_COMPRESSION})

# Add the test for the intra element lossless multi res compression
set(CMC_TEST_MULTI_RES_LOSSLESS_INTRA_ELEM_CODING cmc_test_multi_res_lossless_intra_elem_coding)
add_executable(${CMC_TEST_MULTI_RES_LOSSLESS_INTRA_ELEM_CODING} ${CMAKE_CURRENT_LIST_DIR}/cmc_multi_res_lossless_intra_elem_coding.cxx)
target_link_libraries(${CMC_TEST_MULTI_RES_LOSSLESS_INTRA_ELEM_CODING} PUBLIC cmc)
target_include_directories(${CMC_TEST_MULTI_RES_LOSSLESS_INTRA_ELEM_CODING} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
if(CMC_ENABLE_ALL_WARNINGS)
    target_compile_options(${CMC_TEST_MULTI_RES_LOSSLESS_INTRA_ELEM_CODING} PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME ${CMC_TEST_MULTI_RES_LOSSLESS_INTRA_ELEM_CODING} COMMAND ${CMC_TEST_MULTI_RES_LOSSLESS_INTRA_ELEM_CODING})

# Add the test for the Trixi-related intra element lossy multi res compression with multiple data
set(CMC_TEST_TRIXI_MULTI_RES_LOSSY_MULTI_DATA_COMPRESSION cmc_test_trxi_multi_res_lossy_multi_data_compression)
add_executable(${CMC_TEST_TRIXI_MULTI_RES_LOSSY_MULTI_DATA_COMPRESSION} ${CMAKE_CURRENT_LIST_DIR}/cmc_trixi_multi_res_lossy_compression_multi_data.cxx)
target_link_libraries(${CMC_TEST_TRIXI_MULTI_RES_LOSSY_MULTI_DATA_COMPRESSION} PUBLIC cmc)
target_include_directories(${CMC_TEST_TRIXI_MULTI_RES_LOSSY_MULTI_DATA_COMPRESSION} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_TRIXI_MULTI_RES_LOSSY_MULTI_DATA_COMPRESSION} COMMAND ${CMC_TEST_TRIXI_MULTI_RES_LOSSY_MULTI_DATA_COMPRESSION})

# Add the test for the Trixi-related element lossy multi res compression with single data per element
set(CMC_TEST_TRIXI_MULTI_RES_LOSSY_SINGLE_DATA_COMPRESSION cmc_test_trxi_multi_res_lossy_single_data_compression)
add_executable(${CMC_TEST_TRIXI_MULTI_RES_LOSSY_SINGLE_DATA_COMPRESSION} ${CMAKE_CURRENT_LIST_DIR}/cmc_trixi_multi_res_lossy_compression_single_data.cxx)
target_link_libraries(${CMC_TEST_TRIXI_MULTI_RES_LOSSY_SINGLE_DATA_COMPRESSION} PUBLIC cmc)
target_include_directories(${CMC_TEST_TRIXI_MULTI_RES_LOSSY_SINGLE_DATA_COMPRESSION} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_TRIXI_MULTI_RES_LOSSY_SINGLE_DATA_COMPRESSION} COMMAND ${CMC_TEST_TRIXI_MULTI_RES_LOSSY_SINGLE_DATA_COMPRESSION})

# Add the test for the lossy multi res compression based on RBF interpolation (with single data per element)
set(CMC_TEST_MULTI_RES_RBF_LOSSY_SINGLE_DATA_COMPRESSION cmc_test_multi_res_rbf_lossy_compression)
add_executable(${CMC_TEST_MULTI_RES_RBF_LOSSY_SINGLE_DATA_COMPRESSION} ${CMAKE_CURRENT_LIST_DIR}/cmc_multi_res_rbf_lossy_compression.cxx)
target_link_libraries(${CMC_TEST_MULTI_RES_RBF_LOSSY_SINGLE_DATA_COMPRESSION} PUBLIC cmc)
target_include_directories(${CMC_TEST_MULTI_RES_RBF_LOSSY_SINGLE_DATA_COMPRESSION} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_MULTI_RES_RBF_LOSSY_SINGLE_DATA_COMPRESSION} COMMAND ${CMC_TEST_MULTI_RES_RBF_LOSSY_SINGLE_DATA_COMPRESSION})

# Add the test for the lossy multi res compression based on IDW interpolation (with single data per element) on a hybrid mesh
set(CMC_TEST_MULTI_RES_IDW_LOSSY_SINGLE_DATA_COMPRESSION_HYBRID cmc_test_multi_res_idw_lossy_compression_hybrid_mesh)
add_executable(${CMC_TEST_MULTI_RES_IDW_LOSSY_SINGLE_DATA_COMPRESSION_HYBRID} ${CMAKE_CURRENT_LIST_DIR}/cmc_multi_res_idw_lossy_compression_hybrid_mesh.cxx)
target_link_libraries(${CMC_TEST_MULTI_RES_IDW_LOSSY_SINGLE_DATA_COMPRESSION_HYBRID} PUBLIC cmc)
target_include_directories(${CMC_TEST_MULTI_RES_IDW_LOSSY_SINGLE_DATA_COMPRESSION_HYBRID} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_MULTI_RES_IDW_LOSSY_SINGLE_DATA_COMPRESSION_HYBRID} COMMAND ${CMC_TEST_MULTI_RES_IDW_LOSSY_SINGLE_DATA_COMPRESSION_HYBRID})

endif()
endif()

# Add the test for the bits vector
set(CMC_TEST_BITS_VECTOR cmc_test_bits_vector)
add_executable(${CMC_TEST_BITS_VECTOR} ${CMAKE_CURRENT_LIST_DIR}/cmc_bits_vector.cxx)
target_link_libraries(${CMC_TEST_BITS_VECTOR} PUBLIC cmc)
target_include_directories(${CMC_TEST_BITS_VECTOR} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_BITS_VECTOR} COMMAND ${CMC_TEST_BITS_VECTOR})

# Add the test for the bits vector view
set(CMC_TEST_BITS_VECTOR_VIEW cmc_test_bits_vector_view)
add_executable(${CMC_TEST_BITS_VECTOR_VIEW} ${CMAKE_CURRENT_LIST_DIR}/cmc_bits_vector_view.cxx)
target_link_libraries(${CMC_TEST_BITS_VECTOR_VIEW} PUBLIC cmc)
target_include_directories(${CMC_TEST_BITS_VECTOR_VIEW} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
if(CMC_ENABLE_ALL_WARNINGS)
    target_compile_options(${CMC_TEST_BITS_VECTOR_VIEW} PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME ${CMC_TEST_BITS_VECTOR_VIEW} COMMAND ${CMC_TEST_BITS_VECTOR_VIEW})

# Add the test for the Huffman coder
set(CMC_TEST_HUFFMAN_CODER cmc_test_huffman_coder)
add_executable(${CMC_TEST_HUFFMAN_CODER} ${CMAKE_CURRENT_LIST_DIR}/cmc_huffman_coder.cxx)
target_link_libraries(${CMC_TEST_HUFFMAN_CODER} PUBLIC cmc)
target_include_directories(${CMC_TEST_HUFFMAN_CODER} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
if(CMC_ENABLE_ALL_WARNINGS)
    target_compile_options(${CMC_TEST_HUFFMAN_CODER} PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME ${CMC_TEST_HUFFMAN_CODER} COMMAND ${CMC_TEST_HUFFMAN_CODER})

# Add the test for the lossy multi-res patch-wise IDW compression
set(CMC_TEST_MULTI_RES_LOSSY_IDW_PATCH_COMPRESSION cmc_test_multi_res_idw_lossy_patch_compression)
add_executable(${CMC_TEST_MULTI_RES_LOSSY_IDW_PATCH_COMPRESSION} ${CMAKE_CURRENT_LIST_DIR}/cmc_multi_res_lossy_idw_patch_compression.cxx)
target_link_libraries(${CMC_TEST_MULTI_RES_LOSSY_IDW_PATCH_COMPRESSION} PUBLIC cmc)
target_include_directories(${CMC_TEST_MULTI_RES_LOSSY_IDW_PATCH_COMPRESSION} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
if(CMC_ENABLE_ALL_WARNINGS)
    target_compile_options(${CMC_TEST_MULTI_RES_LOSSY_IDW_PATCH_COMPRESSION} PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME ${CMC_TEST_MULTI_RES_LOSSY_IDW_PATCH_COMPRESSION} COMMAND ${CMC_TEST_MULTI_RES_LOSSY_IDW_PATCH_COMPRESSION})

# Add the test for the lossy multi-res patch-wise RBF compression
set(CMC_TEST_MULTI_RES_LOSSY_RBF_PATCH_COMPRESSION cmc_test_multi_res_rbf_lossy_patch_compression)
add_executable(${CMC_TEST_MULTI_RES_LOSSY_RBF_PATCH_COMPRESSION} ${CMAKE_CURRENT_LIST_DIR}/cmc_multi_res_lossy_rbf_patch_compression.cxx)
target_link_libraries(${CMC_TEST_MULTI_RES_LOSSY_RBF_PATCH_COMPRESSION} PUBLIC cmc)
target_include_directories(${CMC_TEST_MULTI_RES_LOSSY_RBF_PATCH_COMPRESSION} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
if(CMC_ENABLE_ALL_WARNINGS)
    target_compile_options(${CMC_TEST_MULTI_RES_LOSSY_RBF_PATCH_COMPRESSION} PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME ${CMC_TEST_MULTI_RES_LOSSY_RBF_PATCH_COMPRESSION} COMMAND ${CMC_TEST_MULTI_RES_LOSSY_RBF_PATCH_COMPRESSION})

# Add the test for the lossless multi-res patch-wise compression
set(CMC_TEST_MULTI_RES_LOSSLESS_PATCH_COMPRESSION cmc_test_multi_res_lossless_patch_compression)
add_executable(${CMC_TEST_MULTI_RES_LOSSLESS_PATCH_COMPRESSION} ${CMAKE_CURRENT_LIST_DIR}/cmc_multi_res_lossless_patch_compression.cxx)
target_link_libraries(${CMC_TEST_MULTI_RES_LOSSLESS_PATCH_COMPRESSION} PUBLIC cmc)
target_include_directories(${CMC_TEST_MULTI_RES_LOSSLESS_PATCH_COMPRESSION} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
if(CMC_ENABLE_ALL_WARNINGS)
    target_compile_options(${CMC_TEST_MULTI_RES_LOSSLESS_PATCH_COMPRESSION} PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME ${CMC_TEST_MULTI_RES_LOSSLESS_PATCH_COMPRESSION} COMMAND ${CMC_TEST_MULTI_RES_LOSSLESS_PATCH_COMPRESSION})

# Add the test for the fast multi-res patch-wise compression
set(CMC_TEST_MULTI_RES_FAST_PATCH_COMPRESSION cmc_test_multi_res_fast_patch_compression)
add_executable(${CMC_TEST_MULTI_RES_FAST_PATCH_COMPRESSION} ${CMAKE_CURRENT_LIST_DIR}/cmc_multi_res_fast_patch_compression.cxx)
target_link_libraries(${CMC_TEST_MULTI_RES_FAST_PATCH_COMPRESSION} PUBLIC cmc)
target_include_directories(${CMC_TEST_MULTI_RES_FAST_PATCH_COMPRESSION} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
if(CMC_ENABLE_ALL_WARNINGS)
    target_compile_options(${CMC_TEST_MULTI_RES_FAST_PATCH_COMPRESSION} PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME ${CMC_TEST_MULTI_RES_FAST_PATCH_COMPRESSION} COMMAND ${CMC_TEST_MULTI_RES_FAST_PATCH_COMPRESSION})

# Add the test for the lossless multi-res PCP4 patch-wise compression
set(CMC_TEST_MULTI_RES_PCP4_LOSSLESS_PATCH_COMPRESSION cmc_test_multi_res_pcp4_lossless_patch_compression)
add_executable(${CMC_TEST_MULTI_RES_PCP4_LOSSLESS_PATCH_COMPRESSION} ${CMAKE_CURRENT_LIST_DIR}/cmc_multi_res_pcp4_lossless_patch_compression.cxx)
target_link_libraries(${CMC_TEST_MULTI_RES_PCP4_LOSSLESS_PATCH_COMPRESSION} PUBLIC cmc)
target_include_directories(${CMC_TEST_MULTI_RES_PCP4_LOSSLESS_PATCH_COMPRESSION} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
if(CMC_ENABLE_ALL_WARNINGS)
    target_compile_options(${CMC_TEST_MULTI_RES_PCP4_LOSSLESS_PATCH_COMPRESSION} PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME ${CMC_TEST_MULTI_RES_PCP4_LOSSLESS_PATCH_COMPRESSION} COMMAND ${CMC_TEST_MULTI_RES_PCP4_LOSSLESS_PATCH_COMPRESSION})
