if (CMC_WITH_T8CODE)

if (CMC_ENABLE_MPI)
# Add the test for adaptive coarsening
set(CMC_TEST_ADAPTIVE_COARSENING_VARIABLE cmc_test_adaptive_coarsening_variable)
add_executable(${CMC_TEST_ADAPTIVE_COARSENING_VARIABLE} ${CMAKE_CURRENT_LIST_DIR}/cmc_adaptive_coarsening_variable.cxx)
target_link_libraries(${CMC_TEST_ADAPTIVE_COARSENING_VARIABLE} PUBLIC cmc)
target_include_directories(${CMC_TEST_ADAPTIVE_COARSENING_VARIABLE} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_ADAPTIVE_COARSENING_VARIABLE} COMMAND ${CMC_TEST_ADAPTIVE_COARSENING_VARIABLE})

endif()

if (CMC_WITH_NETCDF AND CMC_ENABLE_MPI)

# Add the test for the lossless prefix compression
set(CMC_TEST_PREFIX_COMPRESSION_VARIABLE cmc_test_prefix_compression_variable)
add_executable(${CMC_TEST_PREFIX_COMPRESSION_VARIABLE} ${CMAKE_CURRENT_LIST_DIR}/cmc_prefix_compression_variable.cxx)
target_link_libraries(${CMC_TEST_PREFIX_COMPRESSION_VARIABLE} PUBLIC cmc)
target_include_directories(${CMC_TEST_PREFIX_COMPRESSION_VARIABLE} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_PREFIX_COMPRESSION_VARIABLE} COMMAND ${CMC_TEST_PREFIX_COMPRESSION_VARIABLE})

# Add the test for the lossless multi_res compression
set(CMC_TEST_MULTI_RES_COMPRESSION_VARIABLE cmc_test_multi_res_compression_variable)
add_executable(${CMC_TEST_MULTI_RES_COMPRESSION_VARIABLE} ${CMAKE_CURRENT_LIST_DIR}/cmc_multi_res_compression_variable.cxx)
target_link_libraries(${CMC_TEST_MULTI_RES_COMPRESSION_VARIABLE} PUBLIC cmc)
target_include_directories(${CMC_TEST_MULTI_RES_COMPRESSION_VARIABLE} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_MULTI_RES_COMPRESSION_VARIABLE} COMMAND ${CMC_TEST_MULTI_RES_COMPRESSION_VARIABLE})

endif()

endif()

# Add the test for the byte compression values
set(CMC_TEST_BYTE_COMPRESSION_VALUE cmc_test_byte_compression_value)
add_executable(${CMC_TEST_BYTE_COMPRESSION_VALUE} ${CMAKE_CURRENT_LIST_DIR}/cmc_byte_compression_value.cxx)
target_link_libraries(${CMC_TEST_BYTE_COMPRESSION_VALUE} PUBLIC cmc)
target_include_directories(${CMC_TEST_BYTE_COMPRESSION_VALUE} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_BYTE_COMPRESSION_VALUE} COMMAND ${CMC_TEST_BYTE_COMPRESSION_VALUE})

# Add the test for the bitmap
set(CMC_TEST_BITMAP cmc_test_bitmap)
add_executable(${CMC_TEST_BITMAP} ${CMAKE_CURRENT_LIST_DIR}/cmc_bitmap.cxx)
target_link_libraries(${CMC_TEST_BITMAP} PUBLIC cmc)
target_include_directories(${CMC_TEST_BITMAP} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_BITMAP} COMMAND ${CMC_TEST_BITMAP})

# Add the test for the bitmap-viewer
set(CMC_TEST_BITMAP_VIEW cmc_test_bitmap_view)
add_executable(${CMC_TEST_BITMAP_VIEW} ${CMAKE_CURRENT_LIST_DIR}/cmc_bitmap_view.cxx)
target_link_libraries(${CMC_TEST_BITMAP_VIEW} PUBLIC cmc)
target_include_directories(${CMC_TEST_BITMAP_VIEW} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_BITMAP_VIEW} COMMAND ${CMC_TEST_BITMAP_VIEW})

# Add the test for the bitvector-view
set(CMC_TEST_BITVECTOR_VIEW cmc_test_bitvector_view)
add_executable(${CMC_TEST_BITVECTOR_VIEW} ${CMAKE_CURRENT_LIST_DIR}/cmc_bitvector_view.cxx)
target_link_libraries(${CMC_TEST_BITVECTOR_VIEW} PUBLIC cmc)
target_include_directories(${CMC_TEST_BITVECTOR_VIEW} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_BITVECTOR_VIEW} COMMAND ${CMC_TEST_BITVECTOR_VIEW})

if (CMC_ENABLE_MPI)
# Add the test for the arithmetic encoder of the prefix compression
set(CMC_TEST_PREFIX_COMPRESSION_ARITHMETIC_CODING cmc_test_prefix_compression_arithmetic_coding)
add_executable(${CMC_TEST_PREFIX_COMPRESSION_ARITHMETIC_CODING} ${CMAKE_CURRENT_LIST_DIR}/cmc_prefix_compression_arithmetic_coding.cxx)
target_link_libraries(${CMC_TEST_PREFIX_COMPRESSION_ARITHMETIC_CODING} PUBLIC cmc)
target_include_directories(${CMC_TEST_PREFIX_COMPRESSION_ARITHMETIC_CODING} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_PREFIX_COMPRESSION_ARITHMETIC_CODING} COMMAND ${CMC_TEST_PREFIX_COMPRESSION_ARITHMETIC_CODING})

# Add the test for the arithmetic encoder of the multi_res compression
set(CMC_TEST_MULTI_RES_COMPRESSION_ARITHMETIC_CODING cmc_test_multi_res_compression_arithmetic_coding)
add_executable(${CMC_TEST_MULTI_RES_COMPRESSION_ARITHMETIC_CODING} ${CMAKE_CURRENT_LIST_DIR}/cmc_multi_res_compression_arithmetic_coding.cxx)
target_link_libraries(${CMC_TEST_MULTI_RES_COMPRESSION_ARITHMETIC_CODING} PUBLIC cmc)
target_include_directories(${CMC_TEST_MULTI_RES_COMPRESSION_ARITHMETIC_CODING} PRIVATE ${CMAKE_CURRENT_LIST_DIR} ${CMAKE_CURRENT_LIST_DIR}/..)
add_test(NAME ${CMC_TEST_MULTI_RES_COMPRESSION_ARITHMETIC_CODING} COMMAND ${CMC_TEST_MULTI_RES_COMPRESSION_ARITHMETIC_CODING})

endif()
