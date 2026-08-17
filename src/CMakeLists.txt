target_include_directories(cmc PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
)

install(FILES
    cmc.hxx 
    cmc_config.h 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

target_sources(cmc PRIVATE 
    cmc.cxx 
)

install(FILES
    input/cmc_binary_file_reader.hxx 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/input
)

install(FILES
    utilities/cmc_bits.hxx 
    utilities/cmc_bits_vector.hxx
    utilities/cmc_bits_vector_view.hxx
    utilities/cmc_bits_stream_decoder.hxx 
    utilities/cmc_log_functions.hxx
    utilities/cmc_huffman_coder.hxx
    utilities/cmc_utilities.hxx 
    utilities/cmc_error_domain.hxx 
    utilities/cmc_error_domain_functions.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/utilities
)

install(FILES
    patch/lossless/cmc_multi_res_extraction.hxx 
    patch/lossless/cmc_patch_multi_res_extraction_util.hxx 
    patch/lossless/cmc_multi_res_decompression.hxx 
    patch/lossless/cmc_fast_multi_res_extraction.hxx 
    patch/lossless/cmc_fast_multi_res_decompression.hxx 
    patch/lossless/cmc_multi_res_pcp4_extraction.hxx 
    patch/lossless/cmc_multi_res_pcp4_decompression.hxx 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/patch/lossless
)

install(FILES
    patch/lossy/cmc_multi_res_util.hxx
    patch/lossy/cmc_multi_res_extraction_idw.hxx 
    patch/lossy/cmc_multi_res_decompression_idw.hxx 
    patch/lossy/cmc_multi_res_extraction_rbf.hxx 
    patch/lossy/cmc_multi_res_decompression_rbf.hxx 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/patch/lossless
)

######## WITH MPI
if (CMC_ENABLE_MPI)
# If configured with MPI

target_sources(cmc PRIVATE 
    mpi/cmc_mpi.cxx 
)

install(FILES
    mpi/cmc_mpi.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mpi
)

# End if configured with MPI
endif()

######## WITH T8CODE
if (CMC_WITH_T8CODE)
# If configured with t8code

install(FILES
    t8code/cmc_t8_mesh.hxx 
    t8code/cmc_t8_adaptation_callbacks.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/t8code
)

#End if configured with t8code
endif()

######## WITH T8CODE AND MPI
if (CMC_WITH_T8CODE AND CMC_ENABLE_MPI)
# Lossy AMR compressor needs MPI and t8code
install(FILES
    amr/lossy/cmc_ac_compression_variable.hxx 
    amr/lossy/cmc_adaptive_coarsening.hxx 
    amr/lossy/cmc_simul_ac_compression_variable.hxx 
    amr/lossy/cmc_simul_adaptive_coarsening.hxx 
    amr/lossy/cmc_byte_compression_variable.hxx 
    amr/lossy/cmc_multi_res_trimmed_residuals.hxx 
    amr/lossy/cmc_par_multi_res_extraction_util.hxx 
    amr/lossy/cmc_par_multi_res_extraction.hxx 
    amr/lossy/cmc_par_multi_res_decompression.hxx 
    amr/lossy/cmc_par_multi_res_error_mesh.hxx 
    amr/lossy/cmc_par_multi_res_idw_interpolation_util.hxx 
    amr/lossy/cmc_par_multi_res_extraction_idw.hxx 
    amr/lossy/cmc_par_multi_res_decompression_idw.hxx
    amr/lossy/trixi/cmc_par_multi_res_extraction_idw.hxx 
    amr/lossy/trixi/cmc_par_multi_res_decompression_idw.hxx 
    amr/lossy/trixi/cmc_par_multi_res_util.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/amr/lossy
)

# Lossless AMR compressor needs MPI and t8code 
install(FILES
    amr/lossless/cmc_par_multi_res_extraction_util.hxx 
    amr/lossless/cmc_par_multi_res_extraction.hxx 
    amr/lossless/cmc_par_multi_res_decompression.hxx 
    amr/lossless/cmc_par_fast_multi_res_extraction_util.hxx 
    amr/lossless/cmc_par_fast_multi_res_extraction.hxx 
    amr/lossless/cmc_par_fast_multi_res_decompression.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/amr/lossless
)

#End if configured with t8code and MPI
endif()
