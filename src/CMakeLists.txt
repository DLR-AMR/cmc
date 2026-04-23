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
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/utilities
)

install(FILES
    patch/lossless/cmc_multi_res_extraction.hxx
    patch/lossless/cmc_patch_multi_res_extraction_util.hxx 
    patch/lossless/cmc_multi_res_decompression.hxx 
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
# If configured wit t8code

install(FILES
    t8code/cmc_t8_mesh.hxx 
    t8code/cmc_t8_adaptation_callbacks.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/t8code
)

#End if configured with t8code
endif()


######## WITH T8CODE AND MPI
# Lossless AMR compressor needs MPI and t8code 
if (CMC_WITH_T8CODE AND CMC_ENABLE_MPI)
install(FILES
    amr/lossless/cmc_par_multi_res_extraction_util.hxx 
    amr/lossless/cmc_par_multi_res_extraction.hxx 
    amr/lossless/cmc_par_multi_res_decomrpession.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/amr/lossless
)

#End if configured with t8code and MPI
endif()










######Old 
#
#target_sources(cmc PRIVATE 
#    cmc.cxx 
#    utilities/cmc_utilities.cxx 
#    utilities/cmc_geo_domain.cxx 
#    utilities/cmc_hyperslab.cxx 
#    utilities/cmc_cart_coordinate.cxx 
#    utilities/cmc_geo_utilities.cxx  
#    input/cmc_input_variable.cxx 
#    utilities/cmc_morton.cxx 
#    mpi/cmc_mpi.cxx 
#)
#
#install(FILES
#    cmc.hxx 
#    cmc_config.h 
#    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
#)
#
#install(FILES
#    utilities/cmc_bits.hxx 
#    utilities/cmc_bits_vector.hxx
#    utilities/cmc_bits_vector_view.hxx
#    utilities/cmc_bits_stream_decoder.hxx 
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/utilities
#)
#
#install(FILES
#    utilities/cmc_iface_compression_variable.hxx 
#    utilities/cmc_iface_patch_compression_variable.hxx 
#    utilities/cmc_iface_patch_decompression_variable.hxx 
#    utilities/cmc_log_functions.hxx 
#    utilities/cmc_utilities.hxx 
#    utilities/cmc_compression_settings.hxx 
#    utilities/cmc_dimension_interval.hxx 
#    utilities/cmc_geo_domain.hxx 
#    utilities/cmc_hyperslab.hxx 
#    utilities/cmc_cart_coordinate.hxx  
#    utilities/cmc_huffman_coder.hxx  
#    utilities/cmc_variable_utilities.hxx 
#    utilities/cmc_error_domain.hxx 
#    utilities/cmc_geo_utilities.hxx 
#    utilities/cmc_error_domain_functions.hxx 
#    utilities/cmc_compression_schema.hxx 
#    utilities/cmc_morton.hxx 
#    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/utilities
#)
#
#install(FILES
#    input/cmc_input_variable.hxx 
#    input/cmc_forward_decl_input_variable.hxx 
#    input/cmc_binary_reader.hxx 
#    input/cmc_binary_reader_forward.hxx 
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/input
#)
#
#install(FILES
#    compression_io/cmc_compression_serial_output.hxx 
#    compression_io/cmc_decompression_serial_input.hxx 
#    compression_io/cmc_patch_byte_compression_serial_output.txx 
#    compression_io/cmc_serial_io_util.hxx 
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/compression_io
#)
#
#install(FILES
#    mpi/cmc_mpi.hxx 
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mpi
#)
#
#install(FILES
#    #patch/lossless/cmc_patch_byte_compression_variable.hxx 
#    #patch/lossless/cmc_patch_byte_decompression_variable.hxx 
#    #patch/lossless/cmc_patch_prefix_extraction_plain_suffixes_compression.hxx 
#    #patch/lossless/cmc_patch_prefix_extraction_plain_suffixes_decompression.hxx 
#    #patch/lossless/cmc_patch_multi_res_extraction_compression.hxx 
#    #patch/lossless/cmc_patch_multi_res_extraction_decompression.hxx 
#    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/patch/lossless
#)
#
######### Configuration dependant files ########
######### WITH MPI
#if (CMC_ENABLE_MPI)
## If configured with MPI
#
#target_sources(cmc PRIVATE 
#    mpi/cmc_mpi_io.cxx 
#    mpi/cmc_mpi_data.cxx
#)
#
#install(FILES
#    mpi/cmc_mpi_io.hxx 
#    mpi/cmc_mpi_data.hxx 
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mpi
#)
#
#install(FILES
#    # compression_io/cmc_compression_mpi_io_output.hxx 
#    # compression_io/cmc_decompression_mpi_io_input.hxx 
#    # compression_io/cmc_compression_mpi_io_util.hxx  
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/compression_io
#)
#
## End if configured with MPI
#endif()
######### END WITH MPI
#
#
#
######### WITH NETCDF
#if (CMC_WITH_NETCDF)
## If configured with netCDF
#
#install(FILES
#    input/cmc_netcdf.hxx 
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/input
#)
#
#install(FILES
#    netcdf/cmc_netcdf.hxx 
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/netcdf
#)
#
#target_sources(cmc PRIVATE 
#    input/cmc_netcdf.cxx 
#    netcdf/cmc_netcdf.cxx
#)
#
#if (CMC_ENABLE_MPI AND CMC_WITH_T8CODE)
## If configured with netCDF, MPI and t8code
#
## End If configured with netCDF, MPI and t8code
#endif()
#
## End if configured with netCDF
#endif()
######### END WITH NETCDF
#
#
#
######### WITH T8CODE
#if (CMC_WITH_T8CODE AND CMC_ENABLE_MPI)
## If configured wit t8code
#
#target_sources(cmc PRIVATE 
#    amr/lossy/cmc_adaptive_coarsening.cxx 
#    t8code/cmc_t8_mesh.cxx 
#    t8code/cmc_t8_interpolation.cxx 
#    t8code/cmc_t8_adapt_track_inaccuracy.cxx 
#)
#
#install(FILES 
#    utilities/cmc_iface_compression_adapt_data.hxx 
#    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/utilities
#)
#
#install(FILES
#    t8code/cmc_t8_mesh.hxx 
#    t8code/cmc_t8_interpolation.hxx 
#    t8code/cmc_t8_adapt_track_inaccuracy.hxx 
#    t8code/cmc_t8_adapt_track_inaccuracy_forward.hxx  
#    t8code/cmc_t8_mpi.hxx 
#    t8code/cmc_t8_adaptation_callbacks.hxx 
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/t8code
#)
#
#install(FILES
#    amr/lossy/cmc_ac_compression_variable.hxx 
#    amr/lossy/cmc_adaptive_coarsening.hxx 
#    amr/lossy/cmc_simul_ac_compression_variable.hxx 
#    amr/lossy/cmc_simul_adaptive_coarsening.hxx 
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/amr/lossy
#)
#
#install(FILES
#    amr/lossless/cmc_par_multi_res_extraction_util.hxx 
#    amr/lossless/cmc_par_multi_res_extraction.hxx 
#    amr/lossless/cmc_par_multi_res_decomrpession.hxx 
#    # amr/lossless/cmc_byte_par_compression_variable.hxx 
#    # amr/lossless/cmc_multi_res_par_extraction_compression.hxx 
#    # amr/lossless/cmc_multi_res_par_extraction_util.hxx 
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/amr/lossless
#)
#
#if (CMC_ENABLE_MPI)
## If configured with t8code and MPI
#target_sources(cmc PRIVATE 
#    t8code/cmc_t8_mpi.cxx 
#)
#install(FILES
#    t8code/cmc_t8_mpi.hxx 
#DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/t8code
#)
## End If configured with t8code and MPI
#endif()
#
##End if configured with t8code
#endif()
######### END WITH T8CODE



######## END Configuration dependant files ########

