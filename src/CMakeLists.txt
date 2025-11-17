target_include_directories(cmc PUBLIC
  $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>
)

target_sources(cmc PRIVATE 
    cmc.cxx 
    utilities/cmc_utilities.cxx 
    utilities/cmc_geo_domain.cxx 
    utilities/cmc_hyperslab.cxx 
    utilities/cmc_cart_coordinate.cxx 
    utilities/cmc_geo_utilities.cxx 
    mpi/cmc_mpi.cxx 
    mpi/cmc_mpi_io.cxx 
    mpi/cmc_mpi_data.cxx 
    netcdf/cmc_netcdf.cxx 
    netcdf/cmc_nc_writer.cxx 
    netcdf/cmc_nc_io.cxx 
    netcdf/cmc_nc_reader.cxx 
    t8code/cmc_t8_mesh.cxx 
    t8code/cmc_t8_interpolation.cxx 
    t8code/cmc_t8_adapt_track_inaccuracy.cxx 
    t8code/cmc_t8_morton.cxx 
    t8code/cmc_t8_mpi.cxx 
    input/cmc_netcdf.cxx 
    input/cmc_input_variable.cxx 
    amr/lossy/cmc_adaptive_coarsening.cxx
)

install(FILES
    cmc.hxx 
    cmc_config.h 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(FILES
    utilities/cmc_iface_compression_variable.hxx 
    utilities/cmc_iface_amr_compression_variable.hxx 
    utilities/cmc_iface_embedded_amr_compression_variable.hxx 
    utilities/cmc_log_functions.hxx 
    utilities/cmc_fallbacks.h 
    utilities/cmc_vector_view.hxx 
    utilities/cmc_utilities.hxx 
    utilities/cmc_compression_settings.hxx 
    utilities/cmc_dimension_interval.hxx 
    utilities/cmc_geo_domain.hxx 
    utilities/cmc_hyperslab.hxx 
    utilities/cmc_cart_coordinate.hxx  
    utilities/cmc_bit_vector.hxx 
    utilities/cmc_bit_map.hxx 
    utilities/cmc_huffman.hxx 
    utilities/cmc_arithmetic_encoder.hxx 
    utilities/cmc_ac_model.hxx 
    utilities/cmc_iface_entropy_alphabet.hxx 
    utilities/cmc_variable_utilities.hxx 
    utilities/cmc_error_domain.hxx 
    utilities/cmc_geo_utilities.hxx 
    utilities/cmc_error_domain_functions.hxx 
    utilities/cmc_iface_compression_adapt_data.hxx 
    utilities/cmc_endian.hxx 
    utilities/cmc_byte_iteration.hxx 
    utilities/cmc_byte_value.hxx 
    utilities/cmc_byte_compression_values.hxx 
    utilities/cmc_arithmetic_encoding_frequency_model.hxx 
    utilities/cmc_entropy_coder.hxx 
    utilities/cmc_arithmetic_encoding.hxx 
    utilities/cmc_interpolation_fn.hxx 
    utilities/cmc_compression_schema.hxx 
    utilities/cmc_serialization.hxx 
    utilities/cmc_multi_res_extraction_util.hxx 
    utilities/cmc_byte_compression_arithmetic_encoding_frequency_model.hxx 
    utilities/cmc_byte_compression_arithmetic_encoding.hxx 
    utilities/cmc_multi_res_compression_entropy_alphabet.hxx 
    utilities/cmc_multi_res_entropy_coder.hxx 
    utilities/cmc_prefix_compression_entropy_alphabet.hxx 
    utilities/cmc_prefix_entropy_coder.hxx 
    utilities/cmc_iface_arithmetic_encoding_freq_model.hxx 
    utilities/cmc_embedded_mesh_utilities.hxx 
    utilities/cmc_embedded_variable_attributes.hxx 
    utilities/cmc_trimmed_multi_res_compression_entropy_alphabet.hxx 
    utilities/cmc_trimmed_multi_res_coder.hxx 
    utilities/cmc_residual_quantization_entropy_alphabet.hxx 
    utilities/cmc_residual_quantization_encoder.hxx 
    utilities/cmc_lossy_multi_res_extraction_residual_computation.hxx 
    utilities/cmc_lossless_multi_res_extraction_residual_computation.hxx 
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/utilities
)

install(FILES
    mpi/cmc_mpi_io.hxx 
    mpi/cmc_mpi.hxx 
    mpi/cmc_mpi_data.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mpi
)

install(FILES
    netcdf/cmc_netcdf.hxx 
    netcdf/cmc_nc_io.hxx 
    netcdf/cmc_nc_io_conventions.hxx 
    netcdf/cmc_nc_writer.hxx 
    netcdf/cmc_nc_reader.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/netcdf
)

install(FILES
    input/cmc_netcdf.hxx 
    input/cmc_input_variable.hxx 
    input/cmc_forward_decl_input_variable.hxx 
    input/cmc_binary_reader.hxx 
    input/cmc_binary_reader_forward.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/input
)

install(FILES
    t8code/cmc_t8_mesh.hxx 
    t8code/cmc_t8_interpolation.hxx 
    t8code/cmc_t8_adapt_track_inaccuracy.hxx 
    t8code/cmc_t8_adapt_track_inaccuracy_forward.hxx 
    t8code/cmc_t8_morton.hxx 
    t8code/cmc_t8_mpi.hxx 
    t8code/cmc_t8_adaptation_callbacks.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/t8code
)

install(FILES
    amr/lossy/cmc_ac_compression_variable.hxx 
    amr/lossy/cmc_adaptive_coarsening.hxx 
    amr/lossy/cmc_byte_compression_variable.hxx 
    amr/lossy/cmc_multi_res_trimmed_residuals.hxx 
    amr/lossy/cmc_multi_res_scanline_compression.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/amr/lossy
)

install(FILES
    embedded/lossy/cmc_embedded_byte_compression_variable.hxx 
    embedded/lossy/cmc_embedded_trimmed_multi_res_extraction_compression.hxx 
    embedded/lossy/cmc_embedded_prefix_quantization_compression.hxx 
    embedded/lossy/cmc_embedded_byte_decompression_variable.hxx 
    embedded/lossy/cmc_embedded_prefix_quantization_decompression.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/embedded/lossy
)

install(FILES
    amr/lossless/cmc_byte_compression_variable.hxx 
    amr/lossless/cmc_prefix_extraction_compression.hxx 
    amr/lossless/cmc_prefix_extraction_compression_plain_suffixes.hxx 
    amr/lossless/cmc_multi_res_extraction_compression.hxx 
    amr/lossless/cmc_multi_res_extraction_decompression.hxx 
    amr/lossless/cmc_prefix_extraction_decompression.hxx 
    amr/lossless/cmc_prefix_extraction_decompression_plain_suffixes.hxx 
    amr/lossless/cmc_byte_decompression_variable.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/amr/lossless
)

install(FILES
    embedded/lossless/cmc_embedded_byte_compression_variable.hxx 
    embedded/lossless/cmc_embedded_byte_decompression_variable.hxx 
    embedded/lossless/cmc_embedded_prefix_extraction_compression.hxx 
    embedded/lossless/cmc_embedded_prefix_extraction_decompression.hxx 
    embedded/lossless/cmc_embedded_prefix_extraction_compression_plain_suffixes.hxx 
    embedded/lossless/cmc_embedded_prefix_extraction_decompression_plain_suffixes.hxx 
    embedded/lossless/cmc_embedded_multi_res_extraction_compression.hxx 
    embedded/lossless/cmc_embedded_multi_res_extraction_decompression.hxx 
    embedded/lossless/cmc_test_pcp4_embedded_compression.hxx 
    embedded/lossless/cmc_test_pcp4_embedded_decompression.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/embedded/lossless
)

install(FILES
    mesh_compression/cmc_iface_abstract_mesh_encoder.hxx 
    mesh_compression/cmc_iface_abstract_mesh_decoder.hxx 
    mesh_compression/cmc_iface_mesh_encoder.hxx 
    mesh_compression/cmc_iface_mesh_decoder.hxx 
    mesh_compression/cmc_mesh_encoder.hxx 
    mesh_compression/cmc_mesh_decoder.hxx 
    mesh_compression/cmc_iface_embedded_mesh_encoder.hxx 
    mesh_compression/cmc_iface_embedded_mesh_decoder.hxx 
    mesh_compression/cmc_embedded_mesh_encoder.hxx 
    mesh_compression/cmc_embedded_mesh_decoder.hxx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/mesh_compression
)

install(FILES
    compression_io/cmc_compression_output.hxx 
    compression_io/cmc_compression_attr_names.hxx 
    compression_io/cmc_decompression_input.hxx 
    compression_io/cmc_byte_compression_output.txx 
    compression_io/cmc_embedded_byte_compression_output.txx 
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/compression_io
)

install(TARGETS cmc DESTINATION ${CMAKE_INSTALL_LIBDIR})
install(TARGETS cmc EXPORT CmcTargets)
