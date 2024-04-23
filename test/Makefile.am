#flags_ld_testing = $(AM_LDFLAGS)
#flags_cpp_testing = $(AM_CPPFLAGS)
#flags_cxx_testing = $(AM_CXXFLAGS)

cmc_test_programs = test/cmc_create_compression_settings \
test/cmc_create_input_variables \
test/cmc_transform_input_var_coordinates \
test/cmc_create_error_domains \
test/cmc_create_compression_data \
test/cmc_setup_compression_data \
test/cmc_scaling_and_offset_for_input_var \
test/cmc_split_variables \
test/cmc_lossy_amr_compression \
test/cmc_prefix

test_cmc_create_compression_settings_SOURCES = test/cmc_create_compression_settings.cxx
test_cmc_create_input_variables_SOURCES = test/cmc_create_input_variables.cxx
test_cmc_transform_input_var_coordinates_SOURCES = test/cmc_transform_input_var_coordinates.cxx
test_cmc_create_error_domains_SOURCES = test/cmc_create_error_domains.cxx
test_cmc_create_compression_data_SOURCES = test/cmc_create_compression_data.cxx
test_cmc_setup_compression_data_SOURCES = test/cmc_setup_compression_data.cxx
test_cmc_scaling_and_offset_for_input_var_SOURCES = test/cmc_scaling_and_offset_for_input_var.cxx
test_cmc_split_variables_SOURCES = test/cmc_split_variables.cxx
test_cmc_lossy_amr_compression_SOURCES = test/cmc_lossy_amr_compression.cxx
test_cmc_prefix_SOURCES = test/cmc_prefix.cxx

#test_cmc_create_comrpession_settings_LDFLAGS = $(flags_ld_testing)
#test_cmc_create_comrpession_settings_CXXFLAGS = $(flags_cxx_testing)
#test_cmc_create_comrpession_settings_CPPFLAGS = $(flags_cpp_testing)

check_PROGRAMS += $(cmc_test_programs)
TESTS += $(cmc_test_programs)
