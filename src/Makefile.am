lib_LTLIBRARIES += src/libcmc.la

nodist_include_HEADERS += src/cmc_ac_config.h

include_HEADERS += src/cmc.hxx \
src/cmc_config.h \
src/utilities/cmc_constants_definitions.h \
src/utilities/cmc_log_functions.hxx \
src/utilities/cmc_fallbacks.h \
src/mpi/cmc_mpi_io.h \
src/netcdf/cmc_netcdf.hxx \
src/utilities/cmc_utilities.hxx \
src/mpi/cmc_mpi.hxx \
src/t8code/cmc_t8_data_variables_forward.hxx \
src/t8code/cmc_t8_data_variables.hxx \
src/t8code/cmc_t8_mesh.hxx \
src/t8code/cmc_t8_interpolation.hxx \
src/utilities/cmc_span.hxx \
src/t8code/cmc_t8_adapt_track_inaccuracy.hxx \
src/t8code/cmc_t8_adapt.hxx \
src/lossy/cmc_amr_lossy_compression_settings.hxx \
src/t8code/cmc_t8_data.hxx \
src/utilities/cmc_dimension_interval.hxx \
src/utilities/cmc_geo_domain.hxx \
src/utilities/cmc_hyperslab.hxx \
src/utilities/cmc_cart_coordinate.hxx \
src/utilities/cmc_input_variable_forward.hxx  \
src/utilities/cmc_input_variable.hxx  \
src/t8code/cmc_t8_morton.hxx \
src/t8code/cmc_t8_mpi.hxx \
src/mpi/cmc_mpi_data.hxx \
src/lossy/cmc_amr_lossy_compression.hxx \
src/t8code/cmc_t8_adapt_callbacks.h \
src/utilities/cmc_variable_transformer_forward.hxx \
src/utilities/cmc_variable_transformer.hxx \
src/t8code/cmc_t8_adapt_track_inaccuracy_forward.hxx \
src/utilities/cmc_output_variable_forward.hxx \
src/utilities/cmc_output_variable.hxx \
src/utilities/cmc_geo_utilities.hxx \
src/lossy/cmc_sz_like_compression.hxx \
src/utilities/cmc_prefix.hxx \
src/lossy/cmc_prefix_lossy_compression.hxx \
src/t8code/cmc_t8_byte_variable.hxx \
src/t8code/cmc_t8_variables.hxx \
src/t8code/cmc_t8_prefix_encoding.hxx \
src/netcdf/cmc_nc_decompress.hxx \
src/netcdf/cmc_nc_io.hxx \
src/netcdf/cmc_nc_writer.hxx \
src/netcdf/cmc_nc_reader.hxx \
src/cmc_compression.hxx \
src/netcdf/cmc_nc_io_naming_conventions.hxx

if ENABLE_FORTRAN
include_HEADERS += src/fortran_interfaces/cmcc.h \
src/fortran_interfaces/cmcc_wrapper.h 
endif

src_libcmc_la_SOURCES = src/cmc.cxx \
src/mpi/cmc_mpi.cxx \
src/mpi/cmc_mpi_io.cxx \
src/netcdf/cmc_netcdf.cxx \
src/utilities/cmc_utilities.cxx \
src/t8code/cmc_t8_mesh.cxx \
src/t8code/cmc_t8_data_variables.cxx \
src/t8code/cmc_t8_interpolation.cxx \
src/t8code/cmc_t8_adapt_track_inaccuracy.cxx \
src/t8code/cmc_t8_adapt.cxx \
src/lossy/cmc_amr_lossy_compression_settings.cxx \
src/t8code/cmc_t8_data.cxx \
src/utilities/cmc_geo_domain.cxx \
src/utilities/cmc_hyperslab.cxx \
src/utilities/cmc_input_variable.cxx \
src/t8code/cmc_t8_morton.cxx \
src/t8code/cmc_t8_mpi.cxx \
src/mpi/cmc_mpi_data.cxx \
src/lossy/cmc_amr_lossy_compression.cxx \
src/t8code/cmc_t8_adapt_callbacks.cxx \
src/utilities/cmc_variable_transformer.cxx \
src/utilities/cmc_output_variable.cxx \
src/utilities/cmc_cart_coordinate.cxx \
src/utilities/cmc_geo_utilities.cxx \
src/lossy/cmc_sz_like_compression.cxx \
src/lossy/cmc_prefix_lossy_compression.cxx \
src/t8code/cmc_t8_byte_variable.cxx \
src/netcdf/cmc_nc_decompress.cxx \
src/netcdf/cmc_nc_writer.cxx \
src/netcdf/cmc_nc_io.cxx \
src/netcdf/cmc_nc_reader.cxx \
src/cmc_compression.cxx


if ENABLE_FORTRAN
cmc_module_sources = src/fortran_interfaces/cmc_fortran.f90 \
src/fortran_interfaces/cmc_messy_interface.f90 \
src/fortran_interfaces/cmc_netcdf_interface.f90 \
src/fortran_interfaces/cmc_amr_compression_interface.f90

src_libcmc_la_SOURCES += src/fortran_interfaces/cmcc.c \
src/fortran_interfaces/cmcc_wrapper.c $(cmc_module_sources)

MODSOURCES += $(cmc_module_sources)
endif


src_libcmc_la_CPPFLAGS = $(AM_CPPFLAGS)
src_libcmc_la_LDFLAGS = $(AM_LDFLAGS)
src_libcmc_la_CXXFLAGS = $(AM_CXXFLAGS)
src_libcmc_la_CFLAGS = $(CMC_CFLAGS)
if ENABLE_FORTRAN
src_libcmc_la_FCFLAGS = $(AM_FCFLAGS)
src_libcmc_la_FFLAGS = $(CMC_FFLAGS)
endif
src_libcmc_la_LIBFLAGS = $(NETCDF_LIBFLAGS) $(T8CODE_LIBFLAGS) @CMC_LDADD_T8_SUBMODULE@
src_libcmc_la_LIBADD = @CMC_LDADD_T8_SUBMODULE@
EXTRA_src_libcmc_la_DEPENDENCIES = @CMC_LDADD_T8_SUBMODULE@

LDADD += @top_builddir@/src/libcmc.la @CMC_LDADD_T8_SUBMODULE@

AM_CPPFLAGS += -I@top_srcdir@/src @CMC_T8_CPPFLAGS@
