lib_LTLIBRARIES += src/libcmc.la

nodist_include_HEADERS += src/cmc_ac_config.h

include_HEADERS += src/cmc.h \
src/cmc_config.h \
src/utilities/cmc_constants_definitions.h \
src/utilities/cmc_util.h \
src/utilities/cmc_log_functions.h \
src/utilities/cmc_geo_util.h \
src/utilities/cmc_fallbacks.h \
src/utilities/cmc_container.h \
src/utilities/cmc_container.hxx \
src/mpi/cmc_mpi.h \
src/netcdf/cmc_netcdf.h \
src/netcdf/cmc_netcdf.hxx \
src/t8code/cmc_t8code.h \
src/t8code/cmc_t8code_data.h \
src/t8code/cmc_t8code_data.hxx \
src/t8code/cmc_t8code_geo_mesh.h \
src/t8code/cmc_t8code_geo_data.h \
src/t8code/cmc_t8_adapt_callbacks.h \
src/t8code/cmc_t8_replace_callbacks.h \
src/component_interfaces/cmc_t8_nc.h \
src/component_interfaces/cmc_t8_messy.h \
src/lossy/cmc_amr_compressor.h \
src/messy/cmc_messy.h 

if ENABLE_FORTRAN
include_HEADERS += src/fortran_interfaces/cmcc.h \
src/fortran_interfaces/cmcc_wrapper.h \
cmc_fortran.mod \
cmc_messy_interface.mod \
cmc_netcdf_interface.mod \
cmc_amr_compression_interface.mod
endif

noinst_HEADERS += src/utilities/cmc_container.hxx

src_libcmc_la_SOURCES = src/cmc.cxx \
src/utilities/cmc_util.cxx \
src/utilities/cmc_geo_util.cxx \
src/utilities/cmc_container.cxx \
src/mpi/cmc_mpi.cxx \
src/netcdf/cmc_netcdf.cxx \
src/t8code/cmc_t8code.cxx \
src/t8code/cmc_t8code_data.cxx \
src/t8code/cmc_t8code_geo_mesh.cxx \
src/t8code/cmc_t8code_geo_data.cxx \
src/t8code/cmc_t8_adapt_callbacks.cxx \
src/t8code/cmc_t8_replace_callbacks.cxx \
src/component_interfaces/cmc_t8_nc.cxx \
src/component_interfaces/cmc_t8_messy.cxx \
src/lossy/cmc_amr_compressor.cxx \
src/messy/cmc_messy.cxx 

if ENABLE_FORTRAN
src_libcmc_la_SOURCES += src/fortran_interfaces/cmcc.c \
src/fortran_interfaces/cmcc_wrapper.c \
src/fortran_interfaces/cmc_fortran.f90 \
src/fortran_interfaces/cmc_messy_interface.f90 \
src/fortran_interfaces/cmc_netcdf_interface.f90 \
src/fortran_interfaces/cmc_amr_compression_interface.f90
endif


src_libcmc_la_CPPFLAGS = $(AM_CPPFLAGS)
src_libcmc_la_LDFLAGS = $(AM_LDFLAGS)
src_libcmc_la_LIBFLAGS = $(NETCDF_LIBFLAGS) $(T8CODE_LIBFLAGS)
src_libcmc_la_CXXFLAGS = $(AM_CXXFLAGS)
src_libcmc_la_CFLAGS = $(CMC_CFLAGS)
if ENABLE_FORTRAN
src_libcmc_la_FCFLAGS = $(CMC_FCFLAGS)
src_libcmc_la_FFLAGS = $(CMC_FFLAGS)
endif

LDADD += @top_builddir@/src/libcmc.la $(NETCDF_LIBFLAGS) $(T8CODE_LIBFLAGS)

AM_CPPFLAGS += -I@top_srcdir@/src