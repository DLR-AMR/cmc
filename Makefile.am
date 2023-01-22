ACLOCAL_AMFLAGS = -I config

AM_CPPFLAGS = -I@includedir@ @CPPFLAGS@ 
AM_CXXFLAGS = -std=c++17 -Wall -Werror
AM_LDFLAGS = -L@libdir@ @LDFLAGS@

CLEANFILES = 
DISTCLEANFILES = 
EXTRA_DIST = 
LDADD = 
bin_PROGRAMS =
nodist_include_HEADERS = 
include_HEADERS = 
noinst_HEADERS = 
lib_LTLIBRARIES =  

if ENABLE_FORTRAN
CLEANFILES += *.$(FC_MODEXT)
endif

DISTCLEANFILES += src/cmc_ac_config.h 

EXTRA_DIST += bootstrap

include src/Makefile.am
include example/Makefile.am
