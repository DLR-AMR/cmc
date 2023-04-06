ACLOCAL_AMFLAGS = -I config

AM_CPPFLAGS = -I@includedir@ @CPPFLAGS@ 
AM_CXXFLAGS = -std=c++17 -Wall -Werror
AM_LDFLAGS = -L@libdir@ @LDFLAGS@
AM_FCFLAGS = @FCFLAGS@

CLEANFILES = 
DISTCLEANFILES = 
EXTRA_DIST = 
LDADD = 
bin_PROGRAMS =
nodist_include_HEADERS = 
include_HEADERS = 
noinst_HEADERS = 
lib_LTLIBRARIES =  
MODSOURCES = 

DISTCLEANFILES += src/cmc_ac_config.h 

EXTRA_DIST += bootstrap

include src/Makefile.am
include example/Makefile.am

# Include the Fortran specific makefile
include cmc_fortran_specific.mk
