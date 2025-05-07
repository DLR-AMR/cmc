ACLOCAL_AMFLAGS = -I config @CMC_T8_AMFLAGS@

AM_CPPFLAGS = -I@includedir@ @CPPFLAGS@ 
AM_CXXFLAGS = -std=c++2a -O0
#-Wall -Wextra
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
check_PROGRAMS = 
TESTS = 

SUBDIRS = $(subdirs) $(subdirs_extra)
DIST_SUBDIRS = $(SUBDIRS)

DISTCLEANFILES += src/cmc_ac_config.h 

EXTRA_DIST += bootstrap \
  INSTALL.md README.md

include src/Makefile.am
include example/Makefile.am
include test/Makefile.am
