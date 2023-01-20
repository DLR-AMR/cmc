ACLOCAL_AMFLAGS = -I config

AM_CPPFLAGS = -I@includedir@ @CPPFLAGS@ 
AM_CXXFLAGS = -std=c++17 -Wall -Werror
# @CXXFLAGS@
#AM_CFLAGS = @CFLAGS@
#AM_FCFLAGS = @FC_MODOUT@@includedir@
#AM_FFLAGS = @FC_MODOUT@@includedir@
AM_LDFLAGS = -L@libdir@ @LDFLAGS@
#AM_FCFLAGS = $(FC_MODOUT)$(includedir)
#AM_FFLAGS= $(FC_MODOUT)$(includedir) 
#AM_FCFLAGS_f90=@FC_MODOUT@@includedir@
#fmoddir=/home/niklas/programs/build/cmc_v03/include

CLEANFILES = 
DISTCLEANFILES = 
EXTRA_DIST = 
LDADD = 
bin_PROGRAMS =
nodist_include_HEADERS = 
include_HEADERS = 
noinst_HEADERS = 
lib_LTLIBRARIES =  

CLEANFILES += test -z "$(FC_MODEXT)" || rm -f *.$(FC_MODEXT)
DISTCLEANFILES += src/cmc_ac_config.h 

EXTRA_DIST += bootstrap

include src/Makefile.am
include example/Makefile.am
