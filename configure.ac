AC_INIT([cmc], [0.1.0], [niklas.boeing@dlr.de])
AC_PREREQ(2.61)
AC_CONFIG_SRCDIR([src/cmc.h])
AC_CONFIG_AUX_DIR([build-aux])
AC_CONFIG_MACRO_DIR([config])
AC_PREFIX_DEFAULT([$PWD/local])

# Check for an archiver
AM_PROG_AR

# ==========================================================
# Check whether possible "--enable-..." and "--with-..."-
# options are set by the ./configure command.
# Possible options are:
#   --enable-debug          for a DEBUG-mode configuration
#   --enable-mpi            for compilation with MPI support
#   --enable-fortran        for enabling Fortran code dependent files
#   --with-moddir           the directory path for Fortran module files
#   --with-t8code           for linkage against t8code
#   --with-netcdf           for linkage against netCDF
# ==========================================================

# ==========================================================
# Check if the DEBUG-mode is enabled =======================
AC_ARG_ENABLE([debug],
    [AS_HELP_STRING([--enable-debug], [enable DEBUG mode (additional assertions)])],,
    [enableval=no])

if test "x$enableval" != xno ; then
    AC_DEFINE([ENABLE_DEBUG], 1, [Define if a DEBUG-mode configuration is choosen.])
fi
AM_CONDITIONAL([ENABLE_DEBUG], [test "x$enableval" != xno])
CMC_ENABLE_DEBUG="$enableval"

# ==========================================================
# Check if MPI is enabled ==================================
AC_ARG_ENABLE([mpi],
    [AS_HELP_STRING([--enable-mpi], [compile with MPI support])],,
    [enableval=no])

if test "x$enableval" != xno ; then
    AC_DEFINE([ENABLE_MPI], 1, [Define if MPI is enabled.])
fi
AM_CONDITIONAL([ENABLE_MPI], [test "x$enableval" != xno])
CMC_ENABLE_MPI="$enableval"

# ==========================================================
# Check if FORTRAN Code should be enabled ==================
AC_ARG_ENABLE([fortran],
    [AS_HELP_STRING([--enable-fortran], [activate parts with FORTRAN code])],,
    [enableval=no])

if test "x$enableval" != xno ; then
    AC_DEFINE([ENABLE_FORTRAN], 1, [Define if FORTRAN code should be compiled.])
fi
AM_CONDITIONAL([ENABLE_FORTRAN], [test "x$enableval" != xno])
CMC_ENABLE_FORTRAN="$enableval"
AC_SUBST([CMC_ENABLE_FORTRAN])

# ==========================================================
# Check if a path for the Fortran modules is explicitly 
# specified ====================================================
AC_ARG_WITH([moddir],
    [AS_HELP_STRING([--with-moddir], [Specify the path for Fortran module files])],,
    [withval=no])

if test "x$withval" = xyes ; then
    AC_MSG_ERROR([missing directory path for the module directory])
elif test "x$withval" != xno ; then
    AC_DEFINE([WITH_MODDIR], 1, [Specify the Fortran module files output directory.])
fi
AM_CONDITIONAL([WITH_MODDIR], [test "x$withval" != xno])
CMC_WITH_MODDIR="$withval"

# ==========================================================
# Check if option with-netcdf is given =====================
AC_ARG_WITH([netcdf],
    [AS_HELP_STRING([--with-netcdf], [compile with netCDF support])],,
    [withval=no])

if test "x$withval" != xno ; then
    AC_DEFINE([WITH_NETCDF], 1, [Define if netCDF-linkage is enabled.])
fi
AM_CONDITIONAL([WITH_NETCDF], [test "x$withval" != xno])
CMC_WITH_NETCDF="$withval"

# ==========================================================
# Check if option with-t8code is given =====================
AC_ARG_WITH([t8code],
    [AS_HELP_STRING([--with-t8code], [compile with external t8code support])],,
    [withval=no])

CMC_LIB_ADD_T8_SUBMODULE=
CMC_T8_CONFIG_AMFLAGS=
CMC_T8_SC_P4EST_CPP_FLAGS=
if test "x$withval" != xno ; then
    AC_DEFINE([WITH_T8CODE], 1, [Define if external t8code-linkage is enabled.])
else
    # If the option is not given, we will use t8code as a submodule
    AX_SUBDIRS_CONFIGURE([t8code],
        [],[],[[--with-netcdf=no],[--enable-fortran=no],[--with-moddir=no]],[])
    AC_DEFINE([WITH_T8CODE], 1, [Define if internal t8code-linkage is enabled.])

    CMC_T8_SC_P4EST_CPP_FLAGS+=" -I\$(top_srcdir)/t8code/src -It8code/config -It8code/src"
    CMC_T8_SC_P4EST_CPP_FLAGS+=" -I\$(top_srcdir)/t8code/sc/src -It8code/sc/config -It8code/sc/src"
    CMC_T8_SC_P4EST_CPP_FLAGS+=" -I\$(top_srcdir)/t8code/p4est/src -It8code/p4est/config -It8code/p4est/src"
    CMC_LIB_ADD_T8_SUBMODULE="\$(top_builddir)/t8code/src/libt8.la \$(top_builddir)/t8code/sc/src/libsc.la \$(top_builddir)/t8code/p4est/src/libp4est.la"
    CMC_T8_CONFIG_AMFLAGS="-I \$(top_srcdir)/t8code/config -I \$(top_srcdir)/t8code/sc/config -I \$(top_srcdir)/t8code/p4est/config"
fi
AC_SUBST(CMC_LDADD_T8_SUBMODULE, $CMC_LIB_ADD_T8_SUBMODULE)
AC_SUBST(CMC_T8_AMFLAGS, $CMC_T8_CONFIG_AMFLAGS)
AC_SUBST(CMC_T8_CPPFLAGS, $CMC_T8_SC_P4EST_CPP_FLAGS)

# In either case, we will have an available t8code
AM_CONDITIONAL([WITH_T8CODE], [TRUE])

# We store the value seperately, because in case of an external t8code linkage,
# its access has to be checked.
CMC_WITH_T8CODE="$withval"

# ==========================================================
# Check whether compilers have been set as environment
# variables. If not, try to assign some standard compilers
# (depending on whether MPI should be used or not).
# ==========================================================
# If the option --enable-mpi has been chosen
if test "x$CMC_ENABLE_MPI" != xno; then
    # Check if a MPI C compiler has been set in the environment
    if test "x$CC" = x ; then
        # No MPI C compiler has been set in the environment
        # Define standard mpicxx compiler
        export CC=mpicc
        AC_MSG_NOTICE([CC has been set to $CC])
    fi
    # Check if a MPI C++ compiler has been set in the environment
    if test "x$CXX" = x ; then
        # No MPI C++ compiler has been set in the environment
        # Define standard mpicc compiler
        export CXX=mpicxx
        AC_MSG_NOTICE([CXX has been set to $CXX])
    fi
    if test "x$CMC_ENABLE_FORTRAN" != xno; then
        # If FORTRAN code in the software shall be compiled
        # Check if a MPI F77 compiler has been set in the environment
        if test "x$F77" = x ; then
            # No MPI F77 compiler has been set in the environment
            # Define standard mpif77 compiler
            export F77=mpif77
            AC_MSG_NOTICE([F77 has been set to $F77])
        fi
        # Check if a MPI Fortran compiler has been set in the environment
        if test "x$FC" = x ; then
            # No MPI FC compiler has been set in the environment
            # Define standard mpif90 compiler
            export FC=mpif90
            AC_MSG_NOTICE([FC has been set to $FC])
        fi
    fi
fi
 
# ==========================================================
# Check for compilers and preprocessors
# ==========================================================
AC_PROG_CC
AC_PROG_CC_C_O
AC_PROG_CXX
AC_PROG_CXX_C_O
AC_PROG_CPP
AC_PROG_CXXCPP
if test "x$CMC_ENABLE_FORTRAN" != xno; then
    AC_PROG_F77
    AC_PROG_F77_C_O
    AC_PROG_FC
    AC_PROG_FC_C_O
    AC_F77_WRAPPERS
    AC_FC_WRAPPERS
    AC_FC_MODULE_EXTENSION
    AC_FC_MODULE_FLAG
    AC_FC_MODULE_OUTPUT_FLAG
    if test "x$CMC_WITH_MODDIR" != xno; then
      AC_SUBST(CMC_MODDIR, $CMC_WITH_MODDIR)
    fi
fi

# ==========================================================
# Perform some standard checks
# ==========================================================
AC_C_BIGENDIAN([AC_DEFINE([IS_BIGENDIAN], 1, [Definition to 1 on bigendian machines])])
AC_C_CONST
AC_C_INLINE
AC_TYPE_SIZE_T


# ==========================================================
# Check for some standard C++ headers which are used by this
# software.
# ==========================================================
# Check for C++ standard header files
AC_LANG_PUSH([C++])
AC_CHECK_HEADERS([iostream cstdint string cmath memory array vector cassert algorithm variant limits cstring functional],
    , [AC_MSG_ERROR([standard cpp headers were not found])])
AC_LANG_POP

# ==========================================================
# Check MPI compilation and linkage if MPI was enabled.
# ==========================================================
if test "x$CMC_ENABLE_MPI" != xno; then
# ==========================================================
# Checking for MPI C compiler ==============================
AC_MSG_CHECKING([for compilation and linkage of MPI C program])
AC_LANG_PUSH([C])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
     	#include <stddef.h>
        #include "mpi.h"
    ]],
    [[
        int err = MPI_Init(NULL, NULL);
        if (!(err == MPI_SUCCESS))
        {
            return err;
        }
        err = MPI_Finalize();
        if (!(err == MPI_SUCCESS))
        {
            return err;
        }
    ]]
    )],
    AC_MSG_RESULT([successful]),
    AC_MSG_RESULT([failed])
    AC_MSG_FAILURE([MPI C test has failed]))
AC_LANG_POP

# ==========================================================
# Checking for MPI C++ compiler ============================
AC_MSG_CHECKING([for compilation and linkage of MPI C++ program])
AC_LANG_PUSH([C++])
AC_LINK_IFELSE([AC_LANG_PROGRAM(
    [[
        #include "mpi.h"
    ]],
    [[
        int err = MPI_Init(NULL, NULL);
        if (!(err == MPI_SUCCESS))
        {
            return err;
        }
        err = MPI_Finalize();
        if (!(err == MPI_SUCCESS))
        {
            return err;
        }
    ]]
    )],
    AC_MSG_RESULT([successful]),
    AC_MSG_RESULT([failed])
    AC_MSG_FAILURE([MPI C++ test has failed]))
AC_LANG_POP

    if test "x$CMC_ENABLE_FORTRAN" != xno; then
    # ==========================================================
    # Checking for MPI F90 compiler ============================
    AC_MSG_CHECKING([for compilation and linkage of MPI Fortran program])
    AC_LANG_PUSH([Fortran])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(,
        [[
            use mpi

            INTEGER(KIND=4) ierr

            CALL MPI_Init(ierr)
            CALL MPI_Finalize(ierr)
        ]]
        )],
        AC_MSG_RESULT([successful]),
        AC_MSG_RESULT([failed])
        AC_MSG_FAILURE([MPI Fortran test has failed]))
    AC_LANG_POP
    fi
fi

# ==========================================================
# Check netCDF linkage and whether parallel netCDF routines
# are accesible.
# ==========================================================
cmc_netcdf_use_parallel=no
if test "x$CMC_WITH_NETCDF" != xno ; then
    ac_save_netcdf_CPPFLAGS="$CPPFLAGS"
    ac_save_netcdf_LIBS="$LIBS"
    ac_save_netcdf_LDFLAGS="$LDFLAGS"
    ac_netcdf_LIBS=" -lnetcdf"
    if test "x$CMC_WITH_NETCDF" = xserial ; then
        AC_MSG_NOTICE([A serial netCDF usage was chosen])
    elif test "x$CMC_WITH_NETCDF" != xyes ; then
        ac_netcdf_CPPFLAGS=" -I${CMC_WITH_NETCDF}/include"
        ac_netcdf_LDFLAGS=" -L${CMC_WITH_NETCDF}/lib"
    fi
    LIBS+="$ac_netcdf_LIBS"
    CPPFLAGS+="$ac_netcdf_CPPFLAGS"
    LDFLAGS+="$ac_netcdf_LDFLAGS"
    AC_LANG_PUSH([C])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
        [[
        #include <netcdf.h>
        ]],[[
        int ncid;
        nc_open (NULL, NC_NOWRITE, &ncid);
        ]])],
        [AC_MSG_CHECKING([linkage with netCDF])
        AC_MSG_RESULT([successful])],
        [AC_MSG_ERROR([Unable to link with netCDF library])])
    AC_LANG_POP

    AC_SUBST(NETCDF_CPPFLAGS, $ac_netcdf_CPPFLAGS)
    AC_SUBST(NETCDF_LIBFLAGS, $ac_netcdf_LIBS)
    # Check for parallel netCDF routines
    if test "x$CMC_WITH_NETCDF" != xserial ; then
        if test "x$CMC_ENABLE_MPI" != xno; then
            cmc_netcdf_use_parallel=yes
            AC_CHECK_HEADER([netcdf_par.h], [AC_DEFINE([WITH_NETCDF_PAR], [1], [Define if parallel netCDF routines are accessible])
            cmc_netcdf_use_parallel=yes],
            [AC_DEFINE([WITH_NETCDF_PAR], [0])
            cmc_netcdf_use_parallel=no],
            [#include <netcdf.h>
            ])
        fi
    fi
fi
AM_CONDITIONAL([WITH_NETCDF_PAR], [test "x$cmc_netcdf_use_parallel" != xno])

# ==========================================================
# Check t8code linkage
# ==========================================================
if test "x$CMC_WITH_T8CODE" != xno ; then
    ac_save_t8code_CPPFLAGS="$CPPFLAGS"
    ac_save_t8code_LIBS="$LIBS"
    ac_save_t8code_LDFLAGS="$LDFLAGS"
    ac_t8code_LIBS=" -lt8 -lp4est -lsc -lm -lz"
    if test "x$CMC_WITH_T8CODE" != xyes ; then
        ac_t8code_CPPFLAGS=" -I${CMC_WITH_T8CODE}/include"
        ac_t8code_LDFLAGS=" -L${CMC_WITH_T8CODE}/lib"
    fi
    LIBS+="$ac_t8code_LIBS"
    CPPFLAGS+="$ac_t8code_CPPFLAGS"
    LDFLAGS+="$ac_t8code_LDFLAGS"
    AC_LANG_PUSH([C])
    AC_LINK_IFELSE([AC_LANG_PROGRAM(
        [[
        #include <t8_forest/t8_forest.h>
        ]],[[
        t8_forest_t forest;
        ]])],
        [AC_MSG_CHECKING([linkage with t8code])
        AC_MSG_RESULT([successful])],
        [AC_MSG_ERROR([Unable to link with t8code library])])
    AC_LANG_POP
    AC_SUBST(T8CODE_CPPFLAGS, $ac_t8code_CPPFLAGS)
    AC_SUBST(T8CODE_LIBFLAGS, $ac_t8code_LIBS)
fi

## ==========================================================
## Define supplied/found compilers, preprocessor and Flags
# as variables. 
## ==========================================================
AC_DEFINE_UNQUOTED(CC,          ["${CC}"],          [C Compiler])
AC_DEFINE_UNQUOTED(CFLAGS,      ["${CFLAGS}"],      [C Compiler Flags])
AC_DEFINE_UNQUOTED(CPP,         ["${CPP}"],         [C Preprocessor])
AC_DEFINE_UNQUOTED(CPPFLAGS,    ["${CPPFLAGS}"],    [C Preprocessor Flags])
AC_DEFINE_UNQUOTED(CXX,         ["${CXX}"],         [C++ Compiler])
AC_DEFINE_UNQUOTED(CXXFLAGS,    ["${CXXFLAGS}"],    [C++ Compiler Flags])
AC_DEFINE_UNQUOTED(CXXCPP,      ["${CXXCPP}"],      [C++ Preprocessor Flags])
if test "x$CMC_ENABLE_FORTRAN" != xno; then
AC_DEFINE_UNQUOTED(FC,          ["${FC}"],          [FORTRAN Compiler])
AC_DEFINE_UNQUOTED(F77,         ["${F77}"],         [F77 Compiler])
AC_DEFINE_UNQUOTED(FFLAGS,      ["${FFLAGS}"],      [FORTRAN Compiler Flags])
AC_DEFINE_UNQUOTED(FCFLAGS,     ["${FCFLAGS}"],     [FORTRAN Compiler Flags])
fi
AC_DEFINE_UNQUOTED(LDFLAGS,     ["${LDFLAGS}"],     [Linker Flags])
AC_DEFINE_UNQUOTED(LIBS,        ["${LIBS}"],        [Libraries])

# Check for installer
AC_PROG_INSTALL

# Initialize automake after all compilers and preprocessors have been found
AM_INIT_AUTOMAKE([parallel-tests subdir-objects])

# Initialize Libtool
LT_INIT

# ==========================================================
# Define the output files of the configuration.
# ==========================================================
AC_CONFIG_FILES([
 Makefile
])

# ==========================================================
# Define configuration macro headers and prefix
# ==========================================================
AC_CONFIG_HEADERS([src/config.h])
AX_PREFIX_CONFIG_H([src/cmc_ac_config.h], [CMC])

## ==========================================================
## Give a summary about the supplied/found compilers and pre-
## processors.
## ==========================================================
AS_ECHO([""])
AS_ECHO(["===================================================="])
AS_ECHO(["  Resulting Compilers, Preprocessors and Flags are  "])
AS_ECHO(["===================================================="])
AS_ECHO(["CC            $CC "])
AS_ECHO(["CFLAGS        $CFLAGS "])
AS_ECHO(["CPP           $CPP "])
AS_ECHO(["CPPFLAGS      $CPPFLAGS "])
AS_ECHO([""])
AS_ECHO(["CXX           $CXX "])
AS_ECHO(["CXXFLAGS      $CXXFLAGS "])
AS_ECHO(["CXXCPP        $CXXCPP "])
if test "x$CMC_ENABLE_FORTRAN" != xno; then
    AS_ECHO([""])
    AS_ECHO(["F77           $F77 "])
    AS_ECHO(["FC            $FC "])
    AS_ECHO(["FFLAGS        $FFLAGS "])
    AS_ECHO(["FCFLAGS       $FCFLAGS "])
fi
AS_ECHO([""])
AS_ECHO(["LDFLAGS       $LDFLAGS "])
AS_ECHO(["LIBS          $LIBS "])
AS_ECHO(["===================================================="])
AS_ECHO([""])

AC_OUTPUT