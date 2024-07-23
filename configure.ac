#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT([vartools], [1.41], [jhartman@astro.princeton.edu], [vartools], [http://www.astro.princeton.edu/~jhartman/vartools.html])
AC_CONFIG_SRCDIR([src/doublelinklist.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

SRCDIR=`(cd $srcdir && pwd)`
AC_MSG_RESULT([\"$SRCDIR\"])
AC_SUBST([SRCDIR])

builddir=`pwd`
AC_MSG_CHECKING([location of build directory])
AC_MSG_RESULT([\"$builddir\"])
BUILDDIR=$builddir
AC_SUBST([BUILDDIR])

DATAROOTDIR=`(cd $datarootdir && pwd)`
AC_MSG_RESULT([\"$DATAROOTDIR\"])
AC_SUBST([DATAROOTDIR])

USERLIBSDIR=$datarootdir/vartools/USERLIBS
USERFUNCSDIR=$datarootdir/vartools/USERFUNCS
AC_SUBST([USERLIBSDIR])
AC_SUBST([USERFUNCSDIR])

#AC_CONFIG_AUX_DIR([libltdl/config])

#AC_CONFIG_MACRO_DIR([libltdl/m4])

AM_INIT_AUTOMAKE([subdir-objects])

LT_CONFIG_LTDL_DIR([libltdl])

LT_INIT([dlopen])

LTDL_INIT([recursive])

AC_CONFIG_FILES([libltdl/Makefile])

AC_SUBST([LTDLINCl])
AC_SUBST([LIBLTDL])

# Checks for programs.
AC_PROG_CC
AC_PROG_AWK
AC_PROG_INSTALL
AC_PROG_LN_S
AC_PROG_MAKE_SET
AC_PROG_GREP
AC_PROG_SED
AC_PROG_FC

AC_F77_LIBRARY_LDFLAGS

AC_CHECK_LIB([m], [sqrt])

# Checks for libraries.
AC_ARG_WITH([cfitsio], 
  [AS_HELP_STRING([--with-cfitsio=PATH],
    [prefix directory for the cfitsio library package. Equivalent to --with-cfitsio-include=PATH --with-cfitsio-lib=PATH/libcfitsio.a])],
  [cfitsio_include_path_string=${with_cfitsio}
   cfitsio_lib_path_string=${with_cfitsio}/libcfitsio.a],
  [with_cfitsio=check])

AC_ARG_WITH([cfitsio-include], 
  [AS_HELP_STRING([--with-cfitsio-include=PATH],
    [PATH to the fitsio include file])],
  [if test "x$with_cfitsio_include" != xno ; then
      cfitsio_include_path_string=$with_cfitsio_include ;
   fi],
  [])

AC_ARG_WITH([cfitsio-lib], 
  [AS_HELP_STRING([--with-cfitsio-lib=PATH],
    [Full PATH to the libcfitsio.a library file])],
  [if test "x$with_cfitsio_lib" != xno ; then
      cfitsio_lib_path_string=$with_cfitsio_lib ; 
   fi],
  [])

LIBCFITSIO=
INCLUDECFITSIO=
AS_IF([test "x$cfitsio_include_path_string" != x && test "x$cfitsio_lib_path_string" != x],
  [AS_IF([test -f "${cfitsio_include_path_string}/fitsio.h" && test -f "$cfitsio_lib_path_string"],
       [AC_SUBST([LIBCFITSIO],"$cfitsio_lib_path_string")
        AC_SUBST([INCLUDECFITSIO],"-I ${cfitsio_include_path_string}")
        AC_DEFINE([HAVE_CFITSIO], [1],
                  [Define if you have libcfitsio.a])
       ],
       AC_MSG_WARN([Did not find fitsio.h in directory "${cfitsio_include_path_string}" or did not find file "$cfitsio_lib_path_string"]))],[])


AS_IF([test "x$cfitsio_include_path_string" = x && test "x$cfitsio_lib_path_string" = x && test "x$with_cfitsio" = xcheck],
  [AC_CHECK_LIB([cfitsio], [fits_init_cfitsio],
    [AC_SUBST([LIBCFITSIO], ["-lcfitsio"])
     AC_DEFINE([HAVE_CFITSIO], [1],
               [Define if you have libcfitsio])
    ],
    [if test "x$with_readcfitsio" != xcheck; then
       AC_MSG_FAILURE(
         [--with-cfitsio was given, but test for cfitsio failed])
     fi
    ])],[])


#AC_SEARCH_LIBS([dlopen], [dl dld], 
#  [AC_DEFINE([HAVE_DYNAMICLIB], [1], [Define if you have dlopen])], 
#  [AC_CHECK_HEADER([windows.h],
#   [AC_DEFINE([ISWINDOWS], [1], [Define if the windows header is here])],
#   [AC_MSG_WARN(["Missing library dl and windows.h not found, will not be able to load user-developed libraries."])])])

AC_ARG_WITH([gsl], 
  [AS_HELP_STRING([--with-gsl],
    [support commands using GNU Scientific Library Functions @<:@default=check@:>@])],
  [],
  [with_gsl=check])

LIBGSL=
AS_IF([test "x$with_gsl" != xno],
  [AC_CHECK_LIB([gslcblas], [cblas_dgemm],
    [AC_CHECK_LIB([gsl], [gsl_error],
    [AC_SUBST([LIBGSL], ["-lgsl -lgslcblas"])
     AC_DEFINE([HAVE_GSL], [1],
               [Define if you have libgsl and libgslcblas])
    ],
    [if test "x$with_gsl" != xcheck ; then
       AC_MSG_FAILURE(
         [--with-gsl was given, but test for gsl failed])
     fi],[-lgslcblas])],
    [if test "x$with_gsl" != xcheck ; then
       AC_MSG_FAILURE(
         [--with-gsl was given, but test for gsl failed])
     fi])])

AM_CONDITIONAL([HAVE_GSL], [ test "x$LIBGSL" != x ])

AC_CHECK_LIB([pthread], [pthread_create],
  [AC_SUBST([LIBPTHREAD], ["-lpthread"])
   AC_DEFINE([HAVE_PTHREAD], [1], [Define if you have lpthread])],
  [AC_MSG_WARN(["Missing library pthread, will not be able to do parallel processing."])])

cspice_include_path_string=/usr/include/cspice
cspice_lib_path_string=/usr/lib/cspice.a
AC_ARG_WITH([cspice], 
  [AS_HELP_STRING([--with-cspice=PATH],
    [prefix directory for the cspice library package. Equivalent to --with-cspice-include=PATH/include --with-cspice-lib=PATH/lib/cspice.a])],
  [cspice_include_path_string=${with_cspice}/include
   cspice_lib_path_string=${with_cspice}/lib/cspice.a],
  [])

AC_ARG_WITH([cspice-include], 
  [AS_HELP_STRING([--with-cspice-include=PATH],
    [PATH to the cspice include file])],
  [if test "x$with_cspice_include" != xno ; then
      cspice_include_path_string=$with_cspice_include ;
   fi],
  [])

AC_ARG_WITH([cspice-lib], 
  [AS_HELP_STRING([--with-cspice-lib=PATH],
    [Full PATH to the cspice.a library file])],
  [if test "x$with_cspice_lib" != xno ; then
      cspice_lib_path_string=$with_cspice_lib ; 
   fi],
  [])


LIBCSPICE=
INCLUDECSPICE=
AS_IF([test "x$cspice_include_path_string" != x && test "x$cspice_lib_path_string" != x],
  [AS_IF([test -f "${cspice_include_path_string}/SpiceUsr.h" && test -f "$cspice_lib_path_string"],
       [AC_SUBST([LIBCSPICE],"$cspice_lib_path_string")
        AC_SUBST([INCLUDECSPICE],"-I ${cspice_include_path_string}")
        AC_DEFINE([HAVE_CSPICE], [1],
                  [Define if you have cspice.a])
       ],
       AC_MSG_WARN([Did not find SpiceUsr.h in directory "${cspice_include_path_string}" or did not find file "$cspice_lib_path_string"]))],[])

AC_ARG_WITH([pythonhome], 
  [AS_HELP_STRING([--with-pythonhome=@<:@check,no,PATH@:>@],
    [Hard-code the PYTHONHOME path @<:@default=check@:>@])],
  [],
  [with_pythonhome=check])

AC_ARG_WITH([pythonconfig], 
  [AS_HELP_STRING([--with-pythonconfig=@<:@PATH@:>@],
    [Use python-config to determine the compiler and linker flags @<:@default=no@:>@])],
  [],
  [with_pythonconfig=check])

AC_ARG_WITH([pythonpath], 
  [AS_HELP_STRING([--with-pythonpath=@<:@check,no,PATH@:>@],
    [Hard-code the PYTHONPATH path @<:@default=check@:>@])],
  [],
  [with_pythonpath=check])


AX_PYTHON_DEFAULT( )
AX_PYTHON_ENABLE( "/usr/bin" )
AX_PYTHON_WITH( )
AX_PYTHON_PATH( )
AX_PYTHON_VERSION_ENSURE( [2.6] )


if test x"$ax_python_use" = x"true"; then

   if test x"$with_pythonconfig" = xcheck; then
      AC_MSG_CHECKING([for ${PYTHON}-config])
      if test -x "${PYTHON}-config" ; then
          AC_MSG_RESULT(yes)
          PYTHONCONFIG=${PYTHON}-config
          AC_SUBST([PYTHONCONFIG], [${PYTHONCONFIG}])
          ax_pythonconfig_cflags=`${PYTHONCONFIG} --cflags --embed 2> /dev/null || ${PYTHONCONFIG} --cflags`
          AC_SUBST([PYTHON_CSPEC], [${ax_pythonconfig_cflags}])
          AC_MSG_NOTICE([PYTHON_CSPEC=${ax_pythonconfig_cflags}])
          ax_pythonconfig_lflags=`${PYTHONCONFIG} --ldflags --embed 2> /dev/null || ${PYTHONCONFIG} --ldflags`
          AC_SUBST([PYTHON_LSPEC], [${ax_pythonconfig_lflags}])
          AC_MSG_NOTICE([PYTHON_LSPEC=${ax_pythonconfig_lflags}])
      else      
          AC_MSG_RESULT(no)
          AX_PYTHON_CSPEC( )
          AX_PYTHON_LSPEC( )
      fi
   elif test x"$with_pythonconfig" != xno; then
      ax_pythonconfig_cflags=`(${with_pythonconfig} --cflags --embed 2> /dev/null || ${with_pythonconfig} --cflags) | awk '$1 !~ /Usage/ {print $0}'`
      AC_SUBST([PYTHON_CSPEC], [${ax_pythonconfig_cflags}])
      AC_MSG_NOTICE([PYTHON_CSPEC=${ax_pythonconfig_cflags}])
      ax_pythonconfig_lflags=`(${with_pythonconfig} --ldflags --embed 2> /dev/null || ${with_pythonconfig} --ldflags) | awk '$1 !~ /Usage/ {print $0}'`
      AC_SUBST([PYTHON_LSPEC], [${ax_pythonconfig_lflags}])
      AC_MSG_NOTICE([PYTHON_LSPEC=${ax_pythonconfig_lflags}])
   else
      AX_PYTHON_CSPEC( )
      AX_PYTHON_LSPEC( )
   fi

   AC_MSG_CHECKING([for Numpy])

   CPPFLAGS_SAVE="$CPPFLAGS"
   LDFLAGS_SAVE="$LDFLAGS"

   CPPFLAGS="$PYTHON_CSPEC $CPPFLAGS"
   NUMPY_INCLUDE_DIR=`echo "import numpy; print(numpy.get_include())" | $PYTHON - 2>/dev/null`
   NUMPY_INCLUDE_DIR_FULL="${NUMPY_INCLUDE_DIR}/numpy"
   AC_SUBST(NUMPY_INCLUDE_DIR)
   AC_SUBST(NUMPY_INCLUDE_DIR_FULL)
   AC_CHECK_HEADER([${NUMPY_INCLUDE_DIR}/numpy/arrayobject.h], [NUMPY_HEADER=yes], [AC_MSG_WARN([Numpy extension header not found])], [#include "Python.h"])
   AC_SUBST(NUMPY_HEADER)
   AM_CONDITIONAL(PYTHON_USE, test x$NUMPY_HEADER = xyes)
   if test x"$NUMPY_HEADER" = x"yes"; then
       AC_CHECK_PYTHON_MODULE(numpy, __version__)
       AS_IF([test x"$PYTHON_numpy" = x1],
           [AC_DEFINE_UNQUOTED([PYTHON_NUMPY_VERSION], [${PYTHON_numpy_VERSION}], [NUMPY VERSION])
            AS_VERSION_COMPARE([$PYTHON_numpy_VERSION], [1.6.9], [], 
               [], [AC_DEFINE([HAVE_NUMPY_ASTYPE_COPY], [1], [Define if the copy keyword is available for the astype method in numpy])])], [])
       AX_PYTHON_VERSION_CHECK([3.3], [AC_DEFINE([HAVE_PYTHON_UNICODE], [1], [Define if the Py_UNICODE 32-bit object type is available])], [])
       AC_DEFINE([HAVE_PYTHON], [1], [Define if we are compiling with python])
       AC_DEFINE_UNQUOTED([NUMPY_HEADER_FILE], [<${NUMPY_INCLUDE_DIR}/numpy/arrayobject.h>], [Set to the header file to include for numpy])
       AS_IF([test "x$with_pythonhome" = xcheck],
             [PYTHONBINPATH=`dirname $PYTHON`
              AC_SUBST(PYTHONBINPATH)
              PYTHONHOMEPATH=`dirname $PYTHONBINPATH`
              AC_SUBST(PYTHONHOMEPATH)
              AC_DEFINE_UNQUOTED([PYTHONHOMEPATH], ["${PYTHONHOMEPATH}"], [Set to the PYTHONHOME])],
             [AS_IF([test "x$with_pythonhome" != xno],
                    [AC_DEFINE_UNQUOTED([PYTHONHOMEPATH], ["${with_pythonhome}"], [Set to the PYTHONHOME])], [])])
       AS_IF([test "x$with_pythonpath" = xcheck],
             [PYTHONSEARCHPATH=`echo "import os; import sys; outpath=';'.join(sys.path) if os.name == 'nt' else ':'.join(sys.path); print(outpath)" | $PYTHON - 2>/dev/null`
              AC_SUBST(PYTHONSEARCHPATH)
              AC_DEFINE_UNQUOTED([PYTHONSEARCHPATH], ["${PYTHONSEARCHPATH}"], [Set to the PYTHONPATH])],
             [AS_IF([test "x$with_pythonhome" != xno],
                    [AC_DEFINE_UNQUOTED([PYTHONSEARCHPATH], ["${with_pythonpath}"], [Set to the PYTHONPATH])], [])])
   fi
   LDFLAGS="$PYTHON_LSPEC $LDFLAGS"
   python_lib_base=`echo $PYTHON_LSPEC | awk '{for(i=1; i <= NF; i += 1) {if(\$i ~ /-lpython/) print substr(\$i,3);}}'`
   AC_CHECK_LIB([$python_lib_base],[Py_SetPath],[AC_DEFINE([HAVE_PY_SETPATH],[1],[Define if Py_SetPath is present])],[])
   CPPFLAGS="$CPPFLAGS_SAVE"
   LDFLAGS="$LDFLAGS_SAVE"

fi


AC_ARG_WITH([RHOME], 
  [AS_HELP_STRING([--with-RHOME=PATH],
    [prefix directory for R. The executable for R should be in PATH/bin/R, the libraries in PATH/lib and the header files in PATH/include.])],
  [R_home_string=${with_RHOME}
   LDFLAGS="$LDFLAGS -L${with_RHOME}/lib"
   CPPFLAGS="$CPPFLAGS -I${with_RHOME}/include"],
  [with_RHOME=check])

AC_ARG_WITH([R-include], 
  [AS_HELP_STRING([--with-R-include=PATH],
    [PATH to the R include files])],
  [if test "x$with_R_include" != xno ; then
      CPPFLAGS="$CPPFLAGS -I${with_R_include}/include" ;
   fi],
  [])

AC_ARG_WITH([R-lib], 
  [AS_HELP_STRING([--with-R-lib=PATH],
    [PATH to the directory containing the libR.so library file])],
  [if test "x$with_R_lib" != xno ; then
      LDFLAGS="$LDFLAGS -L${with_R_lib}/lib" ;
   fi],
  [])

AS_IF([test "x$RHOME" = xcheck],
  [AC_CHECK_PROG(ISRBIN,R,yes)
   AS_IF([test "x$ISRBIN" = xyes],
      [RHOME=`R RHOME`],
      [RHOME=]
   )], [])

is_R_found=1

AC_CHECK_LIB([R],[Rf_initEmbeddedR],[],[is_R_found=0; AC_MSG_WARN([Cannot find libR])])

AC_CHECK_HEADERS([Rinternals.h],[],[is_R_found=0; AC_MSG_WARN([Cannot find Rinternals.h])])

AC_CHECK_HEADERS([Rembedded.h],[],[is_R_found=0; AC_MSG_WARN([Cannot find Rembedded.h])])

#AS_IF([test "x$R_include_path_string" != x && test "x$R_lib_path_string" != x],
#  [AS_IF([test -f "${R_include_path_string}/Rinternals.h" && test -f "${R_include_path_string}/Rembedded.h" && test -f "${R_lib_path_string}/libR.so"],
#       [AC_SUBST([LIBR],"-L${R_lib_path_string} -lR")
#        AC_SUBST([INCLUDER],"-I ${R_include_path_string}")
#        AC_DEFINE([HAVE_R], [1],
#                  [Define if you have R])
#        AS_IF([test "x$RHOME" != x],
#            [AC_DEFINE([RHOME], [${RHOME}], [Define the RHOME environment variable])], [])
#       ],
#       AC_MSG_WARN([Did not find Rinternals.h in directory "${R_include_path_string}" or did not find Rembedded.h in directory "${R_include_path_string}" or did not find libR.so in directory "${R_lib_path_string}"]))],[])
#
#
#AS_IF([test "x$R_include_path_string" = x && test "x$R_lib_path_string" = x && test "x$RHOME" != x],
#  [AS_IF([test -f "${RHOME}/include/Rinternals.h" && test -f "${RHOME}/include/Rembedded.h" && test -f "${RHOME}/lib/libR.so"],
#       [AC_SUBST([LIBR],"-L${RHOME}/lib -lR")
#        AC_SUBST([INCLUDER],"-I ${RHOME}/include")
#        AC_DEFINE([HAVE_R], [1],
#                  [Define if you have R])
#        AS_IF([test "x$RHOME" != x],
#            AC_DEFINE_UNQUOTED([RHOME], ["${RHOME}"], [Define the RHOME environment variable]), [])
#       ],
#       AC_MSG_WARN([Did not find Rinternals.h in directory "${RHOME}/include" or did not find Rembedded.h in directory "${RHOME}/include" or did not find libR.so in directory "${RHOME}/lib"]))],[])

AM_CONDITIONAL(R_USE, [test "$is_R_found" -eq 1])


# Checks for header files.
AC_CHECK_HEADERS([float.h malloc.h memory.h stdlib.h string.h strings.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_MALLOC
AC_FUNC_REALLOC
AC_CHECK_FUNCS([alarm floor memset pow sqrt strchr])

AC_F77_WRAPPERS
AC_FC_WRAPPERS

AC_CONFIG_FILES([Makefile
                 src/Makefile
                 USERFUNCS/src/Makefile
                 USERLIBS/src/Makefile])
AC_OUTPUT
