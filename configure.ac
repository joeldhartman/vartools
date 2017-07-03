#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

#AC_PREREQ([2.69])
AC_INIT([vartools], [1.35], [jhartman@astro.princeton.edu], [vartools], [http://www.astro.princeton.edu/~jhartman/vartools.html])
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
  [AS_HELP_STRING([--with-cfitsio],
    [support input and output of fits files @<:@default=check@:>@])],
  [],
  [with_cfitsio=check])

LIBCFITSIO=
AS_IF([test "x$with_cfitsio" != xno],
  [AC_CHECK_LIB([cfitsio], [fits_init_cfitsio],
    [AC_SUBST([LIBCFITSIO], ["-lcfitsio"])
     AC_DEFINE([HAVE_CFITSIO], [1],
               [Define if you have libcfitsio])
    ],
    [if test "x$with_readcfitsio" != xcheck; then
       AC_MSG_FAILURE(
         [--with-cfitsio was given, but test for cfitsio failed])
     fi
    ])])


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
          ax_pythonconfig_cflags=`${PYTHONCONFIG} --cflags`
          AC_SUBST([PYTHON_CSPEC], [${ax_pythonconfig_cflags}])
          AC_MSG_NOTICE([PYTHON_CSPEC=${ax_pythonconfig_cflags}])
          ax_pythonconfig_lflags=`${PYTHONCONFIG} --ldflags`
          AC_SUBST([PYTHON_LSPEC], [${ax_pythonconfig_lflags}])
          AC_MSG_NOTICE([PYTHON_LSPEC=${ax_pythonconfig_lflags}])
      else      
          AC_MSG_RESULT(no)
          AX_PYTHON_CSPEC( )
          AX_PYTHON_LSPEC( )
      fi
   elif test x"$with_pythonconfig" != xno; then
      ax_pythonconfig_cflags=`${with_pythonconfig} --cflags`
      AC_SUBST([PYTHON_CSPEC], [${ax_pythonconfig_cflags}])
      AC_MSG_NOTICE([PYTHON_CSPEC=${ax_pythonconfig_cflags}])
      ax_pythonconfig_lflags=`${with_pythonconfig} --ldflags`
      AC_SUBST([PYTHON_LSPEC], [${ax_pythonconfig_lflags}])
      AC_MSG_NOTICE([PYTHON_LSPEC=${ax_pythonconfig_lflags}])
   else
      AX_PYTHON_CSPEC( )
      AX_PYTHON_LSPEC( )
   fi

   AC_MSG_CHECKING([for Numpy])

   CPPFLAGS_SAVE="$CPPFLAGS"

   CPPFLAGS="$PYTHON_CSPEC $CPPFLAGS"
   NUMPY_INCLUDE_DIR=`echo "import numpy; print numpy.get_include()" | $PYTHON - 2>/dev/null`
   AC_SUBST(NUMPY_INCLUDE_DIR)
   AC_CHECK_HEADER([${NUMPY_INCLUDE_DIR}/numpy/arrayobject.h], [NUMPY_HEADER=yes], [AC_MSG_WARN([Numpy extension header not found])], [#include "Python.h"])
   AC_SUBST(NUMPY_HEADER)
   AM_CONDITIONAL(PYTHON_USE, test x$NUMPY_HEADER = xyes)
   if test x"$NUMPY_HEADER" = x"yes"; then
       AC_CHECK_PYTHON_MODULE(numpy, __version__)
       AS_IF([test x"$PYTHON_numpy" = x1],
           [AC_DEFINE_UNQUOTED([PYTHON_NUMPY_VERSION], [${PYTHON_numpy_VERSION}], [NUMPY VERSION])
            AS_VERSION_COMPARE([$PYTHON_numpy_VERSION], [1.6.9], [], 
               [], [AC_DEFINE([HAVE_NUMPY_ASTYPE_COPY], [1], [Define if the copy keyword is available for the astype method in numpy])])], [])
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
   fi
   CPPFLAGS="$CPPFLAGS_SAVE"
fi



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
