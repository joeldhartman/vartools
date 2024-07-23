userlibsdir = $(USERLIBSDIR)
userlibs_LTLIBRARIES = magadd.la jktebop.la transitbisec.la macula.la medlc.la ftuneven.la
EXTRA_LTLIBRARIES = libjktebop.la libmacula.la

if HAVE_GSL
userlibs_LTLIBRARIES += fastchi2.la splinedetrend.la
endif

ftuneven_la_SOURCES = ftuneven.c ftuneven.h 
ftuneven_la_LIBADD = $(abs_top_srcdir)/src/libvartools.la
ftuneven_la_LDFLAGS = -module

magadd_la_SOURCES = magadd.c magadd.h 
magadd_la_LIBADD = $(abs_top_srcdir)/src/libvartools.la
magadd_la_LDFLAGS = -module

libjktebop_la_SOURCES = jktebop_lib.f
libjktebop_la_LIBADD = $(FLIBS)

nodist_EXTRA_jktebop_la_SOURCES = dummy1.f
jktebop_la_SOURCES = jktebop.c jktebop.h jktebop_lib.h
jktebop_la_LIBADD = libjktebop.la $(abs_top_srcdir)/src/libvartools.la $(FLIBS)
jktebop_la_LDFLAGS = -module

transitbisec_la_SOURCES = transitbisec.c transitbisec.h
transitbisec_la_LIBADD = $(abs_top_srcdir)/src/libvartools.la
transitbisec_la_LDFLAGS = -module

fastchi2_la_SOURCES = fastchi2.c fastchi2.h fastchi2_lib.c fastchi2_lib.h
fastchi2_la_LIBADD = $(abs_top_srcdir)/src/libvartools.la $(GSL_LIBS)
fastchi2_la_CPPFLAGS = $(GSL_CFLAGS)
fastchi2_la_LDFLAGS = -module

splinedetrend_la_SOURCES = splinedetrend.c splinedetrend.h
splinedetrend_la_LIBADD = $(abs_top_srcdir)/src/libvartools.la $(GSL_LIBS)
splinedetrend_la_CPPFLAGS = $(GSL_CFLAGS)
splinedetrend_la_LDFLAGS = -module

libmacula_la_SOURCES = macula.f90
libmacula_la_LIBADD = $(FLIBS)

nodist_EXTRA_macula_la_SOURCES = dummy2.f90
macula_la_SOURCES = macula_c.c macula.h
macula_la_LIBADD = libmacula.la $(abs_top_srcdir)/src/libvartools.la $(FLIBS)
macula_la_LDFLAGS = -module

medlc_la_SOURCES = medlc.c medlc.h
medlc_la_LIBADD = $(abs_top_srcdir)/src/libvartools.la
medlc_la_LDFLAGS = -module
