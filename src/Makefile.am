libdir = $(USERLIBSDIR)
lib_LTLIBRARIES = libvartools.la
libvartools_la_SOURCES = vartools.c

if PYTHON_USE
lib_LTLIBRARIES += libvartoolsrunpython.la
libvartoolsrunpython_la_SOURCES = runpython.c runpython.h
libvartoolsrunpython_la_CFLAGS = @PYTHON_CSPEC@ -O0
libvartoolsrunpython_la_LDFLAGS = @PYTHON_LSPEC@
endif

if R_USE
lib_LTLIBRARIES += libvartoolsrunR.la
libvartoolsrunR_la_SOURCES = runR.c runR.h
#libvartoolsrunR_la_CFLAGS = @INCLUDER@
#libvartoolsrunR_la_LDFLAGS = @LIBR@
endif

# Pipeline library: all vartools sources except main.c, plus pipeline_api.c.
# The vartools binary links against this library instead of re-compiling
# those sources, so every symbol is compiled exactly once.
PIPELINE_SOURCES = alarm.c aov.c aov_harm.c amoeba.c analytic.c binlc.c \
 binarylcio.c chi2rms.c copylc.c decorr.c detrendtfa.c doublelinklist.c \
 error.c fluxtomag.c getampthresh.c ifelse.c initcommands.c inputoutput.c \
 jstet.c killharm.c lombscargle.c mrqmin.c mysort.c nonlinfit.c nr.c \
 parsecommandline.c processcommand.c rescalesig.c runblsfixdurtc.c \
 runblsfixper.c runbls_sn.c sigclip.c starspot.c statistics.c svdcmp.c \
 sysrem.c transitmodel.c util.c transit.c GetBLSsignaltorednoise.c \
 autocorrelation.c harmonicfilter.c injectharm.c dftclean.c injecttransit.c \
 saverestorelc.c saverestorelistdata.c outcolumns.c medianfilter.c \
 findblends.c microlens.c fittrap.c addnoise.c example.c parseinputlist.c \
 parselc.c converttime.c userlib.c linfit.c help.c resample.c \
 restricttimes.c mcmcfit.c wwz.c callrunpython.c callrunR.c fft.c \
 matchcommand.c pipeline_api.c \
 analytic.h functions.h nonlinfit.h userlib.h binarylcio.h ifelse.h \
 outcolumn.h vartools_functionpointers.h commands.h mcmcfit.h programdata.h \
 vartools.h doublelinklist.h mysort.h userfunc.h statistics.h OutText.h

lib_LTLIBRARIES += libvartoolspipeline.la
libvartoolspipeline_la_SOURCES = $(PIPELINE_SOURCES)
libvartoolspipeline_la_CPPFLAGS = $(INCLUDECFITSIO) $(INCLUDECSPICE) \
 $(LTDLINCL) \
 -D VARTOOLSLIB_USERLIBSDIR=\"$(USERLIBSDIR)\" \
 -D VARTOOLSLIB_USERFUNCSDIR=\"$(USERFUNCSDIR)\"
libvartoolspipeline_la_LDFLAGS = -version-info 1:0:0
libvartoolspipeline_la_LIBADD = libvartools.la \
 $(LIBCFITSIO) $(LIBGSL) $(LIBPTHREAD) $(LIBCSPICE) $(LIBLTDL)
libvartoolspipeline_la_DEPENDENCIES = libvartools.la
if PYTHON_USE
libvartoolspipeline_la_LIBADD += libvartoolsrunpython.la
libvartoolspipeline_la_DEPENDENCIES += libvartoolsrunpython.la
endif
if R_USE
libvartoolspipeline_la_LIBADD += libvartoolsrunR.la
libvartoolspipeline_la_DEPENDENCIES += libvartoolsrunR.la
endif

bin_PROGRAMS = vartools
vartools_SOURCES = main.c \
 analytic.h functions.h nonlinfit.h userlib.h binarylcio.h ifelse.h \
 outcolumn.h vartools_functionpointers.h commands.h mcmcfit.h programdata.h \
 vartools.h doublelinklist.h mysort.h userfunc.h statistics.h OutText.h
vartools_CPPFLAGS = $(INCLUDECFITSIO) $(INCLUDECSPICE) $(LTDLINCL) -D VARTOOLSLIB_USERLIBSDIR=\"$(USERLIBSDIR)\" -D VARTOOLSLIB_USERFUNCSDIR=\"$(USERFUNCSDIR)\"
vartools_LDADD = libvartoolspipeline.la libvartools.la
vartools_DEPENDENCIES = $(LTDLDEPS) libvartoolspipeline.la libvartools.la
if PYTHON_USE
vartools_LDADD += libvartoolsrunpython.la
vartools_DEPENDENCIES += libvartoolsrunpython.la
endif
if R_USE
vartools_LDADD += libvartoolsrunR.la
vartools_DEPENDENCIES += libvartoolsrunR.la
endif
vartools_LDADD += $(LIBCFITSIO) $(LIBGSL) $(LIBPTHREAD) $(LIBCSPICE) $(LIBLTDL)
