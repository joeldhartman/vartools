libdir = $(USERLIBSDIR)
lib_LTLIBRARIES = libvartools.la
libvartools_la_SOURCES = vartools.c

if PYTHON_USE
lib_LTLIBRARIES += libvartoolsrunpython.la
libvartoolsrunpython_la_SOURCES = runpython.c runpython.h
libvartoolsrunpython_la_CFLAGS = @PYTHON_CSPEC@
libvartoolsrunpython_la_LDFLAGS = @PYTHON_LSPEC@
endif

bin_PROGRAMS = vartools
vartools_SOURCES = alarm.c aov.c aov_harm.c amoeba.c analytic.c binlc.c binarylcio.c chi2rms.c copylc.c decorr.c detrendtfa.c doublelinklist.c error.c fluxtomag.c getampthresh.c ifelse.c initcommands.c inputoutput.c jstet.c killharm.c lombscargle.c main.c mrqmin.c mysort.c nonlinfit.c nr.c parsecommandline.c processcommand.c rescalesig.c runblsfixdurtc.c runblsfixper.c runbls_sn.c sigclip.c starspot.c statistics.c svdcmp.c sysrem.c transitmodel.c util.c transit.c GetBLSsignaltorednoise.c autocorrelation.c harmonicfilter.c injectharm.c dftclean.c injecttransit.c saverestorelc.c saverestorelistdata.c outcolumns.c medianfilter.c findblends.c microlens.c fittrap.c addnoise.c example.c parseinputlist.c parselc.c converttime.c userlib.c linfit.c help.c resample.c restricttimes.c mcmcfit.c wwz.c callrunpython.c analytic.h functions.h  nonlinfit.h    userlib.h binarylcio.h      ifelse.h     outcolumn.h    vartools_functionpointers.h commands.h        mcmcfit.h    programdata.h  vartools.h doublelinklist.h  mysort.h     userfunc.h statistics.h OutText.h
vartools_CPPFLAGS = $(INCLUDECSPICE) $(LTDLINCL) -D VARTOOLSLIB_USERLIBSDIR=\"$(USERLIBSDIR)\" -D VARTOOLSLIB_USERFUNCSDIR=\"$(USERFUNCSDIR)\"
if PYTHON_USE
vartools_LDADD = libvartools.la libvartoolsrunpython.la $(LIBCFITSIO) $(LIBGSL) $(LIBPTHREAD) $(LIBCSPICE) $(LIBLTDL)
vartools_DEPENDENCIES = $(LTDLDEPS) libvartools.la libvartoolsrunpython.la
else
vartools_LDADD = libvartools.la $(LIBCFITSIO) $(LIBGSL) $(LIBPTHREAD) $(LIBCSPICE) $(LIBLTDL)
vartools_DEPENDENCIES = $(LTDLDEPS) libvartools.la
endif
