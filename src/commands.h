/*     This file is part of VARTOOLS version 1.31                      */
/*                                                                           */
/*     VARTOOLS is free software: you can redistribute it and/or modify      */
/*     it under the terms of the GNU General Public License as published by  */
/*     the Free Software Foundation, either version 3 of the License, or     */
/*     (at your option) any later version.                                   */
/*                                                                           */
/*     This program is distributed in the hope that it will be useful,       */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*     GNU General Public License for more details.                          */
/*                                                                           */
/*     You should have received a copy of the GNU General Public License     */
/*     along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*                                                                           */
/*     Copyright 2007, 2008, 2009  Joel Hartman                              */
/*                                                                           */
/*     This file is part of VARTOOLS version 1.152                      */
/*                                                                           */
/*     VARTOOLS is free software: you can redistribute it and/or modify      */
/*     it under the terms of the GNU General Public License as published by  */
/*     the Free Software Foundation, either version 3 of the License, or     */
/*     (at your option) any later version.                                   */
/*                                                                           */
/*     This program is distributed in the hope that it will be useful,       */
/*     but WITHOUT ANY WARRANTY; without even the implied warranty of        */
/*     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the         */
/*     GNU General Public License for more details.                          */
/*                                                                           */
/*     You should have received a copy of the GNU General Public License     */
/*     along with this program.  If not, see <http://www.gnu.org/licenses/>. */
/*                                                                           */
/*     Copyright 2007, 2008, 2009  Joel Hartman                              */
/*                                                                           */
#include "../config.h"
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "outcolumn.h"
#include "analytic.h"
#include "vt_param_macros.h"
#include "mysort.h"
#include "ifelse.h"

/* This is the general header for the vartools program by J. Hartman */

/*** translate some of the variables set by config.h into variables used
     by vartools ***/
#ifdef HAVE_PTHREAD
#define PARALLEL 1
#endif

#ifdef PACKAGE_VERSION
#define VARTOOLS_VERSION PACKAGE_VERSION
#else
#define VARTOOLS_VERSION "0"
#endif

#ifdef HAVE_CFITSIO
#define USECFITSIO 1
#endif

#ifdef HAVE_GSL
#define _HAVE_GSL 1
#endif

#define _USEBINARY_LC 1

#ifdef HAVE_CSPICE
#define _HAVE_CSPICE 1
#endif

#ifdef HAVE_PYTHON
#define _HAVE_PYTHON 1
#endif

#ifdef HAVE_LIBR
#ifdef HAVE_RINTERNALS_H
#ifdef HAVE_REMBEDDED_H
#define _HAVE_R 1
#endif
#endif
#endif

#define DYNAMICLIB 1

/*
#ifdef HAVE_DYNAMICLIB 
#ifdef HAVE_DLFCN_H 
#define DYNAMICLIB 1
#endif
#endif

#ifndef DYNAMICLIB
#ifdef ISWINDOWS
#define DYNAMICLIB 1
#endif
#endif
*/

#define MAXLEN 2048
#define MINUTESPERDAY 1440.

#define RMSTHYCUT 1.0
#define DEFAULT_JDTOL 0.00001

#define TINY 1.0e-20
#define SIG_CLIP 10.0


#define SQR(A) ((A) * (A))
#define ABS_(A) ((A) > 0 ? (A) : (-(A)))

#define MAX_(A,B) ((A) > (B) ? (A) : (B))
#define MIN_(A,B) ((A) < (B) ? (A) : (B))
#define SIGN(A,B) ((B) >= 0.0 ? fabs(A) : -fabs(A))

#define dmax(A,B) ((A)>=(B)?(A):(B))
#define dmin(A,B) ((A)<=(B)?(A):(B))

#define CNUM_ALARM 0
#define CNUM_AOV 1
#define CNUM_HARMAOV 2
#define CNUM_AUTOCORR 3
#define CNUM_BINLC 4
#define CNUM_BLS 5
#define CNUM_FIXPERBLS 6
#define CNUM_CHANGEERROR 7
#define CNUM_CHI2_NOBIN 8
#define CNUM_CHI2_BIN 9
#define CNUM_CLIP 10
#define CNUM_DECORR 11
#define CNUM_DFTCLEAN 12
#define CNUM_ENSEMBLERESCALESIG 13
#define CNUM_DIFFFLUXTOMAG 14
#define CNUM_GETLSAMPTHRESH 15
#define CNUM_INJECTHARM 16
#define CNUM_INJECTTRANSIT 17
#define CNUM_JSTET 18
#define CNUM_KILLHARM 19
#define CNUM_LS 20
#define CNUM_MANDELAGOLTRANSIT 21
#define CNUM_OUTPUTLCS 22
#define CNUM_PHASE 23
#define CNUM_RESCALESIG 24
#define CNUM_RMS_NOBIN 25
#define CNUM_RMS_BIN 26
#define CNUM_SOFTENEDTRANSIT 27
#define CNUM_STARSPOT 28
#define CNUM_SYSREM 29
#define CNUM_TFA 30
#define CNUM_TFA_SR 31
#define CNUM_SAVELC 32
#define CNUM_RESTORELC 33
#define CNUM_MEDIANFILTER 34
#define CNUM_FINDBLENDS 35
#define CNUM_MICROLENS 36
#define CNUM_FLUXTOMAG 37
#define CNUM_USERCOMMAND 38
#define CNUM_ADDNOISE 39
#define CNUM_CONVERTTIME 40
#define CNUM_EXPRESSION 41
#define CNUM_CHANGEVARIABLE 42
#define CNUM_LINFIT 43
#define CNUM_STATS 44
#define CNUM_IF 45
#define CNUM_RESTRICTTIMES 46
#define CNUM_NONLINFIT 47
#define CNUM_WWZ 48
#define CNUM_COPYLC 49
#define CNUM_RESAMPLE 50
#define CNUM_BLSFIXDURTC 51
#define CNUM_FOURIERFILTER 52
#define CNUM_PYTHON 53
#define CNUM_RESTORETIMES 54
#define CNUM_FFT 55
#define CNUM_R 56
#define CNUM_MATCHCOMMAND 57
#define CNUM_BLSFIXPERDURTC 58
#define CNUM_ADDFITSKEYWORD 59
#define CNUM_SORTLC 60
#define CNUM_PRINT 61
#define CNUM_PDM 62
#define CNUM_FTP 63

#define TOT_CNUMS 61

#define PERTYPE_AOV 0
#define PERTYPE_LS 1
#define PERTYPE_BOTH 2
#define PERTYPE_BLS 3
#define PERTYPE_SPECIFIED 4
#define PERTYPE_FIX 5
#define PERTYPE_UNIFORMRAND 6
#define PERTYPE_LOGRAND 7
#define PERTYPE_UNIFORMRANDFREQ 8
#define PERTYPE_LOGRANDFREQ 9
#define PERTYPE_INJECTHARM 10
#define PERTYPE_FIXCOLUMN 11
#define PERTYPE_AUTOFIND 12
#define PERTYPE_EXPR 13
#define PERTYPE_VAR 14
#define PERTYPE_PDM 15
#define PERTYPE_FTP 16

#define KILLHARM_OUTTYPE_DEFAULT 0
#define KILLHARM_OUTTYPE_AMPPHASE 1
#define KILLHARM_OUTTYPE_AMPRADPHASE 2
#define KILLHARM_OUTTYPE_RPHI 3
#define KILLHARM_OUTTYPE_RRADPHI 4

#define INJECTTR_IDX_PERIOD 0
#define INJECTTR_IDX_RP 1
#define INJECTTR_IDX_MP 2
#define INJECTTR_IDX_PHASE 3
#define INJECTTR_IDX_SINI 4
#define INJECTTR_IDX_E 5
#define INJECTTR_IDX_H 5
#define INJECTTR_IDX_OMEGA 6
#define INJECTTR_IDX_K 6
#define INJECTTR_IDX_MSTAR 7
#define INJECTTR_IDX_RSTAR 8
#define INJECTTR_IDX_DILUTE 9
#define INJECTTR_IDX_LD 10

#define LINEWRAP_LENGTH 80
#define TAB_SPACE_SIZE 4

#define TIMETYPE_MJD 0
#define TIMETYPE_JD 1
#define TIMETYPE_HJD 2
#define TIMETYPE_BJD 3

#define TIMESYSTEM_UTC 0
#define TIMESYSTEM_TDB 1

#define VARTOOLS_STATSTYPE_MEAN 0
#define VARTOOLS_STATSTYPE_WEIGHTEDMEAN 1
#define VARTOOLS_STATSTYPE_MEDIAN 2
#define VARTOOLS_STATSTYPE_STDDEV 3
#define VARTOOLS_STATSTYPE_MEDDEV 4
#define VARTOOLS_STATSTYPE_MEDMEDDEV 5
#define VARTOOLS_STATSTYPE_MAD 6
#define VARTOOLS_STATSTYPE_KURTOSIS 7
#define VARTOOLS_STATSTYPE_SKEWNESS 8
#define VARTOOLS_STATSTYPE_PERCENTILE 9
#define VARTOOLS_STATSTYPE_MAXIMUM 10
#define VARTOOLS_STATSTYPE_MINIMUM 11
#define VARTOOLS_STATSTYPE_SUM 12
#define VARTOOLS_STATSTYPE_MEDIAN_WEIGHT 13
#define VARTOOLS_STATSTYPE_PERCENTILE_WEIGHT 14

#define VARTOOLS_ADDNOISE_WHITE 0
#define VARTOOLS_ADDNOISE_COVAR_SQUAREDEXPONENTIAL 1
#define VARTOOLS_ADDNOISE_COVAR_EXPONENTIAL 2
#define VARTOOLS_ADDNOISE_COVAR_MATERN 3
#define VARTOOLS_ADDNOISE_WAVELET 4

#define VARTOOLS_RESAMPLE_NEAREST 0
#define VARTOOLS_RESAMPLE_LINEAR 1
#define VARTOOLS_RESAMPLE_SPLINE 2
#define VARTOOLS_RESAMPLE_SPLINEMONOTONIC 3
#define VARTOOLS_RESAMPLE_BSPLINE 4
#define VARTOOLS_RESAMPLE_MULTIPLE 5

#define VARTOOLS_FOURIERFILTER_FULLSPEC 0
#define VARTOOLS_FOURIERFILTER_HIGHPASS 1
#define VARTOOLS_FOURIERFILTER_LOWPASS 2
#define VARTOOLS_FOURIERFILTER_BANDPASS 3
#define VARTOOLS_FOURIERFILTER_BANDCUT 4

#define VARTOOLS_FREQSTEPTYPE_FREQ 0
#define VARTOOLS_FREQSTEPTYPE_PERIOD 1
#define VARTOOLS_FREQSTEPTYPE_LOGPERIOD 2

#define VARTOOLS_BINLC_BINTYPE_AVERAGE 0
#define VARTOOLS_BINLC_BINTYPE_MEDIAN 1
#define VARTOOLS_BINLC_BINTYPE_WEIGHTEDAVERAGE 2

#define VARTOOLS_BINLC_TIMETYPE_CENTER 0
#define VARTOOLS_BINLC_TIMETYPE_AVERAGE 1
#define VARTOOLS_BINLC_TIMETYPE_MEDIAN 2
#define VARTOOLS_BINLC_TIMETYPE_NOSHRINK 3



#ifndef _OUTTEXTSTRUCTDEFINE
#include "OutText.h"
#endif

typedef struct {
  double sigclip;
  VT_PARAM_COMPANIONS(sigclip);
  int *Nclip;
  int iter;
  VT_PARAM_COMPANIONS(iter);
  int niter;
  VT_PARAM_COMPANIONS(niter);
  int usemedian;
  int markclip;
  char *clipvarname;
  _Variable *clipvar;
  int noinitmark;
} _Clip;

typedef struct {
  double *rescalefactor;
  double a;
  double b;
  double erssigclip;
  double *chi2_old;
  double *chi2_new;
  int usemask;
  _Variable *maskvar;
} _Ensemblerescalesig;

typedef struct {
  double *rescalefactor;
  double *chi2_old;
  double *chi2_new;
  int usemask;
  _Variable *maskvar;
} _Rescalesig;

typedef struct {
  double *chi2val;
  double *wtave;
  int usemask;
  _Variable *maskvar;
} _Chi2_NoBin;

typedef struct {
  int Nbin;
  double *bintimes;
  double **chi2binvals;
  double **wtavebin;
  int usemask;
  _Variable *maskvar;
} _Chi2_Bin;

typedef struct {
  double *rmsval;
  double *ave;
  double *rmsthy;
  int *ngood;
  int usemask;
  _Variable *maskvar;
} _RMS_NoBin;

typedef struct {
  int Nbin;
  double *bintimes;
  double **rmsbinvals;
  double **rmsthybin;
  int usemask;
  _Variable *maskvar;
} _RMS_Bin;

typedef struct {
  double Jstet_time;
  double *jst;
  double *kur;
  double *lst;
  double wkmax;
  char datesname[MAXLEN];
  int usemask;
  _Variable *maskvar;
} _Jstet;

typedef struct {
  double *alarmvals;
  int usemask;
  _Variable *maskvar;
} _Alarm;

typedef struct {
  double start;
  VT_PARAM_COMPANIONS(start);
  double stop;
  VT_PARAM_COMPANIONS(stop);
  double step;
  VT_PARAM_COMPANIONS(step);
  double errsize;
  char outdir[MAXLEN];
  char suffix[10];
  int usemask;
  _Variable *maskvar;
} _Autocorr;

typedef struct {
  double minp;
  double *minp_vals;
  int minp_source;
  _Variable *minp_var;
  _Expression *minp_expr;
  double maxp;
  double *maxp_vals;
  int maxp_source;
  _Variable *maxp_var;
  _Expression *maxp_expr;
  double subsample;
  double *subsample_vals;
  int subsample_source;
  _Variable *subsample_var;
  _Expression *subsample_expr;
  double finetune;
  double *finetune_vals;
  int finetune_source;
  _Variable *finetune_var;
  _Expression *finetune_expr;
  int Npeaks;
  double *aveaov;
  double *rmsaov;
  double **peakperiods;
  double **peakvalues;
  double **peakSNR;
  double **peakFAP;
  int uselog;
  double clip;
  int clipiter;
  int Nbin;
  int Nbin_source;
  int *Nbin_vals;
  _Variable *Nbin_var;
  _Expression *Nbin_expr;
  int operiodogram;
  char outdir[MAXLEN];
  char suffix[5];
  int whiten;
  double **aveaov_whiten;
  double **rmsaov_whiten;
  int fixperiodSNR;
  int fixperiodSNR_pertype;
  int fixperiodSNR_lastaovindex;
  double fixperiodSNR_fixedperiod;
  double **fixperiodSNR_periods;
  double *fixperiodSNR_peakvalues;
  double *fixperiodSNR_peakSNR;
  double *fixperiodSNR_peakFAP;
  OutColumn *fixperiodSNR_linkedcolumn;
  int usemask;
  _Variable *maskvar;
} _Aov;

typedef struct {
  double minp;
  double *minp_vals;
  int minp_source;
  _Variable *minp_var;
  _Expression *minp_expr;
  double maxp;
  double *maxp_vals;
  int maxp_source;
  _Variable *maxp_var;
  _Expression *maxp_expr;
  double subsample;
  double *subsample_vals;
  int subsample_source;
  _Variable *subsample_var;
  _Expression *subsample_expr;
  double finetune;
  double *finetune_vals;
  int finetune_source;
  _Variable *finetune_var;
  _Expression *finetune_expr;
  int Npeaks;
  double *aveaov;
  double *rmsaov;
  double **peakperiods;
  double **peakvalues;
  double **peakSNR;
  double **peakFAP;
  int **peakNharm;
  int Nharm;
  int *Nharm_vals;
  int Nharm_source;
  _Variable *Nharm_var;
  _Expression *Nharm_expr;
  int operiodogram;
  char outdir[MAXLEN];
  char suffix[10];
  int whiten;
  double **aveaov_whiten;
  double **rmsaov_whiten;
  double clip;
  int clipiter;
  int fixperiodSNR;
  int fixperiodSNR_pertype;
  int fixperiodSNR_lastaovindex;
  double fixperiodSNR_fixedperiod;
  double **fixperiodSNR_periods;
  double *fixperiodSNR_peakvalues;
  double *fixperiodSNR_peakSNR;
  double *fixperiodSNR_peakFAP;
  OutColumn *fixperiodSNR_linkedcolumn;
  int usemask;
  _Variable *maskvar;
} _AovHarm;

typedef struct {
  double minp;
  double *minp_vals;
  _Variable *minp_var;
  _Expression *minp_expr;
  int minp_source;
  double maxp;
  double *maxp_vals;
  _Variable *maxp_var;
  _Expression *maxp_expr;
  int maxp_source;
  double subsample;
  double *subsample_vals;
  _Variable *subsample_var;
  _Expression *subsample_expr;
  int subsample_source;
  int Npeaks;
  double **peakperiods;
  double **peakvalues;
  double **peakFAP;
  double **SNRvalues;
  int operiodogram;
  char outdir[MAXLEN];
  char suffix[4];
  int whiten;
  double clip;
  int clipiter;
  int fixperiodSNR;
  int fixperiodSNR_pertype;
  double **fixperiodSNR_periods;
  double *fixperiodSNR_peakvalues;
  double *fixperiodSNR_FAPvalues;
  double *fixperiodSNR_SNRvalues;
  int fixperiodSNR_lastaovindex;
  double fixperiodSNR_fixedperiod;
  OutColumn *fixperiodSNR_linkedcolumn;
  int use_orig_ls;
  int dobootstrapfap;
  int Nbootstrap;
  int usemask;
  _Variable *maskvar;
} _Ls;

typedef struct {
  int N_globalterms;
  int N_lcterms;
  int N_decorrterms;
  int N_decorrterms_total;
  int N_globaldecorr_JD;
  int correctlc;
  int zeropointterm;
  int subtractfirstterm;
  int size_globaldecorrvector;
  char **global_file_names;
  int *globalfile_order;
  int *lc_order;
  int *lc_columns;
  int *order;
  double ***lcdecorr_terms_in;
  double ***decorr_terms;
  double *globaldecorr_JD;
  char **globaldecorr_stringid;
  int *globaldecorr_stringid_idx;
  double **globaldecorr_terms;
  double *chi2val;
  double **b;
  double **b_err;
  int decorr_vector_size;
  int omodel;
  char modeloutdir[MAXLEN];
  char modelsuffix[14];
  int usemask;
  _Variable *maskvar;
} _Decorr;

typedef struct {
  double gain;
  VT_PARAM_COMPANIONS(gain);
  double SNlimit;
  VT_PARAM_COMPANIONS(SNlimit);
  int nbeam;
  VT_PARAM_COMPANIONS(nbeam);
  double maxfreq;
  VT_PARAM_COMPANIONS(maxfreq);
  int outdspec;
  int finddirtypeaks;
  int outwspec;
  int runclean;
  int outcbeam;
  int outcspec;
  int findcleanpeaks;
  char dirtyspec_outdir[MAXLEN];
  char dirtyspec_suffix[16];
  int Npeaks_dirty;
  char wspec_outdir[MAXLEN];
  char wspec_suffix[16];
  char cbeam_outdir[MAXLEN];
  char cbeam_suffix[16];
  char cspec_outdir[MAXLEN];
  char cspec_suffix[16];
  int Npeaks_clean;
  int clipiter_dirty, clipiter_clean;
  int useampspec, verboseout;
  double clip_dirty;
  VT_PARAM_COMPANIONS(clip_dirty);
  double clip_clean;
  VT_PARAM_COMPANIONS(clip_clean);
  double *aveper_dirty, *stdper_dirty;
  double *aveper_noclip_dirty, *stdper_noclip_dirty;
  double *aveper_clean, *stdper_clean;
  double *aveper_noclip_clean, *stdper_noclip_clean;
  double **peakfreqs_dirty;
  double **peakpows_dirty;
  double **peakfreqs_clean;
  double **peakpows_clean;
  double **SNR_dirty;
  double **SNR_clean;
  int usemask;
  _Variable *maskvar;
} _Dftclean;

typedef struct {
  int pertype;
  int lastlsindex;
  int lastaovindex;
  double **periods;
  double *fixedperiods;
  int *fixedperiods_source;
  _Variable **fixedperiods_var;
  _Expression **fixedperiods_expr;
  int Nper;
  int Nharm;
  int Nsubharm;
  double ***subharmA;
  double ***subharmB;
  double ***harmA;
  double ***harmB;
  double **fundA;
  double **fundB;
  double *mean;
  double **amp;
  int omodel;
  int fitonly;
  int outtype;
  char modeloutdir[MAXLEN];
  /* 32 bytes fits both legacy ".killharm.model" and the longer
     ".harmonicfilter.model" emitted when invoked via the new name. */
  char modelsuffix[32];
  double clip;
  VT_PARAM_COMPANIONS(clip);
  /* Output-column prefix selected by the CLI token used to invoke this
     command: "Killharm" if invoked as -Killharm (legacy), "HarmonicFilter"
     if invoked as -harmonicfilter (preferred).  Column names are built as
     "%s_Mean_Mag_%d", column_prefix, cnum etc. */
  char column_prefix[32];
} _Killharm;

typedef struct {
  int pertype;
  int Nharm;
  int Nsubharm;
  double fixperiod;
  _Variable *fixperiod_var;
  _Expression *fixperiod_expr;
  double minp;
  VT_PARAM_COMPANIONS(minp);
  double maxp;
  VT_PARAM_COMPANIONS(maxp);
  double minf;
  VT_PARAM_COMPANIONS(minf);
  double maxf;
  VT_PARAM_COMPANIONS(maxf);
  double **periods;
  double *periodinject;
  int *harm_amptype;
  int *harm_amprel;
  int *harm_phasetype;
  int *harm_phaserel;
  int *subharm_amptype;
  int *subharm_amprel;
  int *subharm_phasetype;
  int *subharm_phaserel;
  double *harm_ampfix;
  int *harm_ampfix_source;
  _Variable **harm_ampfix_var;
  _Expression **harm_ampfix_expr;
  double ***harm_ampspec;
  double *harm_minamp;
  double *harm_maxamp;
  double *harm_phasefix;
  int *harm_phasefix_source;
  _Variable **harm_phasefix_var;
  _Expression **harm_phasefix_expr;
  double ***harm_phasespec;
  double *subharm_ampfix;
  int *subharm_ampfix_source;
  _Variable **subharm_ampfix_var;
  _Expression **subharm_ampfix_expr;
  double ***subharm_ampspec;
  double *subharm_minamp;
  double *subharm_maxamp;
  double *subharm_phasefix;
  int *subharm_phasefix_source;
  _Variable **subharm_phasefix_var;
  _Expression **subharm_phasefix_expr;
  double ***subharm_phasespec;
  double **harm_amp;
  double **harm_phase;
  double **subharm_amp;
  double **subharm_phase;
  int omodel;
  char modeloutdir[MAXLEN];
  char modelsuffix[18];
} _Injectharm;

typedef struct {
  int Nparam;
  int Nld;
  int paramtype[14];
  double *paraminject[14];
  double **paramspec[14];
  _Expression *paramexpr[14];
  _Variable *paramvar[14];
  double paramfix[14];
  double minp;
  VT_PARAM_COMPANIONS(minp);
  double maxp;
  VT_PARAM_COMPANIONS(maxp);
  double minf;
  VT_PARAM_COMPANIONS(minf);
  double maxf;
  VT_PARAM_COMPANIONS(maxf);
  double minRp;
  VT_PARAM_COMPANIONS(minRp);
  double maxRp;
  VT_PARAM_COMPANIONS(maxRp);
  double minMp;
  VT_PARAM_COMPANIONS(minMp);
  double maxMp;
  VT_PARAM_COMPANIONS(maxMp);
  int eomegatype;
  int ldtype;
  int omodel;
  char modeloutdir[MAXLEN];
  char modelsuffix[21];
} _Injecttransit;

typedef struct {
  int pertype;
  int lastlsindex;
  int lastaovindex;
  double a0;
  VT_PARAM_COMPANIONS(a0);
  double b0;
  VT_PARAM_COMPANIONS(b0);
  double chi0;
  VT_PARAM_COMPANIONS(chi0);
  double inclination0;
  VT_PARAM_COMPANIONS(inclination0);
  double alpha0;
  VT_PARAM_COMPANIONS(alpha0);
  double psi00;
  VT_PARAM_COMPANIONS(psi00);
  double mconst0;
  VT_PARAM_COMPANIONS(mconst0);
  double **period;
  double *a;
  double *b;
  double *chi;
  double *inclination;
  double *alpha;
  double *psi0;
  double *mconst;
  double *chisq;
  int fitP, fita, fitb, fitalpha, fiti, fitchi, fitpsi0, fitmconst;
  int correctlc;
  int omodel;
  char modeloutdir[MAXLEN];
  char modelsuffix[16];
  double fixedperiod;
  _Variable *fixedperiod_var;
  _Expression *fixedperiod_expr;
  OutColumn *linkedcolumn;
} _Starspot;

typedef struct {
  int pertype;
  int lastlsindex;
  int lastaovindex;
  int lastblsindex;
  int lastblsfixperindex;
  int dokillharm;
  int nharm, nsubharm, frombls, fromblsfixper;
  double T00, eta0, cval0, delta0, mconst0, period0, per_harm;
  double *period;
  double *T0;
  double *eta;
  double *cval;
  double *delta;
  double *mconst;
  double *chisq;
  double *per_harm_out;
  double *per_harm_spec;
  double **subharmA, **subharmB;
  double **harmA, **harmB;
  double *fundA, *fundB;
  int fitephem, fiteta, fitcval, fitdelta, fitmconst;
  int correctlc;
  int omodel;
  char modeloutdir[MAXLEN];
  char modelsuffix[22];
} _SoftenedTransit;


typedef struct {
  double **u, **v, *fmin, df, qmin, qmax, rmin, rmax, rho;
  double minexpdurfrac, maxexpdurfrac;
  double **bper, **bt0, **bpow, **sde, **snval, **depth, **qtran;
  double **chisqrplus, *chisqrminus, *bperpos, *meanmagval, **fraconenight, **rednoise, **whitenoise, **sigtopink, **qingress, **OOTmag;
  int **nt, **Nt, **Nbefore, **Nafter;
  double **i1_ph, **i2_ph;
  int nf, *nf2, nbins, Npeak, **i1, **i2, operiodogram, isdf_specified;
  int *sizeuv, rflag;
#ifdef PARALLEL
  double **p;
  int *sizepvec;
#else
  double *p;
  int sizepvec;
#endif
  double minper, maxper, timezone;
  int minper_source, maxper_source, nbins_source;
  _Expression *minper_expr, *maxper_expr, *nbins_expr;
  _Variable *minper_var, *maxper_var, *nbins_var;
  double *minper_val, *maxper_val;
  int *nbins_val;
  double *df_val;
  int *nf_val;
  int rmin_source, rmax_source, qmin_source, qmax_source, rho_source;
  int df_source, nf_source, subsample_source;
  int minexpdurfrac_source, maxexpdurfrac_source;
  _Expression *rmin_expr, *rmax_expr, *qmin_expr, *qmax_expr, *rho_expr;
  _Expression *minexpdurfrac_expr, *maxexpdurfrac_expr;
  _Expression *df_expr, *nf_expr;
  _Expression *subsample_expr;
  _Variable *rmin_var, *rmax_var, *qmin_var, *qmax_var, *rho_var;
  _Variable *minexpdurfrac_var, *maxexpdurfrac_var;
  _Variable *df_var, *nf_var;
  _Variable *subsample_var;
  double *rmin_val, *rmax_val, *rho_val, *minexpdurfrac_val, *maxexpdurfrac_val,
    *qmin_val, *qmax_val, *subsample_val, subsample, *A_val;
  char outdir[MAXLEN], suffix[6];
  int omodel;
  char modeloutdir[MAXLEN];
  char modelsuffix[15];
  int ophcurve;
  char ophcurveoutdir[MAXLEN];
  char ophcurvesuffix[15];
  double phmin, phmax, phstep;
  int ojdcurve;
  char ojdcurveoutdir[MAXLEN];
  char ojdcurvesuffix[15];
  double jdstep;
  int correctlc;
  int extraparams;
  int fittrap;
  int nobinnedrms;
  int freqsteptype;
  int adjust_qmin_mindt;
  int reduce_nb;
  int reportharmonics;
  double **srsum, **ressig, **dipsig, **snrextra, **srshift, **srsig, **dsp, **dspg;
  double **freqlow, **freqhigh;
  double **logprob, **peakarea, **peakmean, **peakdev, **lomblog;
  int **ntv;
  double **gezadsp, **ootsig, **trsig, **ootdftf, **ootdfta;
  double **binsignaltonoise, **maxphasegap, **depth1_2tran, **depth2_2tran, **delchi2_2tran;
  double **sr_sec, **srsum_sec, **q_sec, **epoch_sec, **H_sec, **L_sec, **depth_sec;
  int **nt_sec, **Nt_sec;
  double **sigtopink_sec, **deltachi2transit_sec, **binsignaltonoise_sec;
  double **phaseoffset_sec, **harmmean, **harmA, **harmB, **fundA, **fundB, **harmamp;
  double **harmdeltachi2;
  int usemask;
  _Variable *maskvar;
} _Bls;

typedef struct {
  int pertype;
  int lastlsindex;
  int lastaovindex;
  double perfix;
  _Expression *perexpr;
  double **period;
  double **u, **v, fmin, df, qmin, qmax, rmin, rmax;
  double *bpow, *bt0, *sde, *snval, *depth, *qtran, *i1_ph, *i2_ph;
  double *chisqrplus, *chisqrminus, *bperpos, *meanmagval, *fraconenight, *rednoise, *whitenoise, *sigtopink, *qingress, *OOTmag;
  int *nt, *Nt, *Nbefore, *Nafter;
  int nf, nbins, Npeak, *i1, *i2, operiodogram;
  int *sizeuv, rflag;
  double minper, maxper, timezone;
  char outdir[MAXLEN], suffix[6];
  int omodel;
  char modeloutdir[MAXLEN];
  char modelsuffix[17];
  int correctlc;
  int fittrap;
  OutColumn *linkedcolumn;
  int usemask;
  _Variable *maskvar;
} _BlsFixPer;

typedef struct {
  double **u, **v, *fmin, df;
  int durtype;
  double fixdur;
  _Variable *fixdur_var;
  _Expression *fixdur_expr;
  OutColumn *fixdur_linkedcolumn;
  int TCtype;
  double fixTC;
  _Variable *fixTC_var;
  _Expression *fixTC_expr;
  OutColumn *fixTC_linkedcolumn;
  double *inputTC, *inputdur;
  int fixdepth;
  int depthtype;
  double fixdepthval;
  _Variable *fixdepthval_var;
  _Expression *fixdepthval_expr;
  OutColumn *fixdepth_linkedcolumn;
  int qgresstype;
  double qgressval;
  _Variable *qgressval_var;
  _Expression *qgressval_expr;
  OutColumn *fixqgress_linkedcolumn;
  double *inputdepth, *inputqgress;
  double **bper, **bt0, **bpow, **sde, **snval, **depth, **qtran;
  double **chisqrplus, *chisqrminus, *bperpos, *meanmagval, **fraconenight, **rednoise, **whitenoise, **sigtopink, **qingress, **OOTmag;
  int **nt, **Nt, **Nbefore, **Nafter;
  double **i1_ph, **i2_ph;
  int nf, *nf2, Npeak, **i1, **i2, operiodogram;
  VT_PARAM_COMPANIONS(nf);
  int *sizeuv, rflag;
#ifdef PARALLEL
  double **p;
  int *size_p;
#else
  double *p;
  int size_p;
#endif
  double minper;
  VT_PARAM_COMPANIONS(minper);
  double maxper;
  VT_PARAM_COMPANIONS(maxper);
  double timezone;
  char outdir[MAXLEN], suffix[14];
  int omodel;
  char modeloutdir[MAXLEN];
  char modelsuffix[23];
  int ophcurve;
  char ophcurveoutdir[MAXLEN];
  char ophcurvesuffix[23];
  double phmin, phmax, phstep;
  int ojdcurve;
  char ojdcurveoutdir[MAXLEN];
  char ojdcurvesuffix[23];
  double jdstep;
  int correctlc;
  int fittrap;
  int usemask;
  _Variable *maskvar;
} _BlsFixDurTc;

typedef struct {
  double **u, **v;
  int pertype;
  double fixper;
  _Variable *fixper_var;
  _Expression *fixper_expr;
  OutColumn *fixper_linkedcolumn;
  int durtype;
  double fixdur;
  _Variable *fixdur_var;
  _Expression *fixdur_expr;
  OutColumn *fixdur_linkedcolumn;
  int TCtype;
  double fixTC;
  _Variable *fixTC_var;
  _Expression *fixTC_expr;
  OutColumn *fixTC_linkedcolumn;
  double *inputTC, *inputdur, *inputper;
  int fixdepth;
  int depthtype;
  double fixdepthval;
  _Variable *fixdepthval_var;
  _Expression *fixdepthval_expr;
  OutColumn *fixdepth_linkedcolumn;
  int qgresstype;
  double qgressval;
  _Variable *qgressval_var;
  _Expression *qgressval_expr;
  OutColumn *fixqgress_linkedcolumn;
  double *inputdepth, *inputqgress;
  double *depth, *qtran;
  double *chisqrplus, *meanmagval, *fraconenight, *rednoise, *whitenoise, *sigtopink, *qingress, *OOTmag;
  int *nt, *Nt, *Nbefore, *Nafter;
  int *sizeuv, rflag;
  double timezone;
  int omodel;
  char modeloutdir[MAXLEN];
  /* These suffixes used to be sized 23 — too small for the 23-character
     ".blsfixperdurtc.phcurve" / ".blsfixperdurtc.jdcurve" strings
     written in parsecommandline.c (they need 24 bytes incl. NUL).
     Bumped to 32 to leave headroom and avoid the fortify abort when
     BLSFixPerDurTc is invoked with `ophcurve`/`ojdcurve` set. */
  char modelsuffix[32];
  int ophcurve;
  char ophcurveoutdir[MAXLEN];
  char ophcurvesuffix[32];
  double phmin, phmax, phstep;
  int ojdcurve;
  char ojdcurveoutdir[MAXLEN];
  char ojdcurvesuffix[32];
  double jdstep;
  int correctlc;
  int fittrap;
  int usemask;
  _Variable *maskvar;
} _BlsFixPerDurTc;


typedef struct {
  int pertype;
  int lastlsindex;
  int lastaovindex;
  int lastblsindex;
  int t0type;
  double fixperiod;
  _Variable *fixperiod_var;
  _Expression *fixperiod_expr;
  double fixT0;
  _Variable *fixT0_var;
  _Expression *fixT0_expr;
  double phaseTc;
  double **period;
  double **T0;
  OutColumn *period_linkedcolumn;
  OutColumn *T0_linkedcolumn;
  char *phasevarname;
  _Variable *phasevar;
  double startphase;
} _Phase;

typedef struct {
  int medflag, binsize_Nbins_flag, Nbins, firstbinflag, tflag;
  VT_PARAM_COMPANIONS(Nbins);
  double binsize;
  VT_PARAM_COMPANIONS(binsize);
  double firstbin;
  VT_PARAM_COMPANIONS(firstbin);
  int Nvar;
  int *binstats;
  double *pctval;
  char **binvarnames;
  _Variable **binvars;
  char *binvarstring;
  int only_bin_columns;
  int usemask;
  _Variable *maskvar;
  int T0source;
  double t0fixval;
  _Variable *t0fixval_var;
  OutColumn *t0_linkedcolumn;
  double *t0listval;
  char *t0exprstring;
  _Expression *t0expr;
} _Binlc;

typedef struct {
  int pertype, sizesigfile;
  int lastlsindex;
  int harm_specsigflag;
  double **period;
  int Nsubharm;
  int Nharm;
  double minPer;
  double thresh;
  double *ampthresh_scale;
  double *amp;
  char listfilename[MAXLEN], filename[MAXLEN], *line;
  size_t line_size;
  FILE *listfile, *infile;
  int use_orig_ls;
} _GetLSAmpThresh;

typedef struct {
  int pertype;
  int lastlsindex;
  int lastaovindex;
  int lastblsindex;
  int lastblsfixperindex;
  int frombls, fromblsfixper;
  int refititer;
  double P0;
  VT_PARAM_COMPANIONS(P0);
  double T00;
  VT_PARAM_COMPANIONS(T00);
  double r0;
  VT_PARAM_COMPANIONS(r0);
  double a0;
  VT_PARAM_COMPANIONS(a0);
  double inc0;
  VT_PARAM_COMPANIONS(inc0);
  double bimpact0;
  VT_PARAM_COMPANIONS(bimpact0);
  double sin_i0;
  VT_PARAM_COMPANIONS(sin_i0);
  double e0;
  VT_PARAM_COMPANIONS(e0);
  double omega0;
  VT_PARAM_COMPANIONS(omega0);
  double mconst0;
  VT_PARAM_COMPANIONS(mconst0);
  double ldcoeffs0[4];
  int ldcoeffs0_source[4];
  _Variable *ldcoeffs0_var[4];
  _Expression *ldcoeffs0_expr[4];
  double K0;
  VT_PARAM_COMPANIONS(K0);
  double gamma0;
  VT_PARAM_COMPANIONS(gamma0);
  int type;
  int nldcoeff;
  double *period;
  double *T0;
  double *r;
  double *a;
  double *sin_i;
  double *inc;
  double *bimpact;
  double *e;
  double *omega;
  double *mconst;
  double **ldcoeffs;
  double *K;
  double *gamma;
  double *chisq;
  int fitephem, fitr, fita, fitinclterm, fite, fitomega, fitmconst, fitldcoeffs[4], fitRV, fitK, fitgamma;
  int inputinclterm;
  int correctlc;
  int omodel;
  char modeloutdir[MAXLEN];
  char modelsuffix[25];
  char RVmodeloutfile[MAXLEN];
  char RVinputfile[MAXLEN];
  char *modelvarname;
  _Variable *modelvar;
  int ophcurve;
  char ophcurveoutdir[MAXLEN];
  char ophcurvesuffix[27];
  double phmin, phmax, phstep;
  int ojdcurve;
  char ojdcurveoutdir[MAXLEN];
  char ojdcurvesuffix[27];
  double jdstep;
} _MandelAgolTransit;

typedef struct {
  double **magstar;
  double mag_constant1;
  VT_PARAM_COMPANIONS(mag_constant1);
  double offset;
  VT_PARAM_COMPANIONS(offset);
} _DiffFluxtomag;

typedef struct {
  double mag_constant1;
  VT_PARAM_COMPANIONS(mag_constant1);
  double offset;
  VT_PARAM_COMPANIONS(offset);
} _Fluxtomag;

typedef struct {
  double **trends, *trendx, *trendy, **u, **v, *w1, *JD, clipping, pixelsep, *ave_out, *rms_out, **lcx, **lcy;
  double *trend_prior_means, *trend_prior_stds;
  int *is_trend_prior;
  char **stringid;
  int *stringid_idx;

  double ***u2, ***v2, **w2, **a, **b, **m_out;

  int *Njd_mout;

  int use_trend_coeff_priors;
  int use_lc_errors;
  int weight_by_template_stddev;
  double ave_trend_stddev;
  int Ntrend_priors;
  int Ntrends, Njd, correctlc, ocoeff, omodel, Nskip_trend, JDcol_trend, magcol_trend;
  char trend_list_name[MAXLEN], dates_name[MAXLEN], **trend_names, coeff_outdir[MAXLEN], model_outdir[MAXLEN], coeff_suffix[MAXLEN], model_suffix[MAXLEN];
  char trend_prior_list_name[MAXLEN];
  int jdcol_isfromheader, magcol_isfromheader;
  char jdcol_headername[MAXLEN], magcol_headername[MAXLEN];
  int clippingusemedian;
  int clippinguseMAD;
  int usefitmask;
  char *fitmaskvarname;
  _Variable *fitmaskvar;
  int outputfitmask;
  char *outputfitmaskvarname;
  _Variable *outputfitmaskvar;
} _TFA;

typedef struct {
  double ***u_decorr, ***v_decorr, **w1_decorr, **a_decorr;
  double **trends, *trendx, *trendy, **u, **v, *w1, *JD, clipping, pixelsep, *ave_out, *rms_out, **lcx, **lcy, **periods, fixperiod;
  double **m_out, **signal, ***harmterm, **inputsignal, **bstore, **signal_bin;
  double **b, ***decorr_trends, ***u2, **w2, ***v2, **a;
  char **stringid;
  int *stringid_idx;
  int *Njd_mout, **signal_bin_ids, **signal_bin_N;
  int Ntrends, Njd, correctlc, ocoeff, omodel, Nskip_trend, JDcol_trend, magcol_trend;
  char trend_list_name[MAXLEN], dates_name[MAXLEN], **trend_names, coeff_outdir[MAXLEN], model_outdir[MAXLEN], coeff_suffix[MAXLEN], model_suffix[MAXLEN], signal_listname[MAXLEN], **signalfilenames;
  double ave_trend_stddev;
  int dotfafirst, use_bin, nbins, use_period, maxiter, pertype, lastindex, use_harm, Nharm, Nsubharm;
  double tfathresh;
  int decorrflag;
  int decorr_iterate;
  int decorr_Nlcterms, Ndecorr, Ntfatot;
  int *decorr_lc_order, *decorr_lc_columns;
  double ***lcdecorr_terms_in;
  int jdcol_isfromheader, magcol_isfromheader;
  char jdcol_headername[MAXLEN], magcol_headername[MAXLEN];
  int clippingusemedian;
  int clippinguseMAD;
  int usefitmask;
  char *fitmaskvarname;
  _Variable *fitmaskvar;
  int outputfitmask;
  char *outputfitmaskvarname;
  _Variable *outputfitmaskvar;
} _TFA_SR;

typedef struct {
  double **initial_X, **final_X, **initial_colors, **initial_colors_readin, **final_colors, sigma_clip1, sigma_clip2, saturation, *JD, *mag_ave, *rms_out;
  char dates_name[MAXLEN], model_outdir[MAXLEN], trends_outname[MAXLEN], model_suffix[MAXLEN];
  int Nsysrem_color, Nsysrem_airmass, correctlc, omodel, otrend, subx, suby, sizex, sizey, Njd, Nsysrem_total, useweights;
  char **stringid;
  int *stringid_idx;
} _Sysrem;

typedef struct {
  char outdir[MAXLEN];
  int useformat;
  char format[MAXLEN];
  int usecolumnformat;
  char columnformat[MAXLEN];
  int Nvar;
  int outfits;
  /* Set by the "gzip" / "bzip2" keywords on -o: pipe the output
     light-curve through that compressor.  outgzip implies a trailing
     .gz on the output filename; outbzip2 implies .bz2.  Mutually
     exclusive (parser keeps the last one given). */
  int outgzip;
  int outbzip2;
  int copyheaderfrominput;
  int logcommandline;
  int noclobber;
  _Variable **variables;
  char **printfformats;
  char **varnames;
  char sepchar;
  int useoutnamecommand;
  char *outnamecommand;
  char **descriptions;
  char **units;
  int namefrominlist;
  char **inputlistoutnames;
  /* "changesuffix" keyword: strip a trailing changesuffix_remove from
     the constructed default output basename (if present) and then
     append changesuffix_add.  Applied after the basename is built and
     before any fits/gzip/bzip2 suffixes.  Mutually exclusive with the
     other slot-1 keywords (nameformat / namecommand / namefromlist). */
  int usechangesuffix;
  char changesuffix_remove[MAXLEN];
  char changesuffix_add[MAXLEN];
  /* When set, writelightcurves emits every currently-registered
     VARTOOLS_VECTORTYPE_LC variable using its name as the column name and
     a default printf format chosen per datatype.  Populated in
     CompileAllExpressions once the variable registry is complete. */
  int allcols;
  /* Snapshot of p->NDefinedVariables at the moment the "-o ... allcols"
     argument was parsed, so the allcols expansion only emits variables
     that were defined before this -o command (not ones created by later
     commands on the same command line). */
  int allcols_nvars_snapshot;
  /* When set, treat the "-o" argument as an output *directory* (with the
     per-LC basename derived from p->lcnames) even when running in single-
     file (-i) mode.  Normally the directory-vs-file distinction follows
     listflag (set by -l) vs fileflag (set by -i).  This override exists
     so that the in-process libvartoolspipeline driver, which uses
     fileflag mode, can still produce per-LC output files when called
     with a list of light curves -- and is also available to any -o
     CLI invocation as the "forceoutdirmode" keyword. */
  int force_outdir_mode;
  /* When set, this -o command triggers a memcpy of the current LC
     variables into the matching p->captured[] slot keyed by
     ``capture_id``.  Activated by either the "capture" keyword
     (``-o <id> capture``, no file written) or the "capture_id <id>"
     keyword (``-o <path> ... capture_id <id>``, file written *and*
     captured).  Used by the libvartoolspipeline driver to satisfy
     cmd.o(capture=True). */
  int capture_to_buffer;
  /* When set together with capture_to_buffer, DoOutputLightCurve
     snapshots and returns -- no fopen/fwrite.  Cleared (capture_to_
     buffer still set) means "do the regular file write *and* the
     snapshot".  Determined at parse time: "capture" keyword sets
     skip=1 (the -o argument was an id, not a path); "capture_id"
     keyword sets skip=0 (the -o argument is a real path). */
  int capture_skip_file;
  /* The slot key used to address p->captured[].  Always populated
     when capture_to_buffer is set.  For the "capture" form (skip=1)
     this is a copy of outdir.  For the "capture_id" form (skip=0)
     this is the explicit id argument given after the keyword. */
  char capture_id[MAXLEN];
} _Outputlcs;

typedef struct {
  double *rmsval;
  double *ave;
  int *ngood;
  int usemask;
  _Variable *maskvar;
} _Changeerror;

#define VARTOOLS_CHANGEVAR_TIME 0
#define VARTOOLS_CHANGEVAR_MAG 1
#define VARTOOLS_CHANGEVAR_ERR 2
#define VARTOOLS_CHANGEVAR_ID 3

typedef struct {
  char changevar;
  char newvarname[MAXLEN];
  _Variable *newvar;
} _Changevariable;

typedef struct {
  int savenumber;
  int saveindex;
  int ispartialrestore;
  int Nrestorevars;
  char **restorevarnames;
  _Variable **restorevars;
} _Restorelc;

typedef struct {
  int *runyet;
  int *sizevecs;
  int *sizesvecs;
  int *sizestringid_idxvecs;
  int *NJD;
  int *Ndblterms;
  int *Nsterms;
  int *lclistindx;
  double ***dblterms;
  int **stringid_idx;
  int *Nshterms;
  int *Nlterms;
  int *Niterms;
  int *Nfterms;
  int *Ncterms;
  int maxstring;
  char ***sterms;
  int ***iterms;
  long ***lterms;
  short ***shterms;
  float ***fterms;
  char ***cterms;

} _Savelc;

typedef struct {
  int usemean;
  int replace;
  double time;
  VT_PARAM_COMPANIONS(time);
} _MedianFilter;

typedef struct {
  int sepstarlist;
  char starlistname[MAXLEN];
  double *varx;
  double *vary;
  double **varxyin;
  char **varnames;
  char **varblendnames;
  double *blendamps;
  double matchrad;
  double zeromag;
  int converttoflux;
  double **periods;
  double fixperiod;
  int Nharm;
  int outputmatches;
  int radec;
  int pertype;
  char outmatchesfilename[MAXLEN];
  OutColumn *linkedcolumn;
  OutColumn *linkedcolumn_varname;
} _FindBlends;

typedef struct {
  int fitf0, fitf1, fitu0, fitt0, fittmax;
  int f0_source, f1_source, u0_source, t0_source, tmax_source;
  double **f00, **f10, **u00, **t00, **tmax0;
  double f00_fix, f10_fix, u00_fix, t00_fix, tmax0_fix;
  _Variable *f00_var, *f10_var, *u00_var, *t00_var, *tmax0_var;
  _Expression *f00_expr, *f10_expr, *u00_expr, *t00_expr, *tmax0_expr;
  OutColumn *f0_linkedcolumn, *f1_linkedcolumn, *u0_linkedcolumn, *t0_linkedcolumn, *tmax_linkedcolumn;
  double *f0, *f1, *u0, *t0, *tmax, *chi2_;
  int f0_initialstep, f1_initialstep, u0_initialstep, t0_initialstep, tmax_initialstep;
  double f0_initialstepval, f1_initialstepval, u0_initialstepval, t0_initialstepval, tmax_initialstepval;
  int correctlc, omodel;
  char modeloutdir[MAXLEN];
  char modelsuffix[10];
} _MicroLens;

typedef struct {
  int noise_type;

  double gammaval_fix;
  _Variable *gammaval_var;
  _Expression *gammaval_expr;
  double **gammaval;
  int gammaval_type;
  double sig_r_fix;
  _Variable *sig_r_var;
  _Expression *sig_r_expr;
  double **sig_r;
  int sig_r_type;
  double sig_w_fix;
  _Variable *sig_w_var;
  _Expression *sig_w_expr;
  double **sig_w;
  int sig_w_type;

  double **rho_r;
  int rho_r_type;
  double rho_r_fix;
  _Variable *rho_r_var;
  _Expression *rho_r_expr;

  double **nu_r;
  int nu_r_type;
  double nu_r_fix;
  _Variable *nu_r_var;
  _Expression *nu_r_expr;

  int bintime_type;
  double bintime_fix;
  _Variable *bintime_var;
  _Expression *bintime_expr;
  double **bintime;
} _AddNoise;

typedef struct {
  int inputtimetype;
  int inputsys;
  double inputsubtractval;
  int outputtimetype;
  int outputsys;
  double outputsubtractval;
  int useradec;
  int radec_source;
  _Expression *raval_expr;
  _Expression *decval_expr;
  _Expression *raval_lcexpr;
  _Expression *decval_lcexpr;
  double **ravals;
  double **decvals;
  int raval_lc_col;
  int decval_lc_col;
  double **raval_lcvals;
  double **decval_lcvals;
  double raval_fix;
  double decval_fix;
  double radecepoch;
  int useppm;
  int ppm_source;
  double **ppm_mu_ra_vals;
  double **ppm_mu_dec_vals;
  double ppm_mu_ra_fix;
  double ppm_mu_dec_fix;
  int useinput_radec;
  int inputradec_source;
  _Expression *inputraval_expr;
  _Expression *inputdecval_expr;
  _Expression *inputraval_lcexpr;
  _Expression *inputdecval_lcexpr;
  double **inputravals;
  double **inputdecvals;
  int inputraval_lc_col;
  int inputdecval_lc_col;
  double **inputraval_lcvals;
  double **inputdecval_lcvals;
  double inputraval_fix;
  double inputdecval_fix;
  double inputradecepoch;
  int useinputppm;
  int inputppm_source;
  double **inputppm_mu_ra_vals;
  double **inputppm_mu_dec_vals;
  double inputppm_mu_ra_fix;
  double inputppm_mu_dec_fix;
#ifdef _HAVE_CSPICE
  char ephemfile[MAXLEN];
  char leapsecfile[MAXLEN];
  char planetdatafile[MAXLEN];
  int source_obs_coords;
  double obslat_fixval;
  double obslong_fixval;
  double obsalt_fixval;
  double *obslat_listvals;
  double *obslong_listvals;
  double *obsalt_listvals;
  int obslat_lc_col;
  int obslong_lc_col;
  int obsalt_lc_col;
  double **obslat_lcvals;
  double **obslong_lcvals;
  double **obsalt_lcvals;
  _Expression *obslat_expr;
  _Expression *obslat_lcexpr;
  _Expression *obslong_expr;
  _Expression *obslong_lcexpr;
  _Expression *obsalt_expr;
  _Expression *obsalt_lcexpr;
  int obs_coords_usexyz;
#endif
} _ConvertTime;

typedef struct {
  _Expression *expression;
  _Variable *outputvar;
  char *lhsstring;
  char *rhsstring;
  double **tmpoutvals;
  int lhs_indx_type;
  char *lhsindexstring1;
  char *lhsindexstring2;
  _Expression *lhs_indx_expr1;
  _Expression *lhs_indx_expr2;
  char lhs_indx_range_startmin;
  char lhs_indx_range_stopmax;
  char initialize_output_var;
  int lhs_vectortype_override;
  /* outputcolumn: when set (1) by the optional "outputcolumn" CLI
     keyword, the value computed by this expression is exposed as an
     output column in the result table.  Only valid when
     lhs_vectortype_override is INLIST/SCALAR/CONSTANT (i.e. when the
     user gave "listvar", "scalar", or "const") — for the default
     per-observation LC vectortype an output column would be
     ill-defined, so it is rejected at parse time.  outputcolumn_buffer
     is the per-LC scratch buffer that backs the column; it is
     allocated in initcommands.c and written by RunExpressionCommand. */
  int outputcolumn;
  double *outputcolumn_buffer;
} _ExpressionCommand;

typedef struct {
  char *functionstring;
  char *paramliststring;
  char **paramnames;
  char *modelvarname;
  _Variable *modelvar;
  char *outfile_extension;
  char *outfilename_format;
  int correctlc;
  int omodel;
  char *outdir;
  int Nparams;
  int calcchi2out;
  double *chi2out;
  double **param_outvals;
  double **param_uncertainties;
  int *numrej;
  int *iternum;
  _Variable **params;
  _Expression **expressions;
  _Expression *constantexpression;
  int rejectoutliers;
  double rejsigclip;
  VT_PARAM_COMPANIONS(rejsigclip);
  int rejuseMAD;
  int rejiterate;
  int rejfixnum;
  int rejiternum;
  int usemask;
  char *maskvarname;
  _Variable *maskvar;
} _Linfit;

typedef struct {
  char *functionstring;
  char *paramliststring;
  char **paramnames;
  char **paraminitstrings;
  char **paramerrstrings;
  char *priorliststring;
  char **priorvarnames;
  char **priorstrings;
  char *constraintliststring;
  char **constraintstrings;
  int fittype;
  char *modelvarname;
  _Variable *modelvar;
  int correctlc;
  int omodel;
  char *outdir;
  char *outfile_extension;
  char *outfilename_format;
  int Nparams;
  int Npriors;
  int Nconstraints;
  char *errorstring;
  _Expression *errorexpression;

  int use_covar;
  int covar_type;
  char *Corr_rho_varname;
  char *Corr_rho_exprstring;
  char *Corr_amp_varname;
  char *Corr_amp_exprstring;
  char *Corr_nu_varname;
  char *Corr_nu_exprstring;
  _Expression *Corr_rho_expr;
  _Expression *Corr_amp_expr;
  _Expression *Corr_nu_expr;
  _Variable *Corr_rho_var;
  _Variable *Corr_amp_var;
  _Variable *Corr_nu_var;
  double ***Corr_mat1;
  double ***Corr_mat2;
  int **Corr_Nvec;
  int *Corr_sizemat;
  int *Corr_store_NJD;

  int uselinfit;
  _Linfit *linfit;

  double **param_outvals;
  double **param_uncertainties;
  double *chi2out;
  int *amoeba_isconverged;
  _Variable **params;
  _Expression *functionexpression;
  _Expression **paraminit_expressions;
  _Expression **paramerr_expressions;
  _Expression **prior_expressions;
  _Expression **constraint_expressions;

  double amoeba_tol;
  VT_PARAM_COMPANIONS(amoeba_tol);
  long amoeba_maxsteps;
  VT_PARAM_COMPANIONS(amoeba_maxsteps);

  long mcmc_Naccept;
  VT_PARAM_COMPANIONS(mcmc_Naccept);
  long mcmc_Nlinkstotal;
  VT_PARAM_COMPANIONS(mcmc_Nlinkstotal);
  double mcmc_burninfrac;
  VT_PARAM_COMPANIONS(mcmc_burninfrac);
  double mcmc_eps;
  VT_PARAM_COMPANIONS(mcmc_eps);
  char *mcmc_chain_exprliststring;
  char **mcmc_chain_expr_strings;
  char *mcmc_chain_statsliststring;
  double mcmc_max_mem_store;
  int mcmc_outchains;
  char *mcmc_outchains_dir;
  char *mcmc_outchains_format;
  int mcmc_outchains_print_every;
  int mcmc_skipamoeba;

  int N_mcmc_chain_expressions;
  int N_mcmc_chain_stats;
  int N_mcmc_chain_statstot;

  int *mcmc_statstocalc;
  double *pctval;
  double **mcmc_statsout;

  _Expression **mcmc_chain_stats_expressions;
  int usemask;
  char *maskvarname;
  _Variable *maskvar;
} _Nonlinfit;

#include "userlib.h"
#include "userfunc.h"

#ifdef DYNAMICLIB

typedef struct {
  int datatype;
  int source;
  int Nfixptr;
  size_t size_element;
  int Ncolumns;
  void *dataptr;
  void *inlistdataptr;
  void *inlcdataptr;
  char *outname;
  char **priorname;
  OutColumn **linkedcolumn;
  _Expression *evalexpression;
  void (*initialize_usertype_ptr)(int, void *, void *);
  void *extra_user_data;
  _Variable **existingvariable;
  int expectedvectortype;
} _UserDataPointer;

typedef struct {
  /*_UserLib *lib;*/
  int libnum;
  void *userdata;

  int Nfix;
  int Ninlist;
  int Ninlc;
  int Ncomputed;
  int Nprior;
  int Nptrs;
  int Noutput;
  int Nexpr;

  _UserDataPointer *UserDataPointers;
  _UserDataPointer *FixValues;
  _UserDataPointer *OutputData;

  char **expr_strings;
  _Expression ***UserDataExpressions;

} _UserCommand;

#endif

typedef struct {
  int Nstats;
  int *statstocalc;
  int Nvar;
  char **varnames;
  _Variable **vars;
  int Nstatstot;
  double **statsout;
  double *pctval;
  int usemask;
  char *maskvarname;
  _Variable *maskvar;
} _Stats;

typedef struct {
  int isforward;
  char *inputvarname_real;
  _Variable *inputvar_real;
  char *inputvarname_imag;
  _Variable *inputvar_imag;
  char *outputvarname_real;
  _Variable *outputvar_real;
  char *outputvarname_imag;
  _Variable *outputvar_imag;
} _FFT;

#define VARTOOLS_IFTYPE_IF 0
#define VARTOOLS_IFTYPE_ELIF 1
#define VARTOOLS_IFTYPE_ELSE 2
#define VARTOOLS_IFTYPE_FI 3

typedef struct {
  short iftype;
  int ifindex;
  _IfStruct *ifs;
  char *exprstring;
} _IfCommand;

#define VARTOOLS_RESTRICTTIMES_JDRANGE 0
#define VARTOOLS_RESTRICTTIMES_JDLIST 1
#define VARTOOLS_RESTRICTTIMES_IMAGELIST 2
#define VARTOOLS_RESTRICTTIMES_EXPR 3

typedef struct {
  char restricttype;
  char minJDtype;
  char maxJDtype;
  char exclude; 
  double minJDfixval;
  double maxJDfixval;
  double *minJD;
  double *maxJD;
  OutColumn *minJD_linkedcolumn;
  OutColumn *maxJD_linkedcolumn;
  int N_restrictlist;
  double *JD_restrictlist;
  char **image_restrictlist;
  int *image_restrictlist_indx;
  _Expression *minJDexpr;
  _Expression *maxJDexpr;
  char *minJDexprstring;
  char *maxJDexprstring;
  char *restrictexprstring;
  _Expression *restrictexpr;
  _Savelc *s;
  int saveexcludedpoints;
  int markrestrict;
  char *markvarname;
  _Variable *markvar;
  int noinitmark;
} _RestrictTimes;

typedef struct {
  int restrictnum;
  _RestrictTimes *RestrictTimes;
} _RestoreTimes;

typedef struct {
  double cterm;
  VT_PARAM_COMPANIONS(cterm);
  double maxfreq;
  VT_PARAM_COMPANIONS(maxfreq);
  double freq_sample_factor;
  VT_PARAM_COMPANIONS(freq_sample_factor);
  int auto_tau0;
  int auto_tau1;
  int auto_dtau;
  double tau0;
  VT_PARAM_COMPANIONS(tau0);
  double tau1;
  VT_PARAM_COMPANIONS(tau1);
  double dtau;
  VT_PARAM_COMPANIONS(dtau);
  int outfulltransform;
  int outfulltransform_usefits;
  int outfulltransform_usepm3d;
  char *outfulltransform_dir;
  char *outfulltransform_format;
  int outmaxtransform;
  char *outmaxtransform_dir;
  char *outmaxtransform_format;

  double *max_z;
  double *max_freq;
  double *max_pow;
  double *max_amp;
  double *max_neff;
  double *max_tau;
  double *max_con;
  double *med_z;
  double *med_freq;
  double *med_pow;
  double *med_amp;
  double *med_neff;
  double *med_con;
  int usemask;
  _Variable *maskvar;
} _WWZ;

typedef struct {
  int *lclistindx;
  int *runyet;
  int *Ndblterms;
  int *Nsterms;
  double **dblterms;
  int *Nshterms;
  int *Nlterms;
  int *Niterms;
  int *Nfterms;
  int *Ncterms;
  char **sterms;
  int **iterms;
  long **lterms;
  short **shterms;
  float **fterms;
  char **cterms;
} _SaveListData;

typedef struct {
  int cnum;
  int Ncopies;
  int copycommand_index;
  int priorcopies;
  int *lcid_tothreadid;
  _Savelc *s;
  _SaveListData *SaveListData;
  _IfStack **IfStack;
  int *sizearray_IfStruct_wasfoundtrue_copy;
  char **IfStruct_wasfoundtrue_copy;
} _CopyLC;

typedef struct {
  int resample_method;

  double *tstart;
  double tstart_fix;
  int tstart_source;
  OutColumn *tstart_linkedcolumn;
  _Expression *tstart_expr;
  char *tstart_exprstring;

  double *tstop;
  double tstop_fix;
  int tstop_source;
  OutColumn *tstop_linkedcolumn;
  _Expression *tstop_expr;
  char *tstop_exprstring;

  int *Nresamp;
  int Nresamp_fix;
  int Nresamp_source;
  OutColumn *Nresamp_linkedcolumn;
  _Expression *Nresamp_expr;
  char *Nresamp_exprstring;

  double *delt;
  double delt_fix;
  int delt_source;
  OutColumn *delt_linkedcolumn;
  _Expression *delt_expr;
  char *delt_exprstring;

  char **resample_filenames;
  char *resample_filename_fix;
  int t_column;
  double *t_resamp;
  int N_resamp;
  int resample_filename_source;
  int use_file;
  
  /* Various interpolation parameters */
  double yp1;
  double ypn;
  int bspline_nbreaks;
  int bspline_order;

  /* If a different method is specified for points that are far from an
     observation */
  int use_near_far;
  int resample_method_far;
  double yp1_far;
  double ypn_far;
  int bspline_nbreaks_far;
  int bspline_order_far;
  int minsep_source;
  double *minsep;
  double minsep_fix;
  OutColumn *minsep_linkedcolumn;
  _Expression *minsep_expr;
  char *minsep_exprstring;
  double frac_min_sep_val;
  double frac_med_sep_val;
  double percentile_sep;
  
  /* If a different extrapolation method is specified */
  int use_extrap;
  int resample_method_extrap;
  double yp1_extrap;
  double ypn_extrap;
  int bspline_nbreaks_extrap;
  int bspline_order_extrap;
} _Resample;

typedef struct {
  int filtertype;

  double *maxfreq;
  double maxfreq_fix;
  int maxfreq_source;
  OutColumn *maxfreq_linkedcolumn;
  _Expression *maxfreq_expr;
  char *maxfreq_exprstring;

  double *minfreq;
  double minfreq_fix;
  int minfreq_source;
  OutColumn *minfreq_linkedcolumn;
  _Expression *minfreq_expr;
  char *minfreq_exprstring;
  
  _Expression *filter_expr;
  char *filter_exprstring;
  /* User-visible name of the frequency variable referenced in
     filter_exprstring (default "f").  Substituted at CompileAllExpressions
     time with a unique internal stump name so multiple -fourierfilter
     commands in the same pipeline do not collide and so there is no
     conflict with same-named user variables. */
  char freq_varname[64];

  int calc_full_spec;

  int forcefft;

  /* Optional band-edge taper.  taper_type is one of VARTOOLS_FOURIERFILTER_TAPER_*
     from the enum in fourierfilter.c; NONE means no taper is applied.  The
     taper is centered on each cut edge of the selected filter mode and spans
     +/- taper_deltafreq around that edge.  taper_beta is only used by the
     Kaiser window. */
  int    taper_type;
  double taper_deltafreq;
  double taper_beta;

  /* Optional resample-to-uniform path.  When resample_enabled is set the
     LC is interpolated onto a uniform grid at spacing resample_delta,
     FFT-filtered, IFFT-reconstructed, and the filter result interpolated
     back to the original sample times.  resample_source selects how
     resample_delta is obtained per LC:
       VARTOOLS_SOURCE_FIXED        — resample_delta_fix (scalar on CLI)
       VARTOOLS_SOURCE_INLIST       — resample_delta[lcid] from the list
       VARTOOLS_SOURCE_PRIORCOLUMN  — from an earlier command's output col
       VARTOOLS_SOURCE_EVALEXPRESSION — resample_delta_expr
       -1 (DELMIN)                  — use the minimum dt in this LC
  */
  int    resample_enabled;
  int    resample_source;
  double resample_delta_fix;
  double *resample_delta;
  OutColumn *resample_delta_linkedcolumn;
  _Expression *resample_delta_expr;
  char  *resample_delta_exprstring;

  /* Optional gap-break.  When gapbreak_enabled is set the LC is split at
     any gap >= gap_threshold and each segment is filtered independently;
     segments are concatenated on output.  gapbreak_source selects how
     gap_threshold is obtained:
       VARTOOLS_SOURCE_FIXED        — gapbreak_fix
       VARTOOLS_SOURCE_PRIORCOLUMN  — from an earlier command's output col
       VARTOOLS_SOURCE_INLIST       — gapbreak_threshold[lcid]
       VARTOOLS_SOURCE_EVALEXPRESSION — gapbreak_expr
       -1 (FRAC_MIN_SEP)  — gapbreak_frac_min * min(dt)
       -2 (FRAC_MED_SEP)  — gapbreak_frac_med * median(dt)
       -3 (PERCENTILE_SEP) — percentile(dt, gapbreak_percentile)
  */
  int    gapbreak_enabled;
  int    gapbreak_source;
  double gapbreak_fix;
  double gapbreak_frac_min;
  double gapbreak_frac_med;
  double gapbreak_percentile;
  double *gapbreak_threshold;
  OutColumn *gapbreak_linkedcolumn;
  _Expression *gapbreak_expr;
  char  *gapbreak_exprstring;

  int ofourier;
  char *ofourier_dir;
  int ofourier_formatflag;
  char *ofourier_format;

  /* Optional "nowarn" flag — when set, suppresses all runtime
     per-LC warnings emitted by doFourierFilter (non-uniform without
     resample, within-segment gap vs minfreq, taper overlap, forcefft
     on non-uniform, resample delta <= 0).  Useful in batch pipelines
     where the caller has vetted the data and doesn't need the repeated
     advisories. */
  int nowarn;

  /* Optional edge-padding for the resample path's FFT.  The FFT
     implicitly treats the signal as periodic, so if the first and
     last samples disagree (typical astronomical LCs), the implicit
     wrap-around produces ringing in the filtered output near the
     start/end of the segment.  padmode selects how each segment is
     extended before the FFT and padfrac is the padding length per
     side, as a fraction of the segment length.  The padding is
     discarded before interpolating back to the original times.
       VARTOOLS_FOURIERFILTER_PAD_WRAP     — no padding (default)
       VARTOOLS_FOURIERFILTER_PAD_REFLECT  — mirror at each edge
       VARTOOLS_FOURIERFILTER_PAD_ZERO     — zero-extend each edge
   */
  int    padmode;
  double padfrac;

  _Variable *freq_var;

  /* Per-LC output-column storage.  Allocated in initcommands.c once Nlcs
     is known; written by doFourierFilter at runtime; registered as
     FourierFilter_<name>_N output columns in outcolumns.c. */
  double *mean_mag;          /* DC Fourier coefficient (a_0)                  */
  double *rms_in;            /* RMS of the input light curve                  */
  double *rms_out;           /* RMS after the filter is applied               */
  int    *nfreqcalc;         /* number of Fourier frequency bins computed     */
  int    *nfreqfilt;         /* number of bins inside the filter passband     */

} _FourierFilter;

#ifdef _HAVE_PYTHON
typedef struct {
  int Nvars;
  _Variable **vars;
  int *isvaroutput;
  int Nvars_outonly;
  _Variable **outonlyvars;
  int Nlcvars_nonupdate;
  int **outlcvecs_invars;
  int **outlcvecs_outonlyvars;
  _Variable **lcvars_nonupdate;
  int *IsPythonRunning;
  int **sockets;
  char *progname;
  char *pythoninitializationtext;
  long len_pythoninitializationtextstring;
  char *pythoncommandstring;
  long len_pythoncommandstring;
  char *inputpythonfilename;
  void *pythonobjects;
  int iscontinueprocess;
  void *continueprocesscommandptr;

  int Nchildren;
  void *childcommandptrs;
  int *childcnumvals;

  int cid;

  char *inoutvarliststring;
  char **inoutvarnames;
  int Ninoutvarnames;

  char *invarliststring;
  char *outvarliststring;

  char **invarnames;
  char **outvarnames;

  int Ninvarnames;
  int Noutvarnames;

  char *outcolumnliststring;
  char **outcolumnnames;
  int Noutcolumnvars;
  _Variable **outcolumnvars;

  double **outcolumndata;

  int processallvariables;

  int cnum;

  int RequireReadAll;

  void *FullList;

  int skipfail;

} _PythonCommand;
#else
typedef struct {
  int cnum;
} _PythonCommand;
#endif

#ifdef _HAVE_R
typedef struct {
  int Nvars;
  _Variable **vars;
  int *isvaroutput;
  int Nvars_outonly;
  _Variable **outonlyvars;
  int Nlcvars_nonupdate;
  int **outlcvecs_invars;
  int **outlcvecs_outonlyvars;
  _Variable **lcvars_nonupdate;
  int *IsRRunning;
  int **sockets;
  char *progname;
  char *Rinitializationtext;
  long len_Rinitializationtextstring;
  char *Rcommandstring;
  long len_Rcommandstring;
  char *inputRfilename;
  void *Robjects;
  int iscontinueprocess;
  void *continueprocesscommandptr;

  int Nchildren;
  void *childcommandptrs;
  int *childcnumvals;

  int cid;

  char *inoutvarliststring;
  char **inoutvarnames;
  int Ninoutvarnames;

  char *invarliststring;
  char *outvarliststring;

  char **invarnames;
  char **outvarnames;

  int Ninvarnames;
  int Noutvarnames;

  char *outcolumnliststring;
  char **outcolumnnames;
  int Noutcolumnvars;
  _Variable **outcolumnvars;

  double **outcolumndata;

  int processallvariables;

  int cnum;

  int RequireReadAll;

  void *FullList;

  int verboseR;

} _RCommand;
#else
typedef struct {
  int cnum;
} _RCommand;
#endif

typedef struct {
  int sizevec;
  int Npoints;
  int datatype;
  int incol;
  char *format;
  void *dataptr;
} _MatchData;

typedef struct {
  int cnum;
  char *inputfilename;
  int inlistcolumn;
  int Ninlist;
  char **inputfilenamelist;
  char *opencommand;
  int Nskip;
  int Nskipchar;
  char *skipchars;
  int delimtype;
  char delimchar;
  char *delimstring;
  int matchcolumn;
  char *matchcolumnvarname;
  char *matchcolumn_header_name;
  int getmatchcolumnfromheader;
  _Variable *matchvar;
  int Naddvars;
  _Variable **addvars;
  char **addvar_varnames;
  int *addvar_columns;
  int *addvar_datatypes;
  char **addvar_formats;
  char **addvar_incolumn_header_names;
  int getcolumnsfromheader;
  int missingmethod;
  double missingvalue;
  _MatchData *matchdata;
  int *matchdata_sortindx;
} _MatchCommand;

typedef struct {
  char *keyname;
  int dtype;
  int keyval_source;
  double dbl_fixval;
  int int_fixval;
  long long_fixval;
  char *string_fixval;
  _Variable *keyval_var;
  char *comment_string;
  int hdutouse;
  int updateexisting;
  int combinelckeyword;
  _Variable *lcnumvar;
} _AddFitsKeyword;

typedef struct {
  int issortvar;
  char *sortvarname;
  _Variable *sortvar;
  int isreverse;
  int sortdtype;
} _SortLC;

typedef struct {
  int Nvars;
  _Variable **vars;
  char **varnames;
  int iscolname;
  char **colnames;
  int isformat;
  char **formatstrings;
  void **dataptr;
  int *colindx;
} _PrintCommand;


/* -PDM variant identifiers (shared between parser and pdm.c) */
#define PDM_KIND_STEP        0
#define PDM_KIND_LINTERP     1
#define PDM_KIND_MULTICOVER  2
#define PDM_KIND_TOPHAT      3
#define PDM_KIND_GAUSS       4

typedef struct {
  int kind;                /* PDM_KIND_STEP, PDM_KIND_LINTERP, ... */
  /* Nbin: phase-bin count (binned variants).  May come from a fixed integer,
   * an existing variable, or an analytic expression -- resolved per-LC. */
  int Nbin;
  int *Nbin_vals;
  int Nbin_source;
  _Variable *Nbin_var;
  _Expression *Nbin_expr;
  /* minp / maxp / subsample / finetune: same var/expr/fixed source machinery. */
  double minp;
  double *minp_vals;
  int minp_source;
  _Variable *minp_var;
  _Expression *minp_expr;
  double maxp;
  double *maxp_vals;
  int maxp_source;
  _Variable *maxp_var;
  _Expression *maxp_expr;
  double subsample;
  double *subsample_vals;
  int subsample_source;
  _Variable *subsample_var;
  _Expression *subsample_expr;
  double finetune;
  double *finetune_vals;
  int finetune_source;
  _Variable *finetune_var;
  _Expression *finetune_expr;
  /* Nc: number of phase-shifted bin sets ("covers") for the multicover variant.
   * Nc = 1 for step/linterp; defaults to 2 for multicover. */
  int Nc;
  int *Nc_vals;
  int Nc_source;
  _Variable *Nc_var;
  _Expression *Nc_expr;
  /* dphi: phase-window half-width (tophat) or Gaussian kernel sigma (gauss).
   * Only used by the binless variants; defaults to 0.05 (cuvarbase). */
  double dphi;
  double *dphi_vals;
  int dphi_source;
  _Variable *dphi_var;
  _Expression *dphi_expr;
  /* simple scalars */
  int Npeaks;
  int operiodogram;        /* 0/1: dump periodogram file per LC */
  char outdir[MAXLEN];
  char suffix[8];          /* file suffix, default ".pdm" */
  int useerr;              /* 1: weight by 1/sig^2 (default); 0: noerr keyword */
  double clip;             /* sigma-clip factor for the SNR noise estimate (default 5) */
  int clipiter;            /* 1: iterate clipping until count stable; 0: single pass */
  /* Per-LC outputs */
  double **peakperiods;    /* [Nlcs][Npeaks] */
  double **peakvalues;     /* [Nlcs][Npeaks]  -- theta at each peak */
  double **peakSNR;        /* [Nlcs][Npeaks] */
  double **peakFAP;        /* [Nlcs][Npeaks] */
  double *avetheta;        /* [Nlcs] */
  double *rmstheta;        /* [Nlcs] */
  /* fixperiodSNR: optionally compute theta/SNR/FAP at a specified period
   * (taken from a prior -aov / -ls / -pdm / -Injectharm, or a literal,
   * or a list column, or a fixcolumn back-reference).  Mirrors -aov. */
  int fixperiodSNR;
  int fixperiodSNR_pertype;
  int fixperiodSNR_lastaovindex;
  double fixperiodSNR_fixedperiod;
  double **fixperiodSNR_periods;     /* [Nlcs][1] */
  double *fixperiodSNR_peakvalues;   /* [Nlcs] -- theta at fixed period */
  double *fixperiodSNR_peakSNR;      /* [Nlcs] */
  double *fixperiodSNR_peakFAP;      /* [Nlcs] */
  OutColumn *fixperiodSNR_linkedcolumn;
  /* maskpoints: optional LC vector; points with maskvar > VARTOOLS_MASK_TINY
   * are included, others excluded.  Mirrors -aov's keyword. */
  int usemask;
  _Variable *maskvar;
  /* whiten: iterative pre-whitening between peaks (model = step bin means at
   * the peak period; subtracted from the LC before recomputing the
   * periodogram for the next peak).  Mirrors -aov's keyword. */
  int whiten;
  double **avetheta_whiten;   /* [Nlcs][Npeaks] -- per-cycle periodogram mean */
  double **rmstheta_whiten;   /* [Nlcs][Npeaks] -- per-cycle periodogram rms */
  /* bootstrap: empirical FAP calibration via Nboot shuffled-LC trials.
   * 0 = off; >0 = enabled with that many trials.  When on, the analytic
   * SCz Beta FAP is replaced by an empirical CDF + log-log polynomial
   * tail extrapolation (mirrors -LS's bounded-statistic branch). */
  int bootstrap_Nboot;
} _PDM;


/* ---------------- Fast Template Periodogram (FTP, Hoffman+ 2021) ----------- */

typedef struct {
  /* Template Fourier-series coefficients of length H.
   *   M(phi) = sum_{n=1..H} [c_n cos(n phi) + s_n sin(n phi)]
   *
   * For the "file" and "fitlc" template-source modes, cn/sn are loaded
   * once at parser time and shared across all LCs (inline_mode = 0).
   *
   * For the "inline" mode (inline_mode = 1), cn/sn are per-LC scratch
   * buffers populated at run time in RunFTPCommand from the per-slot
   * (cn_source[k], cn_var[k], cn_expr[k], cn_lit[k]) tuples (same for sn).
   * H = Nharm + 1, following the -harmonicfilter convention. */
  int     H;
  double *cn;             /* size H -- shared (file/fitlc) or per-LC scratch (inline) */
  double *sn;             /* size H */
  char   *template_path;  /* original CLI tag, kept for diagnostics */

  int     inline_mode;    /* 1 when cn/sn are resolved per-LC from var/expr/lit */
  int    *cn_source;      /* size H, VARTOOLS_SOURCE_FIXED/EXISTINGVARIABLE/EVALEXPRESSION */
  int    *sn_source;      /* size H */
  _Variable   **cn_var;   /* size H, non-NULL where source==EXISTINGVARIABLE */
  _Variable   **sn_var;
  _Expression **cn_expr;  /* size H, non-NULL where source==EVALEXPRESSION */
  _Expression **sn_expr;
  double *cn_lit;         /* size H, used where source==FIXED */
  double *sn_lit;

  int     filelist_mode;       /* 1 when each LC has its own template file path */
  char  **template_filenames;  /* [Nlcs], populated by RegisterDataFromInputList */

  int     method;              /* 0=brute (default), 1=poly, 2=verify(both + cmp) */

  /* Negative-amplitude policy: 1 to allow theta_1 < 0 in the best fit
   * (default; flagged via FTP_NegAmp_N_M output column); 0 to reject. */
  int allow_neg_amp;
  int useerr;             /* 1 default; 0 if 'noerr' keyword */

  /* Frequency grid parameters (each has the standard
   *   <"var" name | "expr" expr | fixed_value>
   * source-resolution machinery).  Period units match the LC. */
  double minp, maxp, subsample, finetune;
  int    minp_source, maxp_source, subsample_source, finetune_source;
  _Variable    *minp_var, *maxp_var, *subsample_var, *finetune_var;
  _Expression  *minp_expr, *maxp_expr, *subsample_expr, *finetune_expr;
  double *minp_vals, *maxp_vals, *subsample_vals, *finetune_vals;  /* per-LC */

  int Npeaks;

  /* operiodogram + outdir (mirrors -aov / -PDM). */
  int   operiodogram;
  char  outdir[MAXLEN];
  char  suffix[MAXLEN];

  /* Per-LC outputs */
  double **peakperiods;   /* [Nlcs][Npeaks] */
  double **peakvalues;    /* [Nlcs][Npeaks]  -- P(omega) at each peak in [0,1] */
  double **peakSNR;       /* [Nlcs][Npeaks] */
  int    **peakNegAmp;    /* [Nlcs][Npeaks]  -- 0 / 1 */
  double **peakTheta;     /* [Nlcs][Npeaks]  -- best-fit theta_2 (radians) */
  double *avepower;       /* [Nlcs] -- clipped mean of the periodogram */
  double *rmspower;       /* [Nlcs] -- clipped RMS  of the periodogram */
} _FTP;


typedef struct {
  int cnum;
  int require_sort;
  int require_distinct;
  char *command_outcolumn_suffix;
  _Clip *Clip;
  _Ensemblerescalesig *Ensemblerescalesig;
  _Rescalesig *Rescalesig;
  _Chi2_NoBin *Chi2_NoBin;
  _Chi2_Bin *Chi2_Bin;
  _RMS_NoBin *RMS_NoBin;
  _RMS_Bin *RMS_Bin;
  _Jstet *Jstet;
  _Alarm *Alarm;
  _Aov *Aov;
  _AovHarm *AovHarm;
  _Ls *Ls;
  _Decorr *Decorr;
  _Killharm *Killharm;
  _Injectharm *Injectharm;
  _Injecttransit *Injecttransit;
  _Outputlcs *Outputlcs;
  _Starspot *Starspot;
  _Bls *Bls;
  _BlsFixPer *BlsFixPer;
  _BlsFixDurTc *BlsFixDurTc;
  _BlsFixPerDurTc *BlsFixPerDurTc;
  _Phase *Phase;
  _Binlc *Binlc;
  _SoftenedTransit *SoftenedTransit;
  _GetLSAmpThresh *GetLSAmpThresh;
  _MandelAgolTransit *MandelAgolTransit;
  _DiffFluxtomag *DiffFluxtomag;
  _Fluxtomag *Fluxtomag;
  _TFA *TFA;
  _TFA_SR *TFA_SR;
  _Sysrem *Sysrem;
  _Changeerror *Changeerror;
  _Changevariable *Changevariable;
  _Autocorr *Autocorr;
  _Dftclean *Dftclean;
  _Restorelc *Restorelc;
  _Savelc *Savelc;
  _MedianFilter *MedianFilter;
  _FindBlends *FindBlends;
  _MicroLens *MicroLens;
  _ExpressionCommand *ExpressionCommand;
  _AddNoise *AddNoise;
#ifdef DYNAMICLIB
  _UserCommand *UserCommand;
#endif
  _ConvertTime *ConvertTime;
  _Linfit *Linfit;
  _Nonlinfit *Nonlinfit;
  _Stats *Stats;
  _IfCommand *IfCommand;
  _RestrictTimes *RestrictTimes;
  _RestoreTimes *RestoreTimes;
  _WWZ *WWZ;
  _CopyLC *CopyLC;
  _Resample *Resample;
  _FourierFilter *FourierFilter;
  _PythonCommand *PythonCommand;
  _RCommand *RCommand;
  _FFT *FFT;
  _MatchCommand *MatchCommand;
  _AddFitsKeyword *AddFitsKeyword;
  _SortLC *SortLC;
  _PrintCommand *PrintCommand;
  _PDM *Pdm;
  _FTP *Ftp;

  int N_setparam_expr;
  char **setparam_EvalExprStrings;
  _Expression ***setparam_EvalExpressions;

  int N_prior_vars;
  char *prior_var_datatypes;
  char *prior_var_vectortypes;
  char **prior_var_names;
  _Variable ***prior_vars;
  OutColumn **prior_var_linkedcolumns;
} Command;
