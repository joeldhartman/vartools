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
#define CNUM_HARMONICFILTER 52
#define CNUM_PYTHON 53
#define CNUM_RESTORETIMES 54
#define CNUM_FFT 55

#define TOT_CNUMS 55

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

#define VARTOOLS_HARMONICFILTER_FULLSPEC 0
#define VARTOOLS_HARMONICFILTER_HIGHPASS 1
#define VARTOOLS_HARMONICFILTER_LOWPASS 2
#define VARTOOLS_HARMONICFILTER_BANDPASS 3
#define VARTOOLS_HARMONICFILTER_BANDCUT 4

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
  int *Nclip;
  int iter;
  int niter;
  int usemedian;
} _Clip;

typedef struct {
  double *rescalefactor;
  double a;
  double b;
  double erssigclip;
  double *chi2_old;
  double *chi2_new;
} _Ensemblerescalesig;

typedef struct {
  double *rescalefactor;
  double *chi2_old;
  double *chi2_new;
} _Rescalesig;

typedef struct {
  double *chi2val;
  double *wtave;
} _Chi2_NoBin;

typedef struct {
  int Nbin;
  double *bintimes;
  double **chi2binvals;
  double **wtavebin;
} _Chi2_Bin;

typedef struct {
  double *rmsval;
  double *ave;
  double *rmsthy;
  int *ngood;
} _RMS_NoBin;

typedef struct {
  int Nbin;
  double *bintimes;
  double **rmsbinvals;
  double **rmsthybin;
} _RMS_Bin;

typedef struct {
  double Jstet_time;
  double *jst;
  double *kur;
  double *lst;
  double wkmax;
  char datesname[MAXLEN];
} _Jstet;

typedef struct {
  double *alarmvals;
} _Alarm;

typedef struct {
  double start;
  double stop;
  double step;
  double errsize;
  char outdir[MAXLEN];
  char suffix[10];
} _Autocorr;

typedef struct {
  double minp;
  double maxp;
  double subsample;
  double finetune;
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
} _Aov;

typedef struct {
  double minp;
  double maxp;
  double subsample;
  double finetune;
  int Npeaks;
  double *aveaov;
  double *rmsaov;
  double **peakperiods;
  double **peakvalues;
  double **peakSNR;
  double **peakFAP;
  int **peakNharm;
  int Nharm;
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
} _AovHarm;

typedef struct {
  double minp;
  double maxp;
  double subsample;
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
} _Decorr;

typedef struct {
  double gain;
  double SNlimit;
  int nbeam;
  double maxfreq;
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
  double clip_dirty, clip_clean;
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
} _Dftclean;

typedef struct {
  int pertype;
  int lastlsindex;
  int lastaovindex;
  double **periods;
  double *fixedperiods;
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
  char modelsuffix[16];
  double clip;
} _Killharm;

typedef struct {
  int pertype;
  int Nharm;
  int Nsubharm;
  double fixperiod;
  double minp;
  double maxp;
  double minf;
  double maxf;
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
  double ***harm_ampspec;
  double *harm_minamp;
  double *harm_maxamp;
  double *harm_phasefix;
  double ***harm_phasespec;
  double *subharm_ampfix;
  double ***subharm_ampspec;
  double *subharm_minamp;
  double *subharm_maxamp;
  double *subharm_phasefix;
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
  double paramfix[14];
  double minp;
  double maxp;
  double minf;
  double maxf;
  double minRp;
  double maxRp;
  double minMp;
  double maxMp;
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
  double a0, b0, chi0, inclination0, alpha0, psi00, mconst0;
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
  char modelsuffix[14];
  double fixedperiod;
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
  double **bper, **bt0, **bpow, **sde, **snval, **depth, **qtran;
  double **chisqrplus, *chisqrminus, *bperpos, *meanmagval, **fraconenight, **rednoise, **whitenoise, **sigtopink, **qingress, **OOTmag;
  int **nt, **Nt, **Nbefore, **Nafter;
  double **i1_ph, **i2_ph;
  int nf, *nf2, nbins, Npeak, **i1, **i2, operiodogram;
  int *sizeuv, rflag;
#ifdef PARALLEL
  double **p;
#else
  double *p;
#endif
  double minper, maxper, timezone;
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
  int fittrap;
  int nobinnedrms;
  int freqsteptype;
  int adjust_qmin_mindt;
  int reduce_nb;
  int reportharmonics;
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
} _BlsFixPer;

typedef struct {
  double **u, **v, *fmin, df;
  int durtype;
  double fixdur;
  OutColumn *fixdur_linkedcolumn;
  int TCtype;
  double fixTC;
  OutColumn *fixTC_linkedcolumn;
  double *inputTC, *inputdur;
  int fixdepth;
  int depthtype;
  double fixdepthval;
  OutColumn *fixdepth_linkedcolumn;
  int qgresstype;
  double qgressval;
  OutColumn *fixqgress_linkedcolumn;
  double *inputdepth, *inputqgress;
  double **bper, **bt0, **bpow, **sde, **snval, **depth, **qtran;
  double **chisqrplus, *chisqrminus, *bperpos, *meanmagval, **fraconenight, **rednoise, **whitenoise, **sigtopink, **qingress, **OOTmag;
  int **nt, **Nt, **Nbefore, **Nafter;
  double **i1_ph, **i2_ph;
  int nf, *nf2, Npeak, **i1, **i2, operiodogram;
  int *sizeuv, rflag;
#ifdef PARALLEL
  double **p;
#else
  double *p;
#endif
  double minper, maxper, timezone;
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
} _BlsFixDurTc;

typedef struct {
  int pertype;
  int lastlsindex;
  int lastaovindex;
  int lastblsindex;
  int t0type;
  double fixperiod;
  double fixT0;
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
  double binsize, firstbin;
  int Nvar;
  int *binstats;
  double *pctval;
  char **binvarnames;
  _Variable **binvars;
  char *binvarstring;
  int only_bin_columns;
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
  double P0, T00, r0, a0, inc0, bimpact0, sin_i0, e0, omega0, mconst0, ldcoeffs0[4], K0, gamma0;
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
  double offset;
} _DiffFluxtomag;

typedef struct {
  double mag_constant1;
  double offset;
} _Fluxtomag;

typedef struct {
  double **trends, *trendx, *trendy, **u, **v, *w1, *JD, clipping, pixelsep, *ave_out, *rms_out, **lcx, **lcy;
  char **stringid;
  int *stringid_idx;

  double ***u2, ***v2, **w2, **a, **b, **m_out;

  int *Njd_mout;

  int Ntrends, Njd, correctlc, ocoeff, omodel, Nskip_trend, JDcol_trend, magcol_trend;
  char trend_list_name[MAXLEN], dates_name[MAXLEN], **trend_names, coeff_outdir[MAXLEN], model_outdir[MAXLEN], coeff_suffix[MAXLEN], model_suffix[MAXLEN];
  int jdcol_isfromheader, magcol_isfromheader;
  char jdcol_headername[MAXLEN], magcol_headername[MAXLEN];
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
  int dotfafirst, use_bin, nbins, use_period, maxiter, pertype, lastindex, use_harm, Nharm, Nsubharm;
  double tfathresh;
  int decorrflag;
  int decorr_iterate;
  int decorr_Nlcterms, Ndecorr, Ntfatot;
  int *decorr_lc_order, *decorr_lc_columns;
  double ***lcdecorr_terms_in;
  int jdcol_isfromheader, magcol_isfromheader;
  char jdcol_headername[MAXLEN], magcol_headername[MAXLEN];
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
  int copyheaderfrominput;
  int logcommandline;
  int noclobber;
  _Variable **variables;
  char **printfformats;
  char **varnames;
  char sepchar;
} _Outputlcs;

typedef struct {
  double *rmsval;
  double *ave;
  int *ngood;
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
  double **gammaval;
  int gammaval_type;
  double sig_r_fix;
  double **sig_r;
  int sig_r_type;
  double sig_w_fix;
  double **sig_w;
  int sig_w_type;
  
  double **rho_r;
  int rho_r_type;
  double rho_r_fix;
  
  double **nu_r;
  int nu_r_type;
  double nu_r_fix;

  int bintime_type;
  double bintime_fix;
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
  double **ravals;
  double **decvals;
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
  double **inputravals;
  double **inputdecvals;
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
  _Variable **params;
  _Expression **expressions;
  _Expression *constantexpression;
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
  long amoeba_maxsteps;

  long mcmc_Naccept;
  long mcmc_Nlinkstotal;
  double mcmc_burninfrac;
  double mcmc_eps;
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
} _UserDataPointer;

typedef struct {
  _UserLib *lib;
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
} _RestrictTimes;

typedef struct {
  int restrictnum;
  _RestrictTimes *RestrictTimes;
} _RestoreTimes;

typedef struct {
  double cterm;
  double maxfreq;
  double freq_sample_factor;
  int auto_tau0;
  int auto_tau1;
  int auto_dtau;
  double tau0;
  double tau1;
  double dtau;
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
  int isvarlist;
  char *varliststring;
  _Variable **variables;
  int Nvar;
  int groupmethod;
  char **sourceID_inlist;
  int sourceID_N_split_substr;
  int *sourceID_split_substr;
  char **sourceID_splitstr;
  int *sourceID_split_leftright;
  int *sourceID_substr_startpos;
  int *sourceID_substr_length;

  int mergemethod;
} _Stitch;

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

  int calc_full_spec;

  int forcefft;

  int ofourier;
  char *ofourier_dir;
  int ofourier_formatflag;
  char *ofourier_format;

  _Variable *freq_var;

} _HarmonicFilter;

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

} _PythonCommand;
#else
typedef struct {
  int cnum;
} _PythonCommand;
#endif

typedef struct {
  int cnum;
  int require_sort;
  int require_distinct;
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
  _Stitch *Stitch;
  _HarmonicFilter *HarmonicFilter;
  _PythonCommand *PythonCommand;
  _FFT *FFT;

  int N_setparam_expr;
  char **setparam_EvalExprStrings;
  _Expression ***setparam_EvalExpressions;
} Command;
