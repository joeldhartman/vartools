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
#include "statistics.h"
int main(int, char **);
void example(char *, ProgramData *);
void help(char *, ProgramData *p);
void usage(char *);
int parseone(char *, void *, int);
int parseonedelimstring(char *, void *, int, char *);
int parseonedelimchar(char *, void *, int, char);
int skipone(char *);
void increaseNcommands(ProgramData *p, Command **c);
void dotab(FILE *,int);
void printheader(int, Command *, int, int);
void printresults(int, char *, Command *, int, int, int);
double readdates(char *, double);
void error(int);
void error2(int, char *);
void parsecommandline(int, char **, ProgramData *, Command **);
int ReadAllLightCurves(ProgramData *, Command *);
int ReadSingleLightCurve(ProgramData *, Command *, int, int);
int findX(double *, double, int, int);
void difffluxtomag(double *t, double *mag, double *sig, int N, double mag_star, double mag_constant1, double offset);
void fluxtomag(double *t, double *mag, double *sig, int N, double mag_constant1, double offset);
double binnedchi2(int, double *, double *, double *, double, double *, int *);
double chi2(int, double *, double *, double *, double *, int *);
double binnedrms(int, double *, double *, double *, double, double *, double *, int *);
double rms(int, double *, double *, double *, double *, double *, int *);
void ProcessCommandSingle(ProgramData *, Command *, int, int, int);
void ProcessCommandAll(ProgramData *, Command *, int);
void writelightcurves(ProgramData *p, int threadid, int lcid, char *outname, 
		      int usecolumnformat, int Nvars, _Variable **variables, 
		      char **formats, int noclobber, char sepchar, int logcommandline);
void ReadGlobalDecorr(ProgramData *, Command *);
void DetermineColumns(ProgramData *, Command *);
void Filldecorr_matrix(ProgramData *, Command *, int);
void ReadDatesFiles(ProgramData *, Command *);
void Switchtobasename(ProgramData *, int);
double doalarm(int, double *, double *);
void normalize(int, double *, double *, double *, double *, double *);
//double TestPeriod(int, double *, double *, double, int, _HistType *h);
double TestPeriod_aov_harm(int N, double *t, double *m, double *sig, int Nharm, double testperiod, double *m_noave, double *t_nostart, double *weight, double lcvariance, int *Nharm_used);
void aov_harm(int N, double *t, double *m, double *sig, int Nharm, int Nfreq, double freqmin, double freqstep, double *periodogram, double **out_m_noave, double **out_t_nostart, double **out_weight, double *out_lcvariance);
//void AOVPeriodogram(int, double *, double *, int, double *, double *, int, _HistType *h);
int isDifferentPeriods(double, double, double);
int isDifferentPeriodsDontCheckHarmonics(double, double, double);
void findPeaks_aov(double *t_, double *mag_, double *sig_, int N, double *perpeaks, double *aovpeaks, double *aovSNR, double *aovFAP, int Npeaks, double minP, double maxP, double subsample, double fine_tune, int outflag, char *outname, double *aveaov, double *stddevaov, double *aveaov_whiten, double *stddevaov_whiten, int ascii, int Nbin, int whiten, int uselog, double clip, int clipiter, int fixperiodSNR, double fixperiodSNR_period, double *fixperiodSNR_value, double *fixperiodSNR_SNR, double *fixperiodSNR_FAP);
void findPeaks_aovharm(double *t, double *mag, double *sig, int N, double *perpeaks, double *aovpeaks, double *aovSNR, double *aovFAP, int *Nharm_used, int Npeaks, double minP, double maxP, double subsample, double fine_tune, int outflag, char *outname, double *aveaov, double *stddevaov, double *aveaov_whiten, double *stddevaov_whiten, int ascii, int Nharm, int whiten, double clip, int clipiter, int fixperiodSNR, double fixperiodSNR_period, double *fixperiodSNR_value, double *fixperiodSNR_SNR, double *fixperiodSNR_FAP);
void ludcmp(long double **, int, int *, long double *);
void lubksb(long double **, int, int *, long double *);
void docorr(double *, double *, int, int, double **, int *, double *, double *, double, int);
void magcorr(void *,int,double *, double *, int, int, double **, int *, double *, double *, double *, double,int,char *, int);
void magcorr_chi2only(double *,double *, double *, int, int, double **, int *, double *, double *, double *, double, int, char *, int);
void w_ave(int, double *, double *, double *, double *);
void getJstet(int, double, double, double *, double *, double *, double *, double *, double *, double *);
void dokillharms(int, double *, double *, double *, int, double *, int, int, double **, double **, double **, double **, double *, double *, double *,int,char *, double *, int, int, double);
void doinjectharm(int, double *, double *, double *, int, int, _Injectharm *, char *);
void doinjecttransit(int N, double *t, double *mag, double *sig, int lc, int lcreal, _Injecttransit *c, char *modeloutname);
int gfasper(double *, double *, double *, int, double, double, double *, double *, int, int *, int *);
int fasper(double *, double *, int, double, double, double *, double *, int, int *);
int spread(double, double *, int, double, int);
void realft(double *, int, int);
void four1(double *, int, int);
void avevar(double *, int, double *, double *);
void Lombscargle (int N, double *t, double *mag, double *sig, double minper, double maxper, double subsample, int Npeaks, double *periods, double *peaks, double *probs, double *SNR, int outputflag, char *outfile, int ascii, int whiten, double clip, int clipiter, int fixperiodSNR, double fixperiodSNR_period, double *fixperiodSNR_FAPvalues, double *fixperiodSNR_SNRvalues, double *fixperiodSNR_peakvalues, int use_orig_ls, int dobootstrapfap, int Nbootstrap);
void mysort4(int, double *, double *, double *, double *);
void mysort4_rev(int, double *, double *, double *, double *);
void mysort4ptrint(int, int*, void***, int *, int *);
void mysort3_int(int, double *, double *, int *);
void mysort3(int, double *, double *, double *);
void mysort3_rev(int, double *, double *, double *);
void mysort2(int, double *, double *);
void mysort2_rev(int, double *, double *);
void mysort2dblint(int, double *, int *);
void mysort2ptr(int, int *, double ***);
void mysort3ptrint(int, int *, void ***, int *);
void mysort1int(int, int *);
void mysort1(int, double *);
void mysort1_rev(int, double *);
void mysort2dblint_id(int N, double* data1, int* data2);
void rescalesigma_linear(int, double *, double, double);
void rescalesigma_chi2(int, double *, double);
long double mean(int, double *, double *);
long double mean1(int, double *);
long double mean2(int, double *, double *, double *);
int purge_bad(int *, double *, double *, double *, double, double, double, int);
int sigclip(int, double *, double *, double *, double *, double *, double *, int *, double, int, int, ProgramData *, int, int);
void sigclip_copyterms(int i,int j,ProgramData *p,int lc);
double chisqstarspot(double *, int, int, double *, double *, double *, void *);
int amoeba(double **, double *, int *, int, double, double (*funk)(double *, int, int, double *, double *, double *, void *), int *, int, int, double *, double *, double *, void *);
double amotry(double **, double *, int *, double *, int, int, double (*funk)(double *, int, int, double *, double *, double *, void *), int, double, int, double *, double *, double *, void *);
void starspot(double *, double *, double *, double **, int, int, void *);
void starspot_(int, double *, double *, double, double, double, double, double, double, double);
void fitstarspot_amoeba(int N, double *t, double *mag, double *sig, double *P, double *a, double *b, double *alpha, double *i, double *chi, double *psi0, double *mconst, int fitP, int fita, int fitb, int fitalpha, int fiti, int fitchi, int fitpsi0, int fitmconst, double *chi2_, int correctlc, int omodel, char *modelname);
double getfrac_onenight(int n,double *t,double *u, double *v,double *err,double bper,double depth,double qtran,double bt0,double timezone);
double getclippedsrave(int n, double *sr);
double getclippedstddev(int n, double *pow);
double subtract_binnedrms(int N, double *mag, double bintime, double *aveval, int *ngood, double *binmag, double *binsig);
int eebls(int n, double *t, double *x, double *e, double *u, double *v, int nf, double fmin, double df, int nb, double qmi, double qma, double *p, int Npeak, double *bper, double *bt0, double *bpow, double *sde, double *snval, double *depth, double *qtran, int *in1, int *in2, double *in1_ph, double *in2_ph, double *chisqrplus, double *chisqrminus, double *bperpos, double *meanmagval, double timezone, double *fraconenight, int operiodogram, char *outname, int omodel, char *modelname, int correctlc, int ascii,int *nt, int *Nt, int *Nbefore, int *Nafter, double *rednoise, double *whitenoise, double *sigtopink, int fittrap, double *qingress, double *OOTmag, int ophcurve, char *ophcurvename, double phmin, double phmax, double phstep, int ojdcurve, char *ojdcurvename, double jdstep, int nobinnedrms, int freq_step_type, int adjust_qmin_mindt, int reduce_nb, int reportharmonics);
int eebls_rad(int n, double *t, double *x, double *e, double *u, double *v, int nf, double fmin, double df, int nb, double rmin, double rmax, double *p, int Npeak, double *bper, double *bt0, double *bpow, double *sde, double *snval, double *depth, double *qtran, int *in1, int *in2, double *in1_ph, double *in2_ph, double *chisqrplus, double *chisqrminus, double *bperpos, double *meanmagval, double timezone, double *fraconenight, int operiodogram, char *outname, int omodel, char *modelname, int correctlc, int ascii,int *nt, int *Nt, int *Nbefore, int *Nafter, double *rednoise, double *whitenoise, double *sigtopink, int fittrap, double *qingress, double *OOTmag, int ophcurve, char *ophcurvename, double phmin, double phmax, double phstep, int ojdcurve, char *ojdcurvename, double jdstep, int nobinnedrms, int freq_step_type, int adjust_qmin_mindt, int reduce_nb, int reportharmonics);
void phaselc(int N, double *t, double *mag, double *sig, double period, int is_T0_given, double T0, char *phasevarname, _Variable *phasevar, int threadid, double startphase);
void softened_transit_(int N, double *t, double *mag, double P, double T0, double eta, double c, double delta, double mconst, int Nlin_coeffs, double *lin_coeffs, double **Design_Matrix);
void softened_transit(double *t, double *a_, double *yfit, double **dyda, int ma, int N, int Nlin_coeffs, double *lin_coeffs, double **Design_Matrix, double *y, double *sig, int *varylin_coeffs, void *userparams);
void fitsoftened_transit(int N, double *t, double *mag, double *sig, double *P, double *T0, double *eta, double *cval, double *delta, double *mconst, int fitephem, int fiteta, int fitcval, int fitdelta, int fitmconst, double *chi2_, int correctlc, int omodel, char *modelname, int dokillharm, int nharm, int nsubharm, double perharm, double *subharmA, double *subharmB, double *harmA, double *harmB, double *fundA, double *fundB);
#ifdef PARALLEL
void mrqmin(double *x, double *y, double *sig, int ndata, double *a, int *ia, int ma, double **covar, double **alpha, double *chisq, double *alamda, void (*funcs)(double *, double *, double *, double **, int, int, int, double *, double **, double *, double *, int *, void *), int Nlin_coeff, double **Design_Matrix, double *lin_coeffs, int *varylin_coeffs, int mfit, double *ochisq, double *atry, double *beta, double *da, double **oneda, void *userparams);
#else
void mrqmin(double *x, double *y, double *sig, int ndata, double *a, int *ia, int ma, double **covar, double **alpha, double *chisq, double *alamda, void (*funcs)(double *, double *, double *, double **, int, int, int, double *, double **, double *, double *, int *, void *), int Nlin_coeff, double **Design_Matrix, double *lin_coeffs, int *varylin_coeffs, void *userparams);
#endif
void mrqcof(double *x, double *y, double *sig, int ndata, double *a, int *ia, int ma, double **alpha, double *beta, double *chisq, void (*funcs)(double *, double *, double *, double **, int, int,int,double *, double **, double *, double *, int *, void *), int Nlin_coeff, double **Design_Matrix, double *lin_coeffs, int *varylin_coeffs, void *userparams);
void covsrt(double **covar, int ma, int *ia, int mfit);
int gaussj(double **a, int n, double **b, int m);
void svdcmp(double **a, int m, int n, double *w, double **v);
void svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x);
int eeblsfixper(int n, double *t, double *x, double *e, double *u, double *v, int nb, double qmi, double qma, double *period, double *bt0, double *bpow, double *depth, double *qtran, int *in1, int *in2, double *in1_ph, double *in2_ph, double *chisqrplus, double *chisqrminus, double *meanmagval, double timezone, double *fraconenight, int omodel, char *modelname, int correctlc, int ascii,int *nt, int *Nt, int *Nbefore, int *Nafter, double *rednoise, double *whitenoise, double *sigtopink, int fittrap, double *qingress, double *OOTmag);
int eeblsfixper_rad(int n, double *t, double *x, double *e, double *u, double *v, int nb, double rmin, double rmax, double *period, double *bt0, double *bpow, double *depth, double *qtran, int *in1, int *in2, double *in1_ph, double *in2_ph, double *chisqrplus, double *chisqrminus, double *meanmagval, double timezone, double *fraconenight, int omodel, char *modelname, int correctlc, int ascii,int *nt, int *Nt, int *Nbefore, int *Nafter, double *rednoise, double *whitenoise, double *sigtopink, int fittrap, double *qingress, double *OOTmag);
double ls_oneperiod(double *x, double *y, int n, double period);
double gls_oneperiod(double *x, double *y, double *err, int n, double period);
void getlsampthresh(int N, double *t, double *mag, double *sig, double period, int harm_specsigfile, FILE *signalfile, int Nsubharm, int Nharm, double minPer, double thresh, double *ampthresh_scale, double *amp, int use_orig_ls);
double getlsampthresh_func(double amp_scale, int N, double *t, double *mag_orig, double *sig, double *mag_signal, double *mag_tmp, double period, int Nsubharm, int Nharm, double *subharmA, double *subharmB, double *harmA, double *harmB, double fundA, double fundB, double minPer, double thresh, int use_orig_ls);
double getminsini(double a, double e, double omega, double p);
void fitmandelagoltransit_amoeba(int N, double *t, double *mag, double *sig, double *P, double *T0, double *r, double *a, double *inc, double *bimpact, double *e, double *omega, double *mconst, int type, double *ldcoeffs, int fitephem, int fitr, int fita, int fitinclterm, int fite, int fitomega, int fitmconst, int *fitldcoeffs, double *chi2_, int correctlc, int omodel, char *modelname, int fitRV, char *RVfilename, char *omodelRVcurve, double *K, double *gamma, int fitK, int fitgamma, int refititer, int ophcurve, char *ophcurvename, double phmin, double phmax, double phstep, int ojdcurve, char *ojdcurvename, double jdstep, char *modelvarname, _Variable *modelvar, int threadid);
void initialize_tfa(_TFA *tfa, ProgramData *p);
void detrend_tfa(_TFA *tfa, int N, double *t, double *m, double *e, double lcx, double lcy, char *lc_name, char *coeff_file_name, int coeff_flag, int correctlc, int outlc, char *lc_out_name, double *ave_out, double *rms_out, int matchstringid, char **stringid, int *stringid_idx, int threadid);
void do_sysrem(_Sysrem *Sysrem, int numlc, int *Njd_in, double **t_in, double **mag_in, double **sig_in, char **lcnames, int matchstringid, char ***stringid, int **stringid_idx);
void initialize_sysrem(_Sysrem *Sysrem, int numlcs, int matchstringid);
int binlc_parsevarstring(_Binlc *c);
void binlc(ProgramData *p, _Binlc *c, int lcnum);
double changeerror(int N, double *t, double *mag, double *sig, double *aveval, int *ngood);
#ifdef isinf
#else
int isinf(double);
#endif
double getrms(int N, double *t, double *mag, double *sig, double *aveval, double *rmsthy, int *ngood);
void redwhitenoise(int N, double *t, double *mag, double *sig, double timespan, double *rednoise, double *whitenoise);
void subtractbls(int N, double *t, double *mag, double *sig, double P, double q, double depth, double ph1, int *nt, int *Nt, int *nbefore, int *nafter, double qingress, double OOTmag);
void getsignaltopinknoiseforgivenblsmodel(int N, double *t, double *mag, double *sig, double P, double q, double depth, double in1_ph, int *nt, int *Nt, int *Nbefore, int *Nafter, double *rn, double *wn, double *sigtopink, double qingress, double OOTmag);
void initialize_tfa_sr(_TFA_SR *tfa, int Nlcs, ProgramData *p);
void detrend_tfa_sr(_TFA_SR *tfa, int N, double *t, double *m, double *e, double lcx, double lcy, char *lc_name, char *coeff_file_name, int coeff_flag, int correctlc, int outlc, char *lc_out_name, double *ave_out, double *rms_out, double period, char *signalfilename, int matchstringid, char **stringid, int *stringid_idx, int lcindex, int threadid);
void dodftclean(int N, double *t, double *mag, double *sig, int lc, _Dftclean *c, char *lcbasename, int ascii);
void mysortstringint(int N, int sizestr, char **data1, int *data2);
void autocorrelation(double *t, double *mag, double *sig, int N, double tmin, double tmax, double tstep, char *outname);
void mandelagoltransitmodel(int Npoints, double *phase, double *outlc, int type, double *ldcoeffs, double sin_i, double a, double e, double p, double omega);
void dorestorelc(ProgramData *p, _Savelc *s, _Restorelc *r, int sthreadid, int rthreadid, int lcid);
void dosavelc(ProgramData *p, _Savelc *s, int threadid, int lcid);
void CreateOutputColumns(ProgramData *p, Command *c, int Ncommands);
void printheader_new(ProgramData *p, FILE *outfile);
void printresults_new(ProgramData *p, int lc, int reallc, FILE *outfile);
void AOV_getfixedperiodSNR(double *p, double *pr, int N, double ave, double rms, double Pin, double *peak, double *SNR);
void increaselinkedcols(ProgramData *p, OutColumn **c, char *s, int cmdidx);
void linkcolumns(ProgramData *p);
void medianfilter(int N, double *t, double *mag, double *sig, double timesize, int meanflag, int replace);
void getoutcolumnvalue(OutColumn *c, int lc, int reallc, int outtype, void *outvalue, ...);
void setoutcolumnvalue(OutColumn *c, int lc, int reallc, int intype, void *invalue, ...);
void findblends(int Nvars, int *N, double **t, double **mag, double **sig, _FindBlends *c);
int findXindx(double *x, int *indx, double xval, int i1, int N);
double distxy(double x1, double y1, double x2, double y2);
double distradec(double ra1, double dec1, double ra2, double dec2);
void astr_iarc(double xi, double eta, double rac, double decc, double *ra, double *dec);
void astr_irarc(double xi, double eta, double rac, double decc, double *ra, double *dec);
void microlens(double *t, double *mag, double *sig, int N, int lc, _MicroLens *m, char *outname, double *f0_out, double *f1_out, double *u0_out, double *t0_out, double *tmax_out, double *chi2_);
double chisqtraptransit(double *a_, int ma, int N, double *t, double *mag, double *sig, void *userparams);
void dofittrap_amoeba(int N, double *t, double *mag, double *sig, double P, double *q, double *qingress, double *in1_ph, double *in2_ph, double *depth, double *OOTmag);
void dofittrap_amoeba_fixdur(int N, double *t, double *mag, double *sig, double P, double q, double *qingress, double in1_ph, double in2_ph, double *depth, double *OOTmag);
int sortlcbytime(int size, double *t, int lc, ProgramData *p);
double ran1(void);
double gasdev(void);
void addnoise(ProgramData *p, _AddNoise *c, int threadid, int lcid);
void printtostring_indentwrap(OutText *text, const char *stoadd, int Ntab_indent);
void printtostring(OutText *text, const char *stoadd);
void printtostring_nowrap(OutText *text, const char *stoadd);
void convertUTCtoJD(char *inputUTC, char *UTCformat, int *UTCindex, double *outJD);
void checkUTCFormat(char *UTCformat, int *UTCindex);
void RegisterDataFromInputList(ProgramData *p, void *dataptr, int datatype,
                               int Ncolumns, int cnum, int disjointcolumns,
                               int Nonuniformnames, char *scanformat, ...);
void MemAllocDataFromInputList(ProgramData *p, int Nlc);
void ParseInputList(ProgramData *p, char **inputlines, int Nlcs);
void printinputlistformat(ProgramData *p, FILE *outfile);
void converttime(int N, double *t, int lc, int lcreal, _ConvertTime *c);
double eccentricAnomaly (double M, double e);
void GetOutputFilename(char *lcoutname, char *lcname, char *outdir,
		       char *suffix, char *format, int lc_name_num);
void incrementparameters_foramoeba(int *Nparameters, int *Ntovary, double ***p, int **ia, int varyparam, double initialval, double stepval);
void amoeba_cleanup(int *Nparameters, int *Ntovary, double ***p, int **ia, double **chi2vals);
void amoeba_initializesimplexchi2(int Nparameters, int Ntovary, double **p, double **chi2vals, double (*funk)(double *, int, int, double *, double *, double *, void *), int N, double *t, double *mag, double *err, void * userparam);
void RegisterDataFromLightCurve(ProgramData *p, void *dataptr, int datatype,
				int maxstringlength, int Ncolumns, int cnum, 
				int disjointcolumns,
				int Nonuniformnames, char *scanformat, 
				_Variable *variable, ...);
void InitializeMemAllocDataFromLightCurve(ProgramData *p, int Nthread);
void MemAllocDataFromLightCurve(ProgramData *p, int threadid, int Nterm);
void MemAllocDataFromLightCurveMidProcess(ProgramData *p, int threadid, int Nterm);
void CompileAllExpressions(ProgramData *p, Command *c);
void RunExpressionCommand(int lcindex, int threadindex, 
			  ProgramData *p, _ExpressionCommand *c);
_ExpressionCommand* CreateExpressionCommand(ProgramData *p, char *argv);
_Variable* CreateVariable(ProgramData *p, char *varname, char datatype, char vectortype, void *vptrinput, ...);
double EvaluateExpression(int lcindex, int threadindex, int jdindex, _Expression *expression);
void SetVariable_Value_Double(int lcindex, int threadindex, int jdindex, _Variable *var, double val);
double EvaluateVariable_Double(int lcindex, int threadindex, int jdindex, _Variable *var);
double EvaluateFunctionCall(int lcindex, int threadindex, int jdindex, _FunctionCall *call);
int CheckIsFunctionConstantVariableExpression(char *term, ProgramData *p, char *functionid, double *constval, _Variable **varptr);
_FunctionCall* ParseFunctionCall(char *term, ProgramData *p, char functionid);
_Expression* SplitExpression(char *term, char operatortype, int i1, int sizeterm, ProgramData *p);
_Expression* ParseExpression(char *term, ProgramData *p);
void PrintVartoolsFunctionList(ProgramData *p);
void ParseOutputColumnFormat(_Outputlcs *o);
void SetTimeMagSigPointers(ProgramData *p, int threadid);
void DoLinfit(ProgramData *p, _Linfit *c, int threadid, int lcid);
void InitLinfit(ProgramData *p, _Linfit *c, int cnum);
void SetupLinfitExpression(ProgramData *p, _Linfit *c);
void spline(double *x,double *y,int n,double yp1,double ypn,
            double *y2,double *u);
void splint(double *xa,double *ya,double *y2a,int n,double x,double *y);
void spline_monotonic(int N, double *x, double *y, double *yprime);
double splint_monotonic(int N, double *x, double *y, double *yprime, double xt);
void integratemandelagoltransitmodel(double exptime_phase, int Npoints, double *phase, double *outlc, int type, double *ldcoeffs, double sin_i, double a, double e, double p, double omega, int Nresamp);
double fitpoly(int N, double *x, double *y, double *sig, int order,
	       int subtractfit, double *fitparams, double *paramerrs);
void vsort_generic(int N, int isreverse, int *index, int Nms, mysort_generic_struct *s);
void printinputlcformat(ProgramData *p, FILE *outfile);
void InitCommands(ProgramData *p, Command *c);
int ParseStatsCommand(int *iret, int argc, char **argv, ProgramData *p, _Stats *s);
int ParseLinfitCommand(int *iret, int argc, char **argv, ProgramData *p,
		       _Linfit *c);
void zerooutcolumnvalue(OutColumn *c, int lc, int reallc);
void addcolumn(ProgramData *p, int cnum, int type, int stringsize, void *ptr, char *outputformat, int Ndereference, int usereallc, int lcdereferencecol, ...);
_IfStack *CreateIfStack(void);
char TestIf(_IfStack *stack, ProgramData *p, Command *c, int lcindex, int threadindex);
int ParseIfCommand(int *iret, int argc, char **argv, int cn, ProgramData *p, Command *c);
int check_isspecialchar(char c);
void SetupInListVariable(ProgramData *p, char *varname, int column, int datatype, char *format);
void RunUserCommand(ProgramData *, Command *, int, int);
void RunUserCommand_all_lcs(ProgramData *, Command *);
int listcommands_noexit(char *c, ProgramData *p, OutText *s);
void listcommands(char *c, ProgramData *p);
void usage(char *argv);
void help(char *c, ProgramData *p);
void RestrictTimes_readJDlist(char *filename, double **JDlist, int *Nlist);
void RestrictTimes_readimagelist(char *filename, char ***imagelist, int **imagelist_indx, int *Nlist);
int RestrictTimes_JDrange_apply(int N, double *t,
				int lc, ProgramData *p, _RestrictTimes *c,
				double JDmin, double JDmax, char exclude);
int RestrictTimes_JDlist_apply(int N, double *t,
			       int lc, ProgramData *p, _RestrictTimes *c,
			       double *JDlist, int Nlist, char exclude);
void RestrictTimes_imagelist_apply(int N, char **stringID, int *stringID_indx, 
				   int lc, ProgramData *p, _RestrictTimes *c,
				  char **imagelist, int *imagelist_indx, 
				  int Nlist, char exclude);
void RestrictTimes_ParseExpr(int *iret, int argc, char **argv, ProgramData *p, _RestrictTimes *RestrictTimes, char min_or_max);
void MoveInputListData(ProgramData *p, int isrc, int idest);
void RemoveEmptyLightCurves(ProgramData *p, Command *c);
void DoOutputLightCurve(ProgramData *p, _Outputlcs *c, int lcid, int threadid);
int gnu_getline(char **, size_t *, FILE *);
char * GenerateInternalVariableName(ProgramData *p);
_Variable *SetupScalarVariable(ProgramData *p, char *varname, int datatype);
void RegisterScalarData(ProgramData *p, void *dataptr, int datatype, int Ncolumns);
void MemAllocScalarData(ProgramData *P, int Nthreads);
void MCMC_CleanUp(_MCMC_Chain **);
void MCMC_initialize_chain(_MCMC_Chain *c, double eps, double maxmem, int doprintchain, char *outchainfilename, char *outchainheader, int printevery, int Nauxil, double **auxil_params, double (*funk)(double *, int, int, double *, double *, double *, void *), int N, double *t, double *mag, double *err, void *userparam);
void MCMC_DifferentialEvolution_RunMCMC(int maxNtrials_accepted, int maxNlinks_run, int N_recalc_covar, _MCMC_Chain *c);
int MCMC_GetBestLink(_MCMC_Chain *c);
void MCMC_increment_parameter(_MCMC_Chain **cptr, double value, double delp);
void DoWWZ(ProgramData *p, _WWZ *c, int threadid, int lcid);
int ParseWWZCommand(int *iret, int argc, char **argv, ProgramData *p,
		    _WWZ *c);
void dosaveifstackcopy(ProgramData *p, _CopyLC *c, int threadid);
void dorestoreifstackcopy(ProgramData *p, _CopyLC *c, int sthreadid, int rthreadid);
void turnoffcopies(ProgramData *p, Command *c, int cnum_start, int threadid, int lcid);
void turnoffcopies_onecommand(ProgramData *p, Command *c, int threadid, int lcid);
_CopyLC * CreateCopyLCCommand(ProgramData *p, char *argv, int cnum);
void SetupLCCopies(ProgramData *p, Command *c);
void docopylccommand(ProgramData *p, _CopyLC *c, int threadid, int lcid);
void getlccopy(ProgramData *p, Command *c, int threadid, int lcid);
double gammln(double);
void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp);
double log1minusbetai(double a, double b, double x);
double betai(double a, double b, double x);
void choldc_sparse_neardiag(double **a, int N, int *Nvec, double *p);
void cholsl_sparse_neardiag(double **a, int N, int *Nvec, double *p, double *b, double *x);
void cholmult_sparse_neardiag(double **a, int N, int *Nvec, double *p, double *b, double *x);
void cholsl(double **a, int N, double *p, double *b, double *x);
void choldc(double **a, int N);
void mergeequallctimes(ProgramData *p, int lc);
void DoResample(ProgramData *p, _Resample *c, int threadid, int lcid);
void SetupResampleExpression(ProgramData *p, _Resample *c);
int ParseResampleCommand(int *iret, int argc, char **argv, ProgramData *p,
			 _Resample *c, int cnum);
int eeblsfixdurtc(int n, double *t, double *x, double *e, double *u, double *v, double inputTC, double inputdur, int fixdepth, double inputdepth, double inputqgress, int nf, double fmin, double df, double *p, int Npeak, double *bper, double *bt0, double *bpow, double *sde, double *snval, double *depth, double *qtran, double *chisqrplus, double *chisqrminus, double *bperpos, double *meanmagval, double timezone, double *fraconenight, int operiodogram, char *outname, int omodel, char *modelname, int correctlc, int ascii,int *nt, int *Nt, int *Nbefore, int *Nafter, double *rednoise, double *whitenoise, double *sigtopink, int fittrap, double *qingress, double *OOTmag, int ophcurve, char *ophcurvename, double phmin, double phmax, double phstep, int ojdcurve, char *ojdcurvename, double jdstep);
int load_user_library(char *libname, ProgramData *p, int islib, ...);
void CheckCreateCommandOutputLCVariable(char *varname, _Variable **omodelvar, ProgramData *p);
void occultnl(double rl, double c1, double c2, double c3, double c4, double *b0, double *mulimb0, double **mulimbf, int nb);
void occultquad(double *z0, double u1, double u2, double p, double *muo1, double *mu0, int nz);
void growbuffer(char **buf, int *sizebuf, int newlen);
void printresults_buffer_new(ProgramData *p, int lc, int reallc, char **buf, int *sizebuf);
void emptyresults_buffer(ProgramData *p, FILE *outfile);
_StringBuffer *popFreeBuffer(ProgramData *p);
void pushFullBuffer(ProgramData *p, _StringBuffer *buf);
void pushFreeBuffer(ProgramData *p, _StringBuffer *buf);
void InitializeOutputBufferStacks(ProgramData *p);
int ParseParameterBuiltInCommand(ProgramData *p, int cnum,
				 int *iret, char **argv,
				 int argc, const char *keyword, int Nvec,
				 ...);
void doHarmonicFilter(ProgramData *p, _HarmonicFilter *c, int threadid, int lcid);
int ParseHarmonicFilterCommand(int *iret, int argc, char **argv, ProgramData *p,
			       _HarmonicFilter *c, int cnum);
void ParseDefineAnalyticUserFunction(ProgramData *p, char *argv);
void InitOutTextStruct(OutText *text);
#ifdef _HAVE_PYTHON
//void LoadVartoolsRunPythonLibrary(ProgramData *p);
int ParsePythonCommand(int *inum, int argc, char **argv, ProgramData *p, 
		       _PythonCommand *c, Command *allcommands, int cnum);
_PythonCommand *CreatePythonCommandStruct(ProgramData *p, char *argv0);
void RunPythonCommand(ProgramData *p, int lcindex, int threadindex, int pythreadindex, _PythonCommand *c);
void InitPythonCommand(ProgramData *p, _PythonCommand *c, int Nlcs);
void SetupRunPythonVariables(_PythonCommand *c, ProgramData *p);
void StopRunningPythonCommand(ProgramData *p, int threadindex, _PythonCommand *c);
void KillAllPythonProcesses(ProgramData *p, Command *allcommands);
#endif
#ifdef _HAVE_R
//void LoadVartoolsRunPythonLibrary(ProgramData *p);
int ParseRCommand(int *inum, int argc, char **argv, ProgramData *p, 
		       _RCommand *c, Command *allcommands, int cnum);
_RCommand *CreateRCommandStruct(ProgramData *p, char *argv0);
void RunRCommand(ProgramData *p, int lcindex, int threadindex, int Rthreadindex, _RCommand *c);
void InitRCommand(ProgramData *p, _RCommand *c, int Nlcs);
void SetupRunRVariables(_RCommand *c, ProgramData *p);
void StopRunningRCommand(ProgramData *p, int threadindex, _RCommand *c);
void KillAllRProcesses(ProgramData *p, Command *allcommands);
void StartAllRProcesses(ProgramData *p, Command *allcommands);
#endif
void RestoreTimes(ProgramData *p, _RestoreTimes *RestoreTimes, int sthreadid, int rthreadid);
long randlong(long);
void RunFFTCommand(ProgramData *p, int lcindex, int threadindex, _FFT *f);
int ParseFFTCommand(int *iret, int argc, char **argv, ProgramData *p, _FFT *f);
