/* This file contains code for the -transitbisec command of VARTOOLS
   which is used to calculate the bisectors and BS spans of transits
   that have been fully observed.

Written by J. Hartman. 2012 07 30.

*/

#include "../../src/vartools.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "transitbisec.h"

#ifndef ERR_MEMALLOC
#define ERR_MEMALLOC 3
#endif

double chisq_shiftmandelagoltransit(double *p, int ndim, int N, double *t, 
				    double *mag, double *sig, void *userparams)
/* Function to use with amoeba to fit an integrated Mandel & Agol transit model
   to a transit allowing the time of transit center to be varied. This is used
   in the -transitbisec command */
{
  _TransitBisec_passparams *tb;
  TransitParams *tpar;
  Bisector_Data *bdat;
  double Tc;
  double cos_i, sin_i, magval;
  double *phasemodel, *magmodel, *magresid, chi2val;
  int i, j, Nintegratesamp;

  tb = (_TransitBisec_passparams *) userparams;
  tpar = tb->t;
  bdat = tb->b;

  Tc = p[0];

  cos_i = tpar->bimpact*(1. + tpar->eccen*cos(tpar->omega))/(1. - tpar->eccen*tpar->eccen)/tpar->arstar;

  if(cos_i > 1.)
    cos_i = 1.;
  else if(cos_i < -1.)
    cos_i = -1.;
  sin_i = sqrt(1. - cos_i*cos_i);
  
  Nintegratesamp = ceil(tpar->exptime / TRANSITBISEC_TSAMP_INTEGRATE);
  if(Nintegratesamp < 1) Nintegratesamp = 1;
  

  /* Evalute the integrated transit model at the observed times. */
  if((magmodel = (double *) malloc(N*sizeof(double))) == NULL ||
     (phasemodel = (double *) malloc(N*sizeof(double))) == NULL ||
     (magresid = (double *) malloc(N*sizeof(double))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(i = 0; i < N; i++) {
    phasemodel[i] = (t[i] - Tc)/tpar->P;
    phasemodel[i] = phasemodel[i] - floor(phasemodel[i]);
  }
  VARTOOLS_integratemandelagoltransitmodel(
				  tpar->exptime/tpar->P,
				  N, phasemodel, magmodel,
				  tpar->ldtype, 
				  tpar->ldcoeffs,
				  sin_i, tpar->arstar,
				  tpar->eccen,
				  tpar->rprstar,
				  tpar->omega,
				  Nintegratesamp);
  for(j = 0; j < N; j++) {
    magmodel[j] = -2.5*log(magmodel[j])/M_LN10;
    magresid[j] = mag[j] - magmodel[j];
  }


  /* High-pass filter the residuals using a polynomial fit */
  if(bdat->bisec_resid_fit_poly_order > 0) {
    for(j=0; j < N; j++) {
      phasemodel[j] = t[j] - t[0];
    }
    VARTOOLS_fitpoly(N, phasemodel, magresid, sig,
		     bdat->bisec_resid_fit_poly_order, 1, NULL, NULL);
  }
  else {
    /* Subtract the mean from the residuals if we are not high-pass filtering */
    magval = VARTOOLS_getmean(N,magresid);
    for(j=0; j < N; j++) {
      magresid[j] -= magval;
    }
  }

  chi2val = 0.;
  if(sig != NULL) {
    for(i=0; i < N; i++) {
      chi2val += magresid[i]*magresid[i]/sig[i]/sig[i];
    }
  }
  else {
    for(i=0; i < N; i++) {
      chi2val += magresid[i]*magresid[i];
    }
  }
  free(magmodel);
  free(phasemodel);
  free(magresid);

  return chi2val;
}

void initialize_Bisector_Data(int lc_num, void *userdata,
			      void *userdata2)
/* This function is called when the TransitBisec->bisecdata struct is 
   allocated memory. Here we initialize some of the parameters in the 
   struct */
{
  Bisector_Data *bisec;
  _TransitBisec *TransitBisec;

  bisec = (Bisector_Data *) userdata;
  TransitBisec = (_TransitBisec *) userdata2;

  bisec->Nbisecs = 0;
  bisec->tvals = NULL;
  bisec->Ntvals = 0;
  bisec->sizebisec1 = 0;
  bisec->sizebisec2 = 0;
  bisec->bisec_high_val = TransitBisec->bisec_high_val;
  bisec->bisec_low_val = TransitBisec->bisec_low_val;
  bisec->bisec_step_val = TransitBisec->bisec_step_val;
  bisec->spanvals = NULL;
  bisec->transitnums = NULL;
  bisec->bisec_resid_bin_time = TransitBisec->bisec_resid_bin_time;
  bisec->bisec_resid_fit_poly_order = TransitBisec->bisec_resid_fit_poly_order;
  bisec->bisec_usemoments = TransitBisec->bisec_usemoments;
  bisec->bisec_high_samp_factor = TransitBisec->bisec_high_samp_factor;
}

void initialize_TransitParams(int lc_num, void *userdata,
			      void *userdata2)
/* This function is called when the TransitBisec->transitparams struct is 
   allocated memory. Here we initialize some of the parameters in the 
   struct */
{
  TransitParams *transitparams;
  _TransitBisec *TransitBisec;

  TransitBisec = (_TransitBisec *) userdata2;
  transitparams = (TransitParams *)userdata;
  transitparams->ldcoeffs = malloc(4*sizeof(double));
  transitparams->ldtype = TransitBisec->ldtype;
}

void expand_bisector_data(int Ntransit, Bisector_Data *bisec)
/* This function adds space to the bisec structure to store the bisector for
   another transit */
{
  int i;
  if(!bisec->Ntvals) {
    /* We haven't yet determined the number of time steps we will need to
       store, do that now */
    bisec->Ntvals = ceil(fabs(bisec->bisec_high_val - bisec->bisec_low_val)/bisec->bisec_step_val ) + 1;
  }
  if(bisec->Nbisecs + 1 > bisec->sizebisec1) {
    if(!bisec->sizebisec1) {
      if((bisec->tvals = (double **) malloc((bisec->Nbisecs + 1) * sizeof(double *))) == NULL ||
	 (bisec->transitnums = (int *) malloc((bisec->Nbisecs + 1) * sizeof(int))) == NULL ||
	 (bisec->spanvals = (double *) malloc((bisec->Nbisecs + 1) * sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    } else {
      if((bisec->tvals = (double **) realloc(bisec->tvals, (bisec->Nbisecs + 1)*sizeof(double *))) == NULL ||
	 (bisec->transitnums = (int *) realloc(bisec->transitnums, (bisec->Nbisecs + 1)*sizeof(int))) == NULL ||
	 (bisec->spanvals = (double *) realloc(bisec->spanvals, (bisec->Nbisecs + 1)*sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
    for(i = bisec->sizebisec1; i < bisec->Nbisecs + 1; i++) {
      if((bisec->tvals[i] = (double *) malloc((bisec->Ntvals > bisec->sizebisec2 ? bisec->Ntvals : bisec->sizebisec2)*sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
  }
  if(bisec->Ntvals > bisec->sizebisec2) {
    for(i=0; i < bisec->sizebisec1; i++) {
      if((bisec->tvals[i] = (double *) realloc(bisec->tvals[i], bisec->Ntvals * sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
    }
  }
  if(bisec->Nbisecs + 1 > bisec->sizebisec1)
    bisec->sizebisec1 = bisec->Nbisecs + 1;
  if(bisec->Ntvals > bisec->sizebisec2)
    bisec->sizebisec2 = bisec->Ntvals;
  bisec->transitnums[bisec->Nbisecs] = Ntransit;
  bisec->Nbisecs = bisec->Nbisecs + 1;
}

int get_transit_bisector(int N, double *t, double *mag, int i0, int i1,
			 double BS_highval, double BS_lowval, 
			 int Ntransit,
			 TransitParams *transitparams,
			 Bisector_Data *bisec,
			 int refit_Tc, double *fitted_Tc)
/* This function calculates the bisector for a single transit.
   N - number of points in the light curve.
   t - vector of lc times.
   mag - vecotr of lc mags.
   i0 - index of first point in lc to consider.
   i1 - index of last point in lc to consider.
   BS_highval - fraction of transit depth to use for the high point
                in the bisector span.
   BS_lowval - fraction of transit depth to use for the low point
               in the bisector span.
   Ntransit  - transit number to compute the bisector for.
   transitparams - input struct containing parameters describing the 
                   transit.
   bisec - input struct used to store the bisector. This struct also
           contains parameters defining how the bisector is to be calculated.
   refit_Tc - whether or not to refit the transit center.
*/
{
  double sin_i, cos_i;

  int Nhighsamp, i, j, Nintegratesamp;
  int Nleft, Nright, splitindex;
  double *thighsamp = NULL;
  double *maghighsamp = NULL;
  double *magmodelhighsamp = NULL;
  double *phasehighsamp = NULL;
  double *y2 = NULL;
  double *y2left = NULL;
  double *y2right = NULL;
  double *magmodel = NULL;
  double *magresid = NULL;
  double *phasemodel = NULL;
  double dthighsamp;
  double *leftmagvals = NULL;
  double *lefttvals = NULL;
  double *rightmagvals = NULL;
  double *righttvals = NULL;
  double magval;
  double tleft, tright;
  double BS_thigh, BS_tlow;

  double Tc0;

  double **p_amoeba = NULL;
  int *ia_amoeba = NULL;
  double *chi2vals_amoeba = NULL;
  int Ntovary_amoeba = 0;
  int Npar_amoeba = 0;

  double maghigh, maglow;
  
  double delta_amoeba;

  _TransitBisec_passparams tb;

  int amoeba_val, nfunk_amoeba;
  double ftol_amoeba;

  int *Nbin_moment = NULL;
  double *moment_sum1 = NULL;
  double *moment_sum2 = NULL;

  /* Allocate space for the new bisector */
  expand_bisector_data(Ntransit, bisec);

  /* Plan of action:
     1. evaluate the integrated transit model at each of the observed times.
     2. smooth the residuals if desired.
     3. spline the residuals and add to the non-integrated model evaluated at
        higher time sampling.
     4. Spline each side of the high resolution model and determine the bisector
  */

  /* Refit Tc for this transit if asked to */
  if(refit_Tc) {
    double **p_amoeba = NULL;
    int *ia_amoeba = NULL;
    int Ntovary_amoeba = 0;
    tb.t = transitparams;
    tb.b = bisec;
    delta_amoeba = 0.25*transitparams->P * sqrt((pow((1.0 + transitparams->rprstar),2.0) - transitparams->bimpact*transitparams->bimpact)*(1.0 - transitparams->eccen*transitparams->eccen)/(1.0 + transitparams->eccen*sin(transitparams->omega)))/M_PI/transitparams->arstar;
    VARTOOLS_incrementparameters_foramoeba(&Npar_amoeba, &Ntovary_amoeba,
					   &p_amoeba, &ia_amoeba, 1, 
					   transitparams->Tc0,
					   delta_amoeba);
    /* initialize the chi2vals */
    VARTOOLS_amoeba_initializesimplexchi2(Npar_amoeba, Ntovary_amoeba, p_amoeba,
					  &chi2vals_amoeba, 
					  &chisq_shiftmandelagoltransit,
					  (i1-i0+1), &(t[i0]), &(mag[i0]),
					  NULL, (void *) (&tb));
    nfunk_amoeba = 0;
    ftol_amoeba = TRANSITBISEC_AMOEBA_CONVERGENCELIMIT;
    /* Run the amoeba fit */
    amoeba_val = VARTOOLS_amoeba(p_amoeba, chi2vals_amoeba, ia_amoeba, 
				 Npar_amoeba, ftol_amoeba, 
				 &chisq_shiftmandelagoltransit,
				 &nfunk_amoeba, 0, (i1 - i0 + 1), 
				 &(t[i0]), &(mag[i0]), NULL,
				 (void *) (&tb));
    if(!amoeba_val) {
      if(chi2vals_amoeba[0] < chi2vals_amoeba[i]) {
	Tc0 = p_amoeba[0][0];
      } else {
	Tc0 = p_amoeba[1][0];
      }
    }
    else {
      Tc0 = transitparams->Tc0;
    }
    *fitted_Tc = Tc0;

    VARTOOLS_amoeba_cleanup(&Npar_amoeba, &Ntovary_amoeba, &p_amoeba, &ia_amoeba, &chi2vals_amoeba);
  } 
  else {
    Tc0 = transitparams->Tc0;
  }

    
  /* Re-cast the transit parameters into the form used by the
   *mandelagoltransitmodel functions */
  cos_i = transitparams->bimpact*(1. + transitparams->eccen*cos(transitparams->omega))/(1. - transitparams->eccen*transitparams->eccen)/transitparams->arstar;

  if(cos_i > 1.)
    cos_i = 1.;
  else if(cos_i < -1.)
    cos_i = -1.;
  sin_i = sqrt(1.0 - cos_i*cos_i);

  Nintegratesamp = ceil(transitparams->exptime / TRANSITBISEC_TSAMP_INTEGRATE);
  if(Nintegratesamp < 1) Nintegratesamp = 1;

  /* Evalute the integrated transit model at the observed times. */
  if((magmodel = (double *) malloc((i1 - i0 + 1)*sizeof(double))) == NULL ||
     (phasemodel = (double *) malloc((i1 - i0 + 1)*sizeof(double))) == NULL ||
     (magresid = (double *) malloc((i1 - i0 + 1)*sizeof(double))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(i=i0, j = 0; i <= i1; i++, j++) {
    phasemodel[j] = ((t[i] - (Tc0 + Ntransit*transitparams->P))/(transitparams->P));
    phasemodel[j] = phasemodel[j] - floor(phasemodel[j]);
  }
  VARTOOLS_integratemandelagoltransitmodel(
				  transitparams->exptime/transitparams->P,
				  (i1 - i0 + 1), phasemodel, magmodel,
				  transitparams->ldtype, 
				  transitparams->ldcoeffs,
				  sin_i, transitparams->arstar,
				  transitparams->eccen,
				  transitparams->rprstar,
				  transitparams->omega,
				  Nintegratesamp);
  for(i=i0, j = 0; i <= i1; i++, j++) {
    magmodel[j] = -2.5*log(magmodel[j])/M_LN10;
    magresid[j] = mag[i] - magmodel[j];
  }

  /* Low-pass filter the residuals */
  if(bisec->bisec_resid_bin_time > 0) {
    VARTOOLS_medianfilter((i1-i0+1),&(t[i0]),magresid,NULL,
			  bisec->bisec_resid_bin_time,
			  1, 1);
  }

  /* High-pass filter the residuals using a polynomial fit */
  if(bisec->bisec_resid_fit_poly_order > 0) {
    for(j=0; j < i1-i0+1; j++) {
      phasemodel[j] = t[i0+j] - t[i0];
    }
    VARTOOLS_fitpoly((i1-i0+1), phasemodel, magresid, NULL,
		     bisec->bisec_resid_fit_poly_order, 1, NULL, NULL);
  }
  else {
    /* Subtract the mean from the residuals if we are not high-pass filtering */
    magval = VARTOOLS_getmean((i1-i0+1),magresid);
    for(j=0; j < i1-i0+1; j++) {
      magresid[j] -= magval;
    }
  }

  /* Generate the high time sample model */
  Nhighsamp = bisec->bisec_high_samp_factor*(i1 - i0 + 1);
  if((thighsamp = (double *) malloc(Nhighsamp * sizeof(double))) == NULL ||
     (maghighsamp = (double *) malloc(Nhighsamp * sizeof(double))) == NULL ||
     (magmodelhighsamp = (double *) malloc(Nhighsamp * sizeof(double))) == NULL ||
     (phasehighsamp = (double *) malloc(Nhighsamp * sizeof(double))) == NULL ||
     (y2 = (double *) malloc(Nhighsamp * sizeof(double))) == NULL ||
     (leftmagvals = (double *) malloc(Nhighsamp * sizeof(double))) == NULL ||
     (lefttvals = (double *) malloc(Nhighsamp * sizeof(double))) == NULL ||
     (rightmagvals = (double *) malloc(Nhighsamp * sizeof(double))) == NULL ||
     (righttvals = (double *) malloc(Nhighsamp * sizeof(double))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);

  dthighsamp = (t[i1] - t[i0])/(Nhighsamp - 1);
  for(j=0; j < Nhighsamp; j++)
    {
      thighsamp[j] = t[i0] + j*dthighsamp;
      phasehighsamp[j] = ((thighsamp[j] - (Tc0 + Ntransit*transitparams->P))/(transitparams->P));
      phasehighsamp[j] = phasehighsamp[j] - floor(phasehighsamp[j]);
    }

  VARTOOLS_mandelagoltransitmodel(Nhighsamp,
			 phasehighsamp, magmodelhighsamp,
			 transitparams->ldtype, transitparams->ldcoeffs,
			 sin_i, transitparams->arstar,
			 transitparams->eccen, transitparams->rprstar,
			 transitparams->omega);

  /* Spline the residuals onto the resampled model */
  VARTOOLS_spline_monotonic((i1-i0+1), &(t[i0]), magresid, y2);

  magmodelhighsamp[0] = -2.5*log(magmodelhighsamp[0])/M_LN10;

  maghigh = magmodelhighsamp[0];
  maglow = magmodelhighsamp[0];
  splitindex = 0;
  for(j=0; j < Nhighsamp; j++) {    
    if(j > 0) {
      magmodelhighsamp[j] = -2.5*log(magmodelhighsamp[j])/M_LN10;
    }
    if(magmodelhighsamp[j] > maghigh) {
      splitindex = j;
      maghigh = magmodelhighsamp[j];
    }
    if(magmodelhighsamp[j] < maglow)
      maglow = magmodelhighsamp[j];
    maghighsamp[j] = magmodelhighsamp[j] + 
      VARTOOLS_splint_monotonic((i1-i0+1), &(t[i0]), magresid, y2, thighsamp[j]);
  }

  /* Now compute the bisector */
  
  if(!bisec->bisec_usemoments) {
    for(j=0, Nleft = 0; j <= splitindex; j++) {
      if(Nleft > 0 ? (leftmagvals[Nleft-1] < maghighsamp[j]) : 1) {
	lefttvals[Nleft] = thighsamp[j];
	leftmagvals[Nleft] = maghighsamp[j];
	Nleft++;
      }
    }
    for(j=splitindex, Nright = 0; j < Nhighsamp; j++) {
      if((Nright > 0 ? (rightmagvals[(Nhighsamp - 1) - (Nright - 1)] > maghighsamp[j]) : 1)) {
	righttvals[(Nhighsamp - 1) - Nright] = thighsamp[j];
	rightmagvals[(Nhighsamp - 1) - Nright] = maghighsamp[j];
	Nright++;
      }
    }
    
    if(Nleft <= 1 || Nright <= 1)
      {
	if(thighsamp != NULL) free(thighsamp);
	if(maghighsamp != NULL) free(maghighsamp);
	if(magmodelhighsamp != NULL) free(magmodelhighsamp);
	if(phasehighsamp != NULL) free(phasehighsamp);
	if(y2 != NULL) free(y2);
	if(leftmagvals != NULL) free(leftmagvals);
	if(lefttvals != NULL) free(lefttvals);
	if(rightmagvals != NULL) free(rightmagvals);
	if(righttvals != NULL) free(righttvals);
	if(magmodel != NULL) free(magmodel);
	if(phasemodel != NULL) free(phasemodel);
	if(magresid != NULL) free(magresid);
	return 1;
      }
    if((y2left = (double *) malloc(Nleft * sizeof(double))) == NULL ||
       (y2right = (double *) malloc(Nright * sizeof(double))) == NULL)
      VARTOOLS_error(ERR_MEMALLOC);
    
    VARTOOLS_spline_monotonic(Nleft, leftmagvals, lefttvals, y2left);
    VARTOOLS_spline_monotonic(Nright, &(rightmagvals[Nhighsamp - Nright]),
			      &(righttvals[Nhighsamp - Nright]),
			      y2right);
  
    for(j=0, magval = maglow + bisec->bisec_low_val*(maghigh - maglow);
	j < bisec->Ntvals;
	j++, magval += bisec->bisec_step_val * (maghigh - maglow)) {
      tleft = VARTOOLS_splint_monotonic(Nleft, leftmagvals, lefttvals, y2left, magval);
      tright = VARTOOLS_splint_monotonic(Nright, &(rightmagvals[Nhighsamp - Nright]), &(righttvals[Nhighsamp - Nright]), y2right, magval);
      bisec->tvals[bisec->Nbisecs-1][j] = 0.5*(tleft + tright);
    }
    
    /* Get the bisector span */
    magval = maglow + BS_lowval*(maghigh - maglow);
    tleft = VARTOOLS_splint_monotonic(Nleft, leftmagvals, lefttvals, y2left, magval);
    tright = VARTOOLS_splint_monotonic(Nright, &(rightmagvals[Nhighsamp - Nright]), &(righttvals[Nhighsamp - Nright]), y2right, magval);
    
    BS_tlow = 0.5*(tleft + tright);
    
    magval = maglow + BS_highval*(maghigh - maglow);
    tleft = VARTOOLS_splint_monotonic(Nleft, leftmagvals, lefttvals, y2left, magval);
    tright = VARTOOLS_splint_monotonic(Nright, &(rightmagvals[Nhighsamp - Nright]), &(righttvals[Nhighsamp - Nright]), y2right, magval);
    
    BS_thigh = 0.5*(tleft + tright);

    bisec->spanvals[bisec->Nbisecs-1] = (BS_thigh - BS_tlow);
  } else {
    /* Get the bisector using moments */
    if((Nbin_moment = (int *) malloc(bisec->Ntvals * sizeof(int))) == NULL ||
       (moment_sum1 = (double *) malloc(bisec->Ntvals * sizeof(double))) == NULL ||
       (moment_sum2 = (double *) malloc(bisec->Ntvals * sizeof(double))) == NULL)
      {
	VARTOOLS_error(ERR_MEMALLOC);
      }
    for(j=0; j < bisec->Ntvals; j++) {
      Nbin_moment[j] = 0;
      moment_sum1[j] = 0.;
      moment_sum2[j] = 0.;
    }
    for(j=0; j < Nhighsamp; j++) {
      magval = (maghighsamp[j] - (maglow + (bisec->bisec_low_val - 0.5*bisec->bisec_step_val)*(maghigh - maglow)))/(bisec->bisec_step_val*(maghigh - maglow));
      i = floor(magval);
      if(i >= 0 && i < bisec->Ntvals) {
	moment_sum1[i] += thighsamp[j]*magval;
	moment_sum2[i] += magval;
	Nbin_moment[i] += 1;
      }
    }
    for(i=0; i < bisec->Ntvals; i++) {
      if(moment_sum2[i] > 0.0) {
	bisec->tvals[bisec->Nbisecs-1][i] = moment_sum1[i]/moment_sum2[i];
      } else {
	bisec->tvals[bisec->Nbisecs-1][i] = sqrt(-1.);
      }
    }
    magval = (maglow + BS_lowval*(maghigh - maglow) - (maglow + (bisec->bisec_low_val - 0.5*bisec->bisec_step_val)*(maghigh - maglow)))/(bisec->bisec_step_val*(maghigh - maglow));
    i = floor(magval);
    if(i >= 0 && i < bisec->Ntvals) {
      BS_tlow = bisec->tvals[bisec->Nbisecs-1][i];
    } else {
      BS_tlow = sqrt(-1.);
    }
    magval = (maglow + BS_highval*(maghigh - maglow) - (maglow + (bisec->bisec_low_val - 0.5*bisec->bisec_step_val)*(maghigh - maglow)))/(bisec->bisec_step_val*(maghigh - maglow));
    i = floor(magval);
    if(i >= 0 && i < bisec->Ntvals) {
      BS_thigh = bisec->tvals[bisec->Nbisecs-1][i];
    } else {
      BS_thigh = sqrt(-1.);
    }
    if(!isnan(BS_tlow) && !isnan(BS_thigh)) {
      bisec->spanvals[bisec->Nbisecs-1] = (BS_thigh - BS_tlow);
    } else {
      bisec->spanvals[bisec->Nbisecs-1] = sqrt(-1.);
    
    }
  }

  if(thighsamp != NULL) free(thighsamp);
  if(maghighsamp != NULL) free(maghighsamp);
  if(magmodelhighsamp != NULL) free(magmodelhighsamp);
  if(phasehighsamp != NULL) free(phasehighsamp);
  if(y2 != NULL) free(y2);
  if(leftmagvals != NULL) free(leftmagvals);
  if(lefttvals != NULL) free(lefttvals);
  if(rightmagvals != NULL) free(rightmagvals);
  if(righttvals != NULL) free(righttvals);
  if(magmodel != NULL) free(magmodel);
  if(phasemodel != NULL) free(phasemodel);
  if(magresid != NULL) free(magresid);
  if(y2left != NULL) free(y2left);
  if(y2right != NULL) free(y2right);
  if(Nbin_moment != NULL) free(Nbin_moment);
  if(moment_sum1 != NULL) free(moment_sum1);
  if(moment_sum2 != NULL) free(moment_sum2);

  if(!isnan(bisec->spanvals[bisec->Nbisecs-1]))
    return 0;
  else
    return 1;
}

void transitbisec_getMA_resid_and_model(int N, double *t, double *mag,
					TransitParams *tpar,
					double *outresid,
					double *outmodel)
{
  int i;
  double *phasemodel = NULL;
  double sin_i, cos_i;
  int Nintegratesamp;

  cos_i = tpar->bimpact*(1. + tpar->eccen*cos(tpar->omega))/(1. - tpar->eccen*tpar->eccen)/tpar->arstar;

  if(cos_i > 1.)
    cos_i = 1.;
  else if(cos_i < -1.)
    cos_i = -1.;
  sin_i = sqrt(1. - cos_i*cos_i);

  Nintegratesamp = ceil(tpar->exptime / TRANSITBISEC_TSAMP_INTEGRATE);
  if(Nintegratesamp < 1) Nintegratesamp = 1;
 
  if((phasemodel = (double *) malloc(N*sizeof(double))) == NULL)
    VARTOOLS_error(ERR_MEMALLOC);
  for(i=0; i < N; i++) {
    phasemodel[i] = (t[i] - tpar->Tc0)/tpar->P;
    phasemodel[i] = phasemodel[i] - floor(phasemodel[i]);
  }
  VARTOOLS_integratemandelagoltransitmodel(
				  tpar->exptime/tpar->P,
				  N, phasemodel, outmodel,
				  tpar->ldtype, 
				  tpar->ldcoeffs,
				  sin_i, tpar->arstar,
				  tpar->eccen,
				  tpar->rprstar,
				  tpar->omega,
				  Nintegratesamp);
  for(i=0; i < N; i++) {
    outmodel[i] = -2.5*log(outmodel[i])/M_LN10;
    outresid[i] = mag[i] - outmodel[i];
  }
  free(phasemodel);
}

/* Create a simulated lc with shuffled residuals */
void transitbisec_makemodel_residshuffle(int N, double *magresid, 
					 double *magmodel, double *magsim,
					 double *t,
					 TransitParams *tpar)
{
  int i, j, k;
  double ph, tdurph;

  tdurph = 0.75* sqrt((pow((1.0 + tpar->rprstar),2.0) - tpar->bimpact*tpar->bimpact)*(1.0 - tpar->eccen*tpar->eccen)/(1.0 + tpar->eccen*sin(tpar->omega)))/M_PI/tpar->arstar;


  k = 0;
  j = rand() % N;
  for(k=0, i = j - 1; k < N; k++) {
    do {
      i++;
      if(i >= N) i = 0;
      ph = (t[i] - tpar->Tc0)/tpar->P;
      ph = ph - floor(ph);
    } while (ph < tdurph || ph > (1.0 - tdurph));
    magsim[k] = magmodel[k] + magresid[i];
  }
}


void free_Bisector_Data(Bisector_Data *bdat) {
  int i;
  if(bdat->tvals != NULL) {
    for(i=0; i < bdat->sizebisec1; i++) {
      if(bdat->tvals[i] != NULL)
	free(bdat->tvals[i]);
    }
    free(bdat->tvals);
  }
  if(bdat->spanvals != NULL)
    free(bdat->spanvals);
  if(bdat->transitnums != NULL)
    free(bdat->transitnums);
}
    

int do_transitbisec(int N, double *t, double *mag, 
		    Bisector_Data *bdat,
		    TransitParams *tpar,
		    char *bisecoutname,
		    char *bisecspanoutname,
		    int outputbisec,
		    int outputbisecspan,
		    double max_delta_t,
		    double min_oot_frac,
		    double BS_highval,
		    double BS_lowval,
		    double *out_BSrms,
		    double *out_BSmed,
		    int *out_Ntransits,
		    int geterror,
		    int Nerrorsims,
		    double *out_BSrms_sim,
		    double *out_BSmed_sim,
		    int refit_Tc)
/* This is the main function for the -transitbisec command. It
   determines which transits are present and completely observed in an
   LC, and then calls get_transit_bisector on each of them to
   calculate their bisectors. If requested it also injects transits
   into other parts of the LC, and conducts a prayer-bead analysis on
   the residuals, to determine the uncertainties on the bisector span
   measurements.  */
{
  int i, j, Nt, skiptransit, i0, i1, isim;
  double t0, t1, t_trstart, t_trstop, tdur, magval, Tcval, Tcfit;
  int Ntransit0, Ntransit1;

  FILE *bisec_outfile, *span_outfile;

  double *magmodel = NULL, *magresid = NULL, *magsim = NULL;

  double *BSsim = NULL;

  int havemodelresid = 0, Nsim = 0, sizeBSsim = 0;

  Bisector_Data bdatsim;

  /* Initialize the Bisector_Data structure for this LC */
  bdat->Nbisecs = 0;
  bdat->Ntvals = 0;

  /* Initialize the Bisector_Data structure for the error simulations */
  if(geterror) {
    bdatsim.Nbisecs = 0;
    bdatsim.tvals = NULL;
    bdatsim.Ntvals = 0;
    bdatsim.sizebisec1 = 0;
    bdatsim.sizebisec2 = 0;
    bdatsim.bisec_high_val = bdat->bisec_high_val;
    bdatsim.bisec_low_val = bdat->bisec_low_val;
    bdatsim.bisec_step_val = bdat->bisec_step_val;
    bdatsim.spanvals = NULL;
    bdatsim.transitnums = NULL;
    bdatsim.bisec_resid_bin_time = bdat->bisec_resid_bin_time;
    bdatsim.bisec_resid_fit_poly_order = bdat->bisec_resid_fit_poly_order;
  }

  /* First find transits in the LC that have sufficient coverage to 
     allow bisector measurements */
  
  t0 = t[0];
  t1 = t[N-1];
  Ntransit0 = ceil((t0 - tpar->Tc0)/tpar->P); /* Number of first transit center
						 after first observation */
  Ntransit1 = floor((t1 - tpar->Tc0)/tpar->P); /* Number of last transit center
						  before last observation */
  if(Ntransit1 < Ntransit0 || N < 3) {
    /* No full transit is observed */
    *out_BSrms = sqrt(-1.);
    *out_BSmed = sqrt(-1.);
    *out_Ntransits = 0;
    if(geterror) {
      *out_BSrms_sim = sqrt(-1.);
      *out_BSmed_sim = sqrt(-1.);
    }
    if(magmodel != NULL)
      free(magmodel);
    if(magresid != NULL)
      free(magresid);
    if(magsim != NULL)
      free(magsim);
    if(BSsim != NULL)
      free(BSsim);
    if(geterror)
      free_Bisector_Data(&bdatsim);
    return 1;
  }

  /* The transit duration */
  tdur = tpar->P * sqrt((pow((1.0 + tpar->rprstar),2.0) - tpar->bimpact*tpar->bimpact)*(1.0 - tpar->eccen*tpar->eccen)/(1.0 + tpar->eccen*sin(tpar->omega)))/M_PI/tpar->arstar;

  /* For each transit check if it has sufficient coverage to allow a 
     bisector measurement, and if so, compute it. If not, skip to the
     next transit. Note that we assume the times are sorted. 
     We also require a quarter of a transit duration to be included before
     and after transit */
  *out_Ntransits = 0;
  j = 0;
  for(Nt = Ntransit0; Nt <= Ntransit1; Nt++) {
    t_trstart = tpar->Tc0 + Nt*tpar->P - (0.5 + min_oot_frac)*tdur;
    t_trstop = tpar->Tc0 + Nt*tpar->P + (0.5 + min_oot_frac)*tdur;
    while(j < N ? t[j] < t_trstart : 0) j++;
    if(j >= N) break;
    i0 = j;
    skiptransit = 0;
    if(t[j] - t_trstart > max_delta_t) {
      continue;
    }
    j++;
    while(j < N ? t[j] < t_trstop : 0) {
      if(t[j] - t[j-1] > max_delta_t)
	{
	  skiptransit = 1;
	  break;
	}
      j++;
    }
    if(skiptransit) 
      continue;
    if(fabs(t_trstop - t[j-1]) > max_delta_t)
      continue;
    i1 = j-1;
    if(i1 - i0 < 2)
      continue;

    /* Prepare the data for calculating the error if we're going to
       do that */
    if(geterror && !havemodelresid) {
      if((magmodel = (double *) malloc(N * sizeof(double))) == NULL ||
	 (magresid = (double *) malloc(N * sizeof(double))) == NULL ||
	 (magsim = (double *) malloc(N * sizeof(double))) == NULL)
	VARTOOLS_error(ERR_MEMALLOC);
      transitbisec_getMA_resid_and_model(N, t, mag, tpar, magresid, magmodel);
    }

    /* This transit has good coverage, compute the bisector */
    if(get_transit_bisector(N, t, mag, i0, i1,
			    BS_highval, BS_lowval, 
			    Nt, tpar, bdat, refit_Tc, &Tcfit)) {
      /* There was an error in getting the bisector, skip this transit */
      bdat->Nbisecs -= 1;
      continue;
    }
    (*out_Ntransits) += 1;

    /* Generate simulated BS values by shuffling the residuals and computing
       the BS. These are used to estimate the error on the BS */
    if(geterror) {
      for(isim = 0; isim < Nerrorsims; isim++) {
	transitbisec_makemodel_residshuffle(N, magresid, magmodel, magsim,
					    t, tpar);
	if(!(get_transit_bisector(N, t, magsim, i0, i1, BS_highval,
				  BS_lowval, Nt, tpar, &bdatsim, 0, NULL))) {
	  if(!sizeBSsim) {
	    sizeBSsim = 256;
	    if((BSsim = (double *) malloc(sizeBSsim * sizeof(double))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	  }
	  else if(Nsim >= sizeBSsim) {
	    sizeBSsim *= 2;
	    if((BSsim = (double *) realloc(BSsim, sizeBSsim*sizeof(double))) == NULL)
	      VARTOOLS_error(ERR_MEMALLOC);
	  }
	  BSsim[Nsim] = bdatsim.spanvals[0];
	  Nsim++;
	}
	bdatsim.Nbisecs = 0;
      }
    }
  }

  if((*out_Ntransits) < 1) {
    /* No full transit is observed */
    *out_BSrms = sqrt(-1.);
    *out_BSmed = sqrt(-1.);
    *out_Ntransits = 0;
    if(geterror) {
      *out_BSrms_sim = sqrt(-1.);
      *out_BSmed_sim = sqrt(-1.);
    }
    if(magmodel != NULL)
      free(magmodel);
    if(magresid != NULL)
      free(magresid);
    if(magsim != NULL)
      free(magsim);
    if(BSsim != NULL)
      free(BSsim);
    if(geterror)
      free_Bisector_Data(&bdatsim);
    return 1;
  }
  
  /* Output the bisectors and spans if asked */
  if(outputbisec) {
    if((bisec_outfile = fopen(bisecoutname,"w")) == NULL) {
      fprintf(stderr,"Error: Cannot write to %s\n", bisecoutname);
      exit(1);
    }
    fprintf(bisec_outfile,"# Depth");
    for(i=0; i < bdat->Nbisecs; i++)
      fprintf(bisec_outfile," Bisector_Time_%d", bdat->transitnums[i]);
    fprintf(bisec_outfile,"\n");
    for(j=0, magval = bdat->bisec_low_val; j < bdat->Ntvals; j++, magval += bdat->bisec_step_val) {
      fprintf(bisec_outfile,"%.17g", magval);
      for(i=0; i < bdat->Nbisecs; i++) {
	fprintf(bisec_outfile," %.17g",bdat->tvals[i][j]);
      }
      fprintf(bisec_outfile,"\n");
    }
    fclose(bisec_outfile);
  }
  
  if(outputbisecspan) {
    if((span_outfile = fopen(bisecspanoutname,"w")) == NULL) {
      fprintf(stderr,"Error: Cannot write to %s\n", bisecspanoutname);
      exit(1);
    }
    fprintf(span_outfile,"#Ntr Tc Span\n");
    for(i = 0; i < bdat->Nbisecs; i++) {
      Nt = bdat->transitnums[i];
      Tcval = tpar->Tc0 + Nt*tpar->P;
      fprintf(span_outfile,"%d %.17g %.17g\n", Nt, Tcval, bdat->spanvals[i]);
    }
    fclose(span_outfile);
  }

  /* Compute the RMS and median of the bisector spans */
  if(bdat->Nbisecs > 1) {
    *out_BSrms = VARTOOLS_stddev(bdat->Nbisecs, bdat->spanvals);
  } else {
    *out_BSrms = sqrt(-1.);
  }
  *out_BSmed = VARTOOLS_median(bdat->Nbisecs, bdat->spanvals);

  if(geterror) {
    *out_BSrms_sim = VARTOOLS_stddev(Nsim, BSsim);
    *out_BSmed_sim = VARTOOLS_median(Nsim, BSsim);
  }

  if(magmodel != NULL)
    free(magmodel);
  if(magresid != NULL)
    free(magresid);
  if(magsim != NULL)
    free(magsim);
  if(BSsim != NULL)
    free(BSsim);
  if(geterror)
    free_Bisector_Data(&bdatsim);
  return 0;
    
}


void transitbisec_Initialize(char *commandname, int *RequireReadAll,
			     int *RequireSortLC, int *RequireDistinctTimes,
			     size_t *sizeuserdata)
/* This function defines five things which vartools
   will need to know about the command 

   Every library with the name $LIBNAME.so must contain a function with
   the name $LIBNAME_Initialize to set the five input variables.
*/
{
  /* Set the string used to call the command */
  sprintf(commandname,"-transitbisec");

  /* This is 1 if the command requires all light curves to be read
     at once, otherwise set it to 0. */
  *RequireReadAll = 0;

  /* This is 1 if the command requires input light curves to be sorted by
     time, otherwise set it to 0. */
  *RequireSortLC = 1;
  
  /* This is 1 if the command requires the times of observation in the
     input light curve to be unique */
  *RequireDistinctTimes = 1;

  /* You should define a structure to store the data needed for this
     command (in this example, a vector to hold the values to add to
     the magnitudes of each lc), below you would replace "_TransitBisec"
     with the type name of your structure.

     See TransitBisec.h
*/
  *sizeuserdata = sizeof(_TransitBisec);
}

int transitbisec_ParseCL(ProgramData *p, Command *c,
			   void *userdata, int *iret, char **argv, int argc)
/* Every library with the name $LIBNAME.so must contain a function
   with the name $LIBNAME_ParseCL to check that the user has issued a
   valid syntax, and set-up the command-data accordingly.

   Note that any vectors/arrays used by the command should also be registered
   by this function so that memory will be allocated as needed.

Expected syntax for this example:
   -transitbisec ....

  p = structure containing general program data. In most cases this 
      structure is only passed to other functions (like RegisterDataVector)
      and not used directly by the ParseCL function.

  c = pointer to a generic container for this command. This structure is
      needed by other functions that might be called here (RegisterDataVector)
      but should not be used directly by the ParseCL function.

  userdata = pointer to the structure storing data for this
      command. This should be cast to the appropriate type and the
      data there-in should be modified as needed by the ParseCL
      function. In this example, userdata corresponds to the _TransitBisec
      structure which stores the values to be added to the light
      curves.

  iret = command line argument index. This is initially zero, and on
      return should be set to the index of the last argument parsed by
      this command. Note that index 0 points to the name of the
      command, which does not need to be checked (this function is
      only called if the user issues the command). 

  argv = array of command line arguments. argv[0] is the name of the
      command ("-transitbisec" in this case).

  argc = Number of command line arguments in the argv array (including
      0). It is the user's responsibility to check that i < argc
      before attempting to parse argv[i], failure to do so may lead to
      segmentation violations.
 */
{
  int i = 0;
  int check;
  
  /* Cast the input userdata pointer to the data structure type used
     for this command */
  _TransitBisec *TransitBisec;
  TransitBisec = (_TransitBisec *) userdata;

  /* We'll use i rather than iret to index argv, this is just so we
     don't have to constantly dereference the pointer iret */
  i = *iret;

   /* Our procedure here is to step through each term in the expected
     command-line syntax. When possible we use the 
     VARTOOLS_ParseFixSpecFixcolumn function to parse terms like

     <"fix" value | "list" | "fixcolumn" <colname | colnum> | "expr" expr>

     For other terms we first check to make sure that we haven't
     exceeded the number of command-line arguments (checks like if(i
     >= argc) return 1;), and then use strcmp to test whether the
     command-line term in argv[i] is equal to the expected string (note
     that strcmp returns 0 if the strings are equal). You can use the
     atof() or atoi() functions to parse a string into a double or
     integer respectively. If something on the command-line is not as
     expected return 1 to tell VARTOOLS that there was an error in
     parsing the command. */

  /* Parse the options for each of the transit parameters */
  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"Tc0")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					 VARTOOLS_TYPE_DOUBLE, 
					 (void *) (&TransitBisec->Tc0),
					 0,
					 0,
					 "Tc0");
  if(check) {*iret = i; return 1;}

  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"P")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					 VARTOOLS_TYPE_DOUBLE, 
					 (void *) (&TransitBisec->P),
					 0,
					 0,
					 "Period");
  if(check) {*iret = i; return 1;}

  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"rprstar")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					 VARTOOLS_TYPE_DOUBLE, 
					 (void *) (&TransitBisec->rprstar),
					 0,
					 0,
					 "Rp_over_Rstar");
  if(check) {*iret = i; return 1;}

  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"bimpact")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					 VARTOOLS_TYPE_DOUBLE, 
					 (void *) (&TransitBisec->bimpact),
					 0,
					 0,
					 "bimpact");
  if(check) {*iret = i; return 1;}

  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"arstar")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					 VARTOOLS_TYPE_DOUBLE, 
					 (void *) (&TransitBisec->arstar),
					 0,
					 0,
					 "a_over_Rstar");
  if(check) {*iret = i; return 1;}

  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"eccen")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					 VARTOOLS_TYPE_DOUBLE, 
					 (void *) (&TransitBisec->eccen),
					 0,
					 0,
					 "eccen");
  if(check) {*iret = i; return 1;}

  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"omega")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					 VARTOOLS_TYPE_DOUBLE, 
					 (void *) (&TransitBisec->omega),
					 0,
					 0,
					 "omega");
  if(check) {*iret = i; return 1;}

  if(i >= argc)
    return 1;
  if(strcmp(argv[i],"exptime")) {
    *iret = i; return 1;
  } else {
    i++;
  }
  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					 VARTOOLS_TYPE_DOUBLE, 
					 (void *) (&TransitBisec->exptime),
					 0,
					 0,
					 "exptime");
  if(check) {*iret = i; return 1;}

  if(i >= argc)
    return 1;
  if(!strcmp(argv[i],"quad")) {
    TransitBisec->ldtype = 0;
    TransitBisec->Nld = 2;
  }
  else if(!strcmp(argv[i],"nonlin")) {
    TransitBisec->ldtype = 1;
    TransitBisec->Nld = 4;
  }
  else {
    *iret = i;
    return 1;
  }
  i++;

  check = VARTOOLS_ParseFixSpecFixcolumn(p, c, &i, argv, argc, 1,
					 VARTOOLS_TYPE_DOUBLE,
					 (void *) (&TransitBisec->ldcoeffs),
					 TransitBisec->Nld, 0,
					 "LD_coeff");
  if(check) {*iret = i; return 1;}
  
  if(i < argc) {
    if(!strcmp(argv[i],"bisec_low_value")) {
      i++;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      TransitBisec->bisec_low_val = atof(argv[i]);
      i++;
    } else {
      TransitBisec->bisec_low_val = TRANSITBISEC_DEFAULT_BISEC_LOW_VAL;
    }
  } else {
    TransitBisec->bisec_low_val = TRANSITBISEC_DEFAULT_BISEC_LOW_VAL;
  }
  
  if(i < argc) {
    if(!strcmp(argv[i],"bisec_high_value")) {
      i++;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      TransitBisec->bisec_high_val = atof(argv[i]);
      i++;
    } else {
      TransitBisec->bisec_high_val = TRANSITBISEC_DEFAULT_BISEC_HIGH_VAL;
    }
  } else {
    TransitBisec->bisec_high_val = TRANSITBISEC_DEFAULT_BISEC_HIGH_VAL;
  }

  if(i < argc) {
    if(!strcmp(argv[i],"bisec_stepsize")) {
      i++;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      TransitBisec->bisec_step_val = atof(argv[i]);
      i++;
    } else {
      TransitBisec->bisec_step_val = TRANSITBISEC_DEFAULT_BISEC_STEP;
    }
  } else {
    TransitBisec->bisec_step_val = TRANSITBISEC_DEFAULT_BISEC_STEP;
  }

  if(i < argc) {
    if(!strcmp(argv[i],"span_low_value")) {
      i++;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      TransitBisec->BS_lowval = atof(argv[i]);
      i++;
    } else {
      TransitBisec->BS_lowval = TRANSITBISEC_DEFAULT_SPAN_LOW_VAL;
    }
  } else {
    TransitBisec->BS_lowval = TRANSITBISEC_DEFAULT_SPAN_LOW_VAL;
  }

  if(i < argc) {
    if(!strcmp(argv[i],"span_high_value")) {
      i++;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      TransitBisec->BS_highval = atof(argv[i]);
      i++;
    } else {
      TransitBisec->BS_highval = TRANSITBISEC_DEFAULT_SPAN_HIGH_VAL;
    }
  } else {
    TransitBisec->BS_highval = TRANSITBISEC_DEFAULT_SPAN_HIGH_VAL;
  }

  if(i < argc) {
    if(!strcmp(argv[i],"max_delta_t")) {
      i++;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      TransitBisec->max_delta_t = atof(argv[i]);
      i++;
    } else {
      TransitBisec->max_delta_t = TRANSITBISEC_DEFAULT_MAX_DELTA_T;
    }
  } else {
    TransitBisec->max_delta_t = TRANSITBISEC_DEFAULT_MAX_DELTA_T;
  }

  if(i < argc) {
    if(!strcmp(argv[i],"min_oot_frac")) {
      i++;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      TransitBisec->min_oot_frac = atof(argv[i]);
      i++;
    } else {
      TransitBisec->min_oot_frac = TRANSITBISEC_DEFAULT_MIN_OOT_FRAC;
    }
  } else {
    TransitBisec->min_oot_frac = TRANSITBISEC_DEFAULT_MIN_OOT_FRAC;
  }

  TransitBisec->bisec_usemoments = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"bisec_usemoments")) {
      i++;
      TransitBisec->bisec_usemoments = 1;
    }
  }

  if(i < argc) {
    if(!strcmp(argv[i],"bisec_high_samp_factor")) {
      i++;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      TransitBisec->bisec_high_samp_factor = atoi(argv[i]);
      i++;
    } else {
      TransitBisec->bisec_high_samp_factor = TRANSITBISEC_DEFAULT_BISEC_HIGH_SAMP_FACTOR;
    }
  } else {
    TransitBisec->bisec_high_samp_factor = TRANSITBISEC_DEFAULT_BISEC_HIGH_SAMP_FACTOR;
  }
    

  if(i < argc) {
    if(!strcmp(argv[i],"resid_bin_time_lowpass")) {
      i++;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      TransitBisec->bisec_resid_bin_time = atof(argv[i]);
      i++;
    } else {
      TransitBisec->bisec_resid_bin_time = TRANSITBISEC_DEFAULT_RESID_BIN_TIME;
    }
  } else {
    TransitBisec->bisec_resid_bin_time = TRANSITBISEC_DEFAULT_RESID_BIN_TIME;
  }

  if(i < argc) {
    if(!strcmp(argv[i],"resid_fit_poly_order")) {
      i++;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      TransitBisec->bisec_resid_fit_poly_order = atoi(argv[i]);
      i++;
    } else {
      TransitBisec->bisec_resid_fit_poly_order = TRANSITBISEC_DEFAULT_RESID_FIT_POLY_ORDER;
    }
  } else {
    TransitBisec->bisec_resid_fit_poly_order = TRANSITBISEC_DEFAULT_RESID_FIT_POLY_ORDER;
  }
  
  TransitBisec->refit_Tc = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"refit_Tc")) {
      i++;
      TransitBisec->refit_Tc = 1;
    }
  }

  TransitBisec->geterror = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"geterror")) {
      i++;
      TransitBisec->geterror = 1;
      TransitBisec->Nerrorsims = TRANSITBISEC_DEFAULT_NERRORSIMS;
      if(i < argc) {
	if(!strcmp(argv[i],"Nsims")) {
	  i++;
	  if(i >= argc) {
	    *iret = i;
	    return 1;
	  }
	  TransitBisec->Nerrorsims = atoi(argv[i]);
	  i++;
	}
      }
    }
  }

  TransitBisec->outputbisec = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"obisec")) {
      i++;
      TransitBisec->outputbisec = 1;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      sprintf(TransitBisec->bisecoutdir,"%s",argv[i]);
      i++;
      TransitBisec->outputbisec_useformat = 0;
      if(i < argc) {
	if(!strcmp(argv[i],"nameformat")) {
	  i++;
	  TransitBisec->outputbisec_useformat = 1;
	  if(i >= argc) {
	    *iret = 1;
	    return 1;
	  }
	  sprintf(TransitBisec->outputbisec_format,"%s",argv[i]);
	  i++;
	}
      }
    }
  }
  
  TransitBisec->outputbisecspan = 0;
  if(i < argc) {
    if(!strcmp(argv[i],"ospan")) {
      i++;
      TransitBisec->outputbisecspan = 1;
      if(i >= argc) {
	*iret = i;
	return 1;
      }
      sprintf(TransitBisec->bisecspanoutdir,"%s",argv[i]);
      i++;
      TransitBisec->outputspan_useformat = 0;
      if(i < argc) {
	if(!strcmp(argv[i],"nameformat")) {
	  i++;
	  TransitBisec->outputspan_useformat = 1;
	  if(i >= argc) {
	    *iret = 1;
	    return 1;
	  }
	  sprintf(TransitBisec->outputspan_format,"%s",argv[i]);
	  i++;
	}
      }
    }
  }

  /* Register the transitparams and bisecdata structs. These will be used
     to transit and bisector related parameters and data */

  VARTOOLS_RegisterDataVector(p, c, (void *) (&TransitBisec->bisecdata),
			      VARTOOLS_TYPE_USERDEF, 0,
			      VARTOOLS_SOURCE_COMPUTED,
			      0, NULL, sizeof(Bisector_Data),
			      initialize_Bisector_Data,
			      (void *) TransitBisec);

  VARTOOLS_RegisterDataVector(p, c, (void *) (&TransitBisec->transitparams),
			      VARTOOLS_TYPE_USERDEF, 0,
			      VARTOOLS_SOURCE_COMPUTED,
			      0, NULL, sizeof(TransitParams),
			      initialize_TransitParams,
			      (void *) TransitBisec);

  /* Register the vectors for storing the output data */
  VARTOOLS_RegisterDataVector(p, c, (void *) (&TransitBisec->out_BSrms),
			      VARTOOLS_TYPE_DOUBLE, 0,
			      VARTOOLS_SOURCE_COMPUTED,
			      1, "BS_RMS");
  VARTOOLS_RegisterDataVector(p, c, (void *) (&TransitBisec->out_BSmed),
			      VARTOOLS_TYPE_DOUBLE, 0,
			      VARTOOLS_SOURCE_COMPUTED,
			      1, "BS_median");
  VARTOOLS_RegisterDataVector(p, c, (void *) (&TransitBisec->out_Ntransits),
			      VARTOOLS_TYPE_INT, 0,
			      VARTOOLS_SOURCE_COMPUTED,
			      1, "Ntransits");

  if(TransitBisec->geterror) {
    VARTOOLS_RegisterDataVector(p, c, (void *) (&TransitBisec->out_BSrms_sim),
				VARTOOLS_TYPE_DOUBLE, 0,
				VARTOOLS_SOURCE_COMPUTED,
				1, "BS_RMS_expected");
    
    VARTOOLS_RegisterDataVector(p, c, (void *) (&TransitBisec->out_BSmed_sim),
				VARTOOLS_TYPE_DOUBLE, 0,
				VARTOOLS_SOURCE_COMPUTED,
				1, "BS_median_expected");
  }



  /* If we got this far, the command line was parsed successfully. Set the
     command line index to the new value, and return 0 to indicate success */
  *iret = i;
  return 0;

}

void transitbisec_ShowSyntax(FILE *outfile)
/* Write out the expected syntax for this command to the file outfile */
{
  fprintf(outfile,
	  "-transitbisec\n"
	  "    <\"Tc0\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>>\n"
          "    <\"P\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>>\n"
	  "    <\"rprstar\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>>\n"
	  "    <\"bimpact\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>>\n"
	  "    <\"arstar\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>>\n"
	  "    <\"eccen\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>>\n"
	  "    <\"omega\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>>\n"
	  "    <\"exptime\" <\"fix\" value | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>>\n"
	  "    <<\"quad\" | \"nonlin\">\n"
          "         <\"fix\" v1 v2 [v3 v4] | \"list\" | \"fixcolumn\" <colname | colnum> | \"expr\" expr>>\n"
          "    [\"bisec_low_value\" val] [\"bisec_high_value\" val] [\"bisec_stepsize\" val]\n"
	  "    [\"span_low_value\" val] [\"span_high_value\" val]\n"
	  "    [\"max_delta_t\" val] [\"min_oot_frac\" val]\n"
	  "    [\"bisec_usemoments\"] [\"bisec_high_samp_factor\" val]\n"
	  "    [\"resid_bin_time_lowpass\" val] [\"resid_fit_poly_order\" val]\n"
	  "    [\"refit_Tc\"] [\"geterror\" [\"Nsims\" val]]\n"
          "    [\"obisec\" directory [\"nameformat\" format]]\n"
          "    [\"ospan\" directory [\"nameformat\" format]]\n"
	  );
}


void transitbisec_ShowHelp(FILE *outfile)
/* Output the help information for this command */
{
  fprintf(outfile,
"Compute the bisectors and bisector spans of transits in a light curve. The user specifies the source for parameters used to define the transits. The parameters include:\n"
"\t\"Tc0\" - time of transit center.\n"
"\t\"P\" - transit period.\n"
"\t\"rprstar\" - radius of planet / radius of star.\n"
"\t\"bimpact\" - impact parameter of the transit in units of the stellar radius.\n"
"\t\"arstar\" - semimajor axis / radius of star.\n"
"\t\"eccen\" - orbital eccentricity.\n"
"\t\"omega\" - argument of periastron in degrees.\n"
"\t\"exptime\" - exposure time of observations in days (transit model is\n"
"\t\tintegrated over an exposure.\n"
"\t\"quad\" or \"nonlin\" - The type of limb darkening law. Quadratic or nonlinear\n"
"\t\tlimb darkening is supported. If \"quad\" is used, there are 2 coefficients\n"
"\t\tif \"nonlin\" is used, then there are 4 coefficients.\n"
	  "For each parameter the user specifies the source for the value to use. This can either be the keyword \"fix\" followed by the value to use for all light curves, \"list\" to read the parameter from the light curve list (use the \"column\" keyword folloed by a number to indicate the column number from the list to use, otherwise the next column in the list will be assumed), \"fixcolumn\" followed by the name or number of an output column from a previously executed command, or \"expr\" followed by an analytic expression. Following this are a number of optional parameters used to define how the bisectors are calculated and output. \"bisec_low_value\" is the fraction of the transit depth from the out-of-transit level at which to start computing the bisector, it defaults to 0.05. \"bisec_high_value\" is the fraction of the transit depth from the out-of-transit level at which to stop computing the bisector, it defaults to 0.95. \"bisec_stepsize\" is the step size in fractions of the transit depth at which to determine the bisector (default is 0.005). \"span_low_value\" is the fraction of the transit depth from the out-of-transit level to use for one end of the bisector span, the default is 0.1. \"span_high_value\" is the fraction of the transit depth from the out-of-transit level to use for the other end of the bisector span, the default is 0.9. One can control the maximum allowed gap in time in the sampling across a transit for the bisector of a transit to be computed using the \"max_delta_t\" keyword, the default is 0.03 days = 43 minutes. The minimum amount of out-of-transit before the start of ingress and after the end of egress that must be present for the bisector of a transit to be computed is controlled with the \"min_oot_frac\" keyword. This is expressed in units of the transit duration. The default value is 0.25. By default this command will calculate the bisectors by generating a smooth function of time vs flux and interpolating at the desired flux levels. Alternatively the bisector can be computed by binning the data into flux level bins and computing the temporal moment in each bin, this will be done if \"bisec_usemoments\" keyword is given. The factor by which the model will be resampled can be adjusted with the \"bisec_high_samp_factor\" keyword, the default is 10. The time-scale at which residuals from a Mandel & Agol transit model are smoothed using a low-pass moving mean filter can be adjusted with the \"resid_bin_time_lowpass\" keyword. The default is 0.05 days. Set either value to 0 or a negative number to not apply the filter. The order of a polynomial to fit and subtract from the residuals can be adjusted with the \"resid_fit_poly_order\" keyword. The default is 1. If the keyword \"refit_Tc\" is specified, then the transit centers for each transit will be individually determined. If the keyword \"geterror\" is specified then bisector span uncertainties will be estimated by injecting simulated transits into the light curve (one can specify the number of simulations to run using the optional \"Nsims\" keyword). Finally one may output the bisectors and/or the bisector spans by giving the \"obisec\" and/or \"ospan\" keywords, in each case followed by the directory to output the results to. By default bisectors will be output to the file OBISECDIR/LCNAME.bisec and bisector spans will be output to the file OSPANDIR/LCNAME.span, where LCNAME is the basename of the light curve, and OBISECDIR or OSPANDIR is the specified directory. One can optionally specify a format rule for the output name by giving the \"nameformat\" keyword followed by the formatstring. This uses the same syntax as the \"-o\" command, see \"vartools -help -o\" for an explanation.\n\n");
}

void transitbisec_RunCommand(ProgramData *p, void *userdata, int lc_name_num, int lc_num)
/* This function runs the command on a light curve.  

   p = structure containing various general program data (in
       particular the light curves are contained in this structure).

   userdata = pointer to the structure containing the command specific
       data (including control parameters and vectors to store output
       results).

   lc_name_num = this is the index to use to access the light curve name.

   lc_num = this is the index to use to access the light curve and
   data from Registered vectors (e.g. the addval parameter for this
   command).

*/
{
  int NJD;
  double *t, *mag;
  char bisecoutname[MAXLEN];
  char bisecspanoutname[MAXLEN];
  _TransitBisec *TransitBisec;

  /* Translate pointers and structures input by vartools into
     easier to use forms */
  TransitBisec = (_TransitBisec *)userdata;


  /* The number of points in the light curve, and vectors storing the
     time, magnitudes, and errors are accessed as illustrated
     below. Note that t, mag and err are vectors, t[0] is the time of
     the first observation in the light curve, t[NJD-1] is the time of
     the last observation, etc. */
  NJD=p->NJD[lc_num];
  t = p->t[lc_num];
  mag = p->mag[lc_num];

  /* Get the names of the files to output the bisector and spans to */
  if(TransitBisec->outputbisec) {
    if(!TransitBisec->outputbisec_useformat) {
      VARTOOLS_GetOutputFilename(bisecoutname, p->lcnames[lc_name_num], 
				 TransitBisec->bisecoutdir, "bisec", 
				 NULL, lc_name_num);
    }
    else {
      VARTOOLS_GetOutputFilename(bisecoutname, p->lcnames[lc_name_num], 
				 TransitBisec->bisecoutdir, "bisec", 
				 TransitBisec->outputbisec_format, 
				 lc_name_num);
    }
  }

  if(TransitBisec->outputbisecspan) {
    if(!TransitBisec->outputspan_useformat) {
      VARTOOLS_GetOutputFilename(bisecspanoutname, p->lcnames[lc_name_num], 
				 TransitBisec->bisecspanoutdir, "span", 
				 NULL, lc_name_num);
    }
    else {
      VARTOOLS_GetOutputFilename(bisecspanoutname, p->lcnames[lc_name_num], 
				 TransitBisec->bisecspanoutdir, "span", 
				 TransitBisec->outputspan_format, 
				 lc_name_num);
    }
  }

  /* Set the transit parameters to the appropriate values */
  TransitBisec->transitparams[lc_num].Tc0 = TransitBisec->Tc0[lc_num];
  TransitBisec->transitparams[lc_num].P = TransitBisec->P[lc_num];
  TransitBisec->transitparams[lc_num].rprstar = TransitBisec->rprstar[lc_num];
  TransitBisec->transitparams[lc_num].bimpact = TransitBisec->bimpact[lc_num];
  TransitBisec->transitparams[lc_num].arstar = TransitBisec->arstar[lc_num];
  TransitBisec->transitparams[lc_num].eccen = TransitBisec->eccen[lc_num];
  TransitBisec->transitparams[lc_num].omega = TransitBisec->omega[lc_num];
  TransitBisec->transitparams[lc_num].exptime = TransitBisec->exptime[lc_num];
  memcpy(TransitBisec->transitparams[lc_num].ldcoeffs, TransitBisec->ldcoeffs[lc_num], TransitBisec->Nld*sizeof(double));

  /* Perform the actual routine */
  if(TransitBisec->geterror) {
    do_transitbisec(NJD, t, mag, 
		    &(TransitBisec->bisecdata[lc_num]),
		    &(TransitBisec->transitparams[lc_num]),
		    bisecoutname,
		    bisecspanoutname,
		    TransitBisec->outputbisec,
		    TransitBisec->outputbisecspan,
		    TransitBisec->max_delta_t,
		    TransitBisec->min_oot_frac,
		    TransitBisec->BS_highval,
		    TransitBisec->BS_lowval,
		    &(TransitBisec->out_BSrms[lc_num]),
		    &(TransitBisec->out_BSmed[lc_num]),
		    &(TransitBisec->out_Ntransits[lc_num]),
		    TransitBisec->geterror,
		    TransitBisec->Nerrorsims,
		    &(TransitBisec->out_BSrms_sim[lc_num]),
		    &(TransitBisec->out_BSmed_sim[lc_num]),
		    TransitBisec->refit_Tc);
  }
  else {
    do_transitbisec(NJD, t, mag, 
		    &(TransitBisec->bisecdata[lc_num]),
		    &(TransitBisec->transitparams[lc_num]),
		    bisecoutname,
		    bisecspanoutname,
		    TransitBisec->outputbisec,
		    TransitBisec->outputbisecspan,
		    TransitBisec->max_delta_t,
		    TransitBisec->min_oot_frac,
		    TransitBisec->BS_highval,
		    TransitBisec->BS_lowval,
		    &(TransitBisec->out_BSrms[lc_num]),
		    &(TransitBisec->out_BSmed[lc_num]),
		    &(TransitBisec->out_Ntransits[lc_num]),
		    TransitBisec->geterror,
		    TransitBisec->Nerrorsims,
		    NULL, NULL,
		    TransitBisec->refit_Tc);
  }
}
