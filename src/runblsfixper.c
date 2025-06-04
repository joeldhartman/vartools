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
#include "commands.h"
#include "programdata.h"
#include "functions.h"

/* Default parameters */
#define DEFAULT_MINPER 1.0
#define DEFAULT_MAXPER 10.0
#define DEFAULT_NFREQ 10000
#define DEFAULT_NBINS 400
#define DEFAULT_QMIN 0.01
#define DEFAULT_QMAX 0.10
#define DEFAULT_TIMEZONE -7.
#define MAXSTRINGLENGTH 512

#ifdef MAX_
#undef MAX_
#endif
#define MAX_(A,B) ((A) > (B) ? (A) : (B))

#ifdef MIN_
#undef MIN_
#endif
#define MIN_(A,B) ((A) < (B) ? (A) : (B))

#ifdef ABS_
#undef ABS_
#endif
#define ABS_(A) ((A) > 0 ? (A) : (-(A)))

#define CLIP_FACTOR 3.0
#define BIN_FACTOR 100

/* This version of BLS calculates the average magnitude, delta-chi^2 for the best positive and negative dips (as described in Burke et al. 2006), SDE, SR, the transit period, the depth the phases of transit start and end, assuming phase zero occurs at the first observation */

/* C port of the eebls.f routine -------
c
c------------------------------------------------------------------------
c     >>>>>>>>>>>> This routine computes BLS spectrum <<<<<<<<<<<<<<
c
c         [ see Kovacs, Zucker & Mazeh 2002, A&A, Vol. 391, 369 ]
c
c     This is the slightly modified version of the original BLS routine
c     by considering Edge Effect (EE) as suggested by
c     Peter R. McCullough [ pmcc@stsci.edu ].
c
c     This modification was motivated by considering the cases when
c     the low state (the transit event) happened to be devided between
c     the first and last bins. In these rare cases the original BLS
c     yields lower detection efficiency because of the lower number of
c     data points in the bin(s) covering the low state.
c
c     For further comments/tests see  www.konkoly.hu/staff/kovacs.html
c------------------------------------------------------------------------
c
c     Input parameters:
c     ~~~~~~~~~~~~~~~~~
c
c     n    = number of data points
c     t    = array {t(i)}, containing the time values of the time series
c     x    = array {x(i)}, containing the data values of the time series
c     e    = array {e(i)}, containing the uncertainty values of the time series
c     u    = temporal/work/dummy array, must be dimensioned in the
c            calling program in the same way as  {t(i)}
c     v    = the same as  {u(i)}
c     nb   = number of bins in the folded time series at any test period
c     qmi  = minimum fractional transit length to be tested
c     qma  = maximum fractional transit length to be tested
c
c     Output parameters:
c     ~~~~~~~~~~~~~~~~~~
c
c     period = the fixed period to run BLS at
c     bt0  = the time of the first transit center after the start of the
c            light curve.
c     bpow = value of sr for the best peak
c     depth= depth of the transit
c     qtran= fractional transit length  [ T_transit/bper ]
c     in1  = bin index at the start of the transit [ 0 < in1 < nb+1 ]
c     in2  = bin index at the end   of the transit [ 0 < in2 < nb+1 ]
      chisqrplus = delta_chisqr for the best transit like signal
      chisqrminus = delta_chisqr for the best inverse transit like signal

c -- added sde - the signal detection efficiency
c
c
c     Remarks:
c     ~~~~~~~~
c
c     -- *nb*   MUST be lower than  *nbmax*
c     -- Dimensions of arrays {y(i)} and {ibi(i)} MUST be greater than
c        or equal to  *nbmax*.
c     -- The lowest number of points allowed in a single bin is equal
c        to   MAX(minbin,qmi*N),  where   *qmi*  is the minimum transit
c        length/trial period,   *N*  is the total number of data points,
c        *minbin*  is the preset minimum number of the data points per
c        bin.
c
c========================================================================
c
*/

int eeblsfixper(int n_in, double *t_in, double *x_in, double *e_in, double *u, double *v, int nb, double qmi, double qma, double *period, double *bt0, double *bpow, double *depth, double *qtran, int *in1, int *in2, double *in1_ph, double *in2_ph, double *chisqrplus, double *chisqrminus, double *meanmagval, double timezone, double *fraconenight, int omodel, char *modelname, int correctlc, int ascii,int *nt, int *Nt, int *Nbefore, int *Nafter, double *rednoise, double *whitenoise, double *sigtopink, int fittrap, double *qingress, double *OOTmag, double *srsumout, int lcnum, int lclistnum, int usemask, _Variable *maskvar)
{
  double *y = NULL;
  double *ibi = NULL;
  int minbin = 5;
  int nbmax, nbtot;

  double powerplus, powerminus, bpowminus;
  double sumweights;
  double tot, rnbtot, *weight;
  double rn, s,t1,f0,p0,ph,ph2,phb1,phb2,pow,rn1,rn3,s3,rn4,rn5;
  double kkmi, kk;
  double srsum = 0.0;

  int kmi, kma,nb1,nbkma,i,j,k,jn1,jn2,jnb,nb2,nsr;

  long double sde_sr_ave, sde_srsqr_ave;
  FILE *outfile2;
  int n;
  double *t, *x, *e;
  double *t_mask = NULL, *x_mask = NULL, *e_mask = NULL;

  if(!usemask) {
    n = n_in;
    t = t_in;
    x = x_in;
    e = e_in;
  } else {
    if((t_mask = (double *) malloc(n_in*sizeof(double))) == NULL ||
       (x_mask = (double *) malloc(n_in*sizeof(double))) == NULL ||
       (e_mask = (double *) malloc(n_in*sizeof(double))) == NULL) {
      error(ERR_MEMALLOC);
    }
    n = 0;
    for(i = 0; i < n_in; i++) {
      if(!isnan(x_in[i]) && EvaluateVariable_Double(lclistnum, lcnum, i, maskvar) > VARTOOLS_MASK_TINY) {
	t_mask[n] = t_in[i];
	x_mask[n] = x_in[i];
	e_mask[n] = e_in[i];
	n++;
      }
    }
    t = t_mask;
    x = x_mask;
    e = e_mask;
    if(n <= 1) {
      *bpow = -1.;
      *depth = -1.;
      *bt0 = -1;
      *in1 = -1;
      *in2 = -1;
      *in1_ph = -1.;
      *in2_ph = -1.;
      *qtran = -1.;
      *chisqrplus = -1.;
      *chisqrminus = -1.;
      *meanmagval = -1.;
      *fraconenight = -1.;
      *qingress = -1.;
      *nt = 0;
      *Nt = 0;
      *Nbefore = 0;
      *Nafter = 0;
      *OOTmag = -1.;
      *rednoise = -1.;
      *whitenoise = -1.;
      *sigtopink = -1.;
      if(y != NULL) free(y);
      if(ibi != NULL) free(ibi);
      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
      return(1);
    }
  }

  nbmax = 2*nb;

  /***********************************************************/

  if(nb > nbmax) {
    error(ERR_BLSNBMAX);
  }
  tot = t[n-1] - t[0];

  if((y = (double *) malloc(nbmax * sizeof(double))) == NULL ||
     (ibi = (double *) malloc(nbmax * sizeof(double))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(3);
    }

  /**********************************************************/

  sumweights = 0.;
  weight = (double *) malloc(n * sizeof(double));
  for(i=0;i<n;i++)
    {
      weight[i] = 1./(e[i]*e[i]);
      sumweights += weight[i];
    }
  for(i=0;i<n;i++)
    weight[i] = weight[i] / sumweights;

  rn = (double) n;
  kmi = (int) (qmi*(double)nb);
  if(kmi < 1) kmi = 1;
  kma = ((int) (qma*(double)nb)) + 1;
  kkmi = qmi;
  if(kkmi < (double) minbin / rn) kkmi = (double) minbin / rn;
  //*bpow = 0.;

  bpowminus = 0.;
  sde_sr_ave = 0.;
  sde_srsqr_ave = 0.;

  /**************The following variables are defined for the extension
		 c     of arrays  ibi()  and  y()  [ see below ] ***************/

  nb1 = nb;
  nbkma = nb+kma;

  /*
    c
    c=================================
    c     Set temporal time series
    c=================================
    c
  */

  //sr_ave = 0.;
  //srsqr_ave = 0.;
  nsr = 0;

  s = 0.;
  t1 = t[0];
  for(i=0;i<n;i++)
    {
      u[i]=t[i]-t1;
      s += x[i]*weight[i];
    }
  (*meanmagval) = s;
  //s /= sumweights;
  for(i=0;i<n;i++)
    v[i]=x[i]-s;


  if(period[0] < 0)
    {
      free(weight);
      *bpow = -1.;
      *depth = -1.;
      *bt0 = -1;
      *in1 = -1;
      *in2 = -1;
      *in1_ph = -1.;
      *in2_ph = -1.;
      *qtran = -1.;
      *chisqrplus = -1.;
      *chisqrminus = -1.;
      *meanmagval = -1.;
      *fraconenight = -1.;
      *qingress = -1.;
      *nt = 0;
      *Nt = 0;
      *Nbefore = 0;
      *Nafter = 0;
      *OOTmag = -1.;
      *rednoise = -1.;
      *whitenoise = -1.;
      *sigtopink = -1.;
      if(y != NULL) free(y);
      if(ibi != NULL) free(ibi);
      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
      return(1);
    }


  f0=1./period[0];
  p0=period[0];

      /*
	c
	c======================================================
	c     Compute folded time series with  *p0*  period
	c======================================================
	c
      */


  for(j=0;j<nb;j++)
    {
      y[j] = 0.;
      ibi[j] = 0.;
    }
  nbtot = 0;
  for(i=0;i<n;i++)
    {
      ph = u[i]*f0;
      ph -= (int) ph;
      j = (int) (nb*ph);
      ibi[j] += weight[i];
      nbtot++;
      y[j] += v[i]*weight[i];
    }
  /*      for(i=0;i<nb;i++)
	  y[i] = ibi[i] * y[i] / rn;
  */
  /*
    c
    c-----------------------------------------------
    c     Extend the arrays  ibi()  and  y() beyond
    c     nb   by  wrapping
    c
  */

  for(j=nb1;j<nbkma;j++)
    {
      jnb = j - nb;
      ibi[j] = ibi[jnb];
      nbtot += ibi[j];
      y[j] = y[jnb];
    }
  rnbtot = (double) nbtot;
  /*
    c-----------------------------------------------
    c
    c===============================================
    c     Compute BLS statistics for this period
    c===============================================
    c
  */

  powerplus = 0.;
  powerminus = 0.;
  for(i=0;i<nb;i++)
    {
      s = 0.;
      k = 0;
      kk = 0.;
      nb2 = i+kma;
      for(j=i;j<nb2;j++)
	{
	  k++;
	  kk += ibi[j];
	  s += y[j];
	  if(k >= kmi && kk >= kkmi)
	    {
	      rn1 = (double) kk;
	      rn4 = (double) k;
	      pow = s*s/(rn1*(1. - rn1));
	      if(s > 0. && srsumout != NULL) {
		srsum += sqrt(pow);
	      }
	      if(s > 0. && pow >= powerplus)
		{
		  powerplus = pow;
		  jn1 = i;
		  jn2 = j;
		  rn3 = rn1;
		  rn5 = rn4;
		  s3 = s;
		}
	      else if(s < 0. && pow >= powerminus)
		{
		  powerminus = pow;
		}
	    }
	}
    }
  if(srsumout != NULL) {
    srsum = srsum / (double) kma;
    *srsumout = srsum;
  }
  powerplus = sqrt(powerplus);
  *bpow = powerplus;
  powerminus = sqrt(powerminus);
  bpowminus = powerminus;
  *in1 = jn1;
  *in2 = jn2;
  *in1_ph = ((double) (*in1) / (double) nb);
  *in2_ph = ((double) (*in2) / (double) nb);
  *qtran = (double)(jn2 - jn1 + 1)/(double)nb;
  *depth = powerplus/sqrt(rn3*(1.-rn3));
  if(fittrap) {
    *qingress=0.25;
    *OOTmag=*meanmagval;
    dofittrap_amoeba(n, t, x, e, *period, qtran, qingress, in1_ph, in2_ph, depth, OOTmag);
  } else {
    *qingress = 0.;
    *OOTmag = *meanmagval;
  }
  // Be sure to correct for transits past the edge
  if(*in2 >= nb) *in2 = *in2 - nb;
  *bt0 = t[0] + (0.5*(*qtran)+(*in1_ph))*(*period);


  *chisqrplus = -(*bpow)*(*bpow)*sumweights;
  *chisqrminus = -bpowminus*bpowminus*sumweights;

  *fraconenight = getfrac_onenight(n, t, u, v, e, period[0], *depth, *qtran, (t[0] + ((*in1_ph))*(*period)), timezone);
  getsignaltopinknoiseforgivenblsmodel(n, t, x, e, period[0], *qtran, *depth, *in1_ph, nt, Nt, Nbefore, Nafter, rednoise, whitenoise, sigtopink,*qingress,*OOTmag,NULL);


  //output the model light curve if asked to
  if(omodel)
    {
      if((outfile2 = fopen(modelname,"w")) == NULL)
	error2(ERR_CANNOTWRITE,modelname);

      f0 = 1./(period[0]);

      phb1 = (*qingress)*(*qtran);
      phb2 = (*qtran) - phb1;

      fprintf(outfile2,"#Time  Mag_obs   Mag_model   Error   Phase\n");
      for(i=0;i<n;i++)
	{
	  ph = (u[i] - (*in1_ph)*period[0])*f0;
	  ph -= floor(ph);
	  if(ph >= (*qtran)) {
	    ph2 = ph - 0.5*(*qtran);
	    if(ph2 < 0)
	      ph2 += 1.;
	    fprintf(outfile2,"%f %f %f %f %f\n",t[i], x[i], *OOTmag, e[i], ph2);
	  }
	  else {
	    if(ph >= phb1 && ph <= phb2) {
	      ph2 = ph - 0.5*(*qtran);
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f %f %f\n", t[i], x[i], (*OOTmag)+(*depth), e[i], ph2);
	    }
	    else if(ph < phb1) {
	      ph2 = ph - 0.5*(*qtran);
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f %f %f\n", t[i], x[i], (*OOTmag)+(*depth)*ph/phb1, e[i], ph2);
	    }
	    else {
	      ph2 = ph - 0.5*(*qtran);
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f %f %f\n", t[i], x[i], (*OOTmag)+(*depth)*((*qtran) - ph)/phb1, e[i], ph2);
	    }
	  }
	}
      fclose(outfile2);

    }
  if(correctlc)
    {
      f0 = 1./(period[0]);

      phb1 = (*qingress)*(*qtran);
      phb2 = (*qtran) - phb1;
      if(usemask) {
	n = n_in;
	t = t_in;
	x = x_in;
	e = e_in;
	for(i=0;i<n;i++)
	  {
	    u[i]=t[i]-t1;
	  }
      }
      for(i=0;i<n;i++)
	{
	  ph = (u[i] - (*in1_ph)*period[0])*f0;
	  ph -= floor(ph);
	  if(ph < (*qtran)) {
	    if(ph >= phb1 && ph <= phb2) {
	      x[i] -= (*depth);
	    }
	    else if(ph < phb1) {
	      x[i] -= (*depth)*ph/phb1;
	    }
	    else {
	      x[i] -= (*depth)*((*qtran) - ph)/phb1;
	    }
	  }
	}

    }

  free(weight);

  if(y != NULL) free(y);
  if(ibi != NULL) free(ibi);
  if(t_mask != NULL) free(t_mask);
  if(x_mask != NULL) free(x_mask);
  if(e_mask != NULL) free(e_mask);
  return(0);
}

/* This version adjusts the qmin and qmax according the period using a specified rmin and rmax, it assumes that for P in days and R in solar radii that q is given by:
q = 0.076 * R**(2/3) / P**(2/3)
*/
int eeblsfixper_rad(int n_in, double *t_in, double *x_in, double *e_in, double *u, double *v, int nb, double rmin, double rmax, double *period, double *bt0, double *bpow, double *depth, double *qtran, int *in1, int *in2, double *in1_ph, double *in2_ph, double *chisqrplus, double *chisqrminus, double *meanmagval, double timezone, double *fraconenight, int omodel, char *modelname, int correctlc, int ascii,int *nt, int *Nt, int *Nbefore, int *Nafter, double *rednoise, double *whitenoise, double *sigtopink, int fittrap, double *qingress, double *OOTmag, double *srsumout, int lcnum, int lclistnum, int usemask, _Variable *maskvar)
{
  double *y = NULL;
  double *ibi = NULL;
  double qmi, qma;
  int minbin = 5;
  int nbmax, nbtot;

  double powerplus, powerminus, bpowminus;
  double sumweights;
  double tot, rnbtot, *weight;
  double rn, s,t1,f0,p0,ph,ph2,phb1,phb2,pow_,rn1,rn3,s3,rn4,rn5;
  double kkmi, kk;

  int kmi, kma,nb1,nbkma,i,j,k,jn1,jn2,jnb,nb2,nsr;
  double rminpow, rmaxpow, Ppow;
  long double sde_sr_ave, sde_srsqr_ave;
  FILE *outfile2;
  double srsum = 0.0;
  int n;
  double *t, *x, *e;
  double *t_mask = NULL, *x_mask = NULL, *e_mask = NULL;

  if(!usemask) {
    n = n_in;
    t = t_in;
    x = x_in;
    e = e_in;
  } else {
    if((t_mask = (double *) malloc(n_in*sizeof(double))) == NULL ||
       (x_mask = (double *) malloc(n_in*sizeof(double))) == NULL ||
       (e_mask = (double *) malloc(n_in*sizeof(double))) == NULL) {
      error(ERR_MEMALLOC);
    }
    n = 0;
    for(i = 0; i < n_in; i++) {
      if(!isnan(x_in[i]) && EvaluateVariable_Double(lclistnum, lcnum, i, maskvar) > VARTOOLS_MASK_TINY) {
	t_mask[n] = t_in[i];
	x_mask[n] = x_in[i];
	e_mask[n] = e_in[i];
	n++;
      }
    }
    t = t_mask;
    x = x_mask;
    e = e_mask;
    if(n <= 1) {
      *bpow = -1.;
      *depth = -1.;
      *bt0 = -1;
      *in1 = -1;
      *in2 = -1;
      *in1_ph = -1.;
      *in2_ph = -1.;
      *qtran = -1.;
      *chisqrplus = -1.;
      *chisqrminus = -1.;
      *meanmagval = -1.;
      *fraconenight = -1.;
      *qingress = -1.;
      *nt = 0;
      *Nt = 0;
      *Nbefore = 0;
      *Nafter = 0;
      *OOTmag = -1.;
      *rednoise = -1.;
      *whitenoise = -1.;
      *sigtopink = -1.;
      if(y != NULL) free(y);
      if(ibi != NULL) free(ibi);
      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
      return(1);
    }
  }

  nbmax = 2*nb;

  /***********************************************************/

  if(nb > nbmax) {
    error(ERR_BLSNBMAX);
  }
  tot = t[n-1] - t[0];

  /**********************************************************/

  if((y = (double *) malloc(nbmax * sizeof(double))) == NULL ||
     (ibi = (double *) malloc(nbmax * sizeof(double))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(3);
    }

  rminpow = pow(rmin,0.6666667);
  rmaxpow = pow(rmax,0.6666667);
  Ppow = pow(period[0],0.6666667);
  qmi = 0.076*rminpow/Ppow;
  qma = 0.076*rmaxpow/Ppow;


  sumweights = 0.;
  weight = (double *) malloc(n * sizeof(double));
  for(i=0;i<n;i++)
    {
      weight[i] = 1./(e[i]*e[i]);
      sumweights += weight[i];
    }
  for(i=0;i<n;i++)
    weight[i] = weight[i] / sumweights;

  rn = (double) n;
  kmi = (int) (qmi*(double)nb);
  if(kmi < 1) kmi = 1;
  kma = ((int) (qma*(double)nb)) + 1;
  kkmi = qmi;
  if(kkmi < (double) minbin / rn) kkmi = (double) minbin / rn;
  //*bpow = 0.;

  bpowminus = 0.;
  sde_sr_ave = 0.;
  sde_srsqr_ave = 0.;

  /**************The following variables are defined for the extension
		 c     of arrays  ibi()  and  y()  [ see below ] ***************/

  nb1 = nb;
  nbkma = nb+kma;

  /*
    c
    c=================================
    c     Set temporal time series
    c=================================
    c
  */

  //sr_ave = 0.;
  //srsqr_ave = 0.;
  nsr = 0;

  s = 0.;
  t1 = t[0];
  for(i=0;i<n;i++)
    {
      u[i]=t[i]-t1;
      s += x[i]*weight[i];
    }
  (*meanmagval) = s;
  //s /= sumweights;
  for(i=0;i<n;i++)
    v[i]=x[i]-s;

  if(period[0] < 0)
    {
      free(weight);
      *bpow = -1.;
      *depth = -1.;
      *bt0 = -1.;
      *in1 = -1;
      *in2 = -1;
      *in1_ph = -1.;
      *in2_ph = -1.;
      *qtran = -1.;
      *chisqrplus = -1.;
      *chisqrminus = -1.;
      *meanmagval = -1.;
      *fraconenight = -1.;
      *qingress = -1.;
      *nt = 0;
      *Nt = 0;
      *Nbefore = 0;
      *Nafter = 0;
      *OOTmag = -1.;
      *rednoise = -1.;
      *whitenoise = -1.;
      *sigtopink = -1.;
      if(y != NULL) free(y);
      if(ibi != NULL) free(ibi);
      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
      return(1);
    }


  f0=1./period[0];
  p0=period[0];

      /*
	c
	c======================================================
	c     Compute folded time series with  *p0*  period
	c======================================================
	c
      */


  for(j=0;j<nb;j++)
    {
      y[j] = 0.;
      ibi[j] = 0.;
    }
  nbtot = 0;
  for(i=0;i<n;i++)
    {
      ph = u[i]*f0;
      ph -= (int) ph;
      j = (int) (nb*ph);
      ibi[j] += weight[i];
      nbtot++;
      y[j] += v[i]*weight[i];
    }
  /*      for(i=0;i<nb;i++)
	  y[i] = ibi[i] * y[i] / rn;
  */
  /*
    c
    c-----------------------------------------------
    c     Extend the arrays  ibi()  and  y() beyond
    c     nb   by  wrapping
    c
  */

  for(j=nb1;j<nbkma;j++)
    {
      jnb = j - nb;
      ibi[j] = ibi[jnb];
      nbtot += ibi[j];
      y[j] = y[jnb];
    }
  rnbtot = (double) nbtot;
  /*
    c-----------------------------------------------
    c
    c===============================================
    c     Compute BLS statistics for this period
    c===============================================
    c
  */

  powerplus = 0.;
  powerminus = 0.;
  for(i=0;i<nb;i++)
    {
      s = 0.;
      k = 0;
      kk = 0.;
      nb2 = i+kma;
      for(j=i;j<nb2;j++)
	{
	  k++;
	  kk += ibi[j];
	  s += y[j];
	  if(k >= kmi && kk >= kkmi)
	    {
	      rn1 = (double) kk;
	      rn4 = (double) k;
	      pow_ = s*s/(rn1*(1. - rn1));
	      if(s > 0. && srsumout != NULL) {
		srsum += sqrt(pow_);
	      }
	      if(s > 0. && pow_ >= powerplus)
		{
		  powerplus = pow_;
		  jn1 = i;
		  jn2 = j;
		  rn3 = rn1;
		  rn5 = rn4;
		  s3 = s;
		}
	      else if(s < 0. && pow_ >= powerminus)
		{
		  powerminus = pow_;
		}
	    }
	}
    }
  if(srsumout != NULL) {
    srsum = srsum / (double) kma;
    *srsumout = srsum;
  }
  powerplus = sqrt(powerplus);
  *bpow = powerplus;
  powerminus = sqrt(powerminus);
  bpowminus = powerminus;
  *in1 = jn1;
  *in2 = jn2;
  *in1_ph = ((double) (*in1) / (double) nb);
  *in2_ph = ((double) (*in2) / (double) nb);
  *qtran = (double)(jn2 - jn1 + 1)/(double)nb;
  // Be sure to correct for transits past the edge
  if(*in2 >= nb) *in2 = *in2 - nb;
  *depth = powerplus/sqrt(rn3*(1.-rn3));
  if(fittrap) {
    *qingress=0.25;
    *OOTmag=*meanmagval;
    dofittrap_amoeba(n, t, x, e, *period, qtran, qingress, in1_ph, in2_ph, depth, OOTmag);
  } else {
    *qingress = 0.;
    *OOTmag = *meanmagval;
  }
  *bt0 = t[0] + (0.5*(*qtran)+(*in1_ph))*(*period);

  *chisqrplus = -(*bpow)*(*bpow)*sumweights;
  *chisqrminus = -bpowminus*bpowminus*sumweights;

  *fraconenight = getfrac_onenight(n, t, u, v, e, period[0], *depth, *qtran, (t[0] + ((*in1_ph))*(*period)), timezone);
  getsignaltopinknoiseforgivenblsmodel(n, t, x, e, period[0], *qtran, *depth, *in1_ph, nt, Nt, Nbefore, Nafter, rednoise, whitenoise, sigtopink, *qingress, *OOTmag,NULL);

  //output the model light curve if asked to
  if(omodel)
    {
      if((outfile2 = fopen(modelname,"w")) == NULL)
	error2(ERR_CANNOTWRITE,modelname);

      f0 = 1./(period[0]);

      phb1 = (*qingress)*(*qtran);
      phb2 = (*qtran) - phb1;

      fprintf(outfile2,"#Time  Mag_obs   Mag_model   Error   Phase\n");
      for(i=0;i<n;i++)
	{
	  ph = (u[i] - (*in1_ph)*period[0])*f0;
	  ph -= floor(ph);
	  if(ph >= (*qtran)) {
	    ph2 = ph - 0.5*(*qtran);
	    if(ph2 < 0)
	      ph2 += 1.;
	    fprintf(outfile2,"%f %f %f %f %f\n",t[i], x[i], *OOTmag, e[i], ph2);
	  }
	  else {
	    if(ph >= phb1 && ph <= phb2) {
	      ph2 = ph - 0.5*(*qtran);
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f %f %f\n", t[i], x[i], (*OOTmag)+(*depth), e[i], ph2);
	    }
	    else if(ph < phb1) {
	      ph2 = ph - 0.5*(*qtran);
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f %f %f\n", t[i], x[i], (*OOTmag)+(*depth)*ph/phb1, e[i], ph2);
	    }
	    else {
	      ph2 = ph - 0.5*(*qtran);
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f %f %f\n", t[i], x[i], (*OOTmag)+(*depth)*((*qtran) - ph)/phb1, e[i], ph2);
	    }
	  }
	}
      fclose(outfile2);
    }
  if(correctlc)
    {
      f0 = 1./(period[0]);

      phb1 = (*qingress)*(*qtran);
      phb2 = (*qtran) - phb1;
      if(usemask) {
	n = n_in;
	t = t_in;
	x = x_in;
	e = e_in;
	for(i=0;i<n;i++)
	  {
	    u[i]=t[i]-t1;
	  }
      }
      for(i=0;i<n;i++)
	{
	  ph = (u[i] - (*in1_ph)*period[0])*f0;
	  ph -= floor(ph);
	  if(ph < (*qtran)) {
	    if(ph >= phb1 && ph <= phb2) {
	      x[i] -= (*depth);
	    }
	    else if(ph < phb1) {
	      x[i] -= (*depth)*ph/phb1;
	    }
	    else {
	      x[i] -= (*depth)*((*qtran) - ph)/phb1;
	    }
	  }
	}

    }

  free(weight);

  if(y != NULL) free(y);
  if(ibi != NULL) free(ibi);
  if(t_mask != NULL) free(t_mask);
  if(x_mask != NULL) free(x_mask);
  if(e_mask != NULL) free(e_mask);
  return(0);
}

