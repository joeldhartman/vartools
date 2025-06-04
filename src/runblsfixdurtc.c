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

void getclippedavestddev_blsfixdurtc(int n, double *pow, double *ave_out, double *stddev_out)
{
  int i, n2, n3;
  double ave1, ave2, ave3, ave4;
  ave1 = 0.; ave2 = 0.;
  n3 = 0;
  for(i=0;i<n;i++)
    {
      if(pow[i] > 0.) {
	ave1 += pow[i];
	ave2 += pow[i]*pow[i];
	n3++;
      }
    }
  if(n3 == 0) {
    *ave_out = 0.0;
    *stddev_out = 1.0;
    return;
  }
  ave1 /= n3;
  ave2 = sqrt((ave2 / n3) - (ave1*ave1));
  ave3 = 0.;
  ave4 = 0.;
  n2 = 0;
  for(i=0;i<n;i++)
    {
      if(pow[i] > 0.) {
	if((pow[i] - ave1) < CLIP_FACTOR*ave2)
	  {
	    ave3 += pow[i];
	    ave4 += pow[i]*pow[i];
	    n2++;
	  }
      }
    }
  if(n2 == 0) {
    *ave_out = ave1;
    *stddev_out = ave2;
    return;
  }
  ave3 /= n2;
  ave4 = sqrt((ave4 / n2) - (ave3*ave3));
  *ave_out = ave3;
  *stddev_out = ave4;
  return;
}


int eeblsfixdurtc(int n_in, double *t_in, double *x_in, double *e_in, double *u, double *v, double inputTC, double inputdur, int fixdepth, double inputdepth, double inputqgress, int nf, double fmin, double df, double *p, int Npeak, double *bper, double *bt0, double *bpow, double *sde, double *snval, double *depth, double *qtran, double *chisqrplus, double *chisqrminus, double *bperpos, double *meanmagval, double timezone, double *fraconenight, int operiodogram, char *outname, int omodel, char *modelname, int correctlc, int ascii,int *nt, int *Nt, int *Nbefore, int *Nafter, double *rednoise, double *whitenoise, double *sigtopink, int fittrap, double *qingress, double *OOTmag, int ophcurve, char *ophcurvename, double phmin, double phmax, double phstep, int ojdcurve, char *ojdcurvename, double jdstep, int lcnum, int lclistnum, int usemask, _Variable *maskvar)
{

  double dum1, dum2;
  double y[2000];
  double ibi[2000];
  int minbin = 5;
  int nbmax = 2000, nbtot;
  int nsrvals, nsrvals_minus, test, foundsofar, dumint1;
  double powerplus, powerminus, bpowminus, dumdbl1, dumdbl2, jdtmp;
  double sumweights, phb1, phb2;
  double tot, rnbtot, *weight, sr_minus;
  double rn, s,t1,f0,p0,ph,ph2,pow,rn1,rn3,s3,rn4,rn5;
  double kkmi, kk, allave, allstddev, allave_minus, allstddev_minus, *qtran_array, *depth_array, minbest, qf;
  double *sr_ave, *binned_sr_ave, *binned_sr_sig;
  double in1ph, in2ph;
  int kmi, kma,nb1,nbkma,i,jf,j,k,jn1,jn2,jnb,nb2,nsr,nclippedfreq, *best_id;
  double *p_minus, *bper_array, *sr_ave_minus, *binned_sr_ave_minus, *binned_sr_sig_minus, global_best_sr_ave, global_best_sr_stddev;
  double global_best_sr_ave_inv, global_best_sr_stddev_inv;
  long double sde_sr_ave, sde_srsqr_ave;
  FILE *outfile, *outfile2;

  double inputT0;

  int nb;
  int in1, in2;
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
      for(j=0;j<Npeak;j++)
	{
	  bper[j] = -1.;
          bt0[j] = -1.;
	  snval[j] = -1.;
	  bpow[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	}
      *bperpos = -1.;
      *chisqrminus = 999999.;
      *meanmagval = -1.;
      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
    }
  }

  /***********************************************************/

  if((sr_ave = (double *) malloc(nf * sizeof(double))) == NULL ||
     (binned_sr_ave = (double *) malloc(nf * sizeof(double))) == NULL ||
     (binned_sr_sig = (double *) malloc(nf * sizeof(double))) == NULL ||
     (qtran_array = (double *) malloc(nf * sizeof(double))) == NULL ||
     (depth_array = (double *) malloc(nf * sizeof(double))) == NULL ||
     (bper_array = (double *) malloc(nf * sizeof(double))) == NULL ||
     (sr_ave_minus = (double *) malloc(nf * sizeof(double))) == NULL ||
     (binned_sr_ave_minus = (double *) malloc(nf * sizeof(double))) == NULL ||
     (binned_sr_sig_minus = (double *) malloc(nf * sizeof(double))) == NULL ||
     (p_minus = (double *) malloc(nf * sizeof(double))) == NULL ||
     (best_id = (int *) malloc(Npeak * sizeof(int))) == NULL)

    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(3);
    }

  tot = t[n-1] - t[0];
  if(fmin < 1./tot) {
    error(ERR_BLSFMINTOOSMALL);
  }


  inputT0 = inputTC - inputdur / 2.0;

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
  //*bpow = 0.;

  bpowminus = 0.;
  sde_sr_ave = 0.;
  sde_srsqr_ave = 0.;

  /**************The following variables are defined for the extension
		 c     of arrays  ibi()  and  y()  [ see below ] ***************/


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
      u[i]=t[i]-inputT0;
      s += x[i]*weight[i];
    }
  (*meanmagval) = s;
  //s /= sumweights;
  for(i=0;i<n;i++)
    v[i]=x[i]-s;

  /*
    c
    c******************************
    c     Start period search - we modify this slightly to first compute
the periodogram, and then search it for peaks    *
    c******************************
    c
  */

  for(jf=0;jf<nf;jf++)
    {
      f0=fmin+df*((double)jf);
      p0=1./f0;

      /*
	c
	c======================================================
	c     Compute folded time series with  *p0*  period and
        c     the specified transit duration and T0
	c======================================================
	c
      */

      s = 0.0;
      qf = inputdur/p0;
      k = 0;
      kk = 0.;
      powerplus = 0.;
      powerminus = 0.;
      nsrvals = 0;
      nsrvals_minus = 0;
      if(!fixdepth) {
	for(i=0;i<n;i++)
	  {
	    ph = u[i]*f0;
	    ph -= floor(ph);
	    
	    if(ph < qf) {
	      s += v[i]*weight[i];
	      k++;
	      kk += weight[i];
	    }
	  }
	
	rn1 = (double) kk;
	rn4 = (double) k;
	pow = s*s/(rn1*(1. - rn1));
      } else {
	phb1 = inputqgress*qf;
	phb2 = qf - inputqgress*qf;
	for(i=0;i<n;i++)
	  {
	    ph = u[i]*f0;
	    ph -= floor(ph);
	    if(ph < qf) {
	      if(ph >= phb1 && ph <= phb2) {
		s += weight[i]*((v[i]*v[i]) - (v[i]-inputdepth)*(v[i]-inputdepth));
	      } else if(ph < phb1) {
		s += weight[i]*((v[i]*v[i]) - (v[i]-ph*inputdepth/qf/inputqgress)*(v[i]-ph*inputdepth/qf/inputqgress));
	      }
	      else {
		s += weight[i]*((v[i]*v[i]) - (v[i]-(qf-ph)*inputdepth/qf/inputqgress)*(v[i]-(qf-ph)*inputdepth/qf/inputqgress));
	      }
	      k++;
	      kk += weight[i];
	    }
	  }
	rn1 = (double) kk;
	pow = s;
	if(pow < 0) pow=-pow;
      }
      if(s > 0.)
	{
	  if(pow >= powerplus)
	    {
	      powerplus = pow;
	      rn3 = rn1;
	    }
	}
      else if(s < 0.)
	{
	  if(pow >= powerminus)
	    {
	      powerminus = pow;
	    }
	}
	
      // Find the average value of the srvals
      powerplus = sqrt(powerplus);
      sde_sr_ave += powerplus;
      sde_srsqr_ave += powerplus*powerplus;
      p[jf] = powerplus;
      powerminus = sqrt(powerminus);
      p_minus[jf] = powerminus;
      //sr_ave += powerplus;
      //srsqr_ave += powerplus*powerplus;
      nsr++;
      qtran_array[jf] = qf;
      if(!fixdepth)
	depth_array[jf] = powerplus/sqrt(rn3*(1.-rn3));
      else
	depth_array[jf] = inputdepth;
      bper_array[jf] = p0;
    }


  getclippedavestddev_blsfixdurtc(nf,p,&global_best_sr_ave,&global_best_sr_stddev);
  nclippedfreq = nf;

  /* Now let's find the peaks in the periodogram, first convert the periodogram from SR to SN ratio */

  if(nclippedfreq > Npeak)
    {
      for(i=0; i<nf; i++) {
	p[i] = (p[i] - global_best_sr_ave)/global_best_sr_stddev;
      }
    }
  else
    {
      /* We have no peaks, just put -1. for the values and return to the calling function */
      for(j=0;j<Npeak;j++)
	{
	  bper[j] = -1.;
          bt0[j] = -1.;
	  snval[j] = -1.;
	  bpow[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	}
      *bperpos = -1.;
      *chisqrminus = 999999.;
      *meanmagval = -1.;
      free(weight);
      free(best_id);
      free(sr_ave);
      free(binned_sr_ave);
      free(binned_sr_sig);
      free(qtran_array);
      free(depth_array);
      free(bper_array);
      free(sr_ave_minus);
      free(binned_sr_ave_minus);
      free(binned_sr_sig_minus);
      free(p_minus);
      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
      return 1;
    }

  foundsofar = 0;
  i = 0;
  while(foundsofar < Npeak && i < nf)
    {
      if(p[i] > 0)
	{
	  test = 1;
	  for(j=0;j<foundsofar;j++)
	    {
	      if(!isDifferentPeriods(MIN_(bper[j],bper_array[i]),MAX_(bper[j],bper_array[i]),tot))
		{
		  if(p[i] > snval[j])
		    {
		      bper[j] = bper_array[i];
		      snval[j] = p[i];
		      best_id[j] = i;
		    }
		  test = 0;
		  break;
		}
	    }
	  if(test)
	    {
	      snval[foundsofar] = p[i];
	      bper[foundsofar] = bper_array[i];
	      best_id[foundsofar] = i;
	      foundsofar++;
	    }
	}
      i++;
    }

  if(i < nf)
    {
      mysort3_int(Npeak,snval,bper,best_id);
      minbest = snval[0];
      for(;i<nf;i++)
	{
	  if(p[i] > minbest)
	    {
	      test = 1;
	      for(j=0;j<Npeak;j++)
		{
		  if(!isDifferentPeriods(MIN_(bper[j],bper_array[i]),MAX_(bper[j],bper_array[i]),tot))
		    {
		      if(p[i] > snval[j])
			{
			  snval[j] = p[i];
			  bper[j] = bper_array[i];
			  best_id[j] = i;
			  mysort3_int(Npeak,snval,bper,best_id);
			  minbest = snval[0];
			}
		      test = 0;
		      break;
		    }
		}
	      if(test)
		{
		  snval[0] = p[i];
		  bper[0] = bper_array[i];
		  best_id[0] = i;
		  mysort3_int(Npeak,snval,bper,best_id);
		  minbest = snval[0];
		}
	    }
	}
    }
  else if(foundsofar >= 1)
    {
      /* We have a few peaks, but Npeak of them */
      mysort3_int(foundsofar,snval,bper,best_id);
      for(j=foundsofar;j<Npeak;j++)
	{
	  /* Put -1 for the remaining peaks */
	  bper[j] = -1.;
	  snval[j] = -1.;
	  bpow[j] = -1.;
	  bt0[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	}
    }
  else
    {
      /* We have no peaks, just put -1. for the values and return to the calling function */
      for(j=0;j<Npeak;j++)
	{
	  bper[j] = -1.;
	  snval[j] = -1.;
	  bpow[j] = -1.;
	  bt0[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	}
      *bperpos = -1.;
      *chisqrminus = 999999.;
      *meanmagval = -1.;
      free(weight);
      free(best_id);
      free(sr_ave);
      free(binned_sr_ave);
      free(binned_sr_sig);
      free(qtran_array);
      free(depth_array);
      free(bper_array);
      free(sr_ave_minus);
      free(binned_sr_ave_minus);
      free(binned_sr_sig_minus);
      free(p_minus);
      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
      return 1;
    }
  //fprintf(stderr,"Error Running BLS - no frequencies survive clipping!\n");

  /* invert the snval, bper and best_id vectors */
  for(i = 0, j = foundsofar - 1; i < foundsofar/2 + 1; i++)
    {
      if(i < j)
	{
	  dumdbl1 = snval[j];
	  dumdbl2 = bper[j];
	  dumint1 = best_id[j];
	  snval[j] = snval[i];
	  bper[j] = bper[i];
	  best_id[j] = best_id[i];
	  snval[i] = dumdbl1;
	  bper[i] = dumdbl2;
	  best_id[i] = dumint1;
	}
      j--;
    }

  /* Collect all the output bls parameters for the peaks */
  for(i=0;i<Npeak;i++)
    {
      if(bper[i] > -1)
	{
	  in1ph = (inputT0 - t1)/bper[i] - floor((inputT0 - t1)/bper[i]);
	  in2ph = (inputT0 + inputdur - t1)/bper[i] - floor((inputT0 + inputdur - t1)/bper[i]);
	  nb = 100*ceil(qtran_array[best_id[i]]);
	  bpow[i] = snval[i]*global_best_sr_stddev + global_best_sr_ave;
	  if(fittrap && !fixdepth) {
	    qingress[i]=0.25;
	    OOTmag[i]=*meanmagval;
	    dum1 = in1ph;
	    dum2 = in2ph;
	    dofittrap_amoeba_fixdur(n, t, x, e, bper[i], (qtran_array[best_id[i]]), &(qingress[i]), in1ph, in2ph, &(depth_array[best_id[i]]), &(OOTmag[i]));
	  } else {
	    if(fixdepth)
	      qingress[i] = inputqgress;
	    else
	      qingress[i] = 0.;
	    OOTmag[i] = *meanmagval;
	  }
	  // Be sure to correct for transits past the edge
	  qtran[i] = qtran_array[best_id[i]];
	  bt0[i] = t[0] + (0.5*qtran[i]+in1ph)*bper[i];
	  depth[i] = depth_array[best_id[i]];
	  sde[i] = (bpow[i] - ((double)sde_sr_ave / (double)nsr))/sqrt((double)((sde_srsqr_ave / (long double) nsr) - (sde_sr_ave*sde_sr_ave/((long double)nsr*(long double)nsr))));
	  chisqrplus[i] = -bpow[i]*bpow[i]*sumweights;

	  fraconenight[i] = getfrac_onenight(n, t, u, v, e, bper[i], depth[i], qtran[i], (t[0] + in1ph*bper[i]), timezone);

	  /* Get the signal to pink noise for the peak */
	  getsignaltopinknoiseforgivenblsmodel(n, t, x, e, bper[i], qtran[i], depth[i], in1ph, &nt[i], &Nt[i], &Nbefore[i], &Nafter[i], &rednoise[i], &whitenoise[i], &sigtopink[i], qingress[i], OOTmag[i], NULL);
	}
    }

  /* Now find the maximum inverse transit */
  getclippedavestddev_blsfixdurtc(nf,p_minus,&global_best_sr_ave_inv,&global_best_sr_stddev_inv);
  nclippedfreq = nf;

  if(nclippedfreq > 0.)
    {
      for(i=0;i<nf;i++)
	{
	  p_minus[i] = (p_minus[i] - global_best_sr_ave_inv) / global_best_sr_stddev_inv;
	  if(p_minus[i] > bpowminus)
	    {
	      bpowminus = p_minus[i];
	      sr_minus = p_minus[i]*global_best_sr_stddev_inv + global_best_sr_ave_inv;
	      *bperpos = bper_array[i];
	      *chisqrminus = -sr_minus*sr_minus*sumweights;
	    }
	}
    }
  else
    {
      /* We have no peaks, just put -1. for the values and return to the calling function */
      /*for(j=0;j<Npeak;j++)
	{
	  bper[j] = -1.;
	  snval[j] = -1.;
	  bpow[j] = -1.;
	  in1[j] = -1;
	  in2[j] = -1;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	  }*/
      *bperpos = -1.;
      *chisqrminus = 999999.;
      *meanmagval = -1.;
      free(weight);
      free(best_id);
      free(sr_ave);
      free(binned_sr_ave);
      free(binned_sr_sig);
      free(qtran_array);
      free(depth_array);
      free(bper_array);
      free(sr_ave_minus);
      free(binned_sr_ave_minus);
      free(binned_sr_sig_minus);
      free(p_minus);
      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
      return 1;
    }

  //sde = (*bpow - ((double)sr_ave / (double)nsr))/sqrt((double)((srsqr_ave / (long double) nsr) - (sr_ave*sr_ave/(long double)(nsr*nsr))));

  /*
    c
    c     Edge correction of transit end index
    c
  */

  //output the periodogram if asked to
  if(operiodogram)
    {
      if((outfile = fopen(outname,"w")) == NULL)
	error2(ERR_CANNOTWRITE,outname);
      if(ascii)
	{
	  for(i=0;i<nf;i++)
	    fprintf(outfile,"%f %f\n",bper_array[i],p[i]);
	}
      else
	{
	  fwrite(&nf,4,1,outfile);
	  fwrite(bper_array,8,nf,outfile);
	  fwrite(p,8,nf,outfile);
	}
      fclose(outfile);
    }

  //output the model light curve if asked to
  if(omodel)
    {
      if((outfile2 = fopen(modelname,"w")) == NULL)
	error2(ERR_CANNOTWRITE,modelname);

      f0 = 1./bper[0];
      phb1 = qingress[0]*qtran[0];
      phb2 = qtran[0] - phb1;
      //in1ph = (inputT0 - t1)/bper[0] - floor((inputT0 - t1)/bper[0]);

      fprintf(outfile2,"#Time  Mag_obs   Mag_model   Error   Phase\n");
      for(i=0;i<n;i++)
	{
	  ph = (t[i] - inputT0)*f0;
	  ph -= floor(ph);
	  if(ph >= qtran[0]) {
	    ph2 = ph - 0.5*qtran[0];
	    if(ph2 < 0)
	      ph2 += 1.;
	    fprintf(outfile2,"%f %f %f %f %f\n",t[i], x[i], OOTmag[0], e[i], ph2);
	  }
	  else {
	    if(ph >= phb1 && ph <= phb2) {
	      ph2 = ph - 0.5*qtran[0];
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f %f %f\n", t[i], x[i], OOTmag[0]+depth[0], e[i], ph2);
	    }
	    else if(ph < phb1) {
	      ph2 = ph - 0.5*qtran[0];
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f %f %f\n", t[i], x[i], OOTmag[0]+depth[0]*ph/phb1, e[i], ph2);
	    }
	    else {
	      ph2 = ph - 0.5*qtran[0];
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f %f %f\n", t[i], x[i], OOTmag[0]+depth[0]*(qtran[0] - ph)/phb1, e[i], ph2);
	    }
	  }
	}
      fclose(outfile2);
    }
  // Output the phase curve if asked to.
  if(ophcurve)
    {
      if((outfile2 = fopen(ophcurvename,"w")) == NULL)
	error2(ERR_CANNOTWRITE,ophcurvename);

      fprintf(outfile2,"#Phase Mag_model\n");
      ph2 = phmin;
      phb1 = qingress[0]*qtran[0];
      phb2 = qtran[0] - phb1;
      while(ph2 <= phmax) {
	ph = ph2 + 0.5*qtran[0];
	ph -= floor(ph);
	if(ph >= qtran[0]) {
	  fprintf(outfile2,"%f %f\n",ph2,OOTmag[0]);
	}
	else {
	  if(ph >= phb1 && ph <= phb2) {
	    fprintf(outfile2,"%f %f\n",ph2,OOTmag[0]+depth[0]);
	  }
	  else if(ph < phb1) {
	    fprintf(outfile2,"%f %f\n",ph2,OOTmag[0]+depth[0]*ph/phb1);
	  }
	  else {
	    fprintf(outfile2,"%f %f\n",ph2,OOTmag[0]+depth[0]*(qtran[0] - ph)/phb1);
	  }
	}
	ph2 += phstep;
      }
      fclose(outfile2);
    }

  // Output the JD curve if asked to.
  if(ojdcurve)
    {
      if((outfile2 = fopen(ojdcurvename,"w")) == NULL)
	error2(ERR_CANNOTWRITE,ojdcurvename);

      fprintf(outfile2,"#Time Mag_model Phase\n");
      jdtmp = t[0];
      f0 = 1./bper[0];
      phb1 = qingress[0]*qtran[0];
      phb2 = qtran[0] - phb1;
      while(jdtmp <= t[n-1])
	{
	  ph = (jdtmp - inputT0)*f0;
	  ph -= floor(ph);
	  if(ph >= qtran[0]) {
	    ph2 = ph - 0.5*qtran[0];
	    if(ph2 < 0)
	      ph2 += 1.;
	    fprintf(outfile2,"%f %f %f\n", jdtmp, OOTmag[0], ph2);
	  }
	  else {
	    if(ph >= phb1 && ph <= phb2) {
	      ph2 = ph - 0.5*qtran[0];
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f\n", jdtmp, OOTmag[0]+depth[0], ph2);
	    }
	    else if(ph < phb1) {
	      ph2 = ph - 0.5*qtran[0];
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f\n", jdtmp, OOTmag[0]+depth[0]*ph/phb1, ph2);
	    }
	    else {
	      ph2 = ph - 0.5*qtran[0];
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f\n", jdtmp, OOTmag[0]+depth[0]*(qtran[0] - ph)/phb1, ph2);
	    }
	  }
	  jdtmp += jdstep;
	}
      fclose(outfile2);
    }
  if(correctlc)
    {
      f0 = 1./bper[0];
      phb1 = qingress[0]*qtran[0];
      phb2 = qtran[0] - phb1;
      if(usemask) {
	n = n_in;
	t = t_in;
	x = x_in;
	e = e_in;
      }
      for(i=0;i<n;i++)
	{
	  ph = (t[i] - inputT0)*f0;
	  ph -= floor(ph);
	  if(ph < qtran[0]) {
	    if(ph >= phb1 && ph <= phb2) {
	      x[i] -= depth[0];
	    }
	    else if(ph < phb1) {
	      x[i] -= depth[0]*ph/phb1;
	    }
	    else {
	      x[i] -= depth[0]*(qtran[0] - ph)/phb1;
	    }
	  }
	}

    }

  free(weight);
  free(best_id);
  free(sr_ave);
  free(binned_sr_ave);
  free(binned_sr_sig);
  free(qtran_array);
  free(depth_array);
  free(bper_array);
  free(sr_ave_minus);
  free(binned_sr_ave_minus);
  free(binned_sr_sig_minus);
  free(p_minus);

  if(t_mask != NULL) free(t_mask);
  if(x_mask != NULL) free(x_mask);
  if(e_mask != NULL) free(e_mask);

  return(0);
}


int eeblsfixperdurtc(int n_in, double *t_in, double *x_in, double *e_in, double *u, double *v, double inputper, double inputTC, double inputdur, int fixdepth, double inputdepth, double inputqgress, double *depth, double *qtran, double *chisqrplus, double *meanmagval, double timezone, double *fraconenight, int omodel, char *modelname, int correctlc, int *nt, int *Nt, int *Nbefore, int *Nafter, double *rednoise, double *whitenoise, double *sigtopink, int fittrap, double *qingress, double *OOTmag, int ophcurve, char *ophcurvename, double phmin, double phmax, double phstep, int ojdcurve, char *ojdcurvename, double jdstep, int lcnum, int lclistnum, int usemask, _Variable *maskvar)
{

  double dum1, dum2;
  double y[2000];
  double ibi[2000];
  int minbin = 5;
  int nbmax = 2000, nbtot;
  int nsrvals, nsrvals_minus, test, foundsofar, dumint1;
  double powerplus, powerminus, dumdbl1, dumdbl2, jdtmp;
  double sumweights, phb1, phb2;
  double tot, rnbtot, *weight, sr_minus;
  double rn, s,t1,f0,p0,ph,ph2,pow,rn1,rn3,s3,rn4,rn5;
  double kkmi, kk, allave, allstddev, allave_minus, allstddev_minus, minbest, qf;
  double *sr_ave, *binned_sr_ave, *binned_sr_sig;
  double in1ph, in2ph;
  int kmi, kma,nb1,nbkma,i,jf,j,k,jn1,jn2,jnb,nb2,nsr,nclippedfreq, *best_id;
  double *p_minus, *sr_ave_minus, *binned_sr_ave_minus, *binned_sr_sig_minus, global_best_sr_ave, global_best_sr_stddev;
  double global_best_sr_ave_inv, global_best_sr_stddev_inv;
  double s1in, s1out, s2in, wsumin, wsumout;
  double magin, magout;
  long double sde_sr_ave, sde_srsqr_ave;
  FILE *outfile, *outfile2;

  double inputT0;

  int nb;
  int in1, in2;
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
      *depth = -1.0;
      *chisqrplus = -1.0;
      *meanmagval = -1.0;
      *fraconenight = -1.0;
      *nt = -1;
      *Nt = -1;
      *Nbefore = -1;
      *Nafter = -1;
      *rednoise = -1.0;
      *whitenoise = -1.0;
      *sigtopink = -1.0;
      *qingress = -1.0;
      *OOTmag = -1.0;
      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
      return(1);
    }
  }

  /***********************************************************/


  tot = t[n-1] - t[0];


  inputT0 = inputTC - inputdur / 2.0;

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
  //*bpow = 0.;

  /**************The following variables are defined for the extension
		 c     of arrays  ibi()  and  y()  [ see below ] ***************/


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


  /* Collect all the output bls parameters */
  in1ph = (inputT0 - t1)/inputper - floor((inputT0 - t1)/inputper);
  in2ph = (inputT0 + inputdur - t1)/inputper - floor((inputT0 + inputdur - t1)/inputper);
  nb = 100*ceil(inputdur/inputper);
  /* Get the transit depth */
  s1in = s1out = wsumin = wsumout = 0.0;
  for(i=0; i < n; i++) {
    ph = u[i]/inputper - floor(u[i]/inputper);
    if((in1ph < in2ph && ph >= in1ph && ph <= in2ph) ||
       (in1ph > in2ph && (ph <= in2ph || ph >= in1ph))) {
      s1in += x[i]*weight[i];
      wsumin += weight[i];
    } else {
      s1out += x[i]*weight[i];
      wsumout += weight[i];
    }
  }
  if(wsumin > 0 && wsumout > 0) {
    magin = s1in/wsumin;
    magout = s1out/wsumout;
    *depth = s1in/wsumin - s1out/wsumout;
  } else {
    *depth = -1.0;
    *chisqrplus = -1.0;
    *meanmagval = -1.0;
    *fraconenight = -1.0;
    *nt = -1;
    *Nt = -1;
    *Nbefore = -1;
    *Nafter = -1;
    *rednoise = -1.0;
    *whitenoise = -1.0;
    *sigtopink = -1.0;
    *qingress = -1.0;
    *OOTmag = -1.0;
    free(weight);
    if(t_mask != NULL) free(t_mask);
    if(x_mask != NULL) free(x_mask);
    if(e_mask != NULL) free(e_mask);
    return 1;
  }
    
  if(fittrap && !fixdepth) {
    *qingress=0.25;
    *OOTmag=*meanmagval;
    dum1 = in1ph;
    dum2 = in2ph;
    dofittrap_amoeba_fixdur(n, t, x, e, inputper, (inputdur/inputper), qingress, in1ph, in2ph, depth, OOTmag);
  } else {
    if(fixdepth)
      *qingress = inputqgress;
    else
      *qingress = 0.;
    *OOTmag = *meanmagval;
  }
  // Be sure to correct for transits past the edge
  *qtran = inputdur/inputper;
  
  s1in = s2in = s1out = wsumin = wsumout = 0.0;
  for(i=0; i < n; i++) {
    ph = u[i]/inputper - floor(u[i]/inputper);
    if((in1ph < in2ph && ph >= in1ph && ph <= in2ph) ||
       (in1ph > in2ph && (ph <= in2ph || ph >= in1ph))) {
      s1in += (x[i]-(magin))*(x[i] - (magin))*weight[i];
      s2in += (x[i]-(*meanmagval))*(x[i] - (*meanmagval))*weight[i];
    } else {
      s1in += (x[i]-(magout))*(x[i] - (magout))*weight[i];
      s2in += (x[i]-(*meanmagval))*(x[i] - (*meanmagval))*weight[i];
    }
  }
  *chisqrplus = (s1in - s2in)*sumweights;

  *fraconenight = getfrac_onenight(n, t, u, v, e, inputper, *depth, *qtran, (t[0] + in1ph*inputper), timezone);

  /* Get the signal to pink noise for the peak */
  getsignaltopinknoiseforgivenblsmodel(n, t, x, e, inputper, *qtran, *depth, in1ph, nt, Nt, Nbefore, Nafter, rednoise, whitenoise, sigtopink, *qingress, *OOTmag, NULL);

  //output the model light curve if asked to
  if(omodel)
    {
      if((outfile2 = fopen(modelname,"w")) == NULL)
	error2(ERR_CANNOTWRITE,modelname);

      f0 = 1./inputper;
      phb1 = (*qingress)*(*qtran);
      phb2 = (*qtran) - phb1;
      //in1ph = (inputT0 - t1)/bper[0] - floor((inputT0 - t1)/bper[0]);

      fprintf(outfile2,"#Time  Mag_obs   Mag_model   Error   Phase\n");
      for(i=0;i<n;i++)
	{
	  ph = (t[i] - inputT0)*f0;
	  ph -= floor(ph);
	  if(ph >= *qtran) {
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
  // Output the phase curve if asked to.
  if(ophcurve)
    {
      if((outfile2 = fopen(ophcurvename,"w")) == NULL)
	error2(ERR_CANNOTWRITE,ophcurvename);

      fprintf(outfile2,"#Phase Mag_model\n");
      ph2 = phmin;
      phb1 = (*qingress)*(*qtran);
      phb2 = *qtran - phb1;
      while(ph2 <= phmax) {
	ph = ph2 + 0.5*(*qtran);
	ph -= floor(ph);
	if(ph >= (*qtran)) {
	  fprintf(outfile2,"%f %f\n",ph2,(*OOTmag));
	}
	else {
	  if(ph >= phb1 && ph <= phb2) {
	    fprintf(outfile2,"%f %f\n",ph2,(*OOTmag)+(*depth));
	  }
	  else if(ph < phb1) {
	    fprintf(outfile2,"%f %f\n",ph2,(*OOTmag)+(*depth)*ph/phb1);
	  }
	  else {
	    fprintf(outfile2,"%f %f\n",ph2,(*OOTmag)+(*depth)*((*qtran) - ph)/phb1);
	  }
	}
	ph2 += phstep;
      }
      fclose(outfile2);
    }

  // Output the JD curve if asked to.
  if(ojdcurve)
    {
      if((outfile2 = fopen(ojdcurvename,"w")) == NULL)
	error2(ERR_CANNOTWRITE,ojdcurvename);

      fprintf(outfile2,"#Time Mag_model Phase\n");
      jdtmp = t[0];
      f0 = 1./inputper;
      phb1 = (*qingress)*(*qtran);
      phb2 = (*qtran) - phb1;
      while(jdtmp <= t[n-1])
	{
	  ph = (jdtmp - inputT0)*f0;
	  ph -= floor(ph);
	  if(ph >= (*qtran)) {
	    ph2 = ph - 0.5*(*qtran);
	    if(ph2 < 0)
	      ph2 += 1.;
	    fprintf(outfile2,"%f %f %f\n", jdtmp, (*OOTmag), ph2);
	  }
	  else {
	    if(ph >= phb1 && ph <= phb2) {
	      ph2 = ph - 0.5*(*qtran);
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f\n", jdtmp, (*OOTmag)+(*depth), ph2);
	    }
	    else if(ph < phb1) {
	      ph2 = ph - 0.5*(*qtran);
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f\n", jdtmp, (*OOTmag)+(*depth)*ph/phb1, ph2);
	    }
	    else {
	      ph2 = ph - 0.5*(*qtran);
	      if(ph2 < 0)
		ph2 += 1.;
	      fprintf(outfile2,"%f %f %f\n", jdtmp, (*OOTmag)+(*depth)*((*qtran) - ph)/phb1, ph2);
	    }
	  }
	  jdtmp += jdstep;
	}
      fclose(outfile2);
    }
  if(correctlc)
    {
      f0 = 1./inputper;
      phb1 = (*qingress)*(*qtran);
      phb2 = (*qtran) - phb1;
      if(usemask) {
	n = n_in;
	t = t_in;
	x = x_in;
	e = e_in;
      }
      for(i=0;i<n;i++)
	{
	  ph = (t[i] - inputT0)*f0;
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

  if(t_mask != NULL) free(t_mask);
  if(x_mask != NULL) free(x_mask);
  if(e_mask != NULL) free(e_mask);
  return(0);
}
