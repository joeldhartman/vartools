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
#define SRVALSSIZE 400000

#define BLS_PEAK_EPSILON 2.0
#define BLS_SR_POLY_ORDER 0
#define BLS_Q_OVERSHOOT 0.2
#define BLS_N_OTFREQS 5

#define BLS_DFT_F0  0.125
#define BLS_DFT_F1  8.0

#define G_CGS_DAYS 498.23382528   /* Newton's gravitional constant in cm^3 g^-1 days^-2 */

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

double getfrac_onenight(int n,double *t,double *u, double *v,double *err,double bper,double depth,double qtran,double bt0,double timezone)
{
  int i, j, k, nnights;
  long night, night0, night1;
  double ph, f0, maxfrac;
  double *chisqrnights, chisqrtot, val;
  f0 = 1./bper;

  /* We have to add 0.5 to go from HJD (where days begin at midnight) to roughly JD (where days begin at noon) */
  timezone = timezone/24.0 + 0.5;
  night0 = (long) (t[0] + timezone);
  night1 = (long) (t[n-1] + timezone);
  nnights = (int) ((night1 - night0) + 1);
  if(nnights < 2) nnights = 2;

  chisqrnights = (double *) malloc(nnights * sizeof(double));
  chisqrtot = 0.;
  for(i=0;i<nnights;i++)
    chisqrnights[i] = 0.;

  for(i=0;i<n;i++)
    {
      ph = (t[i]-bt0)*f0;
      ph -= floor(ph);
      if(ph <= qtran)
	{
	  night = (long) (t[i] + timezone);
	  k = (int) (night - night0);
	  if(k < 0) k = 0;
	  if(k >= nnights) k = nnights - 1;
	  val = v[i]*v[i]/(err[i]*err[i]);
	  chisqrnights[k] += val;
	  chisqrtot += val;
	}
    }
  maxfrac = 0.;
  for(k=0;k<nnights;k++)
    {
      val = chisqrnights[k] / chisqrtot;
      if(val > maxfrac) maxfrac = val;
    }
  free(chisqrnights);
  return(maxfrac);
}


double getclippedsrave(int n, double *sr)
{
  int i, n2;
  double ave1, ave2, ave3;
  ave1 = 0.; ave2 = 0.;
  for(i=0;i<n;i++)
    {
      ave1 += sr[i];
      ave2 += sr[i]*sr[i];
    }
  ave3 = 0.;
  n2 = 0;
  ave1 /= n;
  ave2 = sqrt((double)((ave2 / (double) n) - ave1*ave1));
  for(i=0;i<n;i++)
    {
      if((sr[i] - ave1) < CLIP_FACTOR*ave2)
	{
	  ave3 += sr[i];
	  n2++;
	}
    }
  return (ave3 / (double) n2);
}


double getclippedstddev(int n, double *pow)
{
  int i, n2;
  double ave1, ave2, ave3, ave4;
  ave1 = 0.; ave2 = 0.;
  for(i=0;i<n;i++)
    {
      ave1 += pow[i];
      ave2 += pow[i]*pow[i];
    }
  ave1 /= n;
  ave2 = sqrt((ave2 / n) - (ave1*ave1));
  ave3 = 0.;
  ave4 = 0.;
  n2 = 0;
  for(i=0;i<n;i++)
    {
      if((pow[i] - ave1) < CLIP_FACTOR*ave2)
	{
	  ave3 += pow[i];
	  ave4 += pow[i]*pow[i];
	  n2++;
	}
    }
  ave3 /= n;
  ave4 = sqrt((ave4 / n) - (ave3*ave3));
  return(ave4);
}

void getclippedavestddev(int n, double *pow, double *ave_out, double *stddev_out)
{
  int i, n2;
  double ave1, ave2, ave3, ave4;
  ave1 = 0.; ave2 = 0.;
  for(i=0;i<n;i++)
    {
      ave1 += pow[i];
      ave2 += pow[i]*pow[i];
    }
  ave1 /= n;
  ave2 = sqrt((ave2 / n) - (ave1*ave1));
  ave3 = 0.;
  ave4 = 0.;
  n2 = 0;
  for(i=0;i<n;i++)
    {
      if((pow[i] - ave1) < CLIP_FACTOR*ave2)
	{
	  ave3 += pow[i];
	  ave4 += pow[i]*pow[i];
	  n2++;
	}
    }
  ave3 /= n;
  ave4 = sqrt((ave4 / n) - (ave3*ave3));
  *ave_out = ave3;
  *stddev_out = ave4;
}


double subtract_binnedrms(int N, double *mag, double bintime, double *aveval, int *ngood, double *binmag, double *binsig)
{
  int i, n, jmin, jmax, *ngoodpoints, nclippedlast, nclippedthis;
  double avesum1, avesum2, avesum3, rmsval, v;
  double *sumval1, *sumval2, *sumval3, ave;

  *aveval = -1.;
  if(N > 0)
    {
      if((sumval1 = (double *) malloc(N * sizeof(double))) == NULL ||
	 (sumval2 = (double *) malloc(N * sizeof(double))) == NULL ||
	 (sumval3 = (double *) malloc(N * sizeof(double))) == NULL ||
	 (ngoodpoints = (int *) malloc(N * sizeof(int))) == NULL)
	{
	  fprintf(stderr,"Memory Allocation Error\n");
	  exit(2);
	}

      /* First get the clipped average magnitude and rms*/
      nclippedlast = 0;
      nclippedthis = 0;
      *aveval = 0.;
      rmsval = 1000000.;
      do
	{
	  nclippedlast = nclippedthis;
	  nclippedthis = 0;
	  if(!nclippedlast)
	    {
	      avesum1 = mag[0];
	      avesum2 = mag[0]*mag[0];
	      nclippedthis = 1;
	    }
	  else if(ABS_(mag[0] - *aveval) < CLIP_FACTOR*rmsval)
	    {
	      avesum1 = mag[0];
	      avesum2 = mag[0]*mag[0];
	      nclippedthis = 1;
	    }
	  else
	    {
	      avesum1 = 0.;
	      avesum2 = 0.*0.;
	    }
	  for(i=1;i<N;i++)
	    {
	      if(!nclippedlast || ABS_(mag[i] - *aveval) < CLIP_FACTOR*rmsval)
		{
		  avesum1 += mag[i];
		  avesum2 += mag[i]*mag[i];
		  nclippedthis++;
		}
	    }
	  *aveval = avesum1 / nclippedthis;
	  rmsval = sqrt((avesum2 / nclippedthis) - ((*aveval)*(*aveval)));
	} while (nclippedthis > nclippedlast);

      if(ABS_(mag[0] - *aveval) < CLIP_FACTOR*rmsval)
	{
	  ngoodpoints[0] = 1;
	  sumval1[0] = mag[0];
	  sumval2[0] = mag[0]*mag[0];
	}
      else
	{
	  ngoodpoints[0] = 0;
	  sumval1[0] = 0.;
	  sumval2[0] = 0.;
	}
      for(i=1;i<N;i++)
	{
	  if(ABS_(mag[i] - *aveval) < CLIP_FACTOR*rmsval)
	    {
	      sumval1[i] = sumval1[i-1] + mag[i];
	      sumval2[i] = sumval2[i-1] + mag[i]*mag[i];
	      ngoodpoints[i] = ngoodpoints[i-1] + 1;
	    }
	  else
	    {
	      sumval1[i] = sumval1[i-1];
	      sumval2[i] = sumval2[i-1];
	      ngoodpoints[i] = ngoodpoints[i-1];
	    }
	}

      /* Go through the list find the minimum and maximum times to include via bisection and compute the binned average magnitude and error */
      for(i=0;i<N;i++)
	{
	  jmin = MAX_(0,i - bintime);
	  jmax = MIN_(N-1,i + bintime);
	  if(jmin > 0)
	    {
	      if((v = ngoodpoints[jmax] - ngoodpoints[jmin-1]) > 0)
		{
		  binmag[i] = (sumval1[jmax] - sumval1[jmin-1]) / v;
		  binsig[i] = sqrt((sumval2[jmax] - sumval2[jmin-1]) / v - (binmag[i]*binmag[i]));
		}
	      else
		{
		  binmag[i] = 0.;
		  binsig[i] = 0.;
		}
	    }
	  else
	    {
	      if((v = ngoodpoints[jmax]) > 0)
		{
		  binmag[i] = (sumval1[jmax]) / v;
		  binsig[i] = sqrt((sumval2[jmax])/v - (binmag[i]*binmag[i]));
		}
	      else
		{
		  binmag[i] = 0.;
		  binsig[i] = 0.;
		}
	    }
	}
      avesum1 = 0.;
      avesum2 = 0.;
      avesum3 = 0.;
      n = 0;
      for(i=0;i<N;i++)
	{
	  if(binsig[i] > 0.)
	    {
	      avesum1 += binmag[i];
	      avesum2 += (binmag[i] * binmag[i]);
	      n++;
	    }
	}
      if(n > 0)
	{
	  ave = avesum1 / (double) n;
	  *aveval = ave;
	  *ngood = n;
	  rmsval = sqrt((avesum2 / (double) n) - (ave * ave));
	}
      else
	{
	  *ngood = 0;
	  rmsval = -1.;
	}
      free(sumval1);
      free(sumval2);
      free(sumval3);
      free(ngoodpoints);
    }
  else
    {
      *ngood = 0;
      rmsval = -1.;
    }

  return(rmsval);
}

/* This function cuts out the transit from a light curve */

void BLSCalcCutTransit(int *N, double *t, double *mag, double *sig, double *sig2, double *weight, double P, double T0, double q)
{
  double ph;
  double ph0, ph1, sumweight;
  int i, j, k;
  ph0 = q/2.;
  ph1 = 1. - (q/2.);
  sumweight = 0.;
  for(i=0, j = 0; i < (*N); i++)
    {
      ph = (t[i] - T0)/P - floor((t[i] - T0)/P);
      if(ph > ph0 && ph < ph1)
	{
	  if(j != i)
	    {
	      t[j] = t[i];
	      mag[j] = mag[i];
	      sig[j] = sig[i];
	      sig2[j] = sig[i];
	      weight[j] = weight[i];
	      sumweight += weight[j];
	    }
	  j++;
	}
    }
  *N = j;
  for(i=0;i<j;i++)
    weight[i] = weight[i]/sumweight;
}


void BLSCalcFitTwoTransitsHalfFreq(int n, double *t, double *mag, double *sig, double f, double q, double epoch, double Hin, double Lin, double *Hout, double *L1out, double *L2out, double *delchi2)
{
  /* This function takes the frequency, transit duration and epoch 
     from a given bls run and fits a two-transit light curve with half 
     the frequency, phase separations of 0.5 and fixed, equal, durations. 
     It returns the depths of the primary and secondary transits as well as 
     the statistical significance (from delta-chi2) of the improvement in 
     fitting two transits instead of a single transit */
  
  int i, i1_1, i2_1, i1_2, i2_2, j, k;

  double sum1_1, sum2_1, sum1_2, sum2_2, sum1_3, sum2_3, sum3_1, sum3_3, tmin, P, ph, q2, variance, sumweights;

  tmin = t[0];
  sumweights = 1.0/sig[0]/sig[0];

  /* Get the minimum time, don't assume the light curve is sorted */
  for(i=1;i<n;i++)
    {
      if(t[i] < tmin)
	tmin = t[i];
      sumweights += 1.0/sig[i]/sig[i];
    }
  
  P = 1./f;
  while (epoch > tmin)
    epoch -= P;
  
  /* double the period */
  q = q/2.;
  f = f/2.;

  q2 = 0.5 + q;

  sum1_1 = sum2_1 = sum1_2 = sum2_2 = sum1_3 = sum2_3 = sum3_3 = 0.;

  for(i=0;i<n;i++)
    {
      /* Check if this point is in either of the transits */
      ph = (t[i] - epoch)*f;
      ph -= floor(ph);
      
      if(ph < q)
	{
	  /* This point is in the primary transit */
	  sum1_1 += mag[i]/sig[i]/sig[i];
	  sum2_1 += 1.0/sig[i]/sig[i];
	}
      else if(ph >= 0.5 && ph < q2)
	{
	  /* This point is in the secondary transit */
	  sum1_2 += mag[i]/sig[i]/sig[i];
	  sum2_2 += 1.0/sig[i]/sig[i];
	}
      else
	{
	  /* This point is not in transit */
	  sum1_3 += mag[i]/sig[i]/sig[i];
	  sum2_3 += 1.0/sig[i]/sig[i];
	  sum3_3 += mag[i]*mag[i]/sig[i]/sig[i];
	}
    }

  (*L1out) = sum1_1 / sum2_1;
  if(sum2_2 > 0.)
    (*L2out) = sum1_2 / sum2_2;
  else
    (*L2out)= 0.;
  (*Hout) = sum1_3 / sum2_3;
  Hin = (*Hout);
  Lin = (sum1_1 + sum1_2) / (sum2_1 + sum2_2);

  variance = sum3_3/sum2_3 - (*Hout)*(*Hout);
  
  /* Get the chi2 for the previous and new models */
  sum1_1 = sum2_1 = sum1_2 = sum2_2 = sum3_1 = 0.;

  for(i=0;i<n;i++)
    {
      /* Check if this point is in either of the transits */
      ph = (t[i] - epoch)*f;
      ph -= floor(ph);
      
      if(ph < q)
	{
	  /* This point is in the primary transit */
	  sum1_1 += (mag[i] - (*L1out))*(mag[i] - (*L1out))/sig[i]/sig[i];
	  sum2_1 += (mag[i] - Lin)*(mag[i] - Lin)/sig[i]/sig[i];
	}
      else if(ph >= 0.5 && ph < q2)
	{
	  /* This point is in the secondary transit */
	  sum1_1 += (mag[i] - (*L2out))*(mag[i] - (*L2out))/sig[i]/sig[i];
	  sum2_1 += (mag[i] - Lin)*(mag[i] - Lin)/sig[i]/sig[i];
	}
      else
	{
	  /* This point is not in transit */
	  sum1_1 += (mag[i] - (*Hout))*(mag[i] - (*Hout))/sig[i]/sig[i];
	  sum2_1 += (mag[i] - Hin)*(mag[i] - Hin)/sig[i]/sig[i];
	}
    }
  if(variance > 0.)
    (*delchi2) = (sum1_1 - sum2_1);
  else
    (*delchi2) = 0.;
}


void BLSCalcFDFT(double fmin,
		 double fmax, 
		 int fn, 
		 double* t, 
		 double* m, 
		 double* w,
		 int len,
		 double **freqout,
		 double **amplout,
		 double **phaseout,
		 double **specwin)
{
  double f, ph, df, sumr, sumi;
  int i, j;
  double *sph, *cph, *sdph, *cdph, cp,sp;
  
  if(len <= 0 || fn <= 0)
    return;

  df=(fmax-fmin)/(fn-1);

  if((sph = (double *) malloc(len*sizeof(double))) == NULL ||
     (cph = (double *) malloc(len*sizeof(double))) == NULL ||
     (sdph = (double *) malloc(len*sizeof(double))) == NULL ||
     (cdph = (double *) malloc(len*sizeof(double))) == NULL) {
    fprintf(stderr,"Memory Allocation Error\n");
    exit(3);
  }

  for(j=0; j<len; j++){
    sph[j] = sin(fmin*2*M_PI*t[j]);
    cph[j] = cos(fmin*2*M_PI*t[j]);
    sdph[j] = sin(df*2*M_PI*t[j]);
    cdph[j] = cos(df*2*M_PI*t[j]);
  }

  if(amplout != NULL){
    if((*amplout = (double *) malloc(fn*sizeof(double))) == NULL ||
       (*freqout = (double *) malloc(fn*sizeof(double))) == NULL) {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(3);
    }
    if(phaseout != NULL) {
      if((*phaseout = (double *) malloc(fn*sizeof(double))) == NULL) {
	fprintf(stderr,"Memory Allocation Error\n");
	exit(3);
      }
    }
    for(i=0; i<fn; i++){
      f=fmin+df*i;
      (*freqout)[i] = f;
      sumr=0.0;
      sumi=0.0;
      for(j=0; j<len; j++){
	cp=cph[j];
	sp=sph[j];
	sumr+=w[j]*m[j]*cp;
	sumi+=w[j]*m[j]*sp;
	cph[j]=cp*cdph[j]-sp*sdph[j];
	sph[j]=cp*sdph[j]+sp*cdph[j];
      }
      (*amplout)[i] = 2.0*sqrt(sumr*sumr+sumi*sumi);
      if(phaseout != NULL) {
	(*phaseout)[i] = atan2(sumi, sumr);
      }
    }
  }

  if(specwin != NULL){
    if((*specwin = (double *) malloc(fn*sizeof(double))) == NULL) {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(3);
    }
    if(amplout == NULL) {
      if((*freqout = (double *) malloc(fn*sizeof(double))) == NULL) {
	fprintf(stderr,"Memory Allocation Error\n");
	exit(3);
      }
    }
    for(i=0; i<fn; i++){
      f=fmin+df*i;
      if(amplout == NULL) {
	(*freqout)[i] = f;
      }
      f=f*2*M_PI;
      sumr=0.0;
      sumi=0.0;
      for(j=0; j<len; j++){
	sumr+=cos(f);
	sumi+=sin(f);
      }
      (*specwin)[i] = sqrt(sumr*sumr+sumi*sumi);
    }
  }
  free(sph);
  free(cph);
  free(sdph);
  free(cdph);

}


double BLSCalcLomb(double *t, double *m, double *sig, double *freq, int len, int fn)
{
  double wsum, mavg, w2sum, m2sum, wscalefac, wscalefac2;
  double P, Pmax, sig2, tau, om, num, nom, yy;
  double s,c, sumc, sums, expp, f, cdf, zval, expy, effm, outprob; 
  int i, n, ngood = 0;

  for(i=0; i < len; i++) {
    if(!isnan(m[i]) && !isnan(sig[i]) && sig[i] > 0)
      ngood++;
  }

  Pmax = 0.0;
  for(n=0; n<fn; n++){
    P = gls_oneperiod(t, m, sig, len, (1.0/freq[n]));
    if(P > Pmax) Pmax = P;
  }
  effm = log10(fn);
  if(Pmax == 0.0) outprob = 0.0;
  else {
    zval = (ngood - 3)*(Pmax/(1.0 - Pmax))/2.0;
    expy = -(ngood - 3)*log10((1.0 + (2.0*zval)/(ngood - 3)))/2.0;
    outprob = effm + expy;
    if(outprob > -2.0) {
      expy = pow(((1.0 + (2.0+zval)/(ngood - 3))),(-(ngood-3)/2.0));
      outprob = 1.0 - pow(1.0-expy,fn);
      outprob = log10(outprob);
    }
  }
  return(outprob);
}

void BLSCalcProb(double *mag,
		 int Nt,
		 double *SR,
		 int Nf,
		 double *prob)
{
  double meanval;
  double esterr2, tmp, prmax;
  int N, len, i;

  meanval = 0.;
  esterr2 = 0.;
  for(i=0; i < Nt; i++) {
    meanval += mag[i];
  }
  meanval = meanval / ((double) Nt);
  for(i=0; i < Nt; i++) {
    esterr2 += (mag[i] - meanval)*(mag[i] - meanval);
  }
  esterr2 = esterr2 / ((double) (Nt-1));
  for(i=0; i<Nf; i++){
    tmp=0.5*Nt*SR[i]/esterr2; /*We hold the logarithm of probability, cause it
				is easier to handle*/
    if((!i) || (tmp>prmax)) /*We find the maximum logprob to avoid later the overflow 
			      caused by expressions like exp(1000)*/
      prmax=tmp;
    prob[i] = tmp;
  }

  tmp=0.0;
  for(i=0; i<Nf; i++) 
    tmp+=exp(prob[i]-prmax);/*TODO: could be not so simple 
			     when freq distribution is uneven*/
  tmp=log(tmp)+prmax;
  for(i=0; i<Nf; i++) 
    prob[i]=prob[i]-tmp;
}
		     

double * BLSSRPolysub(double *freq,
		      double *SR,
		      int Nf,
		      int order)
{
  double *SRshift;
  double tmpval;
  if(Nf <= 0) return NULL;

  if((SRshift = (double *) malloc(Nf * sizeof(double))) == NULL) {
    fprintf(stderr,"Memory Allocation Error\n");
    exit(3);
  }
  memcpy(SRshift, SR, Nf*sizeof(double));
  tmpval = fitpoly(Nf, freq, SRshift, NULL, order, 1, NULL, NULL);
  return(SRshift);
}


void BLSPeakArea(double *freq,
		 double *prob,
		 int N,
		 int idx,
		 double epsilon,
		 int *idxlow, 
		 int *idxhigh,
		 double *freqlow,
		 double *freqhigh,
		 double *peakarea,
		 double *peakmean,
		 double *peakdeviance)
{
  double sum, df, moment1, moment2, tmp;
  int len, idx1, idx2;
  int tmpidx;

  sum=0.0;
  moment1=0.0;
  moment2=0.0;
  idx1=idx;
  while((idx1>0) && (prob[idx1]>epsilon)){
    if((idx1+1)<N) df=fabs(freq[idx1+1]-freq[idx1]);
    else df=0.0;
    tmp=exp(prob[idx1]);
    sum+=tmp*df;
    moment1+=tmp*freq[idx1]*df;
    moment2+=tmp*freq[idx1]*freq[idx1]*df;
    idx1--;
  }
  idx2=idx+1;
  while((idx2<N) && (prob[idx2]>epsilon)){
    if((idx2+1)<N) df=fabs(freq[idx2+1]-freq[idx2]);
    else df=0.0;
    tmp=exp(prob[idx2]);
    sum+=tmp*df;
    moment1+=tmp*freq[idx2]*df;
    moment2+=tmp*freq[idx2]*freq[idx2]*df;
    idx2++;
  }
  *idxlow=MIN_(idx, idx1+1);
  *idxhigh=MAX_(idx, idx2-1);
  *freqlow=freq[*idxlow];
  *freqhigh=freq[*idxhigh];
  if(*freqlow > *freqhigh) {
    tmp = *freqlow;
    *freqlow = *freqhigh;
    *freqhigh = tmp;
    tmpidx = *idxlow;
    *idxlow = *idxhigh;
    *idxhigh = tmpidx;
  }
  *peakarea=sum/((*freqhigh)-(*freqlow));
  *peakmean=moment1/sum;
  *peakdeviance=sqrt(moment2/sum-(moment1*moment1)/(sum*sum));

}



void GetExtraBLSParameters1(int n, double *mag, int nf,
			    double *srnoshiftvals, double **srshiftvals, 
			    double *probvals,
			    double *freqarray, _Bls *Bls, int lcnum, int Npeak,
			    int *best_id)
{
  int i, indxlow, indxhigh, ll;
  char *delidx;
  double val1, val2, val3;
  double srsig;
  double meanprob = 0.0;
  double probsig = 0.0;

  /* Get the peak areas and then the sr standard deviation after cutting 
     the peaks */
  if((delidx = (char *) malloc(nf * sizeof(char))) == NULL) {
    fprintf(stderr,"Memory Allocation Error\n");
    exit(3);
  }
  for(i=0; i < nf; i++) delidx[i] = 0;
  BLSCalcProb(mag, n, srnoshiftvals, nf, probvals);
  for(i=0; i < nf; i++) meanprob += probvals[i];
  meanprob = meanprob/((double) nf);
  for(i=0; i < nf; i++) probsig += (probvals[i]-meanprob)*(probvals[i]-meanprob);
  probsig = sqrt(probsig/((double) (nf-1)));
  *srshiftvals = BLSSRPolysub(freqarray,
			      srnoshiftvals,
			      nf,
			      BLS_SR_POLY_ORDER);
  for(i=0; i < Npeak; i++) {
    Bls->logprob[lcnum][i] = probvals[best_id[i]];
    BLSPeakArea(freqarray,
		probvals,
		nf,
		best_id[i],
		BLS_PEAK_EPSILON*probsig + meanprob,
		&indxlow,
		&indxhigh,
		&(Bls->freqlow[lcnum][i]),
		&(Bls->freqhigh[lcnum][i]),
		&(Bls->peakarea[lcnum][i]),
		&(Bls->peakmean[lcnum][i]),
		&(Bls->peakdev[lcnum][i]));
    for(ll=MIN_(indxlow,indxhigh);ll<=MAX_(indxlow,indxhigh);ll++){
      delidx[ll] = 1;
    }
  }
  val1 = 0.0;
  val2 = 0.0;
  ll = 0;
  for(i=0; i < nf; i++) {
    if(!delidx[i]) {
      val1 += (double) ((*srshiftvals)[i]);
      ll++;
    }
  }
  val1 = val1/((double) ll);
  for(i=0; i < nf; i++) {
    if(!delidx[i]) {
      val3 = (double) ((*srshiftvals)[i] - val1);
      val2 += val3*val3;
    }
  }
  if(ll > 1) {
    srsig = sqrt((double) (val2/((double) (ll-1))));
  }
  else {
    srsig = 0.0;
  }
  for(i=0; i < Npeak; i++) {
    Bls->srsig[lcnum][i] = srsig;
    Bls->srshift[lcnum][i] = (*srshiftvals)[best_id[i]];
    Bls->snrextra[lcnum][i] = fabs(Bls->srshift[lcnum][i])/srsig;
  }
  free(delidx);
}

void GetExtraBLSParameters2(int n, double *t, double *mag, double *sig, double P, double q, double depth, double in1_ph, double qingress, double OOTmag, _Bls *Bls, int lcnum, int peaknum)
{
  double mavg, w2sum, w, r, r2, s, s2, m2avg, wsum;
  double t0, t1, f0;
  double ph1, ph2, phb1, phb2, tmp, tmp2;
  double *ottime = NULL, *otphase = NULL, *otmagn = NULL, *otweight = NULL;
  double *otsig = NULL;
  double *tmpphase = NULL;
  double otfreqs[BLS_N_OTFREQS];
  int N_ot = 0;

  int i, j;
 
  double Lvar, Hvar, ressig, dipsig, snr, dsp, lombsig;
  double gezadsp;

  double *otfreqsdft = NULL;
  double *dftampl = NULL;

  double maxamp, maxfreq;

  double H_2tran, L1_2tran, L2_2tran, delchi2_2tran;
  double maxph;

  int Ncut;
  double *tcut;
  double *magcut;
  double *sigcut1; 
  double *magcutkillharm;
  double *sigcut2;
  double *weightcut;
  double Tc, harmamp;

  double *harmA, *harmB;
  double fundA, fundB;
  double harmmean;

  int ngood1, ngood2;
  int binj;

  double *u, *v, bt0sec, bpowsec, depthsec, qtransec;
  int in1sec, in2sec;
  double in1_ph_sec, in2_ph_sec, chisqrplussec, chisqrminussec;
  double meanmagvalsec, fraconenighsec;
  int ntsec, Ntsec, Nbeforesec, Naftersec;
  double rednoisesec, whitenoisesec, sigtopinksec, qingresssec, ootmagsec;
  double srsumsec;
  double ph;

  double *binmag = NULL;
  int *Nbin = NULL;
  int Nbins;
  double binval1, binval2;

  if(n <= 0) return;

  harmA = (double *) malloc(sizeof(double));
  harmB = (double *) malloc(sizeof(double));
  
  Nbins = ceil(1. / q) + 1;
  if(Nbins > 0) {
    if((Nbin = (int *) malloc(Nbins * sizeof(int))) == NULL ||
       (binmag = (double *) malloc(Nbins * sizeof(double))) == NULL) {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(3);
    }
    for(i = 0; i < Nbins; i++) {
      Nbin[i] = 0;
      binmag[i] = 0.;
    }
  }

  if((ottime = (double *) malloc(n * sizeof(double))) == NULL ||
     (otphase = (double *) malloc(n * sizeof(double))) == NULL ||
     (otmagn = (double *) malloc(n * sizeof(double))) == NULL ||
     (otweight = (double *) malloc(n * sizeof(double))) == NULL ||
     (otsig = (double *) malloc(n * sizeof(double))) == NULL ||
     (tmpphase = (double *) malloc(n * sizeof(double))) == NULL ||
     (tcut = (double *) malloc(n * sizeof(double))) == NULL ||
     (magcut = (double *) malloc(n * sizeof(double))) == NULL ||
     (magcutkillharm = (double *) malloc(n * sizeof(double))) == NULL ||
     (sigcut1 = (double *) malloc(n * sizeof(double))) == NULL ||
     (sigcut2 = (double *) malloc(n * sizeof(double))) == NULL ||
     (u = (double *) malloc(n * sizeof(double))) == NULL ||
     (v = (double *) malloc(n * sizeof(double))) == NULL ||
     (weightcut = (double *) malloc(n * sizeof(double))) == NULL) {
    fprintf(stderr,"Memory Allocation Error\n");
    exit(3);
  }

  mavg = 0.0;
  w2sum = 0.0;
  wsum = 0.0;
  for(i=0; i < n; i++) {
    w = 1.0/sig[i]/sig[i];
    wsum += w;
    mavg += w*mag[i];
    w2sum += w*w;
  }
  mavg = mavg/wsum;
  w2sum = w2sum/wsum/wsum;

  r = 0.0;
  r2 = 0.0;
  s = 0.0;
  s2 = 0.0;
  m2avg = 0.0;
  ph1 = in1_ph;
  t0 = t[0];
  t1 = t[0] + ph1*P;
  if(ph1 > 1. + q)
    t1 = t1 - P;
  Tc = t1 + P*q/2.0;
  f0 = 1./P;
  ph2 = q;
  phb1 = qingress*q;
  phb2 = q - qingress*q;
  
  Ncut = n;

  for(i=0; i < n; i++) {
    tcut[i] = t[i];
    magcut[i] = mag[i];
    sigcut1[i] = sig[i];
    sigcut2[i] = sig[i]*sig[i];
    w = 1.0/sig[i]/sig[i]/wsum;
    weightcut[i] = w;
    ph = (t[i] - t1)*f0;
    j = (int) ((t[i] - t1)*f0);
    ph = ph - floor(ph);
    binj = (int) (ph / q);
    Nbin[binj]++;
    binmag[binj] += mag[i];
    tmpphase[i] = ph;
    if(ph < q) {
      if(ph >= phb1 && ph <= phb2) {
	tmp = mag[i] - OOTmag - depth;
      } else if(ph < phb1) {
	tmp = mag[i] - OOTmag - ph*depth/q/qingress;
      } else {
	tmp = mag[i] - OOTmag - (q - ph)*depth/q/qingress;
      }
    } else {
      tmp = mag[i] - OOTmag;
    }
    m2avg += w*tmp*tmp;
    if(ph < q) {
      r += w;
      s += w*tmp;
      r2 += w*w;
      s2 += w*tmp*tmp;
    } else if((ph > q*(1.0 + BLS_Q_OVERSHOOT))&&(ph<(1.0-q*BLS_Q_OVERSHOOT))) {
      ottime[N_ot] = t[i];
      otphase[N_ot] = ph;
      otmagn[N_ot] = mag[i];
      otweight[N_ot] = w;
      otsig[N_ot] = sig[i];
      N_ot++;
    }
  }
  Lvar = (r*s2-s*s)/(r*r-r2);
  Hvar=((m2avg-s2)*(1.0-r)-s*s)/((1.0-r)*(1.0-r) - (w2sum-r2));
  ressig = sqrt(Lvar*r+Hvar*(1.0-r));
  dipsig = sqrt(Lvar*r2/(r*r)+Hvar*(w2sum-r2)/((1.0-r)*(1.0-r)));
  dsp=fabs(depth/dipsig);

  for(i=0; i < BLS_N_OTFREQS; i++) {
    otfreqs[i] = f0*(i+1);
  }

  lombsig = BLSCalcLomb(ottime, otmagn, otsig, otfreqs, N_ot, BLS_N_OTFREQS);

  BLSCalcFDFT(BLS_DFT_F0, BLS_DFT_F0, n, otphase, otmagn, otweight, N_ot, 
	      &otfreqsdft, 
	      &dftampl, NULL, NULL);
  
  if(dftampl != NULL && otfreqsdft != NULL) {
    maxamp = dftampl[0];
    maxfreq = otfreqsdft[0];
    for(i=0; i < n; i++) {
      if(dftampl[i] > maxamp) {
	maxamp = dftampl[i];
	maxfreq = otfreqsdft[i];
      }
    }
    gezadsp = fabs(depth)/sqrt(Hvar+maxamp*maxamp);
  } else {
    maxamp = 0.0;
    maxfreq = 0.0;
    gezadsp = 0.0;
  }

  Bls->ressig[lcnum][peaknum] = ressig;
  Bls->dipsig[lcnum][peaknum] = dipsig;
  Bls->dsp[lcnum][peaknum] = dsp;
  Bls->dspg[lcnum][peaknum] = fabs(depth)/sqrt(Hvar);
  Bls->lomblog[lcnum][peaknum] = lombsig;
  Bls->gezadsp[lcnum][peaknum] = gezadsp;
  Bls->ootsig[lcnum][peaknum] = sqrt(Hvar);
  Bls->trsig[lcnum][peaknum] = sqrt(Lvar);
  Bls->ootdftf[lcnum][peaknum] = maxfreq;
  Bls->ootdfta[lcnum][peaknum] = maxamp;

  /* Now calculate the parameters generated by the program blspostprocess */

  /* Get the binned signal to noise */
  binval1 = 0.; binval2 = 0.;
  for(i=1, j=0; i < Nbins; i++)
    {
      if(Nbin[i] > 0) {
	binmag[i] /= Nbin[i];
	binval1 += binmag[i]*binmag[i];
	binval2 += binmag[i];
	j++;
      }
    }

  if(binval1 > 0 && binval1 != binval2*binval2) {
    Bls->binsignaltonoise[lcnum][peaknum] = depth/sqrt(binval1/j - (binval2*binval2/j/j));
  } else {
    Bls->binsignaltonoise[lcnum][peaknum] = -1.;
  }
  
  /* Try fitting 2 transits with double the period */
  BLSCalcFitTwoTransitsHalfFreq(n, t, mag, sig, f0, q, t1, OOTmag, 
				OOTmag+depth, 
				&H_2tran, &L1_2tran, &L2_2tran, &delchi2_2tran);

  Bls->depth1_2tran[lcnum][peaknum] = L1_2tran - H_2tran;
  Bls->depth2_2tran[lcnum][peaknum] = L2_2tran - H_2tran;
  Bls->delchi2_2tran[lcnum][peaknum] = delchi2_2tran;

  /* Find the length in phase of the maximum gap in the phase-folded lc */
  mysort1(n,tmpphase);
  maxph = 1.0 + tmpphase[0] - tmpphase[n-1];
  for(i=1; i < n; i++) {
    tmp = tmpphase[i] - tmpphase[i-1];
    if(tmp > maxph)
      maxph = tmp;
  }

  Bls->maxphasegap[lcnum][peaknum] = maxph;

  /* Cut out the transit, fit a harmonic series, and run a
     fixed-period BLS search */
  BLSCalcCutTransit(&Ncut, tcut, magcut, sigcut1, sigcut2, weightcut, P, Tc, q);

  memcpy(magcutkillharm, magcut, Ncut*sizeof(double));

  dokillharms(Ncut, tcut, magcutkillharm, sigcut1, 1, &P, 0, 1, 
	      NULL, NULL, &harmA, &harmB, &fundA, &fundB, &harmmean, 0, NULL, 
	      &harmamp, 0.0, KILLHARM_OUTTYPE_DEFAULT, -1.0);

  Bls->harmdeltachi2[lcnum][peaknum] =
    chi2(Ncut, tcut, magcutkillharm, sigcut1, &tmp, &ngood1, 0, NULL, 0, 0) -
    chi2(Ncut, tcut, magcut, sigcut1, &tmp2, &ngood2, 0, NULL, 0, 0);

  Bls->harmamp[lcnum][peaknum] = harmamp;

  Bls->harmA[lcnum][peaknum] = harmA[0];
  Bls->harmB[lcnum][peaknum] = harmB[0];
  Bls->fundA[lcnum][peaknum] = fundA;
  Bls->fundB[lcnum][peaknum] = fundB;
  Bls->harmmean[lcnum][peaknum] = harmmean;

  /* Run fix period BLS */
  if(!Bls->rflag) {
    eeblsfixper(Ncut, tcut, magcut, sigcut1, u, v, Bls->nbins_val[lcnum], Bls->qmin_val[lcnum],
		Bls->qmax_val[lcnum], &P, &bt0sec, &bpowsec, &depthsec, &qtransec, &in1sec,
		&in2sec,
		&in1_ph_sec, &in2_ph_sec, &chisqrplussec,
		&chisqrminussec, &meanmagvalsec, Bls->timezone, &fraconenighsec,
		0, NULL, 0, 0, &ntsec, &Ntsec, &Nbeforesec, &Naftersec,
		&rednoisesec, &whitenoisesec, &sigtopinksec, 1,
		&qingresssec, &ootmagsec, &srsumsec, 0, 0, 0, NULL);
  } else {
    eeblsfixper_rad(Ncut, tcut, magcut, sigcut1, u, v, Bls->nbins_val[lcnum], Bls->rmin_val[lcnum],
		Bls->rmax_val[lcnum], &P, &bt0sec, &bpowsec, &depthsec, &qtransec, &in1sec,
		&in2sec,
		&in1_ph_sec, &in2_ph_sec, &chisqrplussec,
		&chisqrminussec, &meanmagvalsec, Bls->timezone, &fraconenighsec,
		0, NULL, 0, 0, &ntsec, &Ntsec, &Nbeforesec, &Naftersec,
		&rednoisesec, &whitenoisesec, &sigtopinksec, 1,
		&qingresssec, &ootmagsec, &srsumsec, 0, 0, 0, NULL);
  }
  Bls->sr_sec[lcnum][peaknum] = bpowsec;
  Bls->srsum_sec[lcnum][peaknum] = srsumsec;
  Bls->q_sec[lcnum][peaknum] = qtransec;
  Bls->epoch_sec[lcnum][peaknum] = bt0sec;
  Bls->H_sec[lcnum][peaknum] = ootmagsec;
  Bls->L_sec[lcnum][peaknum] = ootmagsec + depthsec;
  Bls->depth_sec[lcnum][peaknum] = depthsec;
  Bls->nt_sec[lcnum][peaknum] = ntsec;
  Bls->Nt_sec[lcnum][peaknum] = Ntsec;
  Bls->sigtopink_sec[lcnum][peaknum] = sigtopinksec;
  Bls->deltachi2transit_sec[lcnum][peaknum] = chisqrplussec;
  Bls->phaseoffset_sec[lcnum][peaknum] = (bt0sec - Tc)/P - 
    floor((bt0sec - Tc)/P);

  /* Get the binned signal to noise for the secondary transit */
  if(Nbin != NULL)
    free(Nbin);
  if(binmag != NULL)
    free(binmag);

  Nbins = ceil(1. / qtransec) + 1;
  if(Nbins > 0) {
    if((Nbin = (int *) malloc(Nbins * sizeof(int))) == NULL ||
       (binmag = (double *) malloc(Nbins * sizeof(double))) == NULL) {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(3);
    }
    for(i = 0; i < Nbins; i++) {
      Nbin[i] = 0;
      binmag[i] = 0.;
    }
  }

  t1 = bt0sec - qtransec*P/2.0;
  
  for(i=0; i < Ncut; i++) {
    ph = (tcut[i] - t1)*f0;
    ph = ph - floor(ph);
    binj = (int) (ph / qtransec);
    Nbin[binj]++;
    binmag[binj] += magcut[i];
  }
  
  
  /* Get the binned signal to noise */
  binval1 = 0.; binval2 = 0.;
  for(i=1, j=0; i < Nbins; i++)
    {
      if(Nbin[i] > 0) {
	binmag[i] /= Nbin[i];
	binval1 += binmag[i]*binmag[i];
	binval2 += binmag[i];
	j++;
      }
    }

  if(binval1 > 0 && binval1 != binval2*binval2) {
    Bls->binsignaltonoise_sec[lcnum][peaknum] = depth/sqrt(binval1/j - (binval2*binval2/j/j));
  } else {
    Bls->binsignaltonoise_sec[lcnum][peaknum] = -1.0;
  }



  /* Stopped Working here; need to fill out all the extra BLS parameters, see lc_bls_calc_anal and lc_bls_calc_postproc functions in HATpipe/source/lc/analyse/bls.c */
  if(otfreqsdft != NULL)
    free(otfreqsdft);
  if(dftampl != NULL)
    free(dftampl);
  free(ottime);
  free(otphase);
  free(otmagn);
  free(otweight);
  free(otsig);
  free(tmpphase);
  free(tcut);
  free(magcut);
  free(magcutkillharm);
  free(sigcut1);
  free(sigcut2);
  free(u);
  free(v);
  free(weightcut);
  free(harmA);
  free(harmB);
  if(binmag != NULL)
    free(binmag);
  if(Nbin != NULL)
    free(Nbin);
  return;
}



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
c     nf   = number of frequency points in which the spectrum is computed
c     fmin = minimum frequency (MUST be > 0)
c     df   = frequency step
c     nb   = number of bins in the folded time series at any test period
c     qmi  = minimum fractional transit length to be tested
c     qma  = maximum fractional transit length to be tested
c
c     Output parameters:
c     ~~~~~~~~~~~~~~~~~~
c
c     p    = array {p(i)}, containing the values of the BLS spectrum
c            at the i-th frequency value -- the frequency values are
c            computed as  f = fmin + (i-1)*df
c     bper = period at the highest peak in the frequency spectrum
c     bpow = value of {p(i)} at the highest peak
c     depth= depth of the transit at   *bper*
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
c     -- *fmin* MUST be greater than  *1/total time span*
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

int eebls(int n_in, double *t_in, double *x_in, double *e_in, double *u, double *v, int nf, double fmin, double df, int nb, double qmi, double qma, double *p, int Npeak, double *bper, double *bt0, double *bpow, double *sde, double *snval, double *depth, double *qtran, int *in1, int *in2, double *in1_ph, double *in2_ph, double *chisqrplus, double *chisqrminus, double *bperpos, double *meanmagval, double timezone, double *fraconenight, int operiodogram, char *outname, int omodel, char *modelname, int correctlc,int ascii,int *nt, int *Nt, int *Nbefore, int *Nafter, double *rednoise, double *whitenoise, double *sigtopink, int fittrap, double *qingress, double *OOTmag, int ophcurve, char *ophcurvename, double phmin, double phmax, double phstep, int ojdcurve, char *ojdcurvename, double jdstep, int nobinnedrms, int freq_step_type, int adjust_qmin_mindt, int reduce_nb, int reportharmonics, _Bls *Bls, int lcnum, int lclistnum, int usemask, _Variable *maskvar)
{
  double *y;
  double *ibi;
#ifdef PARALLEL
  int firsttime = 0;
  double *srvals = NULL, *srvals_minus = NULL;
#else
  static int firsttime = 0;
  static double *srvals, *srvals_minus;
#endif
  int minbin = 5;
  int nbmax, nbtot;
  int nsrvals, nsrvals_minus, test, foundsofar, dumint1;
  int indxlow, indxhigh;
  double powerplus, powerminus, bpowminus, dumdbl1, dumdbl2, jdtmp;
  double sumweights, phb1, phb2;
  double tot, rnbtot, *weight, sr_minus;
  double rn, s,t1,f0,p0,ph,ph2,pow,rn1,rn3,s3,rn4,rn5;
  double kkmi, kk, allave, allstddev, allave_minus, allstddev_minus, *qtran_array, *depth_array, minbest;
  double *sr_ave, *binned_sr_ave, *binned_sr_sig;
  int kmi, kma,nb1,nbkma,i,ll,jf,j,k,jn1,jn2,jnb,nb2,nsr,nclippedfreq, *in1_array, *in2_array, *best_id, nbsave;
  char *delidx;
  double *p_minus, *bper_array, *sr_ave_minus, *binned_sr_ave_minus, *binned_sr_sig_minus, global_best_sr_ave, global_best_sr_stddev;
  double *freqarray = NULL;
  double global_best_sr_ave_inv, global_best_sr_stddev_inv;
  double kmisave, kkmisave, kmasave, mindt, qmi_test;
  long double sde_sr_ave, sde_srsqr_ave;
  FILE *outfile, *outfile2;
  long double val1, val2, val3;
  double *probvals = NULL;
  double *srshiftvals = NULL;
  double *srnoshiftvals = NULL;
  int *ntvptr = NULL;
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
	  in1[j] = -1;
	  in2[j] = -1;
	  in1_ph[j] = -1.;
	  in2_ph[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	  if(Bls->extraparams) {
	    Bls->srsum[lcnum][j] = -1.;
	    Bls->ressig[lcnum][j] = -1.;
	    Bls->dipsig[lcnum][j] = -1.;
	    Bls->srshift[lcnum][j] = -1.;
	    Bls->srsig[lcnum][j] = -1.;
	    Bls->snrextra[lcnum][j] = -1.;
	    Bls->dsp[lcnum][j] = -1.;
	    Bls->dspg[lcnum][j] = -1.;
	    Bls->freqlow[lcnum][j] = -1.;
	    Bls->freqhigh[lcnum][j] = -1.;
	    Bls->logprob[lcnum][j] = -1.;
	    Bls->peakarea[lcnum][j] = -1.;
	    Bls->peakmean[lcnum][j] = -1.;
	    Bls->peakdev[lcnum][j] = -1.;
	    Bls->lomblog[lcnum][j] = -1.;
	    Bls->ntv[lcnum][j] = 0;
	    Bls->gezadsp[lcnum][j] = -1.;
	    Bls->ootsig[lcnum][j] = -1.;
	    Bls->trsig[lcnum][j] = -1.;
	    Bls->ootdftf[lcnum][j] = -1.;
	    Bls->ootdfta[lcnum][j] = -1.;
	    Bls->binsignaltonoise[lcnum][j] = -1.;
	    Bls->maxphasegap[lcnum][j] = -1.;
	    Bls->depth1_2tran[lcnum][j] = -1.;
	    Bls->depth2_2tran[lcnum][j] = -1.;
	    Bls->delchi2_2tran[lcnum][j] = -1.;
	    Bls->sr_sec[lcnum][j] = -1.;
	    Bls->srsum_sec[lcnum][j] = -1.;
	    Bls->q_sec[lcnum][j] = -1.;
	    Bls->epoch_sec[lcnum][j] = -1.;
	    Bls->H_sec[lcnum][j] = -1.;
	    Bls->L_sec[lcnum][j] = -1.;
	    Bls->depth_sec[lcnum][j] = -1.;
	    Bls->nt_sec[lcnum][j] = 0;
	    Bls->Nt_sec[lcnum][j] = 0;
	    Bls->sigtopink_sec[lcnum][j] = -1.;
	    Bls->deltachi2transit_sec[lcnum][j] = -1.;
	    Bls->binsignaltonoise_sec[lcnum][j] = -1.;
	    Bls->phaseoffset_sec[lcnum][j] = -1.;
	    Bls->harmmean[lcnum][j] = -1.;
	    Bls->fundA[lcnum][j] = -1.;
	    Bls->fundB[lcnum][j] = -1.;
	    Bls->harmA[lcnum][j] = -1.;
	    Bls->harmB[lcnum][j] = -1.;
	    Bls->harmamp[lcnum][j] = -1.;
	    Bls->harmdeltachi2[lcnum][j] = -1.;
	  }
	}
      *bperpos = -1.;
      *chisqrminus = 999999.;
      *meanmagval = -1.;

      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
      return 1;
    }
  }

  nbmax = 2*nb;

  if(firsttime == 0 && !nobinnedrms)
    {
      if((srvals = (double *) malloc(SRVALSSIZE * sizeof(double))) == NULL ||
	 (srvals_minus = (double *) malloc(SRVALSSIZE * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      firsttime = 1;
    }

  /***********************************************************/

  if((sr_ave = (double *) malloc(nf * sizeof(double))) == NULL ||
     (y = (double *) malloc(nbmax * sizeof(double))) == NULL ||
     (ibi = (double *) malloc(nbmax * sizeof(double))) == NULL ||
     (binned_sr_ave = (double *) malloc(nf * sizeof(double))) == NULL ||
     (binned_sr_sig = (double *) malloc(nf * sizeof(double))) == NULL ||
     (in1_array = (int *) malloc(nf * sizeof(int))) == NULL ||
     (in2_array = (int *) malloc(nf * sizeof(int))) == NULL ||
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
  if(Bls->extraparams) {
    if((freqarray = (double *) malloc(nf * sizeof(double))) == NULL ||
       (probvals = (double *) malloc(nf * sizeof(double))) == NULL ||
       (srnoshiftvals = (double *) malloc(nf * sizeof(double))) == NULL)
      {
	fprintf(stderr,"Memory Allocation Error\n");
	exit(3);
      }
  }

  if(nb > nbmax) {
    error(ERR_BLSNBMAX);
  }
  tot = t[n-1] - t[0];
  if(fmin < 1./tot) {
    error(ERR_BLSFMINTOOSMALL);
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

  if(adjust_qmin_mindt) {
    kmisave = kmi;
    kkmisave = kkmi;
    kmasave = kma;
    mindt = 0;
    if(n > 1) {mindt = t[1] - t[0];}
    for(i=2; i < n; i++) {
      if(t[i]-t[i-1] < mindt) mindt = t[i] - t[i-1];
    }
    if(reduce_nb) {
      nbsave = nb;
    }
  }
    

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
      if(!freq_step_type) {
	f0=fmin+df*((double)jf);
	p0=1./f0;
      } else if(freq_step_type == VARTOOLS_FREQSTEPTYPE_PERIOD) {
	p0 = (1./fmin) - df*((double)jf);
	f0 = 1./p0;
      } else if(freq_step_type == VARTOOLS_FREQSTEPTYPE_LOGPERIOD) {
	f0 = exp(log(fmin) + df*((double)jf));
	p0=1./f0;
      }

      if(adjust_qmin_mindt) {
	qmi_test = mindt*f0;
	if(qmi_test > 1) qmi_test = 1.0;
	if(qmi_test > qmi) {
	  if(reduce_nb) {
	    nb = MIN_(nbsave, ceil(1./(0.5*qmi_test)));
	  }
	  kmi = (int) (qmi_test*(double)nb);
	  if(kmi < 1) kmi = 1;
	  if(qmi_test > qma)
	    kma = ((int) (qmi_test*(double)nb)) + 1;
	  else if(reduce_nb)
	    kma = ((int) (qma*(double)nb)) + 1;
	  nb1 = nb;
	  nbkma = nb+kma;
	  kkmi = qmi;
	  if(kkmi < (double) minbin / rn) kkmi = (double) minbin / rn;
	} else {
	  kmi = kmisave; kkmi = kkmisave;
	  kma = kmasave;
	  if(reduce_nb)
	    nb = nbsave;
	  nb1 = nb;
	  nbkma = nb+kma;
	}
      }


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
      nsrvals = 0;
      nsrvals_minus = 0;
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
		  if(s > 0.)
		    {
		      if(!nobinnedrms) {
			srvals[nsrvals] = sqrt(pow);
			nsrvals++;
		      }
		      if(pow >= powerplus)
			{
			  powerplus = pow;
			  jn1 = i;
			  jn2 = j;
			  rn3 = rn1;
			  rn5 = rn4;
			  s3 = s;
			}
		    }
		  else if(s < 0.)
		    {
		      if(!nobinnedrms) {
			srvals_minus[nsrvals_minus] = sqrt(pow);
			nsrvals_minus++;
		      }
		      if(pow >= powerminus)
			{
			  powerminus = pow;
			}
		    }
		}
	    }
	}
      // Find the average value of the srvals
      if(!nobinnedrms) {
	sr_ave[jf] = getclippedsrave(nsrvals,srvals);
	sr_ave_minus[jf] = getclippedsrave(nsrvals_minus,srvals_minus);
      }
      powerplus = sqrt(powerplus);
      sde_sr_ave += powerplus;
      sde_srsqr_ave += powerplus*powerplus;
      p[jf] = powerplus;
      powerminus = sqrt(powerminus);
      p_minus[jf] = powerminus;
      //sr_ave += powerplus;
      //srsqr_ave += powerplus*powerplus;
      nsr++;
      in1_array[jf] = jn1;
      in2_array[jf] = jn2;
      qtran_array[jf] = (double)(jn2 - jn1 + 1)/(double)nb;
      depth_array[jf] = powerplus/sqrt(rn3*(1.-rn3));
      bper_array[jf] = p0;
      if(freqarray != NULL)
	freqarray[jf] = f0;
      /*      if(powerplus >= *bpow)
	{
	  *bpow = powerplus;
	  *in1 = jn1;
	  *in2 = jn2;
	  *qtran = (double)(jn2 - jn1 + 1)/(double)nb;
	  *depth = powerplus/sqrt(rn3*(1.-rn3));
	  *bper = p0;
	}
      if(powerminus >= bpowminus)
	{
	  bpowminus = powerminus;
	  *bperpos = p0;
	  }*/
    }
  if(Bls->extraparams) {
    /* Save a copy of the uncorrected SR values if we need to compute the
       extra parameters; Note that the code used for the extra parameters
       is derived from the lc/blsanal tool used by HATNet, in which SR is
       treated as the square of the SR values computed by this BLS code.
    */
    memcpy(srnoshiftvals, p, nf*sizeof(double));
    for(jf=0;jf<nf;jf++) {
      srnoshiftvals[jf] = srnoshiftvals[jf]*srnoshiftvals[jf];
    }
  }
  if(!nobinnedrms) {
    allstddev = subtract_binnedrms(nf, sr_ave, BIN_FACTOR, &allave, &nclippedfreq, binned_sr_ave, binned_sr_sig);
  }
  else {
    getclippedavestddev(nf,p,&global_best_sr_ave,&global_best_sr_stddev);
    nclippedfreq = nf;
  }

  /* Now let's find the peaks in the periodogram, first convert the periodogram from SR to SN ratio */

  if(nclippedfreq > Npeak)
    {
      if(!nobinnedrms) {
	for(i=0;i<nf;i++)
	  {
	    if(binned_sr_ave[i] > 0.)
	      {
		p[i] = (p[i] - binned_sr_ave[i]) / allstddev;
	      /*	      if(p[i] > *bpow)
	        {
		  *bpow = p[i];
		  sr_plus = p[i]*allstddev + binned_sr_ave[i];
		  *in1 = in1_array[i];
		  *in2 = in2_array[i];
		  *qtran = qtran_array[i];
		  *depth = depth_array[i];
		  *bper = bper_array[i];
		  }*/
	      }
	    else
	      p[i] = 0.;
	  }
      }
      else {
	for(i=0; i<nf; i++) {
	  p[i] = (p[i] - global_best_sr_ave)/global_best_sr_stddev;
	}
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
	  in1[j] = -1;
	  in2[j] = -1;
	  in1_ph[j] = -1.;
	  in2_ph[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	  if(Bls->extraparams) {
	    Bls->srsum[lcnum][j] = -1.;
	    Bls->ressig[lcnum][j] = -1.;
	    Bls->dipsig[lcnum][j] = -1.;
	    Bls->srshift[lcnum][j] = -1.;
	    Bls->srsig[lcnum][j] = -1.;
	    Bls->snrextra[lcnum][j] = -1.;
	    Bls->dsp[lcnum][j] = -1.;
	    Bls->dspg[lcnum][j] = -1.;
	    Bls->freqlow[lcnum][j] = -1.;
	    Bls->freqhigh[lcnum][j] = -1.;
	    Bls->logprob[lcnum][j] = -1.;
	    Bls->peakarea[lcnum][j] = -1.;
	    Bls->peakmean[lcnum][j] = -1.;
	    Bls->peakdev[lcnum][j] = -1.;
	    Bls->lomblog[lcnum][j] = -1.;
	    Bls->ntv[lcnum][j] = 0;
	    Bls->gezadsp[lcnum][j] = -1.;
	    Bls->ootsig[lcnum][j] = -1.;
	    Bls->trsig[lcnum][j] = -1.;
	    Bls->ootdftf[lcnum][j] = -1.;
	    Bls->ootdfta[lcnum][j] = -1.;
	    Bls->binsignaltonoise[lcnum][j] = -1.;
	    Bls->maxphasegap[lcnum][j] = -1.;
	    Bls->depth1_2tran[lcnum][j] = -1.;
	    Bls->depth2_2tran[lcnum][j] = -1.;
	    Bls->delchi2_2tran[lcnum][j] = -1.;
	    Bls->sr_sec[lcnum][j] = -1.;
	    Bls->srsum_sec[lcnum][j] = -1.;
	    Bls->q_sec[lcnum][j] = -1.;
	    Bls->epoch_sec[lcnum][j] = -1.;
	    Bls->H_sec[lcnum][j] = -1.;
	    Bls->L_sec[lcnum][j] = -1.;
	    Bls->depth_sec[lcnum][j] = -1.;
	    Bls->nt_sec[lcnum][j] = 0;
	    Bls->Nt_sec[lcnum][j] = 0;
	    Bls->sigtopink_sec[lcnum][j] = -1.;
	    Bls->deltachi2transit_sec[lcnum][j] = -1.;
	    Bls->binsignaltonoise_sec[lcnum][j] = -1.;
	    Bls->phaseoffset_sec[lcnum][j] = -1.;
	    Bls->harmmean[lcnum][j] = -1.;
	    Bls->fundA[lcnum][j] = -1.;
	    Bls->fundB[lcnum][j] = -1.;
	    Bls->harmA[lcnum][j] = -1.;
	    Bls->harmB[lcnum][j] = -1.;
	    Bls->harmamp[lcnum][j] = -1.;
	    Bls->harmdeltachi2[lcnum][j] = -1.;
	  }
	}
      *bperpos = -1.;
      *chisqrminus = 999999.;
      *meanmagval = -1.;

      free(weight);
      free(y);
      free(ibi);
      free(best_id);
      free(sr_ave);
      free(binned_sr_ave);
      free(binned_sr_sig);
      free(in1_array);
      free(in2_array);
      free(qtran_array);
      free(depth_array);
      free(bper_array);
      free(sr_ave_minus);
      free(binned_sr_ave_minus);
      free(binned_sr_sig_minus);
      free(p_minus);
      if(freqarray != NULL)
	free(freqarray);
      if(probvals != NULL)
	free(probvals);
      if(srshiftvals != NULL)
	free(srshiftvals);
      if(srnoshiftvals != NULL)
	free(srnoshiftvals);
#ifdef PARALLEL
      if(srvals != NULL) free(srvals);
      if(srvals_minus != NULL) free(srvals_minus);
#endif
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
	      if((!reportharmonics && !isDifferentPeriods(MIN_(bper[j],bper_array[i]),MAX_(bper[j],bper_array[i]),tot)) || (reportharmonics && isDifferentPeriodsDontCheckHarmonics(MIN_(bper[j],bper_array[i]),MAX_(bper[j],bper_array[i]),tot)))
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
		  if((!reportharmonics && !isDifferentPeriods(MIN_(bper[j],bper_array[i]),MAX_(bper[j],bper_array[i]),tot)) || (reportharmonics && !isDifferentPeriodsDontCheckHarmonics(MIN_(bper[j],bper_array[i]),MAX_(bper[j],bper_array[i]),tot)))
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
      /* We have a few peaks, but not Npeak of them */
      mysort3_int(foundsofar,snval,bper,best_id);
      for(j=foundsofar;j<Npeak;j++)
	{
	  /* Put -1 for the remaining peaks */
	  bper[j] = -1.;
	  snval[j] = -1.;
	  bpow[j] = -1.;
	  bt0[j] = -1.;
	  in1[j] = -1;
	  in2[j] = -1;
	  in1_ph[j] = -1.;
	  in2_ph[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	  if(Bls->extraparams) {
	    Bls->srsum[lcnum][j] = -1.;
	    Bls->ressig[lcnum][j] = -1.;
	    Bls->dipsig[lcnum][j] = -1.;
	    Bls->srshift[lcnum][j] = -1.;
	    Bls->srsig[lcnum][j] = -1.;
	    Bls->snrextra[lcnum][j] = -1.;
	    Bls->dsp[lcnum][j] = -1.;
	    Bls->dspg[lcnum][j] = -1.;
	    Bls->freqlow[lcnum][j] = -1.;
	    Bls->freqhigh[lcnum][j] = -1.;
	    Bls->logprob[lcnum][j] = -1.;
	    Bls->peakarea[lcnum][j] = -1.;
	    Bls->peakmean[lcnum][j] = -1.;
	    Bls->peakdev[lcnum][j] = -1.;
	    Bls->lomblog[lcnum][j] = -1.;
	    Bls->ntv[lcnum][j] = 0;
	    Bls->gezadsp[lcnum][j] = -1.;
	    Bls->ootsig[lcnum][j] = -1.;
	    Bls->trsig[lcnum][j] = -1.;
	    Bls->ootdftf[lcnum][j] = -1.;
	    Bls->ootdfta[lcnum][j] = -1.;
	    Bls->binsignaltonoise[lcnum][j] = -1.;
	    Bls->maxphasegap[lcnum][j] = -1.;
	    Bls->depth1_2tran[lcnum][j] = -1.;
	    Bls->depth2_2tran[lcnum][j] = -1.;
	    Bls->delchi2_2tran[lcnum][j] = -1.;
	    Bls->sr_sec[lcnum][j] = -1.;
	    Bls->srsum_sec[lcnum][j] = -1.;
	    Bls->q_sec[lcnum][j] = -1.;
	    Bls->epoch_sec[lcnum][j] = -1.;
	    Bls->H_sec[lcnum][j] = -1.;
	    Bls->L_sec[lcnum][j] = -1.;
	    Bls->depth_sec[lcnum][j] = -1.;
	    Bls->nt_sec[lcnum][j] = 0;
	    Bls->Nt_sec[lcnum][j] = 0;
	    Bls->sigtopink_sec[lcnum][j] = -1.;
	    Bls->deltachi2transit_sec[lcnum][j] = -1.;
	    Bls->binsignaltonoise_sec[lcnum][j] = -1.;
	    Bls->phaseoffset_sec[lcnum][j] = -1.;
	    Bls->harmmean[lcnum][j] = -1.;
	    Bls->fundA[lcnum][j] = -1.;
	    Bls->fundB[lcnum][j] = -1.;
	    Bls->harmA[lcnum][j] = -1.;
	    Bls->harmB[lcnum][j] = -1.;
	    Bls->harmamp[lcnum][j] = -1.;
	    Bls->harmdeltachi2[lcnum][j] = -1.;
	  }
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
	  in1[j] = -1;
	  in2[j] = -1;
	  in1_ph[j] = -1.;
	  in2_ph[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	  if(Bls->extraparams) {
	    Bls->srsum[lcnum][j] = -1.;
	    Bls->ressig[lcnum][j] = -1.;
	    Bls->dipsig[lcnum][j] = -1.;
	    Bls->srshift[lcnum][j] = -1.;
	    Bls->srsig[lcnum][j] = -1.;
	    Bls->snrextra[lcnum][j] = -1.;
	    Bls->dsp[lcnum][j] = -1.;
	    Bls->dspg[lcnum][j] = -1.;
	    Bls->freqlow[lcnum][j] = -1.;
	    Bls->freqhigh[lcnum][j] = -1.;
	    Bls->logprob[lcnum][j] = -1.;
	    Bls->peakarea[lcnum][j] = -1.;
	    Bls->peakmean[lcnum][j] = -1.;
	    Bls->peakdev[lcnum][j] = -1.;
	    Bls->lomblog[lcnum][j] = -1.;
	    Bls->ntv[lcnum][j] = 0;
	    Bls->gezadsp[lcnum][j] = -1.;
	    Bls->ootsig[lcnum][j] = -1.;
	    Bls->trsig[lcnum][j] = -1.;
	    Bls->ootdftf[lcnum][j] = -1.;
	    Bls->ootdfta[lcnum][j] = -1.;
	    Bls->binsignaltonoise[lcnum][j] = -1.;
	    Bls->maxphasegap[lcnum][j] = -1.;
	    Bls->depth1_2tran[lcnum][j] = -1.;
	    Bls->depth2_2tran[lcnum][j] = -1.;
	    Bls->delchi2_2tran[lcnum][j] = -1.;
	    Bls->sr_sec[lcnum][j] = -1.;
	    Bls->srsum_sec[lcnum][j] = -1.;
	    Bls->q_sec[lcnum][j] = -1.;
	    Bls->epoch_sec[lcnum][j] = -1.;
	    Bls->H_sec[lcnum][j] = -1.;
	    Bls->L_sec[lcnum][j] = -1.;
	    Bls->depth_sec[lcnum][j] = -1.;
	    Bls->nt_sec[lcnum][j] = 0;
	    Bls->Nt_sec[lcnum][j] = 0;
	    Bls->sigtopink_sec[lcnum][j] = -1.;
	    Bls->deltachi2transit_sec[lcnum][j] = -1.;
	    Bls->binsignaltonoise_sec[lcnum][j] = -1.;
	    Bls->phaseoffset_sec[lcnum][j] = -1.;
	    Bls->harmmean[lcnum][j] = -1.;
	    Bls->fundA[lcnum][j] = -1.;
	    Bls->fundB[lcnum][j] = -1.;
	    Bls->harmA[lcnum][j] = -1.;
	    Bls->harmB[lcnum][j] = -1.;
	    Bls->harmamp[lcnum][j] = -1.;
	    Bls->harmdeltachi2[lcnum][j] = -1.;
	  }
	}
      *bperpos = -1.;
      *chisqrminus = 999999.;
      *meanmagval = -1.;
      free(weight);
      free(y);
      free(ibi);
      free(best_id);
      free(sr_ave);
      free(binned_sr_ave);
      free(binned_sr_sig);
      free(in1_array);
      free(in2_array);
      free(qtran_array);
      free(depth_array);
      free(bper_array);
      free(sr_ave_minus);
      free(binned_sr_ave_minus);
      free(binned_sr_sig_minus);
      free(p_minus);
      if(freqarray != NULL)
	free(freqarray);
      if(probvals != NULL)
	free(probvals);
      if(srshiftvals != NULL)
	free(srshiftvals);
      if(srnoshiftvals != NULL)
	free(srnoshiftvals);
#ifdef PARALLEL
      if(srvals != NULL) free(srvals);
      if(srvals_minus != NULL) free(srvals_minus);
#endif
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

  if(Bls->extraparams) {

    /* Manually compute srsum for each peak - this is a bit redundant, but
       we want to avoid computing a sqrt at each trial point in the BLS 
       calculation */
    for(jf=0;jf<Npeak;jf++)
      {
	Bls->srsum[lcnum][jf] = 0.0;
	p0 = bper[jf];
	if(p0 <= 0) continue;
	f0 = 1.0 / p0;
	
	if(adjust_qmin_mindt) {
	  qmi_test = mindt*f0;
	  if(qmi_test > 1) qmi_test = 1.0;
	  if(qmi_test > qmi) {
	    if(reduce_nb) {
	      nb = MIN_(nbsave, ceil(1./(0.5*qmi_test)));
	    }
	    kmi = (int) (qmi_test*(double)nb);
	    if(kmi < 1) kmi = 1;
	    if(qmi_test > qma)
	      kma = ((int) (qmi_test*(double)nb)) + 1;
	    else if(reduce_nb)
	      kma = ((int) (qma*(double)nb)) + 1;
	    nb1 = nb;
	    nbkma = nb+kma;
	    kkmi = qmi;
	    if(kkmi < (double) minbin / rn) kkmi = (double) minbin / rn;
	  } else {
	    kmi = kmisave; kkmi = kkmisave;
	    kma = kmasave;
	    if(reduce_nb)
	      nb = nbsave;
	    nb1 = nb;
	    nbkma = nb+kma;
	  }
	}
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
		    if(s > 0.)
		      {
			Bls->srsum[lcnum][jf] += sqrt(pow);
		      }
		  }
	      }
	  }
	if(Bls->extraparams) {
	  Bls->srsum[lcnum][jf] = Bls->srsum[lcnum][jf] / kma;
	}
      }
    /* Call BLS Extra Params 1 */
    GetExtraBLSParameters1(n, x, nf, srnoshiftvals, &srshiftvals, 
			   probvals, freqarray, Bls, lcnum, Npeak,
			   best_id);
  }


  /* Collect all the output bls parameters for the peaks */
  for(i=0;i<Npeak;i++)
    {
      if(bper[i] > -1)
	{
	  if(!nobinnedrms)
	    bpow[i] = snval[i]*allstddev + binned_sr_ave[best_id[i]];
	  else
	    bpow[i] = snval[i]*global_best_sr_stddev + global_best_sr_ave;
	  if(adjust_qmin_mindt && reduce_nb) {
	    qmi_test = mindt/bper[i];
	    if(qmi_test > 1.0) qmi_test = 1.0;
	    if(qmi_test > qmi) {
	      nb = MIN_(nbsave, ceil(1./(0.5*qmi_test)));
	    } else {
	      nb = nbsave;
	    }
	    if(nb != nbsave) {
	      in1[i] = rint(((double) nbsave*in1_array[best_id[i]])/(double) nb);
	      in2[i] = rint(((double) nbsave*in2_array[best_id[i]])/(double) nb);
	    }
	    else {
	      in1[i] = in1_array[best_id[i]];
	      in2[i] = in2_array[best_id[i]];
	    }
	  } else {
	    in1[i] = in1_array[best_id[i]];
	    in2[i] = in2_array[best_id[i]];
	  }
	  in1_ph[i] = ((double) in1_array[best_id[i]]) / ((double) nb);
	  in2_ph[i] = ((double) in2_array[best_id[i]]) / ((double) nb);
	  if(fittrap) {
	    qingress[i]=0.25;
	    OOTmag[i]=*meanmagval;
	    dofittrap_amoeba(n, t, x, e, bper[i], &(qtran_array[best_id[i]]), &(qingress[i]), &(in1_ph[i]), &(in2_ph[i]), &(depth_array[best_id[i]]), &(OOTmag[i]));
	  } else {
	    qingress[i] = 0.;
	    OOTmag[i] = *meanmagval;
	  }
	  // Be sure to correct for transits past the edge
	  if(in2[i] >= nb) in2[i] = in2[i] - nb;
	  qtran[i] = qtran_array[best_id[i]];
	  bt0[i] = t[0] + (0.5*qtran[i]+in1_ph[i])*bper[i];
	  depth[i] = depth_array[best_id[i]];
	  sde[i] = (bpow[i] - ((double)sde_sr_ave / (double)nsr))/sqrt((double)((sde_srsqr_ave / (long double) nsr) - (sde_sr_ave*sde_sr_ave/((long double)nsr*(long double)nsr))));
	  chisqrplus[i] = -bpow[i]*bpow[i]*sumweights;

	  fraconenight[i] = getfrac_onenight(n, t, u, v, e, bper[i], depth[i], qtran[i], (t[0] + in1_ph[i]*bper[i]), timezone);
	  /* Get the signal to pink noise for the peak */
	  if(Bls->extraparams) {
	    ntvptr = &(Bls->ntv[lcnum][i]);
	  }
	  getsignaltopinknoiseforgivenblsmodel(n, t, x, e, bper[i], qtran[i], depth[i], in1_ph[i], &nt[i], &Nt[i], &Nbefore[i], &Nafter[i], &rednoise[i], &whitenoise[i], &sigtopink[i], qingress[i], OOTmag[i], ntvptr);
	  if(Bls->extraparams) {
	    /* Collect the extra BLS parameters */
	    GetExtraBLSParameters2(n, t, x, e, bper[i], qtran[i], depth[i], in1_ph[i], qingress[i], OOTmag[i], Bls, lcnum, i);
	  }
	}
    }

  /* Now find the maximum inverse transit */
  if(!nobinnedrms)
    allstddev_minus = subtract_binnedrms(nf, sr_ave_minus, BIN_FACTOR, &allave_minus, &nclippedfreq, binned_sr_ave_minus, binned_sr_sig_minus);
  else {
    getclippedavestddev(nf,p_minus,&global_best_sr_ave_inv,&global_best_sr_stddev_inv);
    nclippedfreq = nf;
  }

  if(nclippedfreq > 0.)
    {
      if(!nobinnedrms) {
	for(i=0;i<nf;i++)
	  {
	    if(binned_sr_ave[i] > 0.)
	      {
		p_minus[i] = (p_minus[i] - binned_sr_ave_minus[i]) / allstddev_minus;
		if(p_minus[i] > bpowminus)
		  {
		    bpowminus = p_minus[i];
		    sr_minus = p_minus[i]*allstddev_minus + binned_sr_ave_minus[i];
		    *bperpos = bper_array[i];
		    *chisqrminus = -sr_minus*sr_minus*sumweights;
		  }
	      }
	    else
	      p_minus[i] = 0.;
	  }
      }
      else {
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
      free(y);
      free(ibi);
      free(best_id);
      free(sr_ave);
      free(binned_sr_ave);
      free(binned_sr_sig);
      free(in1_array);
      free(in2_array);
      free(qtran_array);
      free(depth_array);
      free(bper_array);
      free(sr_ave_minus);
      free(binned_sr_ave_minus);
      free(binned_sr_sig_minus);
      free(p_minus);
      if(freqarray != NULL)
	free(freqarray);
      if(probvals != NULL)
	free(probvals);
      if(srshiftvals != NULL)
	free(srshiftvals);
      if(srnoshiftvals != NULL)
	free(srnoshiftvals);
#ifdef PARALLEL
      if(srvals != NULL) free(srvals);
      if(srvals_minus != NULL) free(srvals_minus);
#endif
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
	  fprintf(outfile,"#Period  S/N   SR\n");
	  if(!nobinnedrms) {
	    for(i=0;i<nf;i++) {
	      fprintf(outfile,"%.17g %.17g %.17g\n",bper_array[i],p[i],(p[i]*allstddev+binned_sr_ave[i]));
	    }
	  }
	  else {
	    for(i=0;i<nf;i++) {
	      fprintf(outfile,"%.17g %.17g %.17g\n",bper_array[i],p[i],(p[i]*global_best_sr_stddev + global_best_sr_ave));
	    }
	  }
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

      fprintf(outfile2,"#Time  Mag_obs   Mag_model   Error   Phase\n");
      for(i=0;i<n;i++)
	{
	  ph = (u[i] - in1_ph[0]*bper[0])*f0;
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
	  ph = (jdtmp - t[0] - in1_ph[0]*bper[0])*f0;
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
	for(i=0;i<n;i++)
	  {
	    u[i]=t[i]-t1;
	  }
      }

      for(i=0;i<n;i++)
	{
	  ph = (u[i] - in1_ph[0]*bper[0])*f0;
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
  free(y);
  free(ibi);
  free(best_id);
  free(sr_ave);
  free(binned_sr_ave);
  free(binned_sr_sig);
  free(in1_array);
  free(in2_array);
  free(qtran_array);
  free(depth_array);
  free(bper_array);
  free(sr_ave_minus);
  free(binned_sr_ave_minus);
  free(binned_sr_sig_minus);
  free(p_minus);
  if(freqarray != NULL)
    free(freqarray);
  if(probvals != NULL)
    free(probvals);
  if(srshiftvals != NULL)
    free(srshiftvals);
  if(srnoshiftvals != NULL)
    free(srnoshiftvals);
  
#ifdef PARALLEL
  if(srvals != NULL) free(srvals);
  if(srvals_minus != NULL) free(srvals_minus);
#endif
  if(t_mask != NULL) free(t_mask);
  if(x_mask != NULL) free(x_mask);
  if(e_mask != NULL) free(e_mask);
  return(0);
}

/* This version adjusts the qmin and qmax according the period using a specified rmin and rmax, it assumes that for P in days and R in solar radii that q is given by:
q = 0.076 * R**(2/3) / P**(2/3)
*/

int eebls_rad(int n_in, double *t_in, double *x_in, double *e_in, double *u, double *v, int nf, double fmin, double df, int nb, double rmin, double rmax, double *p, int Npeak, double *bper, double *bt0, double *bpow, double *sde, double *snval, double *depth, double *qtran, int *in1, int *in2, double *in1_ph, double *in2_ph, double *chisqrplus, double *chisqrminus, double *bperpos, double *meanmagval, double timezone, double *fraconenight, int operiodogram, char *outname, int omodel, char *modelname, int correctlc, int ascii,int *nt, int *Nt, int *Nbefore, int *Nafter, double *rednoise, double *whitenoise, double *sigtopink, int fittrap, double *qingress, double *OOTmag, int ophcurve, char *ophcurvename, double phmin, double phmax, double phstep, int ojdcurve, char *ojdcurvename, double jdstep, int nobinnedrms, int freq_step_type, int adjust_qmin_mindt, int reduce_nb, int reportharmonics, _Bls *Bls, int lcnum, int lclistnum, int usemask, _Variable *maskvar, int isoptimal, double Aval)
{
  double *y;
  double *ibi;
#ifdef PARALLEL
  int firsttime =0;
  double *srvals = NULL;
  double *srvals_minus = NULL;
#else
  static int firsttime=0;
  static double *srvals;
  static double *srvals_minus;
#endif
  int minbin = 5;
  int nbmax, nbtot;
  int nsrvals, nsrvals_minus, test, foundsofar, dumint1;
  double powerplus, powerminus, bpowminus, dumdbl1, dumdbl2, jdtmp;
  double sumweights, qminP, qmaxP, Ppow,rminpow,rmaxpow;
  double tot, rnbtot, sr_minus, *weight;
  double rn, s,t1,f0,p0,ph,ph2,phb1,phb2,testpow,rn1,rn3,s3,rn4,rn5;
  double kkmi, kk, allave, allstddev, allave_minus, allstddev_minus, *qtran_array, *depth_array, minbest;
  double *sr_ave, *binned_sr_ave, *binned_sr_sig;
  int kmi, kma,nb1,nbkma,i,jf,j,k,jn1,jn2,jnb,nb2,nsr,nclippedfreq, *in1_array, *in2_array, *best_id, nbsave;
  double *p_minus, *bper_array, *sr_ave_minus, *binned_sr_ave_minus, *binned_sr_sig_minus;
  long double sde_sr_ave, sde_srsqr_ave;
  FILE *outfile, *outfile2;
  double global_best_sr_ave, global_best_sr_stddev, global_best_sr_ave_inv, global_best_sr_stddev_inv, qmi_test, mindt;
  double *freqarray = NULL;
  double *probvals = NULL;
  double *srshiftvals = NULL;
  double *srnoshiftvals = NULL;
  int *ntvptr = NULL;
  int n;
  double *t, *x, *e;
  double *t_mask = NULL, *x_mask = NULL, *e_mask = NULL;
  double Cval;
  
  if(isoptimal) {
    Cval = pow(fmin, (1.0/3.0)) - Aval/3.0;
  }
  
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
	  snval[j] = -1.;
	  bpow[j] = -1.;
	  bt0[j] = -1.;
	  in1[j] = -1;
	  in2[j] = -1;
	  in1_ph[j] = -1.;
	  in2_ph[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	  if(Bls->extraparams) {
	    Bls->srsum[lcnum][j] = -1.;
	    Bls->ressig[lcnum][j] = -1.;
	    Bls->dipsig[lcnum][j] = -1.;
	    Bls->srshift[lcnum][j] = -1.;
	    Bls->srsig[lcnum][j] = -1.;
	    Bls->snrextra[lcnum][j] = -1.;
	    Bls->dsp[lcnum][j] = -1.;
	    Bls->dspg[lcnum][j] = -1.;
	    Bls->freqlow[lcnum][j] = -1.;
	    Bls->freqhigh[lcnum][j] = -1.;
	    Bls->logprob[lcnum][j] = -1.;
	    Bls->peakarea[lcnum][j] = -1.;
	    Bls->peakmean[lcnum][j] = -1.;
	    Bls->peakdev[lcnum][j] = -1.;
	    Bls->lomblog[lcnum][j] = -1.;
	    Bls->ntv[lcnum][j] = 0;
	    Bls->gezadsp[lcnum][j] = -1.;
	    Bls->ootsig[lcnum][j] = -1.;
	    Bls->trsig[lcnum][j] = -1.;
	    Bls->ootdftf[lcnum][j] = -1.;
	    Bls->ootdfta[lcnum][j] = -1.;
	    Bls->binsignaltonoise[lcnum][j] = -1.;
	    Bls->maxphasegap[lcnum][j] = -1.;
	    Bls->depth1_2tran[lcnum][j] = -1.;
	    Bls->depth2_2tran[lcnum][j] = -1.;
	    Bls->delchi2_2tran[lcnum][j] = -1.;
	    Bls->sr_sec[lcnum][j] = -1.;
	    Bls->srsum_sec[lcnum][j] = -1.;
	    Bls->q_sec[lcnum][j] = -1.;
	    Bls->epoch_sec[lcnum][j] = -1.;
	    Bls->H_sec[lcnum][j] = -1.;
	    Bls->L_sec[lcnum][j] = -1.;
	    Bls->depth_sec[lcnum][j] = -1.;
	    Bls->nt_sec[lcnum][j] = 0;
	    Bls->Nt_sec[lcnum][j] = 0;
	    Bls->sigtopink_sec[lcnum][j] = -1.;
	    Bls->deltachi2transit_sec[lcnum][j] = -1.;
	    Bls->binsignaltonoise_sec[lcnum][j] = -1.;
	    Bls->phaseoffset_sec[lcnum][j] = -1.;
	    Bls->harmmean[lcnum][j] = -1.;
	    Bls->fundA[lcnum][j] = -1.;
	    Bls->fundB[lcnum][j] = -1.;
	    Bls->harmA[lcnum][j] = -1.;
	    Bls->harmB[lcnum][j] = -1.;
	    Bls->harmamp[lcnum][j] = -1.;
	    Bls->harmdeltachi2[lcnum][j] = -1.;
	  }
	}
      *bperpos = -1.;
      *chisqrminus = 999999.;
      *meanmagval = -1.;
      
      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
      return 1;
    }
  }
  

  if(firsttime == 0 && !nobinnedrms)
    {
      if((srvals = (double *) malloc(SRVALSSIZE * sizeof(double))) == NULL ||
	 (srvals_minus = (double *) malloc(SRVALSSIZE * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      firsttime = 1;
    }

  rminpow = pow(rmin,0.6666667);
  rmaxpow = pow(rmax,0.6666667);
  nbmax = 2*nb;

  /***********************************************************/

  if((sr_ave = (double *) malloc(nf * sizeof(double))) == NULL ||
     (y = (double *) malloc(nbmax * sizeof(double))) == NULL ||
     (ibi = (double *) malloc(nbmax * sizeof(double))) == NULL ||
     (binned_sr_ave = (double *) malloc(nf * sizeof(double))) == NULL ||
     (binned_sr_sig = (double *) malloc(nf * sizeof(double))) == NULL ||
     (in1_array = (int *) malloc(nf * sizeof(int))) == NULL ||
     (in2_array = (int *) malloc(nf * sizeof(int))) == NULL ||
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
  if(Bls->extraparams) {
    if((freqarray = (double *) malloc(nf * sizeof(double))) == NULL ||
       (probvals = (double *) malloc(nf * sizeof(double))) == NULL ||
       (srnoshiftvals = (double *) malloc(nf * sizeof(double))) == NULL)
      {
	fprintf(stderr,"Memory Allocation Error\n");
	exit(3);
      }
  }

  if(nb >= nbmax) {
    error(ERR_BLSNBMAX);
  }
  if(nb < 2) nb = 2;
  tot = t[n-1] - t[0];
  if(fmin < 1./tot) {
    error(ERR_BLSFMINTOOSMALL);
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
  //*bpow = 0.;

  bpowminus = 0.;
  sde_sr_ave = 0.;
  sde_srsqr_ave = 0.;

  /**************The following variables are defined for the extension
		 c     of arrays  ibi()  and  y()  [ see below ] ***************/

  if(adjust_qmin_mindt) {
    mindt = 0;
    if(n > 1) {mindt = t[1] - t[0];}
    for(i=2; i < n; i++) {
      if(t[i]-t[i-1] < mindt) mindt = t[i] - t[i-1];
    }
  }

  nbsave = nb;
  nb1 = nb;

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

  /* Get the minimum time */
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
      if(isoptimal) {
	f0 = ((Aval*(((double)jf)+1.0)/3.0) + Cval);
	f0 = f0*f0*f0;
	p0 = 1./f0;
      } else {
	if(!freq_step_type) {
	  f0=fmin+df*((double)jf);
	  p0=1./f0;
	} else if(freq_step_type == VARTOOLS_FREQSTEPTYPE_PERIOD) {
	  p0 = (1./fmin) - df*((double)jf);
	  f0 = 1./p0;
	} else if(freq_step_type == VARTOOLS_FREQSTEPTYPE_LOGPERIOD) {
	  f0 = exp(log(fmin) + df*((double)jf));
	  p0=1./f0;
	}
      }

      /*
	c
	c======================================================
	c     Compute folded time series with  *p0*  period
	c======================================================
	c
      */

      Ppow = pow(p0,0.6666667);
      qminP = 0.076*rminpow/Ppow;
      qmaxP = 0.076*rmaxpow/Ppow;

      if(qminP > 1.0) qminP = 1.0;
      if(qmaxP > 1.0) qmaxP = 1.0;

      if(adjust_qmin_mindt) {
	qmi_test = mindt*f0;
	if(qmi_test > qminP) {
	  qminP = qmi_test;
	  if(qmi_test > qmaxP)
	    qmaxP = qmi_test;
	}
	if(reduce_nb) {
	  nb = MIN_(nbsave, ceil(1./(0.5*qminP)));
	  if(nb < 2) nb=2;
	  nb1 = nb;
	}
      }

      kmi = (int) (qminP*(double)nb);
      if(kmi < 1) kmi = 1;
      if(kmi > nb-1) kmi = nb-1;
      kma = ((int) (qmaxP*(double)nb)) + 1;
      if(kma > nb-1) kma = nb-1;
      kkmi = qminP;
      if(kkmi < (double) minbin / rn) kkmi = (double) minbin / rn;
      nbkma = nb+kma;


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
      nsrvals = 0;
      nsrvals_minus = 0;
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
		  testpow = s*s/(rn1*(1. - rn1));
		  if(s > 0.)
		    {
		      if(!nobinnedrms) {
			srvals[nsrvals] = sqrt(testpow);
			nsrvals++;
		      }
		      if(testpow >= powerplus)
			{
			  powerplus = testpow;
			  jn1 = i;
			  jn2 = j;
			  rn3 = rn1;
			  rn5 = rn4;
			  s3 = s;
			}
		    }
		  else if(s < 0.)
		    {
		      if(!nobinnedrms) {
			srvals_minus[nsrvals_minus] = sqrt(testpow);
			nsrvals_minus++;
		      }
		      if(testpow >= powerminus)
			{
			  powerminus = testpow;
			}
		    }
		}
	    }
	}
      // Find the average value of the srvals
      if(!nobinnedrms) {
	sr_ave[jf] = getclippedsrave(nsrvals,srvals);
	sr_ave_minus[jf] = getclippedsrave(nsrvals_minus,srvals_minus);
      }
      powerplus = sqrt(powerplus);
      sde_sr_ave += powerplus;
      sde_srsqr_ave += powerplus*powerplus;
      p[jf] = powerplus;
      powerminus = sqrt(powerminus);
      p_minus[jf] = powerminus;
      //sr_ave += powerplus;
      //srsqr_ave += powerplus*powerplus;
      nsr++;
      in1_array[jf] = jn1;
      in2_array[jf] = jn2;
      qtran_array[jf] = (double)(jn2 - jn1 + 1)/(double)nb;
      depth_array[jf] = powerplus/sqrt(rn3*(1.-rn3));
      bper_array[jf] = p0;
      if(freqarray != NULL)
	freqarray[jf] = f0;
      /*      if(powerplus >= *bpow)
	{
	  *bpow = powerplus;
	  *in1 = jn1;
	  *in2 = jn2;
	  *qtran = (double)(jn2 - jn1 + 1)/(double)nb;
	  *depth = powerplus/sqrt(rn3*(1.-rn3));
	  *bper = p0;
	}
      if(powerminus >= bpowminus)
	{
	  bpowminus = powerminus;
	  *bperpos = p0;
	  }*/
    }
  if(Bls->extraparams) {
    /* Save a copy of the uncorrected SR values if we need to compute the
       extra parameters; Note that the code used for the extra parameters
       is derived from the lc/blsanal tool used by HATNet, in which SR is
       treated as the square of the SR values computed by this BLS code.
    */
    memcpy(srnoshiftvals, p, nf*sizeof(double));
    for(jf=0;jf<nf;jf++) {
      srnoshiftvals[jf] = srnoshiftvals[jf]*srnoshiftvals[jf];
    }
  }
  if(!nobinnedrms)
    allstddev = subtract_binnedrms(nf, sr_ave, BIN_FACTOR, &allave, &nclippedfreq, binned_sr_ave, binned_sr_sig);
  else {
    getclippedavestddev(nf,p,&global_best_sr_ave,&global_best_sr_stddev);
    nclippedfreq = nf;
  }

  /* Now let's find the peaks in the periodogram, first convert the periodogram from SR to SN ratio */

  if(nclippedfreq > Npeak)
    {
      if(!nobinnedrms) {
	for(i=0;i<nf;i++)
	  {
	    if(binned_sr_ave[i] > 0.)
	      {
		p[i] = (p[i] - binned_sr_ave[i]) / allstddev;
		/*	      if(p[i] > *bpow)
			      {
			      *bpow = p[i];
			      sr_plus = p[i]*allstddev + binned_sr_ave[i];
			      *in1 = in1_array[i];
			      *in2 = in2_array[i];
			      *qtran = qtran_array[i];
			      *depth = depth_array[i];
			      *bper = bper_array[i];
			      }*/
	      }
	    else
	      p[i] = 0.;
	  }
      }
      else {
	for(i=0; i<nf; i++) {
	  p[i] = (p[i] - global_best_sr_ave)/global_best_sr_stddev;
	}
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
	  in1[j] = -1;
	  in2[j] = -1;
	  in1_ph[j] = -1.;
	  in2_ph[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	  if(Bls->extraparams) {
	    Bls->srsum[lcnum][j] = -1.;
	    Bls->ressig[lcnum][j] = -1.;
	    Bls->dipsig[lcnum][j] = -1.;
	    Bls->srshift[lcnum][j] = -1.;
	    Bls->srsig[lcnum][j] = -1.;
	    Bls->snrextra[lcnum][j] = -1.;
	    Bls->dsp[lcnum][j] = -1.;
	    Bls->dspg[lcnum][j] = -1.;
	    Bls->freqlow[lcnum][j] = -1.;
	    Bls->freqhigh[lcnum][j] = -1.;
	    Bls->logprob[lcnum][j] = -1.;
	    Bls->peakarea[lcnum][j] = -1.;
	    Bls->peakmean[lcnum][j] = -1.;
	    Bls->peakdev[lcnum][j] = -1.;
	    Bls->lomblog[lcnum][j] = -1.;
	    Bls->ntv[lcnum][j] = 0;
	    Bls->gezadsp[lcnum][j] = -1.;
	    Bls->ootsig[lcnum][j] = -1.;
	    Bls->trsig[lcnum][j] = -1.;
	    Bls->ootdftf[lcnum][j] = -1.;
	    Bls->ootdfta[lcnum][j] = -1.;
	    Bls->binsignaltonoise[lcnum][j] = -1.;
	    Bls->maxphasegap[lcnum][j] = -1.;
	    Bls->depth1_2tran[lcnum][j] = -1.;
	    Bls->depth2_2tran[lcnum][j] = -1.;
	    Bls->delchi2_2tran[lcnum][j] = -1.;
	    Bls->sr_sec[lcnum][j] = -1.;
	    Bls->srsum_sec[lcnum][j] = -1.;
	    Bls->q_sec[lcnum][j] = -1.;
	    Bls->epoch_sec[lcnum][j] = -1.;
	    Bls->H_sec[lcnum][j] = -1.;
	    Bls->L_sec[lcnum][j] = -1.;
	    Bls->depth_sec[lcnum][j] = -1.;
	    Bls->nt_sec[lcnum][j] = 0;
	    Bls->Nt_sec[lcnum][j] = 0;
	    Bls->sigtopink_sec[lcnum][j] = -1.;
	    Bls->deltachi2transit_sec[lcnum][j] = -1.;
	    Bls->binsignaltonoise_sec[lcnum][j] = -1.;
	    Bls->phaseoffset_sec[lcnum][j] = -1.;
	    Bls->harmmean[lcnum][j] = -1.;
	    Bls->fundA[lcnum][j] = -1.;
	    Bls->fundB[lcnum][j] = -1.;
	    Bls->harmA[lcnum][j] = -1.;
	    Bls->harmB[lcnum][j] = -1.;
	    Bls->harmamp[lcnum][j] = -1.;
	    Bls->harmdeltachi2[lcnum][j] = -1.;
	  }
	}
      *bperpos = -1.;
      *chisqrminus = 999999.;
      *meanmagval = -1.;

      free(weight);
      free(y);
      free(ibi);
      free(best_id);
      free(sr_ave);
      free(binned_sr_ave);
      free(binned_sr_sig);
      free(in1_array);
      free(in2_array);
      free(qtran_array);
      free(depth_array);
      free(bper_array);
      free(sr_ave_minus);
      free(binned_sr_ave_minus);
      free(binned_sr_sig_minus);
      free(p_minus);
      if(freqarray != NULL)
	free(freqarray);
      if(probvals != NULL)
	free(probvals);
      if(srshiftvals != NULL)
	free(srshiftvals);
      if(srnoshiftvals != NULL)
	free(srnoshiftvals);

#ifdef PARALLEL
      if(srvals != NULL) free(srvals);
      if(srvals_minus != NULL) free(srvals_minus);
#endif
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
	      if((!reportharmonics && !isDifferentPeriods(MIN_(bper[j],bper_array[i]),MAX_(bper[j],bper_array[i]),tot)) || (reportharmonics && !isDifferentPeriodsDontCheckHarmonics(MIN_(bper[j],bper_array[i]),MAX_(bper[j],bper_array[i]),tot)))
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
		  if((!reportharmonics && !isDifferentPeriods(MIN_(bper[j],bper_array[i]),MAX_(bper[j],bper_array[i]),tot)) || (reportharmonics && !isDifferentPeriodsDontCheckHarmonics(MIN_(bper[j],bper_array[i]),MAX_(bper[j],bper_array[i]),tot)))
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
      /* We have a few peaks, but not Npeak of them */
      mysort3_int(foundsofar,snval,bper,best_id);
      for(j=foundsofar;j<Npeak;j++)
	{
	  /* Put -1 for the remaining peaks */
	  bper[j] = -1.;
	  snval[j] = -1.;
	  bpow[j] = -1.;
	  bt0[j] = -1.;
	  in1[j] = -1;
	  in2[j] = -1;
	  in1_ph[j] = -1.;
	  in2_ph[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	  if(Bls->extraparams) {
	    Bls->srsum[lcnum][j] = -1.;
	    Bls->ressig[lcnum][j] = -1.;
	    Bls->dipsig[lcnum][j] = -1.;
	    Bls->srshift[lcnum][j] = -1.;
	    Bls->srsig[lcnum][j] = -1.;
	    Bls->snrextra[lcnum][j] = -1.;
	    Bls->dsp[lcnum][j] = -1.;
	    Bls->dspg[lcnum][j] = -1.;
	    Bls->freqlow[lcnum][j] = -1.;
	    Bls->freqhigh[lcnum][j] = -1.;
	    Bls->logprob[lcnum][j] = -1.;
	    Bls->peakarea[lcnum][j] = -1.;
	    Bls->peakmean[lcnum][j] = -1.;
	    Bls->peakdev[lcnum][j] = -1.;
	    Bls->lomblog[lcnum][j] = -1.;
	    Bls->ntv[lcnum][j] = 0;
	    Bls->gezadsp[lcnum][j] = -1.;
	    Bls->ootsig[lcnum][j] = -1.;
	    Bls->trsig[lcnum][j] = -1.;
	    Bls->ootdftf[lcnum][j] = -1.;
	    Bls->ootdfta[lcnum][j] = -1.;
	    Bls->binsignaltonoise[lcnum][j] = -1.;
	    Bls->maxphasegap[lcnum][j] = -1.;
	    Bls->depth1_2tran[lcnum][j] = -1.;
	    Bls->depth2_2tran[lcnum][j] = -1.;
	    Bls->delchi2_2tran[lcnum][j] = -1.;
	    Bls->sr_sec[lcnum][j] = -1.;
	    Bls->srsum_sec[lcnum][j] = -1.;
	    Bls->q_sec[lcnum][j] = -1.;
	    Bls->epoch_sec[lcnum][j] = -1.;
	    Bls->H_sec[lcnum][j] = -1.;
	    Bls->L_sec[lcnum][j] = -1.;
	    Bls->depth_sec[lcnum][j] = -1.;
	    Bls->nt_sec[lcnum][j] = 0;
	    Bls->Nt_sec[lcnum][j] = 0;
	    Bls->sigtopink_sec[lcnum][j] = -1.;
	    Bls->deltachi2transit_sec[lcnum][j] = -1.;
	    Bls->binsignaltonoise_sec[lcnum][j] = -1.;
	    Bls->phaseoffset_sec[lcnum][j] = -1.;
	    Bls->harmmean[lcnum][j] = -1.;
	    Bls->fundA[lcnum][j] = -1.;
	    Bls->fundB[lcnum][j] = -1.;
	    Bls->harmA[lcnum][j] = -1.;
	    Bls->harmB[lcnum][j] = -1.;
	    Bls->harmamp[lcnum][j] = -1.;
	    Bls->harmdeltachi2[lcnum][j] = -1.;
	  }
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
	  in1[j] = -1;
	  in2[j] = -1;
	  in1_ph[j] = -1.;
	  in2_ph[j] = -1.;
	  qtran[j] = -1.;
	  depth[j] = -1.;
	  sde[j] = -1.;
	  chisqrplus[j] = 999999.;
      	  fraconenight[j] = -1.;
	  if(Bls->extraparams) {
	    Bls->srsum[lcnum][j] = -1.;
	    Bls->ressig[lcnum][j] = -1.;
	    Bls->dipsig[lcnum][j] = -1.;
	    Bls->srshift[lcnum][j] = -1.;
	    Bls->srsig[lcnum][j] = -1.;
	    Bls->snrextra[lcnum][j] = -1.;
	    Bls->dsp[lcnum][j] = -1.;
	    Bls->dspg[lcnum][j] = -1.;
	    Bls->freqlow[lcnum][j] = -1.;
	    Bls->freqhigh[lcnum][j] = -1.;
	    Bls->logprob[lcnum][j] = -1.;
	    Bls->peakarea[lcnum][j] = -1.;
	    Bls->peakmean[lcnum][j] = -1.;
	    Bls->peakdev[lcnum][j] = -1.;
	    Bls->lomblog[lcnum][j] = -1.;
	    Bls->ntv[lcnum][j] = 0;
	    Bls->gezadsp[lcnum][j] = -1.;
	    Bls->ootsig[lcnum][j] = -1.;
	    Bls->trsig[lcnum][j] = -1.;
	    Bls->ootdftf[lcnum][j] = -1.;
	    Bls->ootdfta[lcnum][j] = -1.;
	    Bls->binsignaltonoise[lcnum][j] = -1.;
	    Bls->maxphasegap[lcnum][j] = -1.;
	    Bls->depth1_2tran[lcnum][j] = -1.;
	    Bls->depth2_2tran[lcnum][j] = -1.;
	    Bls->delchi2_2tran[lcnum][j] = -1.;
	    Bls->sr_sec[lcnum][j] = -1.;
	    Bls->srsum_sec[lcnum][j] = -1.;
	    Bls->q_sec[lcnum][j] = -1.;
	    Bls->epoch_sec[lcnum][j] = -1.;
	    Bls->H_sec[lcnum][j] = -1.;
	    Bls->L_sec[lcnum][j] = -1.;
	    Bls->depth_sec[lcnum][j] = -1.;
	    Bls->nt_sec[lcnum][j] = 0;
	    Bls->Nt_sec[lcnum][j] = 0;
	    Bls->sigtopink_sec[lcnum][j] = -1.;
	    Bls->deltachi2transit_sec[lcnum][j] = -1.;
	    Bls->binsignaltonoise_sec[lcnum][j] = -1.;
	    Bls->phaseoffset_sec[lcnum][j] = -1.;
	    Bls->harmmean[lcnum][j] = -1.;
	    Bls->fundA[lcnum][j] = -1.;
	    Bls->fundB[lcnum][j] = -1.;
	    Bls->harmA[lcnum][j] = -1.;
	    Bls->harmB[lcnum][j] = -1.;
	    Bls->harmamp[lcnum][j] = -1.;
	    Bls->harmdeltachi2[lcnum][j] = -1.;
	  }
	}
      *bperpos = -1.;
      *chisqrminus = 999999.;
      *meanmagval = -1.;
      free(weight);
      free(y);
      free(ibi);
      free(best_id);
      free(sr_ave);
      free(binned_sr_ave);
      free(binned_sr_sig);
      free(in1_array);
      free(in2_array);
      free(qtran_array);
      free(depth_array);
      free(bper_array);
      free(sr_ave_minus);
      free(binned_sr_ave_minus);
      free(binned_sr_sig_minus);
      free(p_minus);
      if(freqarray != NULL)
	free(freqarray);
      if(probvals != NULL)
	free(probvals);
      if(srshiftvals != NULL)
	free(srshiftvals);
      if(srnoshiftvals != NULL)
	free(srnoshiftvals);
#ifdef PARALLEL
      if(srvals != NULL) free(srvals);
      if(srvals_minus != NULL) free(srvals_minus);
#endif
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

  if(Bls->extraparams) {

    /* Manually compute srsum for each peak - this is a bit redundant, but
       we want to avoid computing a sqrt at each trial point in the BLS 
       calculation */
    for(jf=0;jf<Npeak;jf++)
      {
	Bls->srsum[lcnum][jf] = 0.0;
	p0 = bper[jf];
	if(p0 <= 0) continue;
	f0 = 1.0 / p0;

	Ppow = pow(p0,0.6666667);
	qminP = 0.076*rminpow/Ppow;
	qmaxP = 0.076*rmaxpow/Ppow;

	if(qminP > 1.0) qminP = 1.0;
	if(qmaxP > 1.0) qmaxP = 1.0;
	
	if(adjust_qmin_mindt) {
	  qmi_test = mindt*f0;
	  if(qmi_test > qminP) {
	    qminP = qmi_test;
	    if(qmi_test > qmaxP)
	      qmaxP = qmi_test;
	  }
	  if(reduce_nb) {
	    nb = MIN_(nbsave, ceil(1./(0.5*qminP)));
	    if(nb < 2) nb=2;
	    nb1 = nb;
	  }
	}

	kmi = (int) (qminP*(double)nb);
	if(kmi < 1) kmi = 1;
	if(kmi > nb-1) kmi = nb-1;
	kma = ((int) (qmaxP*(double)nb)) + 1;
	if(kma > nb-1) kma = nb-1;
	kkmi = qminP;
	if(kkmi < (double) minbin / rn) kkmi = (double) minbin / rn;
	nbkma = nb+kma;

	
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
		    testpow = s*s/(rn1*(1. - rn1));
		    if(s > 0.)
		      {
			Bls->srsum[lcnum][jf] += sqrt(testpow);
		      }
		  }
	      }
	  }
	if(Bls->extraparams) {
	  Bls->srsum[lcnum][jf] = Bls->srsum[lcnum][jf] / kma;
	}
      }
    /* Call BLS Extra Params 1 */
    GetExtraBLSParameters1(n, x, nf, srnoshiftvals, &srshiftvals, 
			   probvals, freqarray, Bls, lcnum, Npeak,
			   best_id);
  }

  /* Collect all the output bls parameters for the peaks */
  for(i=0;i<Npeak;i++)
    {
      if(bper[i] > -1)
	{
	  if(!nobinnedrms)
	    bpow[i] = snval[i]*allstddev + binned_sr_ave[best_id[i]];
	  else
	    bpow[i] = snval[i]*global_best_sr_stddev + global_best_sr_ave;

	  if(adjust_qmin_mindt && reduce_nb) {
	    Ppow = pow(bper[i],0.6666667);
	    qminP = 0.076*rminpow/Ppow;
	    if(qminP > 1.0) qminP = 1.0;
	    qmi_test = mindt/bper[i];
	    if(qmi_test > qminP) {
	      qminP = qmi_test;
	    }
	    nb = MIN_(nbsave, ceil(1./(0.5*qminP)));
	    if(nb != nbsave) {
	      in1[i] = rint(((double) nbsave*in1_array[best_id[i]])/(double) nb);
	      in2[i] = rint(((double) nbsave*in2_array[best_id[i]])/(double) nb);
	    }
	    else {
	      in1[i] = in1_array[best_id[i]];
	      in2[i] = in2_array[best_id[i]];
	    }
	  } else {
	    in1[i] = in1_array[best_id[i]];
	    in2[i] = in2_array[best_id[i]];
	  }
	  in1_ph[i] = ((double) in1_array[best_id[i]]) / ((double) nb);
	  in2_ph[i] = ((double) in2_array[best_id[i]]) / ((double) nb);
	  if(fittrap) {
	    qingress[i]=0.25;
	    OOTmag[i]=*meanmagval;
	    dofittrap_amoeba(n, t, x, e, bper[i], &(qtran_array[best_id[i]]), &(qingress[i]), &(in1_ph[i]), &(in2_ph[i]), &(depth_array[best_id[i]]), &(OOTmag[i]));
	  } else {
	    qingress[i] = 0.;
	    OOTmag[i] = *meanmagval;
	  }
	  // Be sure to correct for transits past the edge
	  if(in2[i] >= nb) in2[i] = in2[i] - nb;
	  qtran[i] = qtran_array[best_id[i]];
	  bt0[i] = t[0] + (0.5*qtran[i]+in1_ph[i])*bper[i];
	  depth[i] = depth_array[best_id[i]];
	  sde[i] = (bpow[i] - ((double)sde_sr_ave / (double)nsr))/sqrt((double)((sde_srsqr_ave / (long double) nsr) - (sde_sr_ave*sde_sr_ave/((long double)nsr*(long double)nsr))));
	  chisqrplus[i] = -bpow[i]*bpow[i]*sumweights;

	  fraconenight[i] = getfrac_onenight(n, t, u, v, e, bper[i], depth[i], qtran[i], (t[0] +in1_ph[i]*bper[i]), timezone);
	  /* Get the signal to pink noise for the peak */
	  if(Bls->extraparams) {
	    ntvptr = &(Bls->ntv[lcnum][i]);
	  }
	  getsignaltopinknoiseforgivenblsmodel(n, t, x, e, bper[i], qtran[i], depth[i], in1_ph[i], &nt[i], &Nt[i], &Nbefore[i], &Nafter[i], &rednoise[i], &whitenoise[i], &sigtopink[i], qingress[i], OOTmag[i],ntvptr);
	  if(Bls->extraparams) {
	    /* Collect the extra BLS parameters */
	    GetExtraBLSParameters2(n, t, x, e, bper[i], qtran[i], depth[i], in1_ph[i], qingress[i], OOTmag[i], Bls, lcnum, i);
	  }
	}

    }

  /* Now find the maximum inverse transit */
  if(!nobinnedrms)
    allstddev_minus = subtract_binnedrms(nf, sr_ave_minus, BIN_FACTOR, &allave_minus, &nclippedfreq, binned_sr_ave_minus, binned_sr_sig_minus);
  else {
    getclippedavestddev(nf,p_minus,&global_best_sr_ave_inv,&global_best_sr_stddev_inv);
    nclippedfreq = nf;
  }
  if(nclippedfreq > 0.)
    {
      if(!nobinnedrms) {
	for(i=0;i<nf;i++)
	  {
	    if(binned_sr_ave[i] > 0.)
	      {
		p_minus[i] = (p_minus[i] - binned_sr_ave_minus[i]) / allstddev_minus;
		if(p_minus[i] > bpowminus)
		  {
		    bpowminus = p_minus[i];
		    sr_minus = p_minus[i]*allstddev_minus + binned_sr_ave_minus[i];
		    *bperpos = bper_array[i];
		    *chisqrminus = -sr_minus*sr_minus*sumweights;
		  }
	      }
	    else
	      p_minus[i] = 0.;
	  }
      }
      else {
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
      free(y);
      free(ibi);
      free(best_id);
      free(sr_ave);
      free(binned_sr_ave);
      free(binned_sr_sig);
      free(in1_array);
      free(in2_array);
      free(qtran_array);
      free(depth_array);
      free(bper_array);
      free(sr_ave_minus);
      free(binned_sr_ave_minus);
      free(binned_sr_sig_minus);
      free(p_minus);
      if(freqarray != NULL)
	free(freqarray);
      if(probvals != NULL)
	free(probvals);
      if(srshiftvals != NULL)
	free(srshiftvals);
      if(srnoshiftvals != NULL)
	free(srnoshiftvals);
#ifdef PARALLEL
      if(srvals != NULL) free(srvals);
      if(srvals_minus != NULL) free(srvals_minus);
#endif
      if(t_mask != NULL) free(t_mask);
      if(x_mask != NULL) free(x_mask);
      if(e_mask != NULL) free(e_mask);
      return 1;
    }


  // *sde = (*bpow - ((double)sr_ave / (double)nsr))/sqrt((double)((srsqr_ave / (long double) nsr) - (sr_ave*sr_ave/(long double)(nsr*nsr))));

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
	  fprintf(outfile,"#Period  S/N   SR\n");
	  if(!nobinnedrms) {
	    for(i=0;i<nf;i++) {
	      fprintf(outfile,"%.17g %.17g %.17g\n",bper_array[i],p[i],(p[i]*allstddev+binned_sr_ave[i]));
	    }
	  }
	  else {
	    for(i=0;i<nf;i++) {
	      fprintf(outfile,"%.17g %.17g %.17g\n",bper_array[i],p[i],(p[i]*global_best_sr_stddev + global_best_sr_ave));
	    }
	  }
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

      fprintf(outfile2,"#Time  Mag_obs   Mag_model   Error   Phase\n");
      for(i=0;i<n;i++)
	{
	  ph = (u[i] - in1_ph[0]*bper[0])*f0;
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
	  ph = (jdtmp - t[0] - in1_ph[0]*bper[0])*f0;
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
	for(i=0;i<n;i++)
	  {
	    u[i]=t[i]-t1;
	  }
      }
      for(i=0;i<n;i++)
	{
	  ph = (u[i] - in1_ph[0]*bper[0])*f0;
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
  free(y);
  free(ibi);
  free(best_id);
  free(sr_ave);
  free(binned_sr_ave);
  free(binned_sr_sig);
  free(in1_array);
  free(in2_array);
  free(qtran_array);
  free(depth_array);
  free(bper_array);
  free(sr_ave_minus);
  free(binned_sr_ave_minus);
  free(binned_sr_sig_minus);
  free(p_minus);
  if(freqarray != NULL)
    free(freqarray);
  if(probvals != NULL)
    free(probvals);
  if(srshiftvals != NULL)
    free(srshiftvals);
  if(srnoshiftvals != NULL)
    free(srnoshiftvals);

#ifdef PARALLEL
  if(srvals != NULL) free(srvals);
  if(srvals_minus != NULL) free(srvals_minus);
#endif
  if(t_mask != NULL) free(t_mask);
  if(x_mask != NULL) free(x_mask);
  if(e_mask != NULL) free(e_mask);
  return(0);
}

void RunBLSCommand(ProgramData *p, _Bls *Bls, int lcnum, int lc_name_num, int thisindex, int threadindex)
{
  int i1, i2;
  char outname[MAXLEN];
  char outname2[MAXLEN];
  char outname3[MAXLEN];
  char outname4[MAXLEN];
  if(p->NJD[lcnum] > 1) {
    if(Bls->omodel)
      {
	i1 = 0;
	i2 = 0;
	while(p->lcnames[lc_name_num][i1] != '\0')
	  {
	    if(p->lcnames[lc_name_num][i1] == '/')
	      i2 = i1 + 1;
	    i1++;
	  }
	sprintf(outname2,"%s/%s%s",Bls->modeloutdir,&p->lcnames[lc_name_num][i2],Bls->modelsuffix);
      }
    if(Bls->ophcurve)
      {
	i1 = 0;
	i2 = 0;
	while(p->lcnames[lc_name_num][i1] != '\0')
	  {
	    if(p->lcnames[lc_name_num][i1] == '/')
	      i2 = i1 + 1;
	    i1++;
	  }
	sprintf(outname3,"%s/%s%s",Bls->ophcurveoutdir,&p->lcnames[lc_name_num][i2],Bls->ophcurvesuffix);
      }
    if(Bls->ojdcurve)
      {
	i1 = 0;
	i2 = 0;
	while(p->lcnames[lc_name_num][i1] != '\0')
	  {
	    if(p->lcnames[lc_name_num][i1] == '/')
	      i2 = i1 + 1;
	    i1++;
	  }
	sprintf(outname4,"%s/%s%s",Bls->ojdcurveoutdir,&p->lcnames[lc_name_num][i2],Bls->ojdcurvesuffix);
      }
    /* First check to see that the u/v vectors are large enough */
    if(Bls->sizeuv[lcnum] == 0)
      {
	Bls->sizeuv[lcnum] = p->NJD[lcnum];
	if((Bls->u[lcnum] = (double *) malloc(Bls->sizeuv[lcnum] * sizeof(double))) == NULL ||
	   (Bls->v[lcnum] = (double *) malloc(Bls->sizeuv[lcnum] * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
      }
    else if(Bls->sizeuv[lcnum] < p->NJD[lcnum])
      {
	Bls->sizeuv[lcnum] = p->NJD[lcnum];
	free(Bls->u[lcnum]);
	free(Bls->v[lcnum]);
	if((Bls->u[lcnum] = (double *) malloc(Bls->sizeuv[lcnum] * sizeof(double))) == NULL ||
	   (Bls->v[lcnum] = (double *) malloc(Bls->sizeuv[lcnum] * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
      }
    
    if(Bls->operiodogram)
      {
	i1 = 0;
	i2 = 0;
	while(p->lcnames[lc_name_num][i1] != '\0')
	  {
	    if(p->lcnames[lc_name_num][i1] == '/')
	      i2 = i1 + 1;
	    i1++;
	  }
	sprintf(outname,"%s/%s%s",Bls->outdir,&p->lcnames[lc_name_num][i2],Bls->suffix);
      }
    if(Bls->minper_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
      Bls->minper_val[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Bls->minper_expr);
    } 
    else if(Bls->minper_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
      Bls->minper_val[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->minper_var);
    }
    else {
      Bls->minper_val[lcnum] = Bls->minper;
    }
    if(Bls->maxper_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
      Bls->maxper_val[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Bls->maxper_expr);
    } 
    else if(Bls->maxper_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
      Bls->maxper_val[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->maxper_var);
    }
    else {
      Bls->maxper_val[lcnum] = Bls->maxper;
    }


    if(!Bls->rflag) {
      if(Bls->qmin_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	Bls->qmin_val[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Bls->qmin_expr);
      } 
      else if(Bls->qmin_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
	Bls->qmin_val[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->qmin_var);
      }
      else {
	Bls->qmin_val[lcnum] = Bls->qmin;
      }
      
      if(Bls->qmax_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	Bls->qmax_val[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Bls->qmax_expr);
      } 
      else if(Bls->qmax_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
	Bls->qmax_val[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->qmax_var);
      }
      else {
	Bls->qmax_val[lcnum] = Bls->qmax;
      }
    }
    else if(Bls->rflag == 1) {
      if(Bls->rmin_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	Bls->rmin_val[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Bls->rmin_expr);
      } 
      else if(Bls->rmin_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
	Bls->rmin_val[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->rmin_var);
      }
      else {
	Bls->rmin_val[lcnum] = Bls->rmin;
      }
      
      if(Bls->rmax_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	Bls->rmax_val[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Bls->rmax_expr);
      } 
      else if(Bls->rmax_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
	Bls->rmax_val[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->rmax_var);
      }
      else {
	Bls->rmax_val[lcnum] = Bls->rmax;
      }
    }
    else if(Bls->rflag == 2) {
      if(Bls->rho_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	Bls->rho_val[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Bls->rho_expr);
      } 
      else if(Bls->rho_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
	Bls->rho_val[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->rho_var);
      }
      else {
	Bls->rho_val[lcnum] = Bls->rho;
      }
      
      if(Bls->minexpdurfrac_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	Bls->minexpdurfrac_val[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Bls->minexpdurfrac_expr);
      } 
      else if(Bls->minexpdurfrac_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
	Bls->minexpdurfrac_val[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->minexpdurfrac_var);
      }
      else {
	Bls->minexpdurfrac_val[lcnum] = Bls->minexpdurfrac;
      }	  
      
      if(Bls->maxexpdurfrac_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	Bls->maxexpdurfrac_val[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Bls->maxexpdurfrac_expr);
      } 
      else if(Bls->maxexpdurfrac_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
	Bls->maxexpdurfrac_val[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->maxexpdurfrac_var);
      }
      else {
	Bls->maxexpdurfrac_val[lcnum] = Bls->minexpdurfrac;
      }
      
      Bls->rmin_val[lcnum] = pow(((0.0848203*Bls->minexpdurfrac_val[lcnum]*pow(Bls->rho_val[lcnum],(-1.0/3.0)))/0.076),1.5);
      Bls->rmax_val[lcnum] = pow(((0.0848203*Bls->maxexpdurfrac_val[lcnum]*pow(Bls->rho_val[lcnum],(-1.0/3.0)))/0.076),1.5);
    }
    
    if(!Bls->isdf_specified) {
      if(Bls->nf_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	Bls->nf_val[lcnum] = ceil(EvaluateExpression(lc_name_num, lcnum, 0, Bls->nf_expr));
      } 
      else if(Bls->nf_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
	Bls->nf_val[lcnum] = ceil(EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->nf_var));
      }
      else {
	Bls->nf_val[lcnum] = Bls->nf;
      }
      Bls->df_val[lcnum] = ((1./Bls->minper_val[lcnum]) - (1./Bls->maxper_val[lcnum])) / (Bls->nf_val[lcnum] - 1);
    } else if(Bls->isdf_specified == 1) {
      if(Bls->df_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	Bls->df_val[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Bls->df_expr);
      } 
      else if(Bls->df_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
	Bls->df_val[lcnum] = ceil(EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->df_var));
      }
      else {
	Bls->df_val[lcnum] = Bls->df;
      }
      Bls->nf_val[lcnum] = ceil(((1./Bls->minper_val[lcnum]) - (1./Bls->maxper_val[lcnum])) / (Bls->df_val[lcnum])) + 1;

    } else if(Bls->isdf_specified == 2) {
      if(Bls->subsample_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	Bls->subsample_val[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Bls->subsample_expr);
      } 
      else if(Bls->subsample_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
	Bls->subsample_val[lcnum] = ceil(EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->subsample_var));
      }
      else {
	Bls->subsample_val[lcnum] = Bls->subsample;
      }
      Bls->A_val[lcnum] = pow((3.0*M_PI/(G_CGS_DAYS * Bls->rho_val[lcnum])),(1.0/3.0))*Bls->subsample_val[lcnum]/(M_PI*(p->t[lcnum][p->NJD[lcnum]-1] - p->t[lcnum][0]));
      Bls->nf_val[lcnum] = ceil((pow((1./Bls->minper_val[lcnum]),(1.0/3.0)) - pow((1./Bls->maxper_val[lcnum]),(1.0/3.0)) + Bls->A_val[lcnum]/3.0)*3.0/Bls->A_val[lcnum]);
    }
    Bls->fmin[lcnum] = dmax((1./(p->t[lcnum][p->NJD[lcnum]-1] - p->t[lcnum][0])),1./Bls->maxper_val[lcnum]);
    if(Bls->isdf_specified == 2) {
      Bls->nf2[lcnum] = Bls->nf_val[lcnum];
    }
    else {
      if(!Bls->freqsteptype) {
	Bls->nf2[lcnum] = floor((((1./Bls->minper_val[lcnum]) - Bls->fmin[lcnum])/Bls->df_val[lcnum])+1.);
      } else if(Bls->freqsteptype == VARTOOLS_FREQSTEPTYPE_PERIOD) {
	Bls->nf2[lcnum] = floor((((1./Bls->fmin[lcnum]) - Bls->minper_val[lcnum])/Bls->df_val[lcnum])+1.);
      } else if(Bls->freqsteptype == VARTOOLS_FREQSTEPTYPE_LOGPERIOD) {
	Bls->nf2[lcnum] = floor(((log(1./Bls->fmin[lcnum]) - log(Bls->minper_val[lcnum]))/Bls->df_val[lcnum])+1.);
      }
    }

    if(Bls->nbins_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
      Bls->nbins_val[lcnum] = ceil(EvaluateExpression(lc_name_num, lcnum, 0, Bls->nbins_expr));
    } 
    else if(Bls->nbins_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
      Bls->nbins_val[lcnum] = ceil(EvaluateVariable_Double(lc_name_num, lcnum, 0, Bls->nbins_var));
    }
    else {
      Bls->nbins_val[lcnum] = Bls->nbins;
    }


    /* Now either run bls using the fixed q range or the fixed stellar radius range */
    if(Bls->nf2[lcnum] > 0 && Bls->nbins_val[lcnum] > 0 && Bls->Npeak > 0) {
#ifdef PARALLEL
      if((Bls->nf2[lcnum]+1) > Bls->sizepvec[threadindex]) {
	if(!(Bls->sizepvec[threadindex])) {
	  Bls->sizepvec[threadindex] = Bls->nf2[lcnum] + 1;
	  if((Bls->p[threadindex] = (double *) malloc(Bls->sizepvec[threadindex]*sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	} else {
	  Bls->sizepvec[threadindex] = Bls->nf2[lcnum] + 1;
	  if((Bls->p[threadindex] = (double *) realloc(Bls->p[threadindex], Bls->sizepvec[threadindex]*sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      }
#else
      if((Bls->nf2[lcnum]+1) > Bls->sizepvec) {
	if(!Bls->sizepvec) {
	  Bls->sizepvec = Bls->nf2[lcnum] + 1;
	  if((Bls->p = (double *) malloc(Bls->sizepvec*sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	} else {
	  Bls->sizepvec = Bls->nf2[lcnum] + 1;
	  if((Bls->p = (double *) realloc(Bls->p, Bls->sizepvec*sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      }
#endif

      if(!Bls->rflag)
	{

	  eebls(p->NJD[lcnum],p->t[lcnum],p->mag[lcnum],p->sig[lcnum],Bls->u[lcnum],Bls->v[lcnum],Bls->nf2[lcnum],Bls->fmin[lcnum],Bls->df_val[lcnum],Bls->nbins_val[lcnum],Bls->qmin_val[lcnum],Bls->qmax_val[lcnum],
#ifdef PARALLEL
		Bls->p[threadindex]
#else
		Bls->p
#endif
		,Bls->Npeak,Bls->bper[lcnum],Bls->bt0[lcnum],Bls->bpow[lcnum],Bls->sde[lcnum],Bls->snval[lcnum],Bls->depth[lcnum],Bls->qtran[lcnum],Bls->i1[lcnum],Bls->i2[lcnum],Bls->i1_ph[lcnum],Bls->i2_ph[lcnum],Bls->chisqrplus[lcnum],&Bls->chisqrminus[lcnum],&Bls->bperpos[lcnum],&Bls->meanmagval[lcnum], Bls->timezone, Bls->fraconenight[lcnum], Bls->operiodogram, outname, Bls->omodel, outname2, Bls->correctlc,p->ascii, Bls->nt[lcnum], Bls->Nt[lcnum], Bls->Nbefore[lcnum], Bls->Nafter[lcnum], Bls->rednoise[lcnum], Bls->whitenoise[lcnum], Bls->sigtopink[lcnum], Bls->fittrap, Bls->qingress[lcnum], Bls->OOTmag[lcnum], Bls->ophcurve, outname3, Bls->phmin, Bls->phmax, Bls->phstep, Bls->ojdcurve, outname4, Bls->jdstep, Bls->nobinnedrms, Bls->freqsteptype, Bls->adjust_qmin_mindt, Bls->reduce_nb, Bls->reportharmonics, Bls, lcnum, lc_name_num, Bls->usemask, Bls->maskvar);
	}
      else
	{
	  
	  eebls_rad(p->NJD[lcnum],p->t[lcnum],p->mag[lcnum],p->sig[lcnum],Bls->u[lcnum],Bls->v[lcnum],Bls->nf2[lcnum],Bls->fmin[lcnum],Bls->df_val[lcnum],Bls->nbins_val[lcnum],Bls->rmin_val[lcnum],Bls->rmax_val[lcnum],
#ifdef PARALLEL
			  Bls->p[threadindex]
#else
			  Bls->p
#endif
			  ,Bls->Npeak,Bls->bper[lcnum],Bls->bt0[lcnum],Bls->bpow[lcnum],Bls->sde[lcnum],Bls->snval[lcnum],Bls->depth[lcnum],Bls->qtran[lcnum],Bls->i1[lcnum],Bls->i2[lcnum],Bls->i1_ph[lcnum],Bls->i2_ph[lcnum],Bls->chisqrplus[lcnum],&Bls->chisqrminus[lcnum],&Bls->bperpos[lcnum],&Bls->meanmagval[lcnum], Bls->timezone, Bls->fraconenight[lcnum], Bls->operiodogram,outname, Bls->omodel, outname2, Bls->correctlc,p->ascii, Bls->nt[lcnum], Bls->Nt[lcnum], Bls->Nbefore[lcnum], Bls->Nafter[lcnum], Bls->rednoise[lcnum], Bls->whitenoise[lcnum], Bls->sigtopink[lcnum], Bls->fittrap, Bls->qingress[lcnum], Bls->OOTmag[lcnum], Bls->ophcurve, outname3, Bls->phmin, Bls->phmax, Bls->phstep, Bls->ojdcurve, outname4, Bls->jdstep, Bls->nobinnedrms,Bls->freqsteptype, Bls->adjust_qmin_mindt, Bls->reduce_nb, Bls->reportharmonics, Bls, lcnum, lc_name_num, Bls->usemask, Bls->maskvar, (Bls->isdf_specified == 2), (Bls->isdf_specified == 2 ? Bls->A_val[lcnum] : 0.0));
	      }
	  } else {
	    if(!p->quiet_mode) {
	      fprintf(stderr,"Warning: skipping -BLS command index %d for light curve number: %d, filename: %s. The light curve is either too short, or an invalid set of parameter options were supplied to BLS.\n", thisindex, lc_name_num, p->lcnames[lc_name_num]);
	    }
	  }
	} else {
	    if(!p->quiet_mode) {
	      fprintf(stderr,"Warning: skipping -BLS command index %d for light curve number: %d, filename: %s. The light curve has too few points for BLS.\n", thisindex, lc_name_num, p->lcnames[lc_name_num]);
	    }
      }
}
