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

/*

This file contains functions for computing the AoV periodogram and period peaks for the program vartools by J. Hartman

The code is adapted from the periodS2 program by J. Devor which is included in his DEBiL package, and is also adapted from the AovPer.c code by A. Schwarzenberg-Czerny which is available on his website: http://users.camk.edu.pl/alex/

*/

/* AOV Period Search Functions */

#define ERROR_SCORE 100000.0  // assumed to be bigger than any valid score (>> 1)
#define MAX_DOUBLE_CHECK_MULTIPLE 19
#define MAX_PERIOD_DIFF_MULTIPLE 5
#define DEFAULT_HIST_SIZE 8

typedef struct {
  double HIST_SIZE;
  int sizehists;
  unsigned long *histN;
  double *histA, *histB, *histC, *histD, *histE;
} _HistType;

  //double HIST_SIZE;
  //unsigned long *histN;
  //double *histA, *histB, *histC, *histD, *histE;

double MIN_PERIOD;

/* Normalize the light curve to have 0 average and unit standard deviation */
void normalize (int size, double *t, double *mag, double *sig, double *ave, double *stddev)
{
  int i, i1, n;
  if(size > 0)
    {
      i=0;
      while(isnan(mag[i]))
	i++;
      *ave = mag[i];
      n = 1;
      i1 = i;
      for(i = i+1; i<size; i++)
	if(!isnan(mag[i]))
	  {
	    *ave += mag[i];
	    n++;
	  }
      *ave /= (double) n;
      *stddev = SQR(mag[i1] - *ave);
      for(i=i1 + 1; i<size; i++)
	if(!isnan(mag[i]))
	  *stddev += SQR(mag[i] - *ave);
      *stddev = sqrt((*stddev) / (double) n);
      for(i = 0; i < size; i++)
	{
	  mag[i] = (mag[i] - *ave) / (*stddev);
	  sig[i] = sig[i] / (*stddev);
	}
    }
}

void whitenlc_aov(int size, double *t, double *mag, double period, _HistType *h)
{
  int i, index;
  double sum, X, Y, s1, L2;
  int N;

  sum = 0.;
  s1 = 0.;
  L2 = 0.;

  /* initialize the histograms */
  memset (h->histN, 0, h->HIST_SIZE * sizeof (unsigned long)) ;
  memset (h->histA, 0, h->HIST_SIZE * sizeof (double)) ;
  memset (h->histB, 0, h->HIST_SIZE * sizeof (double)) ;
  memset (h->histC, 0, h->HIST_SIZE * sizeof (double)) ;
  memset (h->histD, 0, h->HIST_SIZE * sizeof (double)) ;
  memset (h->histE, 0, h->HIST_SIZE * sizeof (double)) ;


  for (i = 0; i < size ; i++)
    {
      if(!isnan(mag[i]))
	{
	  X = fmod(t[i], period) / period ;
	  index = (int)(h->HIST_SIZE * X) ;
	  Y = mag[i] ;
	  h->histN[index]++ ;
	  h->histD[index] += Y;
	}
    }
  for(i = 0; i < h->HIST_SIZE; i++)
    {
      h->histD[i] = h->histD[i] / (double) h->histN[i];
    }
  for (i = 0; i < size ; i++)
    {
      if(!isnan(mag[i]))
	{
	  X = fmod(t[i], period) / period ;
	  index = (int)(h->HIST_SIZE * X) ;
	  mag[i] -= h->histD[index];
	}
    }
  return;
}

double TestPeriod(int size, double *t, double *mag, double period, int uselog, _HistType *h)
{
  int i, index;
  double sum, X, Y, s1, L2;
  int N;

  sum = 0.;
  s1 = 0.;
  L2 = 0.;

  /* initialize the histograms */
  memset (h->histN, 0, h->HIST_SIZE * sizeof (unsigned long)) ;
  memset (h->histA, 0, h->HIST_SIZE * sizeof (double)) ;
  memset (h->histB, 0, h->HIST_SIZE * sizeof (double)) ;
  memset (h->histC, 0, h->HIST_SIZE * sizeof (double)) ;
  memset (h->histD, 0, h->HIST_SIZE * sizeof (double)) ;
  memset (h->histE, 0, h->HIST_SIZE * sizeof (double)) ;

  for (i = 0; i < size ; i++)
    {
      if(!isnan(mag[i]))
	{
	  X = fmod(t[i], period) / period ;
	  index = (int)(h->HIST_SIZE * X) ;
	  Y = mag[i] ;
	  h->histN[index]++ ;
	  h->histA[index] += X ;
	  h->histB[index] += X * X;
	  h->histC[index] += X * Y;
	  h->histD[index] += Y;
	  h->histE[index] += Y * Y;
	}
    }
  for(i = 0; i < h->HIST_SIZE; i++)
    {
      N = h->histN[i];
      if (N <= 1)
	return(ERROR_SCORE);
      Y = h->histD[i];
      X = Y * Y / N;

      sum += Y;
      s1 += X ;
      L2 += (SQR((h->histC[i] * N) - (Y * h->histA[i])) / (N * (SQR(h->histA[i]) - (h->histB[i] * N)))) + h->histE[i] - X;
    }
  s1 -= sum * sum / size ;
  L2 /= size - h->HIST_SIZE;
  if ((s1 <= 0.0) || (L2 <= 0.0))
    return(ERROR_SCORE);
  else if(uselog)
    return(log(L2 / s1));
  else
    return(-s1/L2);
}

#define MINVAR 1.e-32
#define CTMIN 5
/* Compute the AoV Periodogram using A.S-Cz's implementation */
void AOVPeriodogram_ASZ(int size, double *tin, double *mag, int Nperiod, double *periods, double *periodogram, _HistType *h)
{
  int *ncnt, *ind, i, ibin, ip, nbc, ncov;
  long ifr, iflex;
  double *f, *ph, *ave, af, vf, sav;
  double *t, fr, at, dbc, dph;
  int MAXBIN, nbin, nobs;

  ncov = 1;
  MAXBIN = (h->HIST_SIZE + 1)*ncov;
  nbin = h->HIST_SIZE;
  nobs = size;

  ncnt = (int *) malloc((h->HIST_SIZE + 1) * sizeof(int));
  ind = (int *) malloc(size * sizeof(int));
  f = (double *) malloc(size * sizeof(double));
  ph = (double *) malloc(size * sizeof(double));
  ave = (double *) malloc((h->HIST_SIZE + 1) * sizeof(double));
  t = (double *) malloc(size * sizeof(double));

/*   Set variables (incl. mean and variance) */
  nbc = nbin * ncov;
  dbc = (double) nbc;

/*  calculate totals and normalize variables */
  iflex = 0; at = (double) (af = vf = (double)0.);
  for (i = 0; i < nobs; i++) { af += mag[i]; at += tin[i]; }
  af /= (double) nobs; at /= (double) nobs;
  for (i = 0; i < nobs; i++)
    {
      t[i] = tin[i] - at;
      f[i] = (sav = mag[i] - af);
      vf += sav*sav;
    };

/*  Loop over all frequencies  */
  for (ifr = 0; ifr < Nperiod; ifr++)
	{
	  fr = 1. /periods[ifr];
	  /*  Up to two passes over all frequencies, depending on phase coverage  */
	  for (ip = 0; ip < 2 ; ip++)
	    {

	      for (i = 0; i< nbc; i++) { ave[i] = (double)0.; ncnt[i] = 0; };

	      if ( ip == 0) /*   Default fixed size phase bins */
		for (i = 0; i < nobs; i++)
		  {
		    dph=t[i]*fr; /* only dph, t0, fr, must keep TIME precision */
		    ph[i]=(sav=(double)(dph-floor(dph)));
		    ibin=(int)floor(sav*dbc);
		    ave[ibin] += f[i];
		    ++ncnt[ibin];
		  }
	      else /*   For poor phase coverage optional flexible bins  */
		{
		  /* sort index ind using key ph0, any index sort routine would do */
		  ++iflex;
		  for(i=0;i<nobs;i++)
		    ind[i] = i;
		  mysort2dblint(nobs,ph,ind);
		  for (i = 0; i < nobs; i++)
		    {
		      ibin=i*nbc/nobs;
		      ave[ibin] += f[ind[i]];
		      ++ncnt[ibin];
		    }
		}

	      for (i=0; i<ncov; i++) ncnt[i+nbc]=ncnt[i];
	      ibin=0; for (i=ncov+nbc-1;i>=0;i--) ncnt[i]=(ibin+=ncnt[i]);
	      for (i=0;i<nbc;i++) ncnt[i]-=ncnt[i+ncov];
	      for (i = 0; i < nbc ; i++)
		if (ncnt[i] < CTMIN) break;

	      if (i>=nbc) break;
	    }

	  /*    Calculate A.O.V. statistics for a given frequency */
	  for (i=0; i<ncov; i++) ave[i+nbc]=ave[i];
	  sav=(double)0.; for (i=ncov+nbc-1;i>=0;i--) ave[i]=(sav+=ave[i]);
	  for (i=0;i<nbc;i++) ave[i]-=ave[i+ncov];

	  sav=(double)0.; for (i=0;i<nbc;i++)  sav+=(ave[i]*ave[i]/ncnt[i]);
	  sav/=(double)ncov;
	  periodogram[ifr] = -sav/(nbin-1)/MAX_(vf-sav,MINVAR)*(nobs-nbin);
	}

  /* if (iflex > 0) fprintf(stderr,
     "AOV:warning: poor phase coverage at %d frequencies\n",iflex); */
  free(ncnt);
  free(ind);
  free(f);
  free(ph);
  free(ave);
  free(t);

}

/* Return ASZ's AoV at test period 'period'. */
double TestPeriod_ASZ(int size, double *tin, double *mag, double period, _HistType *h)
{
  int *ncnt, *ind, i, ibin, ip, nbc, ncov;
  long ifr, iflex;
  double *f, *ph, *ave, af, vf, sav;
  double *t, fr, at, dbc, dph;
  int MAXBIN, nbin, nobs;
  double periodogram;

  ncov = 1;
  MAXBIN = (h->HIST_SIZE + 1)*ncov;
  nbin = h->HIST_SIZE;
  nobs = size;

  ncnt = (int *) malloc((h->HIST_SIZE + 1) * sizeof(int));
  ind = (int *) malloc(size * sizeof(int));
  f = (double *) malloc(size * sizeof(double));
  ph = (double *) malloc(size * sizeof(double));
  ave = (double *) malloc((h->HIST_SIZE + 1) * sizeof(double));
  t = (double *) malloc(size * sizeof(double));

/*   Set variables (incl. mean and variance) */
  nbc = nbin * ncov;
  dbc = (double) nbc;

/*  calculate totals and normalize variables */
  iflex = 0; at = (double) (af = vf = (double)0.);
  for (i = 0; i < nobs; i++) { af += mag[i]; at += tin[i]; }
  af /= (double) nobs; at /= (double) nobs;
  for (i = 0; i < nobs; i++)
    {
      t[i] = tin[i] - at;
      f[i] = (sav = mag[i] - af);
      vf += sav*sav;
    };

  fr = 1. /period;
  /*  Up to two passes over all frequencies, depending on phase coverage  */

  for (ip = 0; ip < 2 ; ip++)
    {

      for (i = 0; i< nbc; i++) { ave[i] = (double)0.; ncnt[i] = 0; };

      if ( ip == 0) /*   Default fixed size phase bins */
	for (i = 0; i < nobs; i++)
	  {
	    dph=t[i]*fr; /* only dph, t0, fr, must keep TIME precision */
	    ph[i]=(sav=(double)(dph-floor(dph)));
	    ibin=(int)floor(sav*dbc);
	    ave[ibin] += f[i];
	    ++ncnt[ibin];
	  }
      else /*   For poor phase coverage optional flexible bins  */
	{
	  /* sort index ind using key ph0, any index sort routine would do */
	  ++iflex;
	  for(i=0;i<nobs;i++)
	    ind[i] = i;
	  mysort2dblint(nobs,ph,ind);
	  for (i = 0; i < nobs; i++)
	    {
	      ibin=i*nbc/nobs;
	      ave[ibin] += f[ind[i]];
	      ++ncnt[ibin];
	    }
	}

      for (i=0; i<ncov; i++) ncnt[i+nbc]=ncnt[i];
      ibin=0; for (i=ncov+nbc-1;i>=0;i--) ncnt[i]=(ibin+=ncnt[i]);
      for (i=0;i<nbc;i++) ncnt[i]-=ncnt[i+ncov];
      for (i = 0; i < nbc ; i++)
	if (ncnt[i] < CTMIN) break;

      if (i>=nbc) break;
    }

  /*    Calculate A.O.V. statistics for a given frequency */
  for (i=0; i<ncov; i++) ave[i+nbc]=ave[i];
  sav=(double)0.; for (i=ncov+nbc-1;i>=0;i--) ave[i]=(sav+=ave[i]);
  for (i=0;i<nbc;i++) ave[i]-=ave[i+ncov];

  sav=(double)0.; for (i=0;i<nbc;i++)  sav+=(ave[i]*ave[i]/ncnt[i]);
  sav/=(double)ncov;
  periodogram = -sav/(nbin-1)/MAX_(vf-sav,MINVAR)*(nobs-nbin);

  /* if (iflex > 0) fprintf(stderr,
     "AOV:warning: poor phase coverage at %d frequencies\n",iflex); */
  free(ncnt);
  free(ind);
  free(f);
  free(ph);
  free(ave);
  free(t);

  return(periodogram);

}



/* Compute the AOV Periodogram using J. Devor's implementation */
void AOVPeriodogram(int size, double *t, double *mag, int Nperiod, double *periods, double *periodogram, int uselog, _HistType *h)
{
  int i, index, k;
  double sum, X, Y, s1, L2, period;
  int N;

  for(k = 0; k < Nperiod; k++)
    {
      sum = 0.;
      s1 = 0.;
      L2 = 0.;
      period = periods[k];

      /* initialize the histograms */
      memset (h->histN, 0, h->HIST_SIZE * sizeof (unsigned long)) ;
      memset (h->histA, 0, h->HIST_SIZE * sizeof (double)) ;
      memset (h->histB, 0, h->HIST_SIZE * sizeof (double)) ;
      memset (h->histC, 0, h->HIST_SIZE * sizeof (double)) ;
      memset (h->histD, 0, h->HIST_SIZE * sizeof (double)) ;
      memset (h->histE, 0, h->HIST_SIZE * sizeof (double)) ;

      for (i = 0; i < size ; i++)
	{
	  if(!isnan(mag[i]))
	    {
	      X = fmod(t[i], period) / period ;
	      index = (int)(h->HIST_SIZE * X) ;
	      Y = mag[i] ;
	      h->histN[index]++ ;
	      h->histA[index] += X ;
	      h->histB[index] += X * X;
	      h->histC[index] += X * Y;
	      h->histD[index] += Y;
	      h->histE[index] += Y * Y;
	    }
	}
      for(i = 0; i < h->HIST_SIZE; i++)
	{
	  N = h->histN[i];
	  if (N <= 1)
	    periodogram[k] = ERROR_SCORE;
	  Y = h->histD[i];
	  X = Y * Y / N;

	  sum += Y;
	  s1 += X ;
	  L2 += (SQR((h->histC[i] * N) - (Y * h->histA[i])) / (N * (SQR(h->histA[i]) - (h->histB[i] * N)))) + h->histE[i] - X;
	}
      s1 -= sum * sum / size ;
      L2 /= size - h->HIST_SIZE;
      if ((s1 <= 0.0) || (L2 <= 0.0))
	periodogram[k] = ERROR_SCORE;
      else if(uselog)
	periodogram[k] = (log(L2 / s1));
      else
	periodogram[k] = -s1 / L2;
    }
}



/* Determine whether two different periods are the same */
int isDifferentPeriods (double period1, double period2, double T)
{
  int a, b ;
  double period1mul ;

  /*if (T * (period2 - period1) < ((period2 * period2) + (period1 * period1)))
    return (0) ;*/
  if(T * fabs(period2 - period1) < period2*period1)
    return(0);

  for (a = 1 ; a < MAX_PERIOD_DIFF_MULTIPLE ; a++)
    for (b = a + 1 ; b <= MAX_PERIOD_DIFF_MULTIPLE ; b++)   // a < b
      {
	period1mul = period1 * b / a ;
	/*if (T * fabs(period2 - period1mul) < ((period2 * period2) + (period1mul * period1mul)))
	  return (0) ;*/
	if(T * fabs(period2 - period1mul) < period2 * period1mul)
	  return(0);
      }

  return (1) ;
}

/* Determine whether two different periods are the same, dont check multiples */
int isDifferentPeriodsDontCheckHarmonics (double period1, double period2, double T)
{

  /*if (T * (period2 - period1) < ((period2 * period2) + (period1 * period1)))
    return (0) ;*/
  if(T * fabs(period2 - period1) < period2*period1)
    return(0);
  else
    return (1);
}

/* Given a light curve, this function will compute an AOV periodogram and find the top Npeaks peaks */
void findPeaks_aov(double *t_, double *mag_, double *sig_, int N, double *perpeaks, double *aovpeaks, double *aovSNR, double *aovFAP, int Npeaks, double minP, double maxP, double subsample, double fine_tune, int outflag, char *outname, double *aveaov, double *stddevaov, double *aveaov_whiten, double *stddevaov_whiten, int ascii, int Nbin, int whiten, int uselog, double clip, int clipiter, int fixperiodSNR, double fixperiodSNR_period, double *fixperiodSNR_value, double *fixperiodSNR_SNR, double *fixperiodSNR_FAP, int lcnum, int lclistnum, int usemask, _Variable *maskvar)
{
  int i, j, k, peakiter, foundsofar, test, Nperiod, a, b, abest, bbest, ismultiple, Ngood, nclippedthis, nclippedlast, m_eff, sizeHISTvector = 0;

  double *periods = NULL, *periodogram = NULL, **periodogram_whiten = NULL;
  int size_periodogram = 0, size_periodogram_whiten = 0;

  double *aveper_whiten = NULL;
  double *stdper_whiten = NULL;
  int size_aveper_whiten = 0;
  int Norig;

  double *t, *mag, *sig, a_, b_;
  double testperiod, bestscore, lastpoint;
  double T, ave, stddev, freq, minfreq, minbest, freqstep, smallfreqstep, score, aveper, stdper, negln_m_eff;

  long double Sum, Sumsqr;

  FILE *outfile;

  _HistType h;

  if(N <= 1) {
    for(i = 0; i < Npeaks; i++) {
      perpeaks[i] = -1.0;
      aovpeaks[i] = 0.0;
      aovSNR[i] = 0.0;
      aovFAP[i] = 0.0;
    }
    *aveaov = 0.0;
    *stddevaov = 0.0;
    if(whiten) {
      *aveaov_whiten = 0.0;
      *stddevaov_whiten = 0.0;
    }
    if(fixperiodSNR) {
      *fixperiodSNR_value = 0.0;
      *fixperiodSNR_SNR = 0.0;
      *fixperiodSNR_FAP = 0.0;
    }
    return;
  }

  if(!size_aveper_whiten)
    {
      aveper_whiten = (double *) malloc((Npeaks + 1) * sizeof(double));
      stdper_whiten = (double *) malloc((Npeaks + 1) * sizeof(double));
      size_aveper_whiten = Npeaks + 1;
    }

  /* Prepare the HIST_SIZE vectors */
  if(Nbin > sizeHISTvector || (Nbin <= 0 && sizeHISTvector < DEFAULT_HIST_SIZE))
    {
      if(!sizeHISTvector)
	{
	  if(Nbin <= 0)
	    {
	      h.HIST_SIZE = DEFAULT_HIST_SIZE;
	      sizeHISTvector = DEFAULT_HIST_SIZE;
	    }
	  else
	    {
	      h.HIST_SIZE = Nbin;
	      sizeHISTvector = Nbin;
	    }
	  if((h.histN = (unsigned long *) malloc(sizeHISTvector * sizeof(unsigned long))) == NULL ||
	     (h.histA = (double *) malloc(sizeHISTvector * sizeof(double))) == NULL ||
	     (h.histB = (double *) malloc(sizeHISTvector * sizeof(double))) == NULL ||
	     (h.histC = (double *) malloc(sizeHISTvector * sizeof(double))) == NULL ||
	     (h.histD = (double *) malloc(sizeHISTvector * sizeof(double))) == NULL ||
	     (h.histE = (double *) malloc(sizeHISTvector * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      else
	{
	  if(Nbin <= 0)
	    {
	      h.HIST_SIZE = DEFAULT_HIST_SIZE;
	      sizeHISTvector = DEFAULT_HIST_SIZE;
	    }
	  else
	    {
	      h.HIST_SIZE = Nbin;
	      sizeHISTvector = Nbin;
	    }
	  if((h.histN = (unsigned long *) realloc(h.histN, sizeHISTvector * sizeof(unsigned long))) == NULL ||
	     (h.histA = (double *) realloc(h.histA, sizeHISTvector * sizeof(double))) == NULL ||
	     (h.histB = (double *) realloc(h.histB, sizeHISTvector * sizeof(double))) == NULL ||
	     (h.histC = (double *) realloc(h.histC, sizeHISTvector * sizeof(double))) == NULL ||
	     (h.histD = (double *) realloc(h.histD, sizeHISTvector * sizeof(double))) == NULL ||
	     (h.histE = (double *) realloc(h.histE, sizeHISTvector * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
    }
  else
    {
      if(Nbin <= 0)
	h.HIST_SIZE = DEFAULT_HIST_SIZE;
      else
	h.HIST_SIZE = Nbin;
    }

  /* copy the t mag and sig vectors */

  if((t = (double *) malloc(N * sizeof(double))) == NULL ||
     (mag = (double *) malloc(N * sizeof(double))) == NULL ||
     (sig = (double *) malloc(N * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  if(!usemask) {
    for(i=0;i<N;i++)
      {
	t[i] = t_[i];
	mag[i] = mag_[i];
	sig[i] = sig_[i];
      }
  } else {
    Norig = N;
    N = 0;
    for(i = 0; i < Norig; i++) {
      if(!isnan(mag_[i]) && sig_[i] > 0. && EvaluateVariable_Double(lclistnum, lcnum, i, maskvar) > VARTOOLS_MASK_TINY) {
	t[N] = t_[i];
	mag[N] = mag_[i];
	sig[N] = sig_[i];
	N++;
      }
    }
    if(N <= 1) {
      for(i = 0; i < Npeaks; i++) {
	perpeaks[i] = -1.0;
	aovpeaks[i] = 0.0;
	aovSNR[i] = 0.0;
	aovFAP[i] = 0.0;
      }
      *aveaov = 0.0;
      *stddevaov = 0.0;
      if(whiten) {
	*aveaov_whiten = 0.0;
	*stddevaov_whiten = 0.0;
      }
      if(fixperiodSNR) {
	*fixperiodSNR_value = 0.0;
	*fixperiodSNR_SNR = 0.0;
	*fixperiodSNR_FAP = 0.0;
      }
      free(t);
      free(mag);
      free(sig);
      if(sizeHISTvector > 0) {
	free(h.histN);
	free(h.histA);
	free(h.histB);
	free(h.histC);
	free(h.histD);
	free(h.histE);
      }
      
      if(periods != NULL)
	free(periods);
      
      if(periodogram != NULL)
	free(periodogram);
      
      if(periodogram_whiten != NULL) {
	for(i=0; i < Npeaks + 1; i++)
	  free(periodogram_whiten[i]);
	free(periodogram_whiten);
      }
      
      if(aveper_whiten != NULL)
	free(aveper_whiten);
      
      if(stdper_whiten != NULL)
	free(stdper_whiten);
      return;
    }
  }

  /* initialize some of the period search variables*/
  T = t[N - 1] - t[0];

  /* normalize the light curve to have 0 average and unit standard deviation */
  normalize(N, t, mag, sig, &ave, &stddev);

  // Initialize the periodogram if it hasn't already been initialized

  Nperiod = 0;
  freq = 1./minP;
  minfreq = 1./maxP;
  freqstep = subsample / T;
  while(freq >= minfreq)
    {
      freq -= freqstep;
      Nperiod++;
    }

  /* Set the estimate for the effective number of independent frequency samples */
  m_eff = (int) ((1./minP - 1./maxP) * T);
  if(subsample > 1.)
    m_eff = (int) ((double) m_eff / subsample);

  negln_m_eff = -log((double) m_eff);

  if(!size_periodogram && !whiten)
    {
      if(!size_periodogram && !size_periodogram_whiten)
	{
	  if((periods = (double *) malloc(Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      else if(Nperiod > MAX_(size_periodogram_whiten,size_periodogram))
	{
	  if((periods = (double *) realloc(periods, Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      size_periodogram = Nperiod;
      if((periodogram = (double *) malloc((Nperiod) * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  else if(size_periodogram < Nperiod && !whiten)
    {
      if(Nperiod > MAX_(size_periodogram_whiten,size_periodogram))
	{
	  if((periods = (double *) realloc(periods, Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      size_periodogram = Nperiod;
      if((periodogram = (double *) realloc(periodogram, Nperiod * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  if(!size_periodogram_whiten && whiten)
    {
      if(!size_periodogram && !size_periodogram_whiten)
	{
	  if((periods = (double *) malloc(Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      else if(Nperiod > MAX_(size_periodogram_whiten,size_periodogram))
	{
	  if((periods = (double *) realloc(periods, Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      size_periodogram_whiten = Nperiod;
      if((periodogram_whiten = (double **) malloc((Npeaks + 1) * sizeof(double *))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0;i<Npeaks+1;i++)
	{
	  if((periodogram_whiten[i] = (double *) malloc(Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
    }
  else if(size_periodogram_whiten < Nperiod && whiten)
    {
      if(Nperiod > MAX_(size_periodogram_whiten,size_periodogram))
	{
	  if((periods = (double *) realloc(periods, Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
      size_periodogram_whiten = Nperiod;
      for(i=0;i<Npeaks+1;i++)
	{
	  if((periodogram_whiten[i] = (double *) realloc(periodogram_whiten[i], Nperiod * sizeof(double))) == NULL)
	    error(ERR_MEMALLOC);
	}
    }
  freq = 1./minP;
  minfreq = 1./maxP;
  i = 0;
  while(freq >= minfreq)
    {
      periods[i] = 1./freq;
      freq -= freqstep;
      i++;
    }

  /* Here we do things differently depending on whether or not we're going to whiten the signal between finding peaks */
  if(!whiten)
    {
      /* Compute the periodogram */
      if(uselog)
	AOVPeriodogram(N, t, mag, Nperiod, periods, periodogram, uselog, &h);
      else
	AOVPeriodogram_ASZ(N, t, mag, Nperiod, periods, periodogram, &h);

      /* Compute the average and std-dev of the periodogram */
      Sum = 0.;
      Sumsqr = 0.;
      Ngood = 0;
      for(i=0;i<Nperiod;i++)
	{
	  if(periodogram[i] < ERROR_SCORE && (periodogram[i]*0.0 == 0.0))
	    {
	      Sum += (long double) periodogram[i];
	      Sumsqr += (long double) (periodogram[i] * periodogram[i]);
	      Ngood++;
	    }
	}

      Sum /= Ngood;
      Sumsqr /= Ngood;

      aveper = (double) Sum;
      stdper = sqrt((double)(Sumsqr - (Sum*Sum)));

      nclippedthis = 0;
      do {
	nclippedlast = nclippedthis;
	nclippedthis = 0;

	Sum = 0.;
	Sumsqr = 0.;
	Ngood = 0;
	for(i=0;i<Nperiod;i++)
	  {
	    if(periodogram[i] < ERROR_SCORE && (periodogram[i]*0.0 == 0.0))
	      {
		if(periodogram[i] > aveper - clip*stdper)
		  {
		    Sum += (long double) periodogram[i];
		    Sumsqr += (long double) (periodogram[i] * periodogram[i]);
		    Ngood++;
		  }
		else
		  nclippedthis++;
	      }
	  }

	Sum /= Ngood;
	Sumsqr /= Ngood;

	aveper = (double) Sum;
	stdper = sqrt((double)(Sumsqr - (Sum*Sum)));
      } while(clipiter && nclippedthis > nclippedlast);

      *aveaov = aveper;
      *stddevaov = stdper;

      /* Get the peak value and SNR for the fixed period if we're doing that */
      if(fixperiodSNR)
	{
	  AOV_getfixedperiodSNR(periods,periodogram,Nperiod,aveper,stdper,fixperiodSNR_period,fixperiodSNR_value,fixperiodSNR_SNR);
	  if(uselog)
	    *fixperiodSNR_value = (aveper - *fixperiodSNR_value) / stdper;
	  else
	    {
	      *fixperiodSNR_SNR = (aveper - *fixperiodSNR_value) / stdper;
	      *fixperiodSNR_value = -(*fixperiodSNR_value);
	      a_ = 0.5*((double) (h.HIST_SIZE - 1));
	      b_ = 0.5*((double) (N - h.HIST_SIZE));
	      *fixperiodSNR_FAP = -log1minusbetai(a_, b_, ((2. * a_ * (*fixperiodSNR_value))/(2. * (b_ + a_*(*fixperiodSNR_value)))));
	    }
	}

      /* Write out the periodogram if asked to */
      if(outflag)
	{
	  if((outfile = fopen(outname,"w")) == NULL)
	    {
	      fprintf(stderr,"Cannot Write periodogram to %s\n",outname);
	      exit(3);
	    }
	  if(ascii)
	    {
	      fprintf(outfile,"#Period");
	      if(uselog) {
		fprintf(outfile," AOV_LOGSNR\n");
	      }
	      else
		fprintf(outfile," AOV\n");
	      for(i=0;i<Nperiod;i++)
		{
		  if(periodogram[i] < ERROR_SCORE && 0.0*periodogram[i] == 0.0)
		    {
		      if(uselog)
			fprintf(outfile,"%f %f\n",periods[i], ((aveper - periodogram[i]) / stdper));
		      else
			fprintf(outfile,"%f %f\n",periods[i], -periodogram[i]);
		    }
		}
	    }
	  else
	    {
	      fwrite(&Nperiod,4,1,outfile);
	      fwrite(&aveper,8,1,outfile);
	      fwrite(&stdper,8,1,outfile);
	      fwrite(periods,8,Nperiod,outfile);
	      fwrite(periodogram,8,Nperiod,outfile);
	    }

	  fclose(outfile);
	}

      /* Replace the periodogram with only points that are local minima */
      lastpoint = periodogram[0] - 1.;
      for(i=0,k=0;k<Nperiod-1;k++)
	{
	  if(periodogram[k] < lastpoint && periodogram[k] < periodogram[k + 1])
	    {
	      lastpoint = periodogram[k];
	      periodogram[i] = periodogram[k];
	      periods[i] = periods[k];
	      i++;
	    }
	  else
	    {
	      lastpoint = periodogram[k];
	    }
	}
      if(periodogram[k] < lastpoint)
	{
	  periodogram[i] = periodogram[k];
	  periods[i] = periods[k];
	  i++;
	}
      Nperiod = i;

      /* Search through the periodogram to identify the best Npeaks periods */

      foundsofar = 0;
      i = 0;
      while(foundsofar < Npeaks && i < Nperiod)
	{
	  if(periodogram[i] < ERROR_SCORE && (periodogram[i]*0.0 == 0.0))
	    {
	      test = 1;
	      for(j=0;j<foundsofar;j++)
		{
		  if(!isDifferentPeriods(perpeaks[j],periods[i],T))
		    {
		      if(periodogram[i] < aovpeaks[j])
			{
			  perpeaks[j] = periods[i];
			  aovpeaks[j] = periodogram[i];
			}
		      test = 0;
		      break;
		    }
		}
	      if(test)
		{
		  perpeaks[foundsofar] = periods[i];
		  aovpeaks[foundsofar] = periodogram[i];
		  foundsofar++;
		}
	    }
	  i++;
	}

      if(foundsofar < Npeaks)
	{
	  for(k=foundsofar;k<Npeaks;k++)
	    {
	      perpeaks[k] = ERROR_SCORE + 1.;
	      aovpeaks[k] = ERROR_SCORE + 1.;
	    }
	}

      mysort2(Npeaks,aovpeaks,perpeaks);

      minbest = aovpeaks[Npeaks - 1];

      for(;i<Nperiod;i++)
	{
	  if(periodogram[i] < ERROR_SCORE && periodogram[i]*0.0 == 0.0)
	    {
	      if(periodogram[i] < minbest)
		{
		  test = 1;
		  for(j=0;j<Npeaks;j++)
		    {
		      if(!isDifferentPeriods(periods[i],perpeaks[j],T))
			{
			  if(periodogram[i] < aovpeaks[j])
			    {
			      aovpeaks[j] = periodogram[i];
			      perpeaks[j] = periods[i];
			      mysort2(Npeaks,aovpeaks,perpeaks);
			      minbest = aovpeaks[Npeaks - 1];
			    }
			  test = 0;
			  break;
			}
		    }
		  if(test)
		    {
		      perpeaks[Npeaks - 1] = periods[i];
		      aovpeaks[Npeaks - 1] = periodogram[i];
		      mysort2(Npeaks,aovpeaks,perpeaks);
		      minbest = aovpeaks[Npeaks - 1];
		    }
		}
	    }
	}

      /* Now perform the high-resolution period scan on the peaks */

      smallfreqstep = fine_tune/T;

      for(j=0;j<Npeaks;j++)
	{
	  freq = dmin((1./perpeaks[j]) + freqstep,(1./minP));
	  minfreq = dmax((1./perpeaks[j]) - freqstep,(1./maxP));
	  while (freq >= minfreq)
	    {
	      testperiod = 1./freq;
	      if(uselog)
		score = TestPeriod(N, t, mag,testperiod, uselog,&h);
	      else
		score = TestPeriod_ASZ(N, t, mag, testperiod,&h);
	      if(score < ERROR_SCORE && (score * 0.0 == 0.0))
		{
		  if(score < aovpeaks[j])
		    {
		      aovpeaks[j] = score;
		      perpeaks[j] = testperiod;
		    }
		}
	      freq -= smallfreqstep;
	    }
	  // double-check period multiples
	  bestscore = aovpeaks[j];
	  ismultiple = 0;
	  for(a=1 ; a <= MAX_DOUBLE_CHECK_MULTIPLE ; a++)
	    for(b=1 ; b <= MAX_DOUBLE_CHECK_MULTIPLE ; b++)
	      if(a != b)
		{
		  testperiod = perpeaks[j] * a / b;
		  if((testperiod > minP) && (testperiod < maxP))
		    {
		      if(uselog)
			score = TestPeriod(N, t, mag, testperiod, uselog,&h);
		      else
			score = TestPeriod_ASZ(N, t, mag, testperiod,&h);
		      if(score < ERROR_SCORE && (score * 0.0 == 0.0))
			{
			  if (score < bestscore)
			    {
			      ismultiple = 1;
			      abest = a;
			      bbest = b;
			      bestscore = score;
			    }
			}
		    }
		}
	  if(ismultiple)
	    {
	      perpeaks[j] = perpeaks[j] * abest / bbest;
	      aovpeaks[j] = bestscore;
	    }

	}

      mysort2(Npeaks,aovpeaks,perpeaks);


      for(j=0;j<Npeaks;j++)
	{
	  if(aovpeaks[j] < ERROR_SCORE)
	    {
	      if(uselog)
		aovpeaks[j] = (aveper - aovpeaks[j]) / stdper;
	      else
		{
		  aovSNR[j] = (aveper - aovpeaks[j]) / stdper;
		  aovpeaks[j] = -aovpeaks[j];
		  a_ = 0.5*((double) (h.HIST_SIZE - 1));
		  b_ = 0.5*((double) (N - h.HIST_SIZE));
		  aovFAP[j] = -log1minusbetai(a_, b_, ((2. * a_ * aovpeaks[j])/(2. * (b_ + a_*aovpeaks[j])))) + negln_m_eff;
		}
	    }
	  else
	    {
	      aovpeaks[j] = -ERROR_SCORE - 1.;
	      perpeaks[j] = 1.0;
	    }
	}
    }
  else
    {
      /* We need to whiten the light curve and re-compute the periodogram after finding each peak */
      /* Compute the periodogram */
      if(uselog)
	AOVPeriodogram(N, t, mag, Nperiod, periods, periodogram_whiten[0], uselog, &h);
      else
	AOVPeriodogram_ASZ(N, t, mag, Nperiod, periods, periodogram_whiten[0], &h);
      for(peakiter=0;peakiter < Npeaks; peakiter++)
	{
	  /* Find the peak period */
	  perpeaks[peakiter] = 1.0;
	  aovpeaks[peakiter] = ERROR_SCORE;
	  for(i=0; i < Nperiod; i++)
	    {
	      if(periodogram_whiten[peakiter][i] < ERROR_SCORE && periodogram_whiten[peakiter][i]*0.0 == 0.0)
		{
		  if(periodogram_whiten[peakiter][i] < aovpeaks[peakiter])
		    {
		      test = 1;
		      for(j=0;j<peakiter;j++)
			{
			  if(!isDifferentPeriods(periods[i],perpeaks[j],T))
			    {
			      test = 0;
			      break;
			    }
			}
		      if(test)
			{
			  perpeaks[peakiter] = periods[i];
			  aovpeaks[peakiter] = periodogram_whiten[peakiter][i];
			}
		    }
		}
	    }

	  /* Now perform the high-resolution period scan on the peak */

	  smallfreqstep = fine_tune/T;

	  j = peakiter;
	  freq = dmin((1./perpeaks[j]) + freqstep,(1./minP));
	  minfreq = dmax((1./perpeaks[j]) - freqstep,(1./maxP));
	  while (freq >= minfreq)
	    {
	      testperiod = 1./freq;
	      if(uselog)
		score = TestPeriod(N, t, mag,testperiod, uselog, &h);
	      else
		score = TestPeriod_ASZ(N, t, mag,testperiod, &h);
	      if(score < ERROR_SCORE && (score * 0.0 == 0.0))
		{
		  if(score < aovpeaks[j])
		    {
		      aovpeaks[j] = score;
		      perpeaks[j] = testperiod;
		    }
		}
	      freq -= smallfreqstep;
	    }
	  // double-check period multiples
	  bestscore = aovpeaks[j];
	  ismultiple = 0;
	  for(a=1 ; a <= MAX_DOUBLE_CHECK_MULTIPLE ; a++)
	    for(b=1 ; b <= MAX_DOUBLE_CHECK_MULTIPLE ; b++)
	      if(a != b)
		{
		  testperiod = perpeaks[j] * a / b;
		  if((testperiod > minP) && (testperiod < maxP))
		    {
		      if(uselog)
			score = TestPeriod(N, t, mag, testperiod, uselog, &h);
		      else
			score = TestPeriod_ASZ(N, t, mag, testperiod, &h);
		      if(score < ERROR_SCORE && (score * 0.0 == 0.0))
			{
			  if (score < bestscore)
			    {
			      ismultiple = 1;
			      abest = a;
			      bbest = b;
			      bestscore = score;
			    }
			}
		    }
		}
	  if(ismultiple)
	    {
	      perpeaks[j] = perpeaks[j] * abest / bbest;
	      aovpeaks[j] = bestscore;
	    }

	  if(aovpeaks[peakiter] < ERROR_SCORE)
	    {
	      /* Whiten the light curve at this period */
	      whitenlc_aov(N, t, mag, perpeaks[peakiter],&h);

	      /* Normalize the whitened light curve */
	      normalize(N, t, mag, sig, &ave, &stddev);

	      /* Compute the periodogram of the whitened light curve */
	      if(uselog)
		AOVPeriodogram(N, t, mag, Nperiod, periods, periodogram_whiten[peakiter + 1], uselog, &h);
	      else
		AOVPeriodogram_ASZ(N, t, mag, Nperiod, periods, periodogram_whiten[peakiter + 1], &h);

	      /* Get the mean and standard deviation of this periodogram */
	      Sum = 0.;
	      Sumsqr = 0.;
	      Ngood = 0;
	      for(i=0;i<Nperiod;i++)
		{
		  if(periodogram_whiten[peakiter + 1][i] < ERROR_SCORE && (periodogram_whiten[peakiter + 1][i]*0.0 == 0.0))
		    {
		      Sum += (long double) periodogram_whiten[peakiter + 1][i];
		      Sumsqr += (long double) (periodogram_whiten[peakiter + 1][i] * periodogram_whiten[peakiter + 1][i]);
		      Ngood++;
		    }
		}

	      Sum /= Ngood;
	      Sumsqr /= Ngood;

	      aveper_whiten[peakiter] = (double) Sum;
	      stdper_whiten[peakiter] = sqrt((double)(Sumsqr - (Sum*Sum)));

	      nclippedthis = 0;
	      do {
		nclippedlast = nclippedthis;
		nclippedthis = 0;

		Sum = 0.;
		Sumsqr = 0.;
		Ngood = 0;
		for(i=0;i<Nperiod;i++)
		  {
		    if(periodogram_whiten[peakiter + 1][i] < ERROR_SCORE && (periodogram_whiten[peakiter + 1][i]*0.0 == 0.0))
		      {
			if(periodogram_whiten[peakiter + 1][i] > aveper_whiten[peakiter] - clip*stdper_whiten[peakiter])
			  {
			    Sum += (long double) periodogram_whiten[peakiter+1][i];
			    Sumsqr += (long double) (periodogram_whiten[peakiter+1][i] * periodogram_whiten[peakiter+1][i]);
			    Ngood++;
			  }
			else
			  nclippedthis++;
		      }
		  }

		Sum /= Ngood;
		Sumsqr /= Ngood;

		aveper_whiten[peakiter] = (double) Sum;
		stdper_whiten[peakiter] = sqrt((double)(Sumsqr - (Sum*Sum)));
	      } while(clipiter && nclippedthis > nclippedlast);

	      aveaov_whiten[peakiter] = aveper_whiten[peakiter];
	      stddevaov_whiten[peakiter] = stdper_whiten[peakiter];
	      /* Get the peak value and SNR for the fixed period if we're doing that */
	      if(fixperiodSNR && !peakiter)
		{
		  AOV_getfixedperiodSNR(periods,periodogram_whiten[peakiter],Nperiod,aveper_whiten[peakiter],stdper_whiten[peakiter],fixperiodSNR_period,fixperiodSNR_value,fixperiodSNR_SNR);
		  if(uselog)
		    *fixperiodSNR_value = (aveper_whiten[peakiter] - *fixperiodSNR_value) / stdper_whiten[peakiter];
		  else
		    {
		      *fixperiodSNR_SNR = (aveper_whiten[peakiter] - *fixperiodSNR_value) / stdper_whiten[peakiter];
		      *fixperiodSNR_value = -(*fixperiodSNR_value);
		      a_ = 0.5*((double) (h.HIST_SIZE - 1));
		      b_ = 0.5*((double) (N - h.HIST_SIZE));
		      *fixperiodSNR_FAP = -log1minusbetai(a_, b_, ((2. * a_ * (*fixperiodSNR_value))/(2. * (b_ + a_*(*fixperiodSNR_value)))));
		    }
		}
	    }
	  else
	    {
	      for(;peakiter < Npeaks;peakiter++)
		{
		  aovpeaks[peakiter] = ERROR_SCORE;
		  perpeaks[peakiter] = 1.;
		  aveaov_whiten[peakiter] = 0.;
		  stddevaov_whiten[peakiter] = 0.;
		}
	    }
	}


      /* Write out the periodogram if asked to */
      if(outflag)
	{
	  if((outfile = fopen(outname,"w")) == NULL)
	    {
	      fprintf(stderr,"Cannot Write periodogram to %s\n",outname);
	      exit(3);
	    }
	  if(ascii)
	    {
	      fprintf(outfile,"#Period");
	      for(j=0;j<Npeaks;j++) {
		if(uselog) {
		  fprintf(outfile," AOV_LOGSNR_WhitenCycle_%d", j);
		}
		else
		  fprintf(outfile," AOV_WhitenCycle_%d", j);
	      }
	      fprintf(outfile,"\n");
	      for(i=0;i<Nperiod;i++)
		{
		  if(periodogram_whiten[0][i] < ERROR_SCORE && 0.0*periodogram_whiten[0][i] == 0.0)
		    {
		      fprintf(outfile,"%f",periods[i]);
		      for(j=0;j<Npeaks;j++)
			{
			  if(uselog)
			    fprintf(outfile," %f",((aveper_whiten[j] - periodogram_whiten[j][i]) / stdper_whiten[j]));
			  else
			    fprintf(outfile," %f",-periodogram_whiten[j][i]);
			}
		      fprintf(outfile,"\n");
		    }
		}
	    }
	  else
	    {
	      fwrite(&Nperiod,4,1,outfile);
	      fwrite(periods,8,Nperiod,outfile);
	      for(j=0;j<Npeaks;j++)
		{
		  fwrite(&(aveper_whiten[j]),8,1,outfile);
		  fwrite(&(stdper_whiten[j]),8,1,outfile);
		  fwrite(periodogram_whiten[j],8,Nperiod,outfile);
		}

	    }
	  fclose(outfile);
	}

      for(j=0;j<Npeaks;j++)
	{
	  if(aovpeaks[j] < ERROR_SCORE)
	    {
	      if(uselog)
		aovpeaks[j] = (aveper_whiten[j] - aovpeaks[j]) / stdper_whiten[j];
	      else
		{
		  aovSNR[j] = (aveper_whiten[j] - aovpeaks[j]) / stdper_whiten[j];
		  aovpeaks[j] = -aovpeaks[j];
		  a_ = 0.5*((double) h.HIST_SIZE - 1);
		  b_ = 0.5*((double) (N - h.HIST_SIZE));
		  aovFAP[j] = -log1minusbetai(a_, b_, ((2. * a_ * aovpeaks[j])/(2. * (b_ + a_*aovpeaks[j])))) + negln_m_eff;
		}
	    }
	  else
	    {
	      aovpeaks[j] = -ERROR_SCORE - 1.;
	      perpeaks[j] = 1.0;
	    }
	}

    }

  free(t);
  free(mag);
  free(sig);
  if(sizeHISTvector > 0) {
    free(h.histN);
    free(h.histA);
    free(h.histB);
    free(h.histC);
    free(h.histD);
    free(h.histE);
  }

  if(periods != NULL)
    free(periods);

  if(periodogram != NULL)
    free(periodogram);

  if(periodogram_whiten != NULL) {
    for(i=0; i < Npeaks + 1; i++)
      free(periodogram_whiten[i]);
    free(periodogram_whiten);
  }

  if(aveper_whiten != NULL)
    free(aveper_whiten);

  if(stdper_whiten != NULL)
    free(stdper_whiten);

}

void RunAOVCommand(ProgramData *p, Command *c, _Aov *Aov, int lcnum, int lc_name_num, int thisindex)
{
  int i1, i2;
  char outname[MAXLEN];
  double d1;
  double *d1ptr, *d2ptr, *d3ptr;
  /* Calculate the AoV with phase binning */
  if(Aov->operiodogram)
    {
      i1 = 0;
      i2 = 0;
      while(p->lcnames[lc_name_num][i1] != '\0')
	{
	  if(p->lcnames[lc_name_num][i1] == '/')
	    i2 = i1 + 1;
	  i1++;
	}
      sprintf(outname,"%s/%s%s",Aov->outdir,&p->lcnames[lc_name_num][i2],Aov->suffix);
    }

  if(Aov->minp_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    Aov->minp_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Aov->minp_expr);
  }
  else if(Aov->minp_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    Aov->minp_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Aov->minp_var);
  }
  else {
    Aov->minp_vals[lcnum] = Aov->minp;
  }

  if(Aov->maxp_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    Aov->maxp_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Aov->maxp_expr);
  }
  else if(Aov->maxp_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    Aov->maxp_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Aov->maxp_var);
  }
  else {
    Aov->maxp_vals[lcnum] = Aov->maxp;
  }

  if(Aov->subsample_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    Aov->subsample_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Aov->subsample_expr);
  }
  else if(Aov->subsample_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    Aov->subsample_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Aov->subsample_var);
  }
  else {
    Aov->subsample_vals[lcnum] = Aov->subsample;
  }

  if(Aov->finetune_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    Aov->finetune_vals[lcnum] = EvaluateExpression(lc_name_num, lcnum, 0, Aov->finetune_expr);
  }
  else if(Aov->finetune_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    Aov->finetune_vals[lcnum] = EvaluateVariable_Double(lc_name_num, lcnum, 0, Aov->finetune_var);
  }
  else {
    Aov->finetune_vals[lcnum] = Aov->finetune;
  }

  if(Aov->Nbin_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    Aov->Nbin_vals[lcnum] = ceil(EvaluateExpression(lc_name_num, lcnum, 0, Aov->Nbin_expr));
  }
  else if(Aov->Nbin_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
    Aov->Nbin_vals[lcnum] = ceil(EvaluateVariable_Double(lc_name_num, lcnum, 0, Aov->Nbin_var));
  }
  else {
    Aov->Nbin_vals[lcnum] = Aov->Nbin;
  }


  if(Aov->fixperiodSNR)
    {
      if(Aov->fixperiodSNR_pertype == PERTYPE_AOV)
	{
	  i1=Aov->fixperiodSNR_lastaovindex;
	  if(c[i1-thisindex].cnum == CNUM_AOV)
	    Aov->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].Aov->peakperiods[lcnum][0];
	  else if(c[i1-thisindex].cnum == CNUM_HARMAOV)
	    Aov->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].AovHarm->peakperiods[lcnum][0];
	  
	}
      else if(Aov->fixperiodSNR_pertype == PERTYPE_LS)
	{
	  i1 = Aov->fixperiodSNR_lastaovindex;
	  Aov->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].Aov->peakperiods[lcnum][0];
	}
      else if(Aov->fixperiodSNR_pertype == PERTYPE_INJECTHARM)
	{
	  i1 = Aov->fixperiodSNR_lastaovindex;
	  Aov->fixperiodSNR_periods[lcnum][0] = c[i1-thisindex].Injectharm->periodinject[lcnum];
	}
      else if(Aov->fixperiodSNR_pertype == PERTYPE_FIX)
	{
	  Aov->fixperiodSNR_periods[lcnum][0] = Aov->fixperiodSNR_fixedperiod;
	}
      else if(Aov->fixperiodSNR_pertype == PERTYPE_FIXCOLUMN)
	{
	  getoutcolumnvalue(Aov->fixperiodSNR_linkedcolumn, lcnum, lc_name_num, VARTOOLS_TYPE_DOUBLE, &(Aov->fixperiodSNR_periods[lcnum][0]));
	}
      d1 = Aov->fixperiodSNR_periods[lcnum][0];
      d1ptr = &(Aov->fixperiodSNR_peakvalues[lcnum]);
      d2ptr = &(Aov->fixperiodSNR_peakSNR[lcnum]);
      d3ptr = &(Aov->fixperiodSNR_peakFAP[lcnum]);
    }
  else
    {
      d1 = 1.;
      d1ptr = NULL;
      d2ptr = NULL;
      d3ptr = NULL;
    }
  findPeaks_aov(p->t[lcnum], p->mag[lcnum], p->sig[lcnum], p->NJD[lcnum], Aov->peakperiods[lcnum], Aov->peakvalues[lcnum], Aov->peakSNR[lcnum], Aov->peakFAP[lcnum], Aov->Npeaks, Aov->minp_vals[lcnum], Aov->maxp_vals[lcnum], Aov->subsample_vals[lcnum], Aov->finetune_vals[lcnum], Aov->operiodogram, outname, &Aov->aveaov[lcnum], &Aov->rmsaov[lcnum], Aov->aveaov_whiten[lcnum], Aov->rmsaov_whiten[lcnum], p->ascii, Aov->Nbin_vals[lcnum], Aov->whiten, Aov->uselog, Aov->clip, Aov->clipiter, Aov->fixperiodSNR, d1, d1ptr, d2ptr, d3ptr, lcnum, lc_name_num, Aov->usemask, Aov->maskvar);

}

