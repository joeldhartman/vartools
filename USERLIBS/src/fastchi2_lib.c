/* This file fastchi2_lib.c is a copy of the file fastchi2.c taken from
   David Palmer's FastChi2 program version 1.03. Here this file is used
   as a library to implement the fastchi2 algorithm in the 
   VARTOOLS program. The original reference for the fastchi2 algorithm is
     Palmer 2009, ApJ, 695, 496.
   Below is the copyright notice included in the original version of
   fastchi2.h.

   The following modifications have been made relative to the original 
   version from Palmer:

   The functions findBestGridPoint and findBestFrequency have been modified 
       to find multiple peaks.

*/

/* fastchi2.c */
/* routines to implement the fast chi^2 algorithm */
/*
 Copyright 2007.  Los Alamos National Security, LLC. This material was 
 produced under U.S. Government contract DE-AC52-06NA25396 for 
 Los Alamos National Laboratory (LANL), which is operated by 
 Los Alamos National Security, LLC for the U.S. Department of Energy. 
 The U.S. Government has rights to use, reproduce, and distribute this software.  
 NEITHER THE GOVERNMENT NOR LOS ALAMOS NATIONAL SECURITY, LLC MAKES ANY 
 WARRANTY, EXPRESS OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF 
 THIS SOFTWARE.  If software is modified to produce derivative works, 
 such modified software should be clearly marked, so as not to confuse 
 it with the version available from LANL.
 
 
 Additionally, this program is free software; you can redistribute it 
 and/or modify it under the terms of the GNU General Public License 
 as published by the Free Software Foundation; version 2.0 of the License. 
 Accordingly, this program is distributed in the hope that it will be useful, 
 but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY 
 or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
 for more details
*/

#include <stdlib.h>
#include <memory.h>
#include <assert.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_fft_halfcomplex_float.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_cdf.h>

#include "fastchi2_lib.h"

#include <math.h>

char* version_fastchi2 = "1.00";

#define VERBOSE 0
#define PRINTMATRIX 0

/* The line below was added to this modified version of this file by
   J. Hartman */
int VARTOOLS_isDifferentPeriods (double period1, double period2, double TimeSpan);


/*********************** Prototypes and inlines **************/

int my_vector_fprintf (FILE * stream, const gsl_vector * v,
                  const char *format);
int my_matrix_fprintf (FILE * stream, const gsl_matrix * v,
                  const char *format);

inline double dmin(double d1, double d2) { return d1 < d2 ? d1 : d2;}
inline double dmax(double d1, double d2) { return d1 > d2 ? d1 : d2;}

            /* Produce a vector of chisquared values as a function of frequency */
void    fastChi(const double *sampletimes, const float *vals, const float *stderrs, 
                long ndata, long nharmonics, double deltat, 
                long nrealpoints, long nchivalues, double tstart,
                float *dovstorage, float *oovstorage,
                float *chireductionlist)
{
    int d;
    int i;
    int h;
    float invdeltat = 1/deltat;

    gsl_matrix *alpha = gsl_matrix_alloc(2 * nharmonics + 1,2 * nharmonics + 1);
    gsl_vector *beta =  gsl_vector_alloc(2 * nharmonics + 1);
    gsl_vector *x    =  gsl_vector_alloc(2 * nharmonics + 1);

    float *cosdata = malloc((nharmonics + 1) * sizeof(float));
    float *sindata = malloc((nharmonics + 1) * sizeof(float));
    float *cosvar = malloc((2*nharmonics + 1) * sizeof(float));
    float *sinvar = malloc((2*nharmonics + 1) * sizeof(float));

        /* if the pointers point somewhere, assume that they point to the right place */
    assert(sampletimes && vals && stderrs && dovstorage && oovstorage && chireductionlist);
    assert(nchivalues * 2 * nharmonics <= nrealpoints); /* FIXME: change 2 to 1 when */
                                                        /* using aliasing */
    assert(cosdata && sindata && cosvar && sinvar);
    
    memset(dovstorage, 0, nrealpoints * sizeof(float));
    memset(oovstorage, 0, nrealpoints * sizeof(float));
    
    for (d = 0 ; d < ndata ; d++)
    {
        long tindex = (sampletimes[d] - tstart) * invdeltat;
        float inv_variance = 1/(stderrs[d] * stderrs[d]); 
        assert(tindex >= 0 && tindex < nrealpoints);
                /* adding rather than setting averages measurments that */
                /* fall in same timebin.  It does, however, reduce the d.o.f. */
                /* and subsequently improperly reduces the chisquared.  */
                /* This is a minor point. */
            /* dov: _d_ata _o_ver _v_ariance  oov: _o_ne  _o_ver _v_ariance */
        dovstorage[tindex] += vals[d] * inv_variance;
        oovstorage[tindex] += inv_variance;
    }
    gsl_fft_real_float_radix2_transform(dovstorage, 1, nrealpoints);
    gsl_fft_real_float_radix2_transform(oovstorage, 1, nrealpoints);

    cosdata[0] = dovstorage[0];     /* The constant (0f) terms are the same for all f */
    cosvar[0] = oovstorage[0];

    for (i = 1 ; i < nchivalues ; i++)
    {
        for (h = 1 ; h <= nharmonics ; h++)
        {
            cosdata[h] = dovstorage[i*h];
            sindata[h] = dovstorage[nrealpoints - (i * h)] * FFTSIGNCONVENTION;
        }
        for (h = 1 ; h <= 2*nharmonics ; h++)
        {                               /* put alias wrapping in here */
            cosvar[h] = oovstorage[i*h];
            sinvar[h] = oovstorage[nrealpoints - (i * h)] * FFTSIGNCONVENTION;
        }
        chireductionlist[i] 
                = calcChiReduction(nharmonics,
                                cosdata, sindata, cosvar, sinvar,
                                alpha, beta, x);
    }

        
    gsl_matrix_free(alpha);
    gsl_vector_free(beta);
    gsl_vector_free(x);

    free(cosdata);
    free(sindata);
    free(cosvar);
    free(sinvar);

}


                
            /* Calculate the chi-squared reduction at a single frequency, returning the harmonic */
            /* coefficients as well */
double  singleChi(const double *sampletimes, const float *vals, const float *stderrs, 
                long ndata, long nharmonics, double tstart,
                double frequency, float *harmoniccoeffs)
{
#define MAXHARM 17
    int i, j, d;
    int result;
    double  alpha[2*MAXHARM+1][2*MAXHARM+1];
    double  beta[2*MAXHARM+1];
    double chisqred;

    gsl_matrix *alpha_gsl = gsl_matrix_alloc(2 * nharmonics + 1,2 * nharmonics + 1);
    gsl_vector *beta_gsl =  gsl_vector_alloc(2 * nharmonics + 1);
    gsl_vector *x_gsl    =  gsl_vector_alloc(2 * nharmonics + 1);

    assert(nharmonics <= MAXHARM);

    for (i = 0 ; i < 2*nharmonics+1 ; i++)
    {
        beta[i] = 0;
        for (j = 0 ; j < 2*nharmonics+1 ; j++)
            alpha[i][j] = 0;
    }

    for (d = 0 ; d < ndata ; d++)
    {
        for (i = 0 ; i < 2*nharmonics+1 ; i++)
        {
            double invvar = 1/(stderrs[d] * stderrs[d]);
            beta[i] += vals[d] 
                    * invvar 
                    * Mfunc(i, (sampletimes[d] - tstart) * frequency * 2 * M_PI);
            for (j = 0 ; j < 2*nharmonics+1 ; j++)
                alpha[i][j] += invvar 
                            * Mfunc(i, (sampletimes[d] - tstart) * frequency * 2 * M_PI)
                            * Mfunc(j, (sampletimes[d] - tstart) * frequency * 2 * M_PI);
        }
    }    

    for (i = 0 ; i < 2*nharmonics+1 ; i++)
    {
        gsl_vector_set(beta_gsl, i, beta[i]);
        for (j = 0 ; j < 2*nharmonics+1 ; j++)
            gsl_matrix_set(alpha_gsl, i, j, alpha[i][j]);
    }

    result = gsl_linalg_HH_solve(alpha_gsl, beta_gsl, x_gsl);
            /* x = alpha^-1 beta */
    assert(result == 0);
    gsl_blas_ddot(x_gsl, beta_gsl, &chisqred);

    if (harmoniccoeffs)
    {
        for (i = 0 ; i < 2*nharmonics+1 ; i++)
        {
            harmoniccoeffs[i] = gsl_vector_get(x_gsl, i);
        }
    }
    
    gsl_matrix_free(alpha_gsl);
    gsl_vector_free(beta_gsl);
    gsl_vector_free(x_gsl);

    return chisqred;
}


/* Optimize the chi-squared near a single frequency +/- fhwidth,   */
/* returns the chisquare improvement, the harmonic coeffs, and the  */
/* optimal frequency */
double	localOptimizeChi(const double *sampletimes, const float *vals, const float *stderrs, 
                         long ndata, long nharmonics, double tstart,
                         double frequency, double fhwidth, 
                         float *harmoniccoeffs, double *bestf)
{
    int i;
    double chisqvalues[3];  /* Array for parabolic fit */
    double firstderiv, secondderiv, peakloc;
    double bestchireduction;
    
    for (i = 0 ; i < 3 ; i++)
    {
        chisqvalues[i] = singleChi(sampletimes, vals, stderrs, 
                                   ndata, nharmonics, tstart,
                                   frequency + fhwidth*(i-1.0), harmoniccoeffs);
    }
    
    firstderiv = (chisqvalues[2] - chisqvalues[0])/2;
    secondderiv = (chisqvalues[2] + chisqvalues[0]) -  2 * chisqvalues[1];
    peakloc = (-firstderiv)/secondderiv;
    if (peakloc > 1.0 || peakloc < -1.0)
    {
        /* This can happen because, although the gridded fast-chi-squared is */
        /* a local minimum compared to that +/- one bin, the ungridded values */
        /* might cause enough of a change that this need not be the case. */
        /* fprintf(stderr, "# Frequency peak found outside of allowed range\n"); */
    }
    *bestf = frequency + peakloc * fhwidth;
    
    bestchireduction = singleChi(sampletimes, vals, stderrs, 
                                 ndata, nharmonics, tstart,
                                 *bestf, harmoniccoeffs);
    if (bestchireduction < chisqvalues[1])
    {
        /* This can also happen if the chi-squared surface is iinsufficiently */
        /* parabolic in this vicinity.  Just use the center frequency as Good Enough */
        /* fprintf(stderr, "# Parabolic fit does not improve chi-squared\n"); */
        bestchireduction = chisqvalues[1];
        *bestf = frequency;
    }
    return bestchireduction;
}



/* Calculate the height of the parabolic maximum near p0, if any, */
/* and zero if not peak.  This is not a general routine: */
/* it is only useful if the values are all positive and you don't */
/* care what the answer is except near the maximum */
inline float parabolic_array_peak(float* p0)
{
    float vm1 = *(p0-1);    /* Value at [i-1] */
    float v0 = *p0;
    float vp1 = *(p0+1);    /* Trust compiler to not do this access if not necessary (p0 < vm1) */
    if (v0 < vm1 || v0 <= vp1) 
    {
        return 0;
    }
    else
    {   /* v_peak = v0 - v'^2/2v'' */
        float firstderiv = (vp1 - vm1)/2;
        float secondderiv = (vp1 + vm1) -  2 * v0;
        return v0 - 0.5 * firstderiv * firstderiv/secondderiv;
    }
}


/* Return the index of the best chi-squared in the array, where 'best' may */
/* be adjusted by the number of trials. */
void    findBestGridPoint(float *chireductionvalues, 
                         int findexmin, int findexmax,
                         int nharmonics,
                         int fAdjusted,
			 int *findexout,
			 int Ntofind)
{
    int findex;
    float chireductiontest = 0.0;
    float chireduction_parab;
    int i, j, havefound = 0;

    double *pbest;
    double ndof = 2*nharmonics;
    
    pbest = (double *) malloc(Ntofind * sizeof(double));

    for(i=0; i < Ntofind; i++) {
      pbest[i] = findexmax;
      findexout[i] = 0;
    }

    /* Do a fast check to find the chi-squared or adjusted-probability minimum */
    /* among the gridded points, without spending the CPU time to optimize */
    /* each potential minimum-so-far.  */
    for (findex = findexmin ; findex < findexmax ; findex++)
    {
        /* Check if it is the best chisquared value so far.  */
        /* Also check to make sure that it is a local best */
      chireduction_parab = parabolic_array_peak(&chireductionvalues[findex]);
      if(chireduction_parab > chireductiontest) {
	double p = Pmetric(chireduction_parab, ndof, findex, fAdjusted);
	if(p >= pbest[Ntofind-1])
	  continue;
	for(i=0; i < Ntofind; i++) {
	  if(p < pbest[i]) {
	    for(j=Ntofind-1; j > i; j--) {
	      pbest[j] = pbest[j-1];
	      findexout[j] = findexout[j-1];
	    }
	    pbest[i] = p;
	    findexout[i] = findex;
	    havefound++;
	    break;
	  }
	}
	if(havefound >= Ntofind)
	  chireductiontest = Pmetric_inverse(pbest[Ntofind-1], ndof, findex, fAdjusted);
      }
    }
    free(pbest);
}

void findBestFrequency(int fAdjusted,
                         const double *sampletimes, const float *vals, const float *stderrs, 
                         long ndata, int nharmonics,
                         double freqmin, double deltaf, double t0, float chimargin,
                         float *chireductionvalues, int nchivalues, 
                         float *pbest_peakup, double *bestfreq,
			 int Ntofind)
{
    float harmweights[2*MAXORDER + 1];
    double ndof = 2*nharmonics;
    
    int findex;
    double freqlocalbest;
    
    double  pbest;
    int *findexbest_vals;
    
    double chireduction;
    float chireduction_parab;
    float  chireductiontest;   /* The gridded peak must be better than this to have a chance to */
                                /* be near the best fit point */
    
    double best_peakup = 0;      /* Maximum improvement of chisquared from searching off the gridpoints */ 
    float usemargin;            /* If the improvement is better FOR THIS STAR at some frequency, use that margin for the rest of the star */
    
    /* Do not look at the endpoint values since you can't tell whether they are local optima */
    double timerange = sampletimes[ndata - 1] - sampletimes[0];
    int findexmin = dmax(1, dmax(freqmin, 1/timerange)/deltaf);
    int findexmax = nchivalues-2;
    int findexmin_finesearch;
    int findexmax_finesearch;
    int findexbest;
    int i;

    findexbest_vals = (int *) malloc(Ntofind * sizeof(int));

    findBestGridPoint(chireductionvalues, 
                                   findexmin, findexmax,
                                   nharmonics,
                                   fAdjusted,
				   findexbest_vals,
				   Ntofind);
    

    for(i=0; i < Ntofind; i++) {
      /* Now set up the initial test values */
      findexbest = findexbest_vals[i];

      bestfreq[i] = 0.0;
      chireduction = localOptimizeChi(sampletimes, vals, stderrs, 
				      ndata, nharmonics, t0,
				      findexbest * deltaf, deltaf,
				      harmweights, &(bestfreq[i]));
    
      /* The margin you should use for gridded vs. ungridded is at least twice as */
      /* large as the one you found for the best gridded point.   (The twice is   */
      /* so that the margin is relative to the grid point rather than the         */
      /* optimized point.)                                                        */
      best_peakup = chireduction - parabolic_array_peak(&chireductionvalues[findexbest]);
      usemargin = dmax(chimargin, 2*best_peakup);
      pbest = Pmetric(chireduction, ndof, findexbest, fAdjusted);
    
    /* The chi reduction will not be a best point if it is not at least this good */
      chireductiontest = Pmetric_inverse(pbest, ndof, findexmin, fAdjusted) - usemargin;
    
      /* to be better than the best of the gridded values around it. */
    
      findexmin_finesearch = (findexbest - 4 >= findexmin ? findexbest - 4 : findexmin);
      findexmax_finesearch = (findexbest + 5 <= findexmax ? findexbest + 5 : findexmax);
      for (findex = findexmin_finesearch ; findex < findexmax_finesearch ; findex++)
	{
	  /* Check if it is a possible best chisquared.  */
	  /* Also check to make sure that it is a local best */
	  if ((chireduction_parab = parabolic_array_peak(&chireductionvalues[findex]))
	      > chireductiontest)
	    {
	      double pnew;
	      double peakup;
            
	      /* More carefully calculate the chisquared reduction */ 
	      chireduction = localOptimizeChi(sampletimes, vals, stderrs, 
					      ndata, nharmonics, t0,
					      findex * deltaf, deltaf,
					      harmweights, &freqlocalbest);
            
	      peakup = chireduction - chireduction_parab;
	      if (peakup > best_peakup)
		{   /* keep track of how much the ungridded data improves chi squared */
		  best_peakup = peakup;
		  usemargin = dmax(usemargin, best_peakup);
		}
	      
	      pnew = Pmetric(chireduction, ndof, findexbest, fAdjusted);
	      if (pnew < pbest)
		{
		  pbest = pnew;
		  findexbest = findex;
		  bestfreq[i] = freqlocalbest;
		  usemargin = dmax(usemargin, 2*peakup);
		}
	      /* This test value can become more stringent even if the Pmetric */
	      /* isn't improved, due to the increased number of trials */
	      /* if fAdjusted is true, or less stringent if usemargin increases. */
	      chireductiontest = Pmetric_inverse(pbest, ndof, findex, fAdjusted) - usemargin;
	      
	      /* Strictly speaking, if the usemargin increases, the entire array */
	      /* should be rescanned.  For simplicity, this is not done in this */
	      /* implementation. Instead, warnings are sent to stderr if -m     */
	      /* is too small.  */
	    }
	}
      pbest_peakup[i] = best_peakup;
    }
    free(findexbest_vals);
    return;
}

            /* given the set of FFT values for data/variance and 1/variance     */
            /* calculate the best chisquared reduction.  The gsl_{matrix,vector}s must be */
            /* pre-allocated for efficiency, and are used to solve: */
            /*      alpha beta = x      */
            /* cos*[0] is the mean value (f=0), sin*[0] is ignored, */
            /* *data[0..nharmonics], *var[0..2*nharmonics]  */
double  calcChiReduction(int nharmonics, 
                const float *cosdata, const float *sindata,
                const float *cosvar, const float *sinvar,
                gsl_matrix *alpha, gsl_vector *beta, gsl_vector *x)
{
    double chisqred;


    calcAlphaBeta(nharmonics, 
                cosdata, sindata,
                cosvar, sinvar,
                alpha, beta);

        /* alpha and beta are now fully set.  Time to solve */
        /* Stick in your thumb and pull out a plum: solve by Householder */
    gsl_linalg_HH_solve(alpha, beta, x);
            /* x = alpha^-1 beta */
    
    gsl_blas_ddot(x, beta, &chisqred);
    return chisqred;
}                


void  calcAlphaBeta(int nharmonics, 
                const float *cosdata, const float *sindata,
                const float *cosvar, const float *sinvar,
                gsl_matrix *alpha, gsl_vector *beta)
{
    int i, j;
    int f1, f2;
    
    gsl_matrix_set(alpha, 0, 0, cosvar[0]);
    gsl_vector_set(beta, 0, cosdata[0]);
    
                /* calculate the lower triangular matrix of alpha */
    for (f1 = 1 ; f1 <= nharmonics ; f1++)
    {    
                /* Beta */
        gsl_vector_set(beta, 2*f1-1, sindata[f1]);
        gsl_vector_set(beta, 2*f1,   cosdata[f1]);
                
                /* left column of alpha */
        gsl_matrix_set(alpha, 2*f1-1, 0,      sinvar[f1]);
        gsl_matrix_set(alpha, 2*f1,   0,      cosvar[f1]);
        
                /* diagonal block of alpha */
        gsl_matrix_set(alpha, 2*f1-1, 2*f1-1, (cosvar[0] - cosvar[2*f1]) * 0.5);
        gsl_matrix_set(alpha, 2*f1,   2*f1,   (cosvar[0] + cosvar[2*f1]) * 0.5);
        gsl_matrix_set(alpha, 2*f1,   2*f1-1,   sinvar[2*f1] * 0.5);
                            /* [2*f1-1,2*f1] set below when symmetrizing */
                             
                /* lower blocks of alpha */
        for (f2 = 1 ; f2 < f1 ; f2++)
        {
            gsl_matrix_set(alpha, 2*f1-1, 2*f2-1, (cosvar[f1-f2] - cosvar[f1+f2])*0.5);
            gsl_matrix_set(alpha, 2*f1-1, 2*f2,   (sinvar[f1+f2] + sinvar[f1-f2])*0.5);
            gsl_matrix_set(alpha, 2*f1,   2*f2-1, (sinvar[f1+f2] - sinvar[f1-f2])*0.5);
            gsl_matrix_set(alpha, 2*f1,   2*f2,   (cosvar[f1-f2] + cosvar[f1+f2])*0.5);            
        } 
 
    }

   
            /* fill in the rest of alpha by symmetry */
            
    for (i = 0 ; i < 2 * nharmonics + 1 ; i++)
        for (j = 0 ; j < i ; j++)
            gsl_matrix_set(alpha, j, i, gsl_matrix_get(alpha, i, j));
}


double  Mfunc(int m_index, double phi)
{           /* cos(0 phi), sin(1 phi), cos(1 phi), sin(2 phi), cos(2 phi)... */
    if (m_index % 2 == 0)
    {
        return cos((phi - (2 * M_PI * (long)(phi/(2 * M_PI)))) * (m_index/2));
    }
    else
    {
        return sin((phi - (2 * M_PI * (long)(phi/(2 * M_PI)))) * ((m_index + 1)/2));
    }
}

int my_vector_fprintf (FILE * stream, const gsl_vector * v,
                  const char *format)
{
    unsigned int i;
    for (i = 0 ; i < v->size ; i++)
        fprintf(stream, format, v->data[i * v->stride]);
    
    fprintf(stderr, "\n");
    
    return 0;
}

int my_matrix_fprintf (FILE * stream, const gsl_matrix * v,
                  const char *format)
{
    unsigned int i, j;
    for (i = 0 ; i < v->size1 ; i++)
    {
        for (j = 0 ; j < v->size2 ; j++)
            fprintf(stream, format, v->data[i*v->tda + j]);
        fprintf(stream, "\n");
    }
    return 0;
}

/* Return a value which is monotonic in the probability of achieving */
/* the given chi-squared reduction by chance, adjusted for the number */
/* of trials if desired.  This allows relative comparison of two chisquared */
/* probabilities (which must have the same ndof and fAdjusted) without having */
/* to worry about e.g. Prob(chisquared) being zero due to underflow. */
#define MAX_CHISQ_USEFUL 1000
/* if the chisquared reduction is greater than 1000, then the probabilities */
/* are essentially meaningless, and the calculation would underflow anyway. */

double Pmetric(double chisq_reduction, double ndof, int ntrials, int fAdjusted)
{
    if (fAdjusted && chisq_reduction < MAX_CHISQ_USEFUL)
    {
        double p = gsl_cdf_chisq_Q(chisq_reduction, ndof) * ntrials;
        if (p > 0)
        {
            return p;
        }
        else
        {   /* underflow (chisquared of order 1000 or more) */
            return -chisq_reduction;
        }
    }
    else
    {
        return -chisq_reduction;
    }
}

/* return the chiqsuared reduction required for this pmetric value */
double Pmetric_inverse(double p, double ndof, int ntrials, int fAdjusted)
{
    if (p <= 0)
    {
        return -p;
    }
    else
    {
        assert(fAdjusted);
        if (p >= ntrials)
        {
            return 0;
        }
        else
        {
            return gsl_cdf_chisq_Qinv(p/ntrials, ndof);
        }
    }
}
