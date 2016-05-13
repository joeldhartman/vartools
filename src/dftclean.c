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

#define TWOPI 6.28318530717958647692528676656

void DFT(double *x, double *y, int N, double df, double Nf, double *f_r, double *f_i)
/* This computes the slow DFT of a signal with uneven time sampling. The number of points in x is N, the number of output frequencies is 2*Nf + 1 with a step-size of df. The real part of the DFT is stored in f_r, the imaginary part in f_i, Note the routine will subract the average x and y values.*/
{
  int i, j, Ntot;
  double avx, avy, cosval, sinval, fval, val;

  /* First get the average x and y values */
  for(i=0, avx = 0., avy = 0.;i<N;i++)
    {
      avx += x[i];
      avy += y[i];
    }
  avx = avx/N;
  avy = avy/N;

  Ntot = 2*Nf + 1;
  for(i=0,fval = -Nf*df;i<Ntot;i++, fval += df)
    {
      for(j=0,f_r[i] = 0., f_i[i] = 0.;j<N;j++)
	{
	  val = -TWOPI * fval * (x[j] - avx);
	  cosval = cos(val);
	  sinval = sin(val);
	  f_r[i] += (y[j] - avy)*cosval;
	  f_i[i] += (y[j] - avy)*sinval;
	}
      f_r[i] /= (double) N;
      f_i[i] /= (double) N;
    }
}

void DFT_Wfunc(double *x, int N, double df, double Nf, double *W_r, double *W_i)
/* Same as above, this time for computing the window spectrum (no y) */
{
  int i, j, Ntot;
  double avx, cosval, sinval, fval, val;

  /* First get the average x and y values */
  for(i=0, avx = 0.;i<N;i++)
    {
      avx += x[i];
    }
  avx = avx/N;

  Ntot = 2*Nf + 1;
  for(i=0,fval = -Nf*df;i<Ntot;i++, fval += df)
    {
      for(j=0,W_r[i] = 0., W_i[i] = 0.;j<N;j++)
	{
	  val = -TWOPI * fval * (x[j] - avx);
	  cosval = cos(val);
	  sinval = sin(val);
	  W_r[i] += cosval;
	  W_i[i] += sinval;
	}
      W_r[i] /= (double) N;
      W_i[i] /= (double) N;
    }
}

void FDFT(double *x, double *y, int N, double df, double Nf, double *f_r, double *f_i)
/* This computes the fast DFT of a signal with uneven time sampling using the algorithm from Kurtz D.W. 1985, MNRAS, 213, 773. The number of points in x is N, the number of output frequencies is 2*Nf + 1 with a step-size of df. The real part of the DFT is stored in f_r, the imaginary part in f_i, Note the routine will subract the average x and y values.*/
{
  int i, j, Ntot;
  double avx, avy, c, s, fval, val, yval;
  double xf0, xdf, sdf, cdf, s0, c0, xt, f0;

  /* First get the average x and y values */
  for(i=0, avx = 0., avy = 0.;i<N;i++)
    {
      avx += x[i];
      avy += y[i];
    }
  avx = avx/N;
  avy = avy/N;

  Ntot = 2*Nf + 1;
  for(i=0;i<Ntot;i++)
    {
      f_r[i] = 0.;
      f_i[i] = 0.;
    }

  for(j=0;j<N;j++)
    {
      f0 = (-Nf - 1)*df;
      xt = TWOPI * (x[j] - avx);
      xf0 = xt*f0;
      xdf = xt*df;
      s0 = sin(xf0);
      c0 = cos(xf0);
      sdf = sin(xdf);
      cdf = cos(xdf);
      yval = y[j] - avy;
      for(i=0;i<Ntot;i++)
	{
	  c = c0*cdf - s0*sdf;
	  s = s0*cdf + c0*sdf;
	  f_r[i] += yval*c;
	  f_i[i] += yval*s;
	  c0 = c;
	  s0 = s;
	}
    }
  for(i=0;i<Ntot;i++)
    {
      f_r[i] /= (double) N;
      f_i[i] /= (double) N;
    }
}

void FDFT_Wfunc(double *x, int N, double df, double Nf, double *W_r, double *W_i)
/* Same as above, this time for computing the window spectrum (no y) */
{
  int i, j, Ntot;
  double avx, avy, c, s, fval, val, yval;
  double xf0, xdf, sdf, cdf, s0, c0, xt, f0;

  /* First get the average x and y values */
  for(i=0, avx = 0.;i<N;i++)
    {
      avx += x[i];
    }
  avx = avx/N;

  Ntot = 2*Nf + 1;
  for(i=0;i<Ntot;i++)
    {
      W_r[i] = 0.;
      W_i[i] = 0.;
    }

  for(j=0;j<N;j++)
    {
      f0 = (-Nf - 1)*df;
      xt = TWOPI * (x[j] - avx);
      xf0 = xt*f0;
      xdf = xt*df;
      s0 = sin(xf0);
      c0 = cos(xf0);
      sdf = sin(xdf);
      cdf = cos(xdf);
      for(i=0;i<Ntot;i++)
	{
	  c = c0*cdf - s0*sdf;
	  s = s0*cdf + c0*sdf;
	  W_r[i] += c;
	  W_i[i] += s;
	  c0 = c;
	  s0 = s;
	}
    }
  for(i=0;i<Ntot;i++)
    {
      W_r[i] /= (double) N;
      W_i[i] /= (double) N;
    }
}

#define DFTTINY 1.0e-12

void getamplitude(double *f_r, double *f_i, double *W_r, double *W_i, int Nf, int peakfindx, double *a_r, double *a_i)
/* Computes equation 25 of Roberts et al. 1987, the real part of the amplitude will be stored in a_r, the imaginary part in a_i, Nf is the number of frequency points in the f spectrum */
{
  double denom;
  int i_f, i_w;

  i_f = peakfindx;
  i_w = 2*peakfindx;

  denom = 1. - W_r[i_w]*W_r[i_w] - W_i[i_w]*W_i[i_w];
  if(denom < DFTTINY)
    {
      *a_r = 0.;
      *a_i = 0.;
    }
  else
    {
      *a_r = (f_r[i_f] - f_r[i_f]*W_r[i_w] - f_i[i_f]*W_i[i_w])/denom;
      *a_i = (f_i[i_f] + f_i[i_f]*W_r[i_w] - f_r[i_f]*W_i[i_w])/denom;
    }
}

void doclean(double *R_r, double *R_i, double *W_r, double *W_i, int Nf, double *C_r, double *C_i, double gainval, double SNlimit, double df)
/* This runs the clean algorithm on the dirty spectrum D with window spectrum W. The output clean spectrum is in C, with R being replaced by the final residual spectrum. Note that we do not find the clean beam or add the residual spectrum to the clean spectrum at this stage - that will be done in a separate function. The routine will iterate until the peak in the residual power spectrum is less than SNlimit* the RMS of the residual */
{
  int i, j, Nftot;
  int peaki, iter;
  double peakval, rmsval, ave1, ave2, val, a_r, a_i, lastsum, thissum;


  Nftot = 2*Nf + 1;
  /* Initialize the clean spectrum to 0. */
  for(i=0;i<Nftot;i++)
    {
      C_r[i] = 0.;
      C_i[i] = 0.;
    }

  lastsum = 0.;
  iter = 0;
  while(1)
    {
      iter++;
      /* First find the peak in the residual power spectrum (searching only positive frequencies), also compute the RMS of the power spectrum */
      for(i=Nf, peaki=-1, peakval = -1., ave1 = 0., ave2 = 0.;i<Nftot;i++)
	{
	  val = R_r[i]*R_r[i] + R_i[i]*R_i[i];
	  if(val > peakval)
	    {
	      peakval = val;
	      peaki = i;
	    }
	  ave1 += val;
	  ave2 += val*val;
	}
      ave1 = ave1 / (double) (Nf - 1);
      ave2 = ave2 / (double) (Nf - 1);
      rmsval = sqrt(ave2 - (ave1*ave1));

      /*printf("Iteration %d: %d %f %f %f %f\n",iter,peaki,df*((double) (peaki - Nf)),peakval, rmsval, (peakval - ave1)/rmsval);*/

      /* Stop the iteration if the peak is less than SNlimit times the noise */
      if(peakval - ave1 < SNlimit * rmsval)
	break;

      /* Get the amplitude of the peak */
      getamplitude(R_r, R_i, W_r, W_i, Nf, peaki, &a_r, &a_i);

      peaki -= Nf;
      /* Calculate the new residual spectrum */
      for(j=0;j<Nftot;j++)
	{
	  R_r[j] -= gainval*(a_r*W_r[Nf+j-peaki] - a_i*W_i[Nf+j-peaki] + a_r*W_r[Nf+j+peaki] + a_i*W_i[Nf+j+peaki]);
	  R_i[j] -= gainval*(a_i*W_r[Nf+j-peaki] + a_r*W_i[Nf+j-peaki] - a_i*W_r[Nf+j+peaki] + a_r*W_i[Nf+j+peaki]);
	}

      thissum = lastsum;
      /* Update the clean spectrum */
      thissum -= (C_r[Nf + peaki]*C_r[Nf + peaki] + C_i[Nf + peaki]*C_i[Nf + peaki]);
      C_r[Nf + peaki] += gainval*a_r;
      C_i[Nf + peaki] += gainval*a_i;
      C_r[Nf - peaki] += gainval*a_r;
      C_i[Nf - peaki] -= gainval*a_i;
      thissum += (C_r[Nf + peaki]*C_r[Nf + peaki] + C_i[Nf + peaki]*C_i[Nf + peaki]);

      /* Break if the clean spectrum isn't changing appreciably any more */
      if(thissum < DFTTINY && lastsum < DFTTINY)
	break;
      if(thissum - lastsum < DFTTINY)
	break;
      if(lastsum > 0 ? (thissum < DFTTINY*lastsum) : 0)
	break;
      lastsum=thissum;
    }
}

#define BEAMTOL 1.0e-5

void GetCleanBeam(double *W_r, double *W_i, int Nf, double **B_r, int *Nb)
{
  /* This function fits a gaussian to the central peak of the window function, it will allocate memory for the beam and automatically determine an appropriate size for it */
  int i, j, Nfit, Nftot, x, nincrease;
  double centval, lastval, val, sigval, x2, f2, f3, val1, val2;

  /* First determine an appropriate size for fitting the beam */
  centval = sqrt(W_r[2*Nf]*W_r[2*Nf] + W_i[2*Nf]*W_i[2*Nf]);
  Nftot = 4*Nf + 1;
  lastval = centval;
  nincrease = 0;
  for(i=2*Nf+1,j=0;i<Nftot;i++,j++)
    {
      val = sqrt(W_r[i]*W_r[i] + W_i[i]*W_i[i]);
      if(val < BEAMTOL * centval || (val > lastval && nincrease > 2))
	break;
      if(val > lastval)
	nincrease++;
      else
	nincrease = 0;
      lastval = val;
    }
  Nfit = j;

  /*printf("%d\n",Nfit);*/

  val1 = 0.;
  val2 = 0.;
  for(i=0,j=2*Nf-(Nfit),x=-(Nfit);i<(2*(Nfit) + 1);i++,j++,x++)
    {
      x2 = (double) x*x;
      val = sqrt(W_r[j]*W_r[j] + W_i[j]*W_i[j]);
      f2 = val*val;
      f3 = log(val)*f2;
      val1 -= x2*f3;
      val2 += x2*x2*f2;
    }
  sigval = val1 / val2;

  /*printf("Clean beam sigma: %f\n",sigval);*/

  (*Nb) = ceil(sqrt(11.5/sigval));

  if(((*B_r) = (double *) malloc((2*(*Nb) + 1)* sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  for(i=0,x=-(*Nb);i<(2*(*Nb) + 1);i++,x++)
    {
      x2 = (double) -x*x*sigval;
      (*B_r)[i] = exp(x2);
    }
}

#define MIN_(A,B) ((A) < (B) ? (A) : (B))
#define MAX_(A,B) ((A) > (B) ? (A) : (B))

void Convolve_CleanBeam(double *C_r, double *C_i, double *R_r, double *R_i, double *S_r, double *S_i, int Nf, int Nb, double *B_r)
{
  int j, k,l,Nftot,Nbtot;
  Nftot = 2*Nf + 1;
  Nbtot = 2*Nb + 1;
  for(j=0;j<Nftot;j++)
    {
      S_r[j] = R_r[j];
      S_i[j] = R_i[j];
      for(k=MAX_((j-Nb),0);k<MIN_((j+Nb),Nftot);k++, l++)
	{
	  l = Nb + j - k;
	  S_r[j] += C_r[k]*B_r[l];
	  S_i[j] += C_i[k]*B_r[l];
	}
    }
}

void finddftpeaks(int Nf, double df, double *pow, int Npeaks, double *peaks, double *peakvals, double clip, int clipiter, double *SNR, int useampspec, double *stdper_out, double *aveper_out, double *stdper_noclip_out, double *aveper_noclip_out)
{
  /* This searches a power spectrum for the Npeaks highest peaks */
  double *pow_cpy, f, maxf, maxpow, lastpow;
  int Nftot, i, j, k, maxi;

  int nclippedlast, nclippedthis, Ngood;

  long double Sum, Sumsqr;
  double aveper, stdper;

  Nftot = 2*Nf + 1;

  if((pow_cpy = (double *) malloc(Nftot * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  memcpy(pow_cpy, pow, Nftot * sizeof(double));

  if(useampspec) {
    for(i=0; i < Nftot; i++)
      pow_cpy[i] = sqrt(pow_cpy[i]);
  }

  /* Compute the average and std-dev of the amplitude-spectrum */
  Sum = 0.;
  Sumsqr = 0.;
  Ngood = 0;
  for(i=Nf;i<Nftot;i++)
    {
      Sum += (long double) (pow_cpy[i]);
      Sumsqr += (long double) (pow_cpy[i]*pow_cpy[i]);
      Ngood++;
    }
  Sum /= Ngood;
  Sumsqr /= Ngood;
  aveper = (double) Sum;
  stdper = sqrt((double)(Sumsqr - (Sum*Sum)));
  *stdper_noclip_out = stdper;
  *aveper_noclip_out = aveper;
  nclippedthis = 0;
  do {
    nclippedlast = nclippedthis;
    nclippedthis = 0;
    Sum = 0.;
    Sumsqr = 0.;
    Ngood = 0;
    for(i=Nf;i<Nftot;i++)
      {
	if((pow_cpy[i]) < aveper + clip*stdper)
	  {
	    Sum += (long double) (pow_cpy[i]);
	    Sumsqr += (long double) (pow_cpy[i]*pow_cpy[i]);
	    Ngood++;
	  }
	else
	  nclippedthis++;
      }
    Sum /= Ngood;
    Sumsqr /= Ngood;
    aveper = (double) Sum;
    stdper = sqrt((double)(Sumsqr - (Sum*Sum)));
  } while(clipiter && nclippedthis > nclippedlast);

  *stdper_out = stdper;
  *aveper_out = aveper;

  for(k=0;k<Npeaks;k++)
    {
      maxpow = 0.;
      maxi = -1;
      maxf = -1.;
      for(i=Nf,f = 0.;i<Nftot;i++, f += df)
	{
	  if((pow_cpy[i]) > maxpow)
	    {
	      maxpow = (pow_cpy[i]);
	      maxf = f;
	      maxi = i;
	    }
	}
      if(maxi > 0)
	{
	  lastpow = maxpow;
	  for(j=maxi-1; j > -1; j--)
	    {
	      if((pow_cpy[j]) > lastpow)
		{
		  break;
		}
	      else
		{
		  lastpow = (pow_cpy[j]);
		  pow_cpy[j] = 0.;
		}
	    }
	  lastpow = maxpow;
	  for(j=maxi + 1; j < Nftot; j++)
	    {
	      if((pow_cpy[j]) > lastpow)
		{
		  break;
		}
	      else
		{
		  lastpow = (pow_cpy[j]);
		  pow_cpy[j] = 0.;
		}
	    }
	  pow_cpy[maxi] = 0.;
	}
      peaks[k] = maxf;
      peakvals[k] = maxpow;
      SNR[k] = (maxpow - aveper)/stdper;
    }
  free(pow_cpy);
}

void dodftclean(int N, double *t, double *mag, double *sig, int lc, _Dftclean *c, char *lcbasename, int ascii)
{
  double df, gain, SNlimit, maxfreq, T, delmin, fval;
  int Nf, nb, Nb, i, Nftot, Nwtot, Nbtot;
  FILE *outfile;
  double *R_r, *R_i, *R_pow, *C_r, *C_i, *C_pow, *W_r, *W_i, *B_r, *S_r, *S_i;
  char outname[MAXLEN];

  gain = c->gain;
  SNlimit = c->SNlimit;
  nb = c->nbeam;
  maxfreq = c->maxfreq;

  T = t[N-1] - t[0];
  /* Determine the frequency spacing and number of frequencies */
  if(maxfreq < 0)
    {
      if(t[1] != t[0])
	delmin = t[1] - t[0];
      else
	delmin = T;
      for(i=2;i<N;i++)
	{
	  if(t[i] - t[i-1] < delmin && t[i] - t[i-1] > 0.)
	    delmin = t[i] - t[i-1];
	}
      maxfreq = 1./(2.0*delmin);
    }
  Nf = floor(((double) nb) * T * maxfreq);
  df = maxfreq / (double) Nf;

  /* Allocate memory for the spectra */
  Nftot = 2*Nf + 1;
  Nwtot = 4*Nf + 1;
  if((R_r = (double *) malloc(Nftot * sizeof(double))) == NULL ||
     (R_i = (double *) malloc(Nftot * sizeof(double))) == NULL ||
     (R_pow = (double *) malloc(Nftot * sizeof(double))) == NULL ||
     (C_r = (double *) malloc(Nftot * sizeof(double))) == NULL ||
     (C_i = (double *) malloc(Nftot * sizeof(double))) == NULL ||
     (C_pow = (double *) malloc(Nftot * sizeof(double))) == NULL ||
     (W_r = (double *) malloc(Nwtot * sizeof(double))) == NULL ||
     (W_i = (double *) malloc(Nwtot * sizeof(double))) == NULL ||
     (S_r = (double *) malloc(Nftot * sizeof(double))) == NULL ||
     (S_i = (double *) malloc(Nftot * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  /* Compute the dirty spectrum */
  FDFT(t, mag, N, df, Nf, R_r, R_i);
  for(i=0;i<Nftot;i++)
    R_pow[i] = R_r[i]*R_r[i] + R_i[i]*R_i[i];

  /* Output the dirty spectrum if we are asked to */
  if(c->outdspec)
    {
      sprintf(outname,"%s/%s%s",c->dirtyspec_outdir,lcbasename,c->dirtyspec_suffix);
      if((outfile = fopen(outname,"w")) == NULL)
	error2(ERR_CANNOTWRITE,outname);
      if(ascii)
	{
	  fprintf(outfile,"#Freq Pow DFT_real DFT_imag\n");
	  for(i=0,fval=-Nf*df;i<Nftot;i++, fval += df)
	    {
	      fprintf(outfile,"%f %g %g %g\n",fval, R_pow[i], R_r[i], R_i[i]);
	    }
	}
      else
	{
	  fwrite(&Nftot,4,1,outfile);
	  fwrite(&df,8,1,outfile);
	  fwrite(R_pow,8,Nftot,outfile);
	  fwrite(R_r,8,Nftot,outfile);
	  fwrite(R_i,8,Nftot,outfile);
	}
      fclose(outfile);
    }

  /* Find the peaks in the dirty spectrum if we are asked to */
  if(c->finddirtypeaks)
    finddftpeaks(Nf, df, R_pow, c->Npeaks_dirty, c->peakfreqs_dirty[lc], c->peakpows_dirty[lc],c->clip_dirty,c->clipiter_dirty,c->SNR_dirty[lc],c->useampspec,&(c->stdper_dirty[lc]), &(c->aveper_dirty[lc]), &(c->stdper_noclip_dirty[lc]), &(c->aveper_noclip_dirty[lc]));

  /* Compute the window function if we're asked to output it, or if we are going to run clean */
  if(c->outwspec || c->runclean)
    {
      FDFT_Wfunc(t, N, df, 2*Nf, W_r, W_i);

      /* Output the window function if we are asked to */
      if(c->outwspec)
	{
	  sprintf(outname,"%s/%s%s",c->wspec_outdir,lcbasename,c->wspec_suffix);
	  if((outfile = fopen(outname,"w")) == NULL)
	    error2(ERR_CANNOTWRITE,outname);
	  if(ascii)
	    {
	      fprintf(outfile,"#Freq Wspec_real Wspec_imag\n");
	      for(i=0,fval=-2*Nf*df;i<Nwtot;i++, fval += df)
		{
		  fprintf(outfile,"%f %g %g\n",fval, W_r[i], W_i[i]);
		}
	    }
	  else
	    {
	      fwrite(&Nwtot,4,1,outfile);
	      fwrite(&df,8,1,outfile);
	      fwrite(W_r,8,Nwtot,outfile);
	      fwrite(W_i,8,Nwtot,outfile);
	    }
	  fclose(outfile);
	}
    }

  /* Do clean if we are asked to */
  if(c->runclean)
    {
      doclean(R_r, R_i, W_r, W_i, Nf, C_r, C_i, gain, SNlimit, df);
      GetCleanBeam(W_r, W_i, Nf, &B_r, &Nb);
      Convolve_CleanBeam(C_r, C_i, R_r, R_i, S_r, S_i, Nf, Nb, B_r);

      for(i=0;i<Nftot;i++)
	C_pow[i] = S_r[i]*S_r[i] + S_i[i]*S_i[i];

      /* Output the clean beam if asked to */
      if(c->outcbeam)
	{
	  sprintf(outname,"%s/%s%s",c->cbeam_outdir,lcbasename,c->cbeam_suffix);
	  if((outfile = fopen(outname,"w")) == NULL)
	    error2(ERR_CANNOTWRITE,outname);
	  if(ascii)
	    {
	      fprintf(outfile,"#Freq C_beam\n");
	      for(i=0,fval=-Nb*df;i<(2*Nb + 1);i++, fval += df)
		{
		  fprintf(outfile,"%f %g\n",fval, B_r[i]);
		}
	    }
	  else
	    {
	      Nbtot = 2*Nb + 1;
	      fwrite(&Nbtot,4,1,outfile);
	      fwrite(&df,8,1,outfile);
	      fwrite(B_r,8,Nbtot,outfile);
	    }
	  fclose(outfile);
	}

      /* Output the clean spectrum if asked to */
      if(c->outcspec)
	{
	  sprintf(outname,"%s/%s%s",c->cspec_outdir,lcbasename,c->cspec_suffix);
	  if((outfile = fopen(outname,"w")) == NULL)
	    error2(ERR_CANNOTWRITE,outname);
	  if(ascii)
	    {
	      fprintf(outfile,"#Freq Pow_clean DFT_real_clean DFT_imag_clean DFT_real_clean_deconv DFT_imag_clean_deconv\n");
	      for(i=0,fval=-Nf*df;i<Nftot;i++, fval += df)
		{
		  fprintf(outfile,"%f %g %g %g %g %g\n",fval, C_pow[i], S_r[i], S_i[i], C_r[i], C_i[i]);
		}
	    }
	  else
	    {
	      fwrite(&Nftot,4,1,outfile);
	      fwrite(&df,8,1,outfile);
	      fwrite(C_pow,8,Nftot,outfile);
	      fwrite(S_r,8,Nftot,outfile);
	      fwrite(S_i,8,Nftot,outfile);
	      fwrite(C_r,8,Nftot,outfile);
	      fwrite(C_i,8,Nftot,outfile);
	    }
	  fclose(outfile);
	}

      /* Get the peaks in the clean spectrum if asked to */
      if(c->findcleanpeaks)
	finddftpeaks(Nf, df, C_pow, c->Npeaks_clean, c->peakfreqs_clean[lc], c->peakpows_clean[lc],c->clip_clean,c->clipiter_clean,c->SNR_clean[lc],c->useampspec,&(c->aveper_clean[lc]), &(c->stdper_clean[lc]), &(c->aveper_noclip_clean[lc]), &(c->stdper_noclip_clean[lc]));
    }

  /* Free all the allocated vectors */

  free(R_r);
  free(R_i);
  free(R_pow);
  free(C_r);
  free(C_i);
  free(C_pow);
  free(W_r);
  free(W_i);
  free(S_r);
  free(S_i);

}
