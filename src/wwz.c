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
/* This file includes functions used by the -wwz command in vartools
   This command is an implementation of the Wavelet Transform for 
   non-uniformly sampled data as defined by Foster 1996, AJ, 112, 1709.
*/

#include "commands.h"
#include "programdata.h"
#include "functions.h"

#ifdef USECFITSIO
#include "fitsio.h"
#endif

#ifndef TWOPI
#define TWOPI 6.28318530717958647692528676656
#endif

#define WWZ_TINY 1.0e-9

void invertsym3x3matrix(double *A, double *Ainv)
/* Computes inverse of a symmetric matrix A, returns result in Ainv; Assumes
   A[i][j] = i*3 + j, and A[i][j] = A[j][i];
*/
{
  double Adet;

  Ainv[0] = A[8]*A[4] - A[5]*A[5];
  Ainv[1] = A[2]*A[5] - A[8]*A[1];
  Ainv[2] = A[1]*A[5] - A[2]*A[4];
  Adet = Ainv[0]*A[0] + Ainv[1]*A[1] + Ainv[2]*A[2];
  Ainv[4] = A[0]*A[8] - A[2]*A[2];
  Ainv[5] = A[2]*A[1] - A[0]*A[5];
  Ainv[8] = A[0]*A[4] - A[1]*A[1];
  Ainv[0] /= Adet;
  Ainv[1] /= Adet;
  Ainv[2] /= Adet;
  Ainv[4] /= Adet;
  Ainv[5] /= Adet;
  Ainv[8] /= Adet;
  Ainv[3] = Ainv[1];
  Ainv[6] = Ainv[2];
  Ainv[7] = Ainv[5];
}


void WWZA(double *x, double *y, int N, double c, double f0, double df, int Nf, double tau0, double dtau, int Ntau, int savefull, double **ya, double *WWZ_out, double *WWA_out, double *neff_out, double *mz_out, double *mfreq_out, double *mpow_out, double *mamp_out, double *mcon_out, double *mneff_out)
/* This computes the WWZ and WWA transforms of a signal with uneven time
   sampling. The method is that of Foster 1996, AJ, 112, 1709.

   The function is based on the FORTRAN implementation of WWZ by M. Templeton
   version 1.1 wwz11.f distributed by the AAVSO, which in turn is based on the
   original BASIC code due to Foster.

   x - input times of observation
   y - input measurements (magnitudes)
   N - number of observations
   c - c parameter in wavelet transform
   f0 - minimum frequency for transform (cycles per day)
   df - frequency step size (cycles per day)
   Nf - number of frequencies
   tau0 - minimum time shift (days)
   dtau - time shift step size (days)
   Ntau - number of time shifts
   savefull - whether or not to save the full wavelet in ya, WWZ_out,
              WWA_out, and neff_out. If this is 1 then it will be saved,
              if it is 0 then it will not be saved. ya, WWZ_out, WWA_out and
              neff_out may be NULL if savefull == 0.
   ya - coefficients for the wavelets.  This matrix should have
        dimensions [3][Ntau*Nf]. And is indexed such that
        ya[k][i] is the coefficient for trial function k (0 = 1;
           1 = cos(2*pi*f(t-tau)); 2 = sin(2*pi*f(t-tau)))
        and i = tau_i*Nf + f_i where f_i indexes the frequency and 
            tau_i indexes the time shift.
   WWZ_out - the WWZ transform. This vector should have dimension
        Ntau*Nf; It is indexed in the same manner as the second
        component of ya.
   WWA_out - the WWA transform. This vector should have dimension
        Ntau*Nf; It is indexed in the same manner as the second
        component of ya.
   neff_out - the effective number of data points. This vector should have 
        dimension Ntau*Nf; It is indexed in the same manner as the second
        component of ya.
   mz_out - the maximum value of WWZ as a function of tau; This should have
        dimension Ntau;
   mfreq_out - the frequency at which WWZ is maximum as a function tau; This
        should have dimension Ntau;
   mpow_out - the power at which WWZ is maximum as a function of tau; This
        should have dimension Ntau;
   mamp_out - The amplitude (WWA) at which WWZ is maximum as a function of
        tau; This should have dimension Ntau;
   mcon_out - The signal mean value as a function of tau; This should have
        dimension Ntau;
   mneff_out - The effective number of datapoints at which WWZ is maximum as
        a function of tau; This should have dimension Ntau;

 */
{
  int i, j, k, l, ll;
  
  double tau, f, phi1, phi2, phi3, coeff[3], S[9], Sinv[9], V[3], mz, z, weight, weight2;
  double cc, cw, ss, sw, xw, varw = 0.0, damp, dpower, neff, avew, powz;

  int indx, kstart;

  f0 = TWOPI * f0;
  df = TWOPI * df;

  phi1 = 1.0;

  for(i=0,indx=0,tau=tau0;i<Ntau;i++,tau += dtau)
    {
      kstart = 0;
      mz = 0.0;
      for(j=0,f=f0;j<Nf;j++,f+=df,indx++)
	{
	  for(l=0; l < 3; l++)
	    V[l] = 0.;
	  for(l=0; l < 9; l++)
	    S[l] = 0.;
	  weight2 = 0.;
	  varw = 0.0;
	  for(k=kstart; k < N; k++) {
	    z = f*(x[k]-tau);
	    weight = exp(-1.0*c*z*z);
	    if(weight > WWZ_TINY) {
	      cc = cos(z);
	      cw = weight*cc;
	      ss = sin(z);
	      sw = weight*ss;
	      S[0] += weight;
	      weight2 += weight*weight;
	      S[1] += cw;
	      S[2] += sw;
	      S[4] += cw*cc;
	      S[5] += cw*ss;
	      S[8] += sw*ss;
	      xw = weight*y[k];
	      V[0] += xw;
	      varw += xw*y[k];
	      V[1] += cw*y[k];
	      V[2] += sw*y[k];
	    } else {
	      if(z > 0.)
		break;
	      else
		kstart++;
	    }
	  }
	  S[3] = S[1];
	  S[6] = S[2];
	  S[7] = S[5];
	  dpower = 0.;
	  damp = 0.;
	  for(l=0; l < 3; l++)
	    coeff[l] = 0.;
	  if(weight2 > 0.)
	    neff=(S[0]*S[0])/weight2;
	  else
	    neff = 0.;
	  if(neff > 3.0) {
	    for(l=0; l < 3; l++) {
	      V[l] = V[l]/S[0];
	    }
	    for(l=1; l < 9; l++)
	      S[l] /= S[0];
	    if(S[0] > 0.0) 
	      varw = varw / S[0];
	    else
	      varw = 0.;
	    S[0] = 1.0;
	    avew = V[0];
	    varw = varw - avew*avew;
	    if(varw <= 0.0) varw = 1.0e-12;
	    invertsym3x3matrix(S, Sinv);
	    for(l=0; l < 3; l++) {
	      for(ll=0; ll < 3; ll++) {
		coeff[l] += Sinv[3*l+ll]*V[ll];
	      }
	      dpower += coeff[l]*V[l];
	    }
	    dpower = dpower - avew*avew;
	    powz = (neff - 3.0)*dpower/(varw-dpower)/2.0;
	    dpower = (neff - 1.0)*dpower/varw/2.0;
	    damp = sqrt((coeff[1]*coeff[1])+(coeff[2]*coeff[2]));
	  } else {
	    powz = 0.;
	    dpower = 0.;
	    damp = 0.;
	    if(neff < 1.0e-9) neff=0.0;
	  }
	  //if(damp < 1.0e-9) damp = 0.0;
	  //if(dpower < 1.0e-9) dpower = 0.0;
	  //if(powz < 1.0e-9) powz = 0.0;
	  if(savefull) {
	    WWZ_out[indx] = powz;
	    WWA_out[indx] = damp;
	    neff_out[indx] = neff;
	    for(l=0; l < 3; l++)
	      ya[l][indx] = coeff[l];
	  }
	  if(powz > mz) {
	    mz = powz;
	    mz_out[i] = powz;
	    mfreq_out[i] = f/TWOPI;
	    mpow_out[i] = dpower;
	    mamp_out[i] = damp;
	    mcon_out[i] = coeff[0];
	    mneff_out[i] = neff;
	  }
	}
    }
}

#ifdef USECFITSIO
void WriteFitsImage_WWZ(char *outname, double f0, double df, int Nf,
			double tau0, double dtau, int Ntau,
			double **ya, double *WWZ_out, double *WWA_out,
			double *neff_out) {
  fitsfile *fptr;
  int status = 0, ii, jj, k;
  long fpixel, nelements, exposure;
  int bitpix = DOUBLE_IMG;
  int Nextension = 6;
  long naxis = 2;
  long naxes[2];
  long tmp;
  double dtmp;
  double **extensions;
  char strtmp[64];

  naxes[0] = Nf;
  naxes[1] = Ntau;
  
  remove(outname); /* delete old file if it already exists */
  status = 0;

  extensions = (double **) malloc(6*sizeof(double *));
  extensions[0] = WWZ_out;
  extensions[1] = WWA_out;
  extensions[2] = neff_out;
  extensions[3] = ya[0];
  extensions[4] = ya[1];
  extensions[5] = ya[2];

  if (fits_create_file(&fptr, outname, &status)) {
    error2(ERR_CANNOTWRITE,outname);
  }
  
  for(k=0; k < Nextension; k++) {

    if(!k) {
      if ( fits_create_img(fptr, bitpix, naxis, naxes, &status)) {
	error2(ERR_CANNOTWRITE,outname);
      }
    } else {
      if ( fits_insert_img(fptr, bitpix, naxis, naxes, &status)) {
	error2(ERR_CANNOTWRITE,outname);
      }
    }
  
    fpixel = 1;
    nelements = naxes[0] * naxes[1];

    if( fits_write_img(fptr, TDOUBLE, fpixel, nelements, extensions[k], &status) )
      error2(ERR_CANNOTWRITE,outname);
  
    tmp = 2;
    fits_update_key(fptr, TLONG, "WCSAXES", &tmp, "", &status);
    dtmp = f0;
    fits_update_key(fptr, TDOUBLE, "CRVAL1", &dtmp, "", &status);
    dtmp = 1.0;
    fits_update_key(fptr, TDOUBLE, "CRPIX1", &dtmp, "", &status);
    dtmp = df;
    fits_update_key(fptr, TDOUBLE, "CDELT1", &dtmp, "", &status);
    sprintf(strtmp,"Frequency");
    fits_update_key(fptr, TSTRING, "CTYPE1", strtmp, "", &status);
    sprintf(strtmp,"Cyles Per Day");
    fits_update_key(fptr, TSTRING, "CUNIT1", strtmp, "", &status);
    dtmp = tau0;
    fits_update_key(fptr, TDOUBLE, "CRVAL2", &dtmp, "", &status);
    dtmp = 1.0;
    fits_update_key(fptr, TDOUBLE, "CRPIX2", &dtmp, "", &status);
    dtmp = dtau;
    fits_update_key(fptr, TDOUBLE, "CDELT2", &dtmp, "", &status);
    sprintf(strtmp,"Time Shift");
    fits_update_key(fptr, TSTRING, "CTYPE2", strtmp, "", &status);
    sprintf(strtmp,"Days");
    fits_update_key(fptr, TSTRING, "CUNIT2", strtmp, "", &status);
    /*
    dtmp = 1.0;
    fits_update_key(fptr, TDOUBLE, "PC1_1", &dtmp, "", &status);
    dtmp = 1.0;
    fits_update_key(fptr, TDOUBLE, "PC2_2", &dtmp, "", &status);
    dtmp = 0.0;
    fits_update_key(fptr, TDOUBLE, "PC1_2", &dtmp, "", &status);
    dtmp = 0.0;
    fits_update_key(fptr, TDOUBLE, "PC2_1", &dtmp, "", &status);
    dtmp = 0.0;
    fits_update_key(fptr, TDOUBLE, "CD1_1", &dtmp, "", &status);
    dtmp = 0.0;
    fits_update_key(fptr, TDOUBLE, "CD1_2", &dtmp, "", &status);
    dtmp = 0.0;
    fits_update_key(fptr, TDOUBLE, "CD2_1", &dtmp, "", &status);
    dtmp = 0.0;
    fits_update_key(fptr, TDOUBLE, "CD2_2", &dtmp, "", &status);
    */
    switch(k) {
    case 0:
      sprintf(strtmp,"WWZ");
      break;
    case 1:
      sprintf(strtmp,"WWA");
      break;
    case 2:
      sprintf(strtmp,"Neff");
      break;
    case 3:
      sprintf(strtmp,"Constant Coeff");
      break;
    case 4:
      sprintf(strtmp,"Cosine Coeff");
      break;
    case 5:
      sprintf(strtmp,"Sine Coeff");
      break;
    }
    fits_update_key(fptr, TSTRING, "ARRAY", strtmp, "", &status);
  }

  if (fits_close_file(fptr, &status)) {
    error2(ERR_CANNOTWRITE,outname);
  }

  return;

}
#endif

void WriteAsciiData_WWZ(char *outfullfilename, double f0, double df, int Nf, 
			double tau0, double dtau, int Ntau,
			double **ya, double *WWZ_out, double *WWA_out, 
			double *neff_out, int ispm3d) {
  FILE *outfile;
  int indx, i, j;
  double f, tau;

  if((outfile = fopen(outfullfilename,"w")) == NULL)
    error2(ERR_CANNOTWRITE,outfullfilename);

  fprintf(outfile,"# Time_Shift   Frequency   WWZ   WWA   Neff   Constant_Coeff   Cosine_Coeff  Sine_Coeff\n");
  fprintf(outfile,"# [1]          [2]         [3]   [4]   [5]    [6]              [7]           [8]\n");

  for(i=0,indx=0,tau=tau0;i<Ntau;i++,tau += dtau)
    {
      for(j=0,f=f0;j<Nf;j++,f+=df,indx++)
	{
	  fprintf(outfile,
		  "%.17g %.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
		  tau, f, WWZ_out[indx], WWA_out[indx], neff_out[indx],
		  ya[0][indx], ya[1][indx], ya[2][indx]);
	}
      if(ispm3d && i < (Ntau-1)) {
	fprintf(outfile,"\n");
      }
    }
  fclose(outfile);
}


void WriteAsciiDataMax_WWZ(char *outmaxfilename, double tau0, double dtau, 
			   int Ntau, double *mz_out, double *mfreq_out, 
			   double *mpow_out, double *mamp_out, 
			   double *mcon_out, double *mneff_out) {
  FILE *outfile;
  int i;
  double tau;
  
  if((outfile = fopen(outmaxfilename,"w")) == NULL)
    error2(ERR_CANNOTWRITE,outmaxfilename);
  
  fprintf(outfile,"# Time_Shift   WWZ   Frequency   Power   Amplitude   Local_Average   Neff\n");
  fprintf(outfile,"# [1]          [2]   [3]         [4]     [5]         [6]             [7]\n");
  
  for(i=0,tau=tau0;i<Ntau;i++,tau += dtau) {
    fprintf(outfile,"%.17g %.17g %.17g %.17g %.17g %.17g %.17g\n",
	    tau, mz_out[i], mfreq_out[i], mpow_out[i], mamp_out[i],
	    mcon_out[i], mneff_out[i]);
  }
  fclose(outfile);
}

	  
void DoWWZ(ProgramData *p, _WWZ *c, int threadid, int lcid) {
  /* This function executes the -WWZ command */
  int savefull;
  char outfullfilename[MAXLEN];
  char outmaxfilename[MAXLEN];
  int NJD, i; 
  double *t;
  double *mag, T, delmin, maxfreq, FreqSampFact;
  double cterm, f0, df, tau0, dtau, tau1, 
    **ya = NULL, *WWZ_out = NULL, *WWA_out = NULL, *neff_out = NULL,
    *mz_out = NULL, *mfreq_out = NULL, *mpow_out = NULL, *mamp_out = NULL, 
    *mcon_out = NULL, *mneff_out = NULL;
  int Nf, Ntau, Nterms;
  double *t_mask = NULL;
  double *mag_mask = NULL;

  if(p->NJD[threadid] <= 1) {
    c->max_z[threadid] = 0.;
    c->max_freq[threadid] = 0.;
    c->max_pow[threadid] = 0.;
    c->max_amp[threadid] = 0.;
    c->max_neff[threadid] = 0.;
    c->max_tau[threadid] = 0.;
    c->max_con[threadid] = 0.;
    c->med_z[threadid] = 0.;
    c->med_freq[threadid] = 0.;
    c->med_pow[threadid] = 0.;
    c->med_amp[threadid] = 0.;
    c->med_neff[threadid] = 0.;
    c->med_con[threadid] = 0.;
    return;
  }


  if(!c->usemask) {
    NJD = p->NJD[threadid];
    t = p->t[threadid];
    mag = p->mag[threadid];
  } else {
    if((t_mask = (double *) malloc(p->NJD[threadid]*sizeof(double))) == NULL ||
       (mag_mask = (double *) malloc(p->NJD[threadid]*sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    NJD = 0;
    for(i=0; i < p->NJD[threadid]; i++) {
      if(!isnan(p->mag[threadid][i]) && EvaluateVariable_Double(lcid, threadid, i, c->maskvar) > VARTOOLS_MASK_TINY) {
	t_mask[NJD] = p->t[threadid][i];
	mag_mask[NJD] = p->mag[threadid][i];
	NJD++;
      }
    }
    t = t_mask;
    mag = mag_mask;
    if(NJD <= 1) {
      c->max_z[threadid] = 0.;
      c->max_freq[threadid] = 0.;
      c->max_pow[threadid] = 0.;
      c->max_amp[threadid] = 0.;
      c->max_neff[threadid] = 0.;
      c->max_tau[threadid] = 0.;
      c->max_con[threadid] = 0.;
      c->med_z[threadid] = 0.;
      c->med_freq[threadid] = 0.;
      c->med_pow[threadid] = 0.;
      c->med_amp[threadid] = 0.;
      c->med_neff[threadid] = 0.;
      c->med_con[threadid] = 0.;
      if(t_mask != NULL) free(t_mask);
      if(mag_mask != NULL) free(mag_mask);
      return;
    }
  }
      

  cterm = c->cterm;
  maxfreq = c->maxfreq;
  FreqSampFact = c->freq_sample_factor;
  T = (t[NJD-1] - t[0]);
  

  delmin = -1.0;

  /* Determine the maximum frequency if it is automated */
  if(maxfreq < 0)
    {
      if(t[1] != t[0])
	delmin = t[1] - t[0];
      else
	delmin = T;
      for(i=2; i < NJD; i++) {
	if(t[i] - t[i-1] < delmin && t[i] - t[i-1] > 0.)
	  delmin = t[i] - t[i-1];
      }
      maxfreq = 1.0/(2.0*delmin);
    }

  Nf = floor((T / FreqSampFact) * maxfreq);
  df = maxfreq / (double) Nf;
  f0 = df;

  if(!c->auto_tau0) {
    tau0 = c->tau0;
  } else {
    tau0 = t[0];
  }
  if(!c->auto_tau1) {
    tau1 = c->tau1;
  } else {
    tau1 = t[NJD - 1];
  }
  if(!c->auto_dtau) {
    dtau = c->dtau;
  } else {
    if(delmin <= 0.0) {
      if(t[1] != t[0])
	delmin = t[1] - t[0];
      else
	delmin = T;
      for(i=2; i < NJD; i++) {
	if(t[i] - t[i-1] < delmin && t[i] - t[i-1] > 0.)
	  delmin = t[i] - t[i-1];
      }
    }
    dtau = delmin;
  }
  Ntau = floor((tau1 - tau0)/dtau);

  if(Ntau <= 0 || Nf <= 0) {
    c->max_z[threadid] = 0.;
    c->max_freq[threadid] = 0.;
    c->max_pow[threadid] = 0.;
    c->max_amp[threadid] = 0.;
    c->max_neff[threadid] = 0.;
    c->max_tau[threadid] = 0.;
    c->max_con[threadid] = 0.;
    c->med_z[threadid] = 0.;
    c->med_freq[threadid] = 0.;
    c->med_pow[threadid] = 0.;
    c->med_amp[threadid] = 0.;
    c->med_neff[threadid] = 0.;
    c->med_con[threadid] = 0.;
    if(t_mask != NULL) free(t_mask);
    if(mag_mask != NULL) free(mag_mask);
    return;
  }


  if(c->outfulltransform) {
#ifdef USECFITSIO
    if(c->outfulltransform_usefits) {
      GetOutputFilename(outfullfilename,p->lcnames[lcid],
			c->outfulltransform_dir,
			"wwz.fits", c->outfulltransform_format, lcid);
    } else {
#endif
      GetOutputFilename(outfullfilename,p->lcnames[lcid],
			c->outfulltransform_dir,
			"wwz", c->outfulltransform_format, lcid);
#ifdef USECFITSIO
    }
#endif
    savefull = 1;
  }
  else {
    outfullfilename[0] = '\0';
    savefull = 0;
  }
  if(c->outmaxtransform) {
    GetOutputFilename(outmaxfilename,p->lcnames[lcid],
		      c->outmaxtransform_dir,
		      "mwwz", c->outmaxtransform_format, lcid);
  }
  else {
    outmaxfilename[0] = '\0';
  }

  Nterms = Ntau*Nf;
  if(savefull) {
    if((ya = (double **) malloc(3*sizeof(double *))) == NULL ||
       (WWZ_out = (double *) malloc(Nterms*sizeof(double))) == NULL ||
       (WWA_out = (double *) malloc(Nterms*sizeof(double))) == NULL ||
       (neff_out = (double *) malloc(Nterms*sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < 3; i++) {
      if((ya[i] = (double *) malloc(Nterms*sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  }

  if((mz_out = (double  *) malloc(Ntau*sizeof(double))) == NULL ||
     (mfreq_out = (double *) malloc(Ntau*sizeof(double))) == NULL ||
     (mpow_out = (double *) malloc(Ntau*sizeof(double))) == NULL ||
     (mamp_out = (double *) malloc(Ntau*sizeof(double))) == NULL ||
     (mcon_out = (double *) malloc(Ntau*sizeof(double))) == NULL ||
     (mneff_out = (double *) malloc(Ntau*sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  /* Run the transform */
  WWZA(t, mag, NJD, cterm, f0, df, Nf, tau0, dtau, Ntau, savefull, ya, WWZ_out,
       WWA_out, neff_out, mz_out, mfreq_out, mpow_out, mamp_out, mcon_out,
       mneff_out);

  /* Write out the full transform */
  if(savefull) {
#ifdef USECFITSIO
    if(c->outfulltransform_usefits) {
      WriteFitsImage_WWZ(outfullfilename, f0, df, Nf, tau0, dtau, Ntau,
			 ya, WWZ_out, WWA_out, neff_out);
    } else {
#endif
      WriteAsciiData_WWZ(outfullfilename, f0, df, Nf, tau0, dtau, Ntau,
			 ya, WWZ_out, WWA_out, neff_out,
			 c->outfulltransform_usepm3d);
#ifdef USECFITSIO
    }
#endif
  }

  if(c->outmaxtransform) {
    WriteAsciiDataMax_WWZ(outmaxfilename, tau0, dtau, Ntau, mz_out, 
			  mfreq_out, mpow_out, mamp_out, mcon_out, mneff_out);
  }

  /* Get the maximum over tau */
  c->max_z[threadid] = 0.;
  for(i=0; i < Ntau; i++) {
    if(mz_out[i] > c->max_z[threadid]) {
      c->max_z[threadid] = mz_out[i];
      c->max_freq[threadid] = mfreq_out[i];
      c->max_pow[threadid] = mpow_out[i];
      c->max_amp[threadid] = mamp_out[i];
      c->max_con[threadid] = mcon_out[i];
      c->max_neff[threadid] = mneff_out[i];
      c->max_tau[threadid] = tau0 + dtau*i;
    }
  }

  c->med_z[threadid] = median(Ntau, mz_out);
  c->med_freq[threadid] = median(Ntau, mfreq_out);
  c->med_pow[threadid] = median(Ntau, mpow_out);
  c->med_amp[threadid] = median(Ntau, mamp_out);
  c->med_neff[threadid] = median(Ntau, mneff_out);
  c->med_con[threadid] = median(Ntau, mcon_out);

  if(ya != NULL) {
    for(i=0; i < 3; i++) {
      if(ya[i] != NULL)
	free(ya[i]);
    }
    free(ya);
  }
  if(WWZ_out != NULL)
    free(WWZ_out);
  if(WWA_out != NULL)
    free(WWA_out);
  if(neff_out != NULL)
    free(neff_out);
  if(mz_out != NULL)
    free(mz_out);
  if(mfreq_out != NULL)
    free(mfreq_out);
  if(mpow_out != NULL)
    free(mpow_out);
  if(mamp_out != NULL)
    free(mamp_out);
  if(mcon_out != NULL)
    free(mcon_out);
  if(mneff_out != NULL)
    free(mneff_out);
  if(t_mask != NULL) free(t_mask);
  if(mag_mask != NULL) free(mag_mask);

  return;

}

int ParseWWZCommand(int *iret, int argc, char **argv, ProgramData *p,
		    _WWZ *c, Command *cs)
{
  int i;

  i = *iret;

  c->cterm = 0.0126651479552922;
  c->maxfreq = -1.;
  c->freq_sample_factor = 0.25;
  c->auto_tau0 = 1;
  c->auto_tau1 = 1;
  c->auto_dtau = 1;
  c->outfulltransform = 0;
  c->outfulltransform_usefits = 0;
  c->outfulltransform_usepm3d = 0;
  c->outmaxtransform = 0;

  i++;
  if(i >= argc) {*iret = i; return 1;}
  if(!strcmp(argv[i],"maxfreq")) {
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if(!strcmp(argv[i],"auto")) {
	c->maxfreq = -1.;
      } else {
	c->maxfreq = atof(argv[i]);
	if(c->maxfreq <= 0.0)
	  error2(ERR_INVALID_PARAMETERVALUE,"-wwz, maxfreq must be > 0");
      }
  } else {*iret = i; return 1;} 

  i++;
  if(i >= argc) {*iret = i; return 1;}
  if(!strcmp(argv[i],"freqsamp")) {
    i++;
    if(i >= argc) {*iret = i; return 1;}
    c->freq_sample_factor = atof(argv[i]);
    if(c->freq_sample_factor <= 0) {
      error2(ERR_INVALID_PARAMETERVALUE,"-wwz, freqsamp must be > 0");
    }
  } else {*iret = i; return 1;}
  
  i++;
  if(i >= argc) {*iret = i; return 1;}
  if(!strcmp(argv[i],"tau0")) {
    i++;
    if(i >= argc) {*iret = i; return 1;}
    if(!strcmp(argv[i],"auto")) {
      c->auto_tau0 = 1;
    } else {
      c->auto_tau0 = 0;
      c->tau0 = atof(argv[i]);
    }
  } else {*iret = i; return 1;}

  i++;
  if(i >= argc) {*iret = i; return 1;}
  if(!strcmp(argv[i],"tau1")) {
    i++;
    if(i >= argc) {*iret = i; return 1;}
    if(!strcmp(argv[i],"auto")) {
      c->auto_tau1 = 1;
    } else {
      c->auto_tau1 = 0;
      c->tau1 = atof(argv[i]);
      if(!c->auto_tau0) {
	if(c->tau1 <= c->tau0) {
	  error2(ERR_INVALID_PARAMETERVALUE,"-wwz, tau0 must be < tau1");
	}
      }
    }
  } else {*iret = i; return 1;}

  i++;
  if(i >= argc) {*iret = i; return 1;}
  if(!strcmp(argv[i],"dtau")) {
    i++;
    if(i >= argc) {*iret = i; return 1;}
    if(!strcmp(argv[i],"auto")) {
      c->auto_dtau = 1;
    } else {
      c->auto_dtau = 0;
      c->dtau = atof(argv[i]);
      if(c->dtau <= 0) {
	error2(ERR_INVALID_PARAMETERVALUE,"-wwz, dtau must be > 0");
      }
    }
  } else {*iret = i; return 1;}

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"c")) {
      i++;
      if(i >= argc) {*iret = i; return 1;}
      c->cterm = atof(argv[i]);
      if(c->cterm <= 0) {
	error2(ERR_INVALID_PARAMETERVALUE,"-wwz, c must be > 0");
      }
    } else
      i--;
  } else 
    i--;


  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"outfulltransform")) {
      c->outfulltransform = 1;
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->outfulltransform_dir = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->outfulltransform_dir,"%s",argv[i]);
      
      if((c->outfulltransform_format = (char *) malloc(1)) == NULL)
	error(ERR_MEMALLOC);
      c->outfulltransform_format[0] = '\0';
      
      /* Check if the user gave the "fits" keyword */
      c->outfulltransform_usefits = 0;
#ifdef USECFITSIO
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"fits")) {
	  c->outfulltransform_usefits = 1;
	} else 
	  i--;
      } else 
	i--;
#endif      

      c->outfulltransform_usepm3d = 0;
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"pm3d")) {
	  c->outfulltransform_usepm3d = 1;
	} else 
	  i--;
      } else 
	i--;
      
      if(c->outfulltransform_usefits && 
	 c->outfulltransform_usepm3d) {
	error2(ERR_INVALID_PARAMETERVALUE,"-wwz, cannot use both fits and pm3d keywords");
      }

      /* Check if the user gave the "format" keyword */
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"format")) {
	  i++;
	  if(i < argc) {
	    c->outfulltransform_format = (char *) realloc(c->outfulltransform_format, (strlen(argv[i])+1)*sizeof(char));
	    sprintf(c->outfulltransform_format,"%s",argv[i]);
	    i++;
	  } else {
	    *iret = i; return 1;
	  }
	} else
	  i--;
      }
      else
	i--;
    } else
      i--;
  } else
    i--;
  
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"outmaxtransform")) {
      c->outmaxtransform = 1;
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if((c->outmaxtransform_dir = (char *) malloc((strlen(argv[i])+1))) == NULL)
	error(ERR_MEMALLOC);
      sprintf(c->outmaxtransform_dir,"%s",argv[i]);
      
      if((c->outmaxtransform_format = (char *) malloc(1)) == NULL)
	error(ERR_MEMALLOC);
      c->outmaxtransform_format[0] = '\0';
      
      /* Check if the user gave the "format" keyword */
      i++;
      if(i < argc) {
	if(!strcmp(argv[i],"format")) {
	  i++;
	  if(i < argc) {
	    c->outmaxtransform_format = (char *) realloc(c->outmaxtransform_format, (strlen(argv[i])+1)*sizeof(char));
	    sprintf(c->outmaxtransform_format,"%s",argv[i]);
	    i++;
	  } else {
	    *iret = i; return 1;
	  }
	} else
	  i--;
      }
      else
	i--;
    } else
      i--;
  } else
    i--;

  c->usemask = 0;
  c->maskvar = NULL;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"maskpoints")) {
      c->usemask = 1;
      i++;
      if(i >= argc) {*iret = i; return 1;}
      parse_setparam_existingvariable(cs, argv[i], &(c->maskvar), VARTOOLS_VECTORTYPE_LC, VARTOOLS_TYPE_NUMERIC);
    } else
      i--;
  } else 
    i--;

  RegisterScalarData(p, (void *)(&(c->max_z)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->max_freq)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->max_pow)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->max_amp)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->max_neff)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->max_tau)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->max_con)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->med_z)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->med_freq)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->med_pow)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->med_amp)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->med_neff)), VARTOOLS_TYPE_DOUBLE, 0);
  RegisterScalarData(p, (void *)(&(c->med_con)), VARTOOLS_TYPE_DOUBLE, 0);

  *iret = i;
  return 0;
  
}
