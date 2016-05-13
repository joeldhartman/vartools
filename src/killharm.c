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

/* This file contains functions to remove periodic signals from light curves for the program vartools by J. Hartman */

#define TWOPI 6.28318530717958647692528676656

#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f)

#define TOLERANCE 1.0E-9
#define ITMAX 100
#define ZEPS 1.0e-10

double f_tomin(double x, int minmaxflag, int Nsubharm, int Nharm, double **subharmA, double **subharmB, double **harmA, double **harmB, double *fundA, double *fundB, double period, int p_per)
{
  int p, k, freqnum;
  double val, st, ct, ct_, st_, d1, d2;
  val = 0.;
  p = p_per;
  for(k=0;k<Nsubharm;k++)
    {
      freqnum = (Nsubharm - k) + 1;
      st = sin(x * TWOPI / (double) freqnum);
      ct = cos(x * TWOPI / (double) freqnum);
      val += st * subharmA[p][k];
      val += ct * subharmB[p][k];
    }
  st = sin(x * TWOPI);
  ct = cos(x * TWOPI);
  st_ = st;
  ct_ = ct;
  val += st_ * fundA[p];
  val += ct_ * fundB[p];
  for(k=0;k<Nharm;k++)
    {
      d1 = ct_*st + st_*ct;
      d2 = ct_*ct - st_*st;
      st_ = d1;
      ct_ = d2;
      val += st_ * harmA[p][k];
      val += ct_ * harmB[p][k];
    }
  val = val*minmaxflag;
  return(val);
}

double dfdx_tomin(double x, int minmaxflag, int Nsubharm, int Nharm, double **subharmA, double **subharmB, double **harmA, double **harmB, double *fundA, double *fundB, double period, int p_per)
{
  int p, freqnum, k;
  double val, st, ct, st_, ct_, d1, d2;
  val = 0.;
  p = p_per;
  for(k=0;k<Nsubharm;k++)
    {
      freqnum = (Nsubharm - k) + 1;
      st = sin(x * TWOPI / (double) freqnum);
      ct = cos(x * TWOPI / (double) freqnum);
      val += ct * subharmA[p][k] * TWOPI / (double) freqnum;
      val -= st * subharmB[p][k] * TWOPI / (double) freqnum;
    }
  st = sin(x * TWOPI);
  ct = cos(x * TWOPI);
  st_ = st;
  ct_ = ct;
  val += ct_ * fundA[p] * TWOPI;
  val -= st_ * fundB[p] * TWOPI;
  for(k=0;k<Nharm;k++)
    {
      d1 = ct_*st + st_*ct;
      d2 = ct_*ct - st_*st;
      st_ = d1;
      ct_ = d2;
      val += ct_ * harmA[p][k] * TWOPI * ((double) (k + 2));
      val -= st_ * harmB[p][k] * TWOPI * ((double) (k + 2));
    }
  val = val*minmaxflag;
  return(val);
}

double dbrent(double am, double bm, double cm, int minmaxflag, double tol, double *mmin, int *iter, int Nsubharm, int Nharm, double **subharmA, double **subharmB, double **harmA, double **harmB, double *fundA, double *fundB, double period, int p_per)
{
  /* Uses dbrent to minimize or maximize the function */

  int ok1, ok2;
  double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
  double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

  a=(am < cm ? am : cm);
  b=(am > cm ? am : cm);
  x=w=v=bm;

  fw=fv=fx=f_tomin(x,minmaxflag, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, period, p_per);
  dw=dv=dx=dfdx_tomin(x,minmaxflag, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, period, p_per);
  for ((*iter) = 1;(*iter)<=ITMAX;(*iter)++)
    {
      xm=0.5*(a+b);
      tol1=tol*fabs(x)+ZEPS;
      tol2=2.0*tol1;
      if(fabs(x-xm) <= (tol2-0.5*(b-a))) {
	*mmin = x;
	return fx;
      }
      if (fabs(e) > tol1) {
	d1=2.0*(b-a);
	d2=d1;
	if (dw != dx) d1=(w-x)*dx/(dx-dw);
	if (dv != dx) d2=(v-x)*dx/(dx-dv);
	u1=x+d1;
	u2=x+d2;
	ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
	ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
	olde=e;
	e=d;
	if (ok1 || ok2) {
	  if (ok1 && ok2)
	    d=(fabs(d1) < fabs(d2) ? d1 : d2);
	  else if(ok1)
	    d=d1;
	  else
	    d=d2;
	  if (fabs(d) <= fabs(0.5*olde)) {
	    u=x+d;
	    if (u-a < tol2 || b-u < tol2)
	      d=SIGN(tol1,xm-x);
	  } else {
	    d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
	  }
	} else {
	  d = 0.5*(e=(dx >= 0.0 ? a-x : b-x));
	}
      } else {
	d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
      }
      if (fabs(d) >= tol1) {
	u=x+d;
	fu=f_tomin(u,minmaxflag, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, period, p_per);
      } else {
	u=x+SIGN(tol1,d);
	fu=f_tomin(u,minmaxflag, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, period, p_per);
	if (fu > fx) {
	  *mmin = x;
	  return fx;
	}
      }
      du=dfdx_tomin(u,minmaxflag, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, period, p_per);
      if (fu <= fx) {
	if (u >= x) a=x; else b=x;
	MOV3(v,fv,dv, w,fw,dw);
	MOV3(w,fw,dw, x,fx,dx);
	MOV3(x,fx,dx, u,fu,du);
      } else {
	if (u < x) a=u; else b=u;
	if (fu <= fw || w == x) {
	  MOV3(v,fv,dv, w,fw,dw);
	  MOV3(w,fw,dw, u,fu,du);
	} else if(fu < fv || v == x || v == w) {
	  MOV3(v,fv,dv, u,fu,du);
	}
      }
    }
  *mmin = -999.;
  return -1.;
}

void dokillharms(int N, double *t, double *mag, double *sig, int Nper, double *periods, int Nsubharm, int Nharm, double **subharmA, double **subharmB, double **harmA, double **harmB, double *fundA, double *fundB, double *mean, int omodel, char *modelname, double *amp, int fitonly, int outtype, double clip)
{
  int i, j, k, ncomp, p, u, l, ngood, test, iter;
  double st, st_, d1, d2;
  double ct, ct_;
  double freq, sum1, sum2, mag_ave, modelterm, rms_resid;
  int freqnum;
  double *magmodel = NULL;
  double **Amatrix, **DesignMatrix, *Avector, *Bvector, *w, **v, wmax, thresh;
  FILE *outfile = NULL;
  double am, bm, cm, step, minval, maxval, x, try, dum;
  int mini, maxi, nstep;
  int *isclip = NULL;
  int Nclip = 0;
  ncomp = 1;
  test = 0;
  for(i=0; i<Nper;i++)
    {
      if(periods[i] <= 0.)
	{
	  test = 1;
	}
      }
  if(test)
    {
      *mean = 9999999.;
      for(i=0;i<Nper;i++)
	{
	  fundA[i] = -1.;
	  fundB[i] = -1.;
	  for(j=0;j<Nharm;j++)
	    {
	      harmA[i][j] = -1.;
	      harmB[i][j] = -1.;
	    }
	  for(j=0;j<Nsubharm;j++)
	    {
	      subharmA[i][j] = -1.;
	      subharmB[i][j] = -1.;
	    }
	}
      return;
    }



  for(i = 0; i<Nper;i++)
    {
      ncomp += 2 + (2*Nharm + 2*Nsubharm);
    }
  if((Amatrix = (double **) malloc(N * sizeof(double *))) == NULL ||
     (DesignMatrix = (double **) malloc(N * sizeof(double *))) == NULL ||
     (Avector = (double *) malloc(ncomp * sizeof(double))) == NULL ||
     (Bvector = (double *) malloc(N * sizeof(double))) == NULL ||
     (w = (double *) malloc(ncomp * sizeof(double))) == NULL ||
     (v = (double **) malloc(ncomp * sizeof(double *))) == NULL ||
     (isclip = (int *) malloc(N * sizeof(int))) == NULL ||
     (magmodel = (double *) malloc(N * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0;i<N;i++) {
    isclip[i] = 0;
    if((Amatrix[i] = (double *) malloc(ncomp * sizeof(double))) == NULL ||
       (DesignMatrix[i] = (double *) malloc(ncomp * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  }
  for(i=0;i<ncomp;i++)
    if((v[i] = (double *) malloc(ncomp * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  /*  if((a = (long double **) malloc(ncomp * sizeof(long double *))) == NULL ||
     (indx = (int *) malloc(ncomp * sizeof(int))) == NULL ||
     (b = (long double *) malloc(ncomp * sizeof(long double))) == NULL)
    error(ERR_MEMALLOC);

  for(i=0;i<ncomp;i++)
    if((a[i] = (long double *) malloc(ncomp * sizeof(long double))) == NULL)
      error(ERR_MEMALLOC);

  for(i=0;i<ncomp;i++)
    {
      b[i] = 0.;
      for(j=0;j<ncomp;j++)
	a[i][j] = 0.;
    }
  */
  /* Prepare the design matrix */
  do {
    for(i=0, l=0;i<N;i++)
      {
	if(!isnan(mag[i]) && sig[i] > 0. && !isclip[i])
	  {
	    Bvector[l] = mag[i] / sig[i];
	    l++;
	  }
      }
    ngood = l;

    for(i=0,l=0;i<N;i++)
      {
	if(!isnan(mag[i]) && sig[i] > 0. && !isclip[i])
	  {
	    DesignMatrix[l][0] = 1.;
	    Amatrix[l][0] = DesignMatrix[l][0]/sig[i];
	    u = 1;
	    for(p=0;p<Nper;p++)
	      {
		freq = TWOPI / periods[p];
		for(k=0;k<Nsubharm;k++)
		  {
		    freqnum = (Nsubharm - k) + 1;
		    st = sin(t[i] * (freq / (double) freqnum));
		    ct = cos(t[i] * (freq / (double) freqnum));
		    DesignMatrix[l][u] = st;
		    Amatrix[l][u] = st/sig[i];
		    //mag[i] -= st * (double) b[u];
		    u++;
		    DesignMatrix[l][u] = ct;
		    Amatrix[l][u] = ct/sig[i];
		    //mag[i] -= ct * (double) b[u];
		    u++;
		  }
		st = sin(t[i] * freq);
		ct = cos(t[i] * freq);
		st_ = st;
		ct_ = ct;
		DesignMatrix[l][u] = st_;
		Amatrix[l][u] = st_/sig[i];
		//mag[i] -= st_ * (double) b[u];
		u++;
		DesignMatrix[l][u] = ct_;
		Amatrix[l][u] = ct_/sig[i];
		//mag[i] -= ct_ * (double) b[u];
		u++;
		for(k=0;k<Nharm;k++)
		  {
		    d1 = ct_*st + st_*ct;
		    d2 = ct_*ct - st_*st;
		    st_ = d1;
		    ct_ = d2;
		    DesignMatrix[l][u] = st_;
		    Amatrix[l][u] = st_/sig[i];
		    //mag[i] -= st_ * (double) b[u];
		    u++;
		    DesignMatrix[l][u] = ct_;
		    Amatrix[l][u] = ct_/sig[i];
		    //mag[i] -= ct_ * (double) b[u];
		    u++;
		  }
	      }
	    l++;
	  }
      }

    /* Solve the matrix via Singular Value Decomposition */
    svdcmp(Amatrix,ngood,ncomp,w,v);
    wmax = 0.;
    for(j=0;j<ncomp;j++)
      if(w[j] > wmax) wmax = w[j];
    thresh = TOLERANCE * wmax;
    for(j=0;j<ncomp;j++)
      if(w[j] < thresh) w[j] = 0.0;
    svbksb(Amatrix,w,v,ngood,ncomp,Bvector,Avector);

    /* Correct the light curves and store the fitting parameters in the input matrices */

    // First find the average
    sum1 = 0.;
    sum2 = 0.;
    for(i=0;i<N;i++)
      if(!isnan(mag[i]) && sig[i] > 0. && !isclip[i])
	{
	  sum1 += mag[i] / (sig[i]*sig[i]);
	  sum2 += 1. / (sig[i]*sig[i]);
	}
    mag_ave = (sum1 / sum2);
    *mean = Avector[0];

    sum1 = 0.;
    sum2 = 0.;

    for(i=0,l=0;i<N;i++)
      if(!isnan(mag[i]) && sig[i] > 0.)
	{
	  if(!isclip[i]) {
	    magmodel[i] = (double) Avector[0];
	    //	mag[i] -= (double) b[0];
	    //	mag[i] += mag_ave;
	    u = 1;
	    for(p=0;p<Nper;p++)
	      {
		freq = TWOPI / periods[p];
		for(k=0;k<Nsubharm;k++)
		  {
		    //freqnum = (Nsubharm - k) + 1;
		    //st = sin(t[i] * (freq / (double) freqnum));
		    //ct = cos(t[i] * (freq / (double) freqnum));
		    //modelterm += st * (double) b[u];
		    magmodel[i] += DesignMatrix[l][u] * Avector[u];
		    //mag[i] -= st * (double) b[u];
		    subharmA[p][k] = Avector[u];
		    u++;
		    magmodel[i] += DesignMatrix[l][u] * Avector[u];
		    //modelterm += ct * (double) b[u];
		    //mag[i] -= ct * (double) b[u];
		    subharmB[p][k] = Avector[u];
		    u++;
		  }
		//st = sin(t[i] * freq);
		//ct = cos(t[i] * freq);
		//st_ = st;
		//ct_ = ct;
		magmodel[i] += DesignMatrix[l][u] * Avector[u];
		//modelterm += st_ * (double) b[u];
		//mag[i] -= st_ * (double) b[u];
		fundA[p] = Avector[u];
		u++;
		magmodel[i] += DesignMatrix[l][u] * Avector[u];
		//modelterm += ct_ * (double) b[u];
		//mag[i] -= ct_ * (double) b[u];
		fundB[p] = Avector[u];
		u++;
		for(k=0;k<Nharm;k++)
		  {
		    d1 = ct_*st + st_*ct;
		    d2 = ct_*ct - st_*st;
		    st_ = d1;
		    ct_ = d2;
		    magmodel[i] += DesignMatrix[l][u] * Avector[u];
		    //modelterm += st_ * (double) b[u];
		    //mag[i] -= st_ * (double) b[u];
		    harmA[p][k] = Avector[u];
		    u++;
		    magmodel[i] += DesignMatrix[l][u] * Avector[u];
		    //modelterm += ct_ * (double) b[u];
		    //mag[i] -= ct_ * (double) b[u];
		    harmB[p][k] = Avector[u];
		    u++;
		  }
	      }
	    sum1 += (mag[i] - magmodel[i]);
	    sum2 += (mag[i] - magmodel[i])*(mag[i] - magmodel[i]);
	    l++;
	  }
	  else {
	    magmodel[i] = (double) Avector[0];
	    //	mag[i] -= (double) b[0];
	    //	mag[i] += mag_ave;
	    u = 1;
	    for(p=0;p<Nper;p++)
	      {
		freq = TWOPI / periods[p];
		for(k=0;k<Nsubharm;k++)
		  {
		    freqnum = (Nsubharm - k) + 1;
		    st = sin(t[i] * (freq / (double) freqnum));
		    ct = cos(t[i] * (freq / (double) freqnum));
		    magmodel[i] += st * Avector[u];
		    subharmA[p][k] = Avector[u];
		    u++;
		    magmodel[i] += ct * Avector[u];
		    subharmB[p][k] = Avector[u];
		    u++;
		  }
		st = sin(t[i] * freq);
		ct = cos(t[i] * freq);
		st_ = st;
		ct_ = ct;
		magmodel[i] += st_ * Avector[u];
		fundA[p] = Avector[u];
		u++;
		magmodel[i] += ct_ * Avector[u];
		fundB[p] = Avector[u];
		u++;
		for(k=0;k<Nharm;k++)
		  {
		    d1 = ct_*st + st_*ct;
		    d2 = ct_*ct - st_*st;
		    st_ = d1;
		    ct_ = d2;
		    magmodel[i] += st_ * Avector[u];
		    harmA[p][k] = Avector[u];
		    u++;
		    magmodel[i] += ct_ * Avector[u];
		    harmB[p][k] = Avector[u];
		    u++;
		  }
	      }
	  }
	}
    if(clip > 0 && Nclip == 0) {
      /* Find the points to clip */
      rms_resid = sqrt((sum2 / (double) l) - (sum1 * sum1 ) / ((double) l*l));
      Nclip = 0;
      for(i=0; i < N; i++) {
	if(!isnan(mag[i]) && sig[i] > 0. && !isclip[i])
	  {
	    if(fabs(mag[i] - magmodel[i]) > clip*rms_resid) {
	      isclip[i] = 1;
	      Nclip++;
	    }
	  }
      }
    }
    else Nclip = 0;
  } while (Nclip > 0); /* Redo the fit if points have been clipped */

  if(omodel)
    {
      if((outfile = fopen(modelname,"w")) == NULL)
	error2(ERR_CANNOTWRITE,modelname);
    }
  if(omodel || !fitonly) {
    for(i=0; i < N; i++) {
      if(!isnan(mag[i]) && sig[i] > 0.) {
	if(omodel)
	  fprintf(outfile,"%f %f %f %f\n",t[i],mag[i],magmodel[i],sig[i]);
	if(!fitonly)
	  mag[i] = mag[i] - magmodel[i] + mag_ave;
      }
    }
  }
  if(omodel)
    fclose(outfile);

  /* Get the amplitudes of the model fourier series */
  nstep = 100*(Nharm + Nsubharm + 1);
  step = 1./(double) nstep;
  for(p=0;p<Nper;p++)
    {
      mini = maxi = 0;
      minval = f_tomin(0., 1, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, periods[p], p);
      maxval = minval;
      x = step;
      for(j=1;j<nstep;j++, x += step)
	{
	  try = f_tomin(x, 1, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, periods[p], p);
	  if(try < minval)
	    {
	      minval = try;
	      mini = j;
	    }
	  if(try > maxval)
	    {
	      maxval = try;
	      maxi = j;
	    }
	}
      bm = mini*step;
      am = bm - step;
      cm = bm + step;
      minval = dbrent(am, bm, cm, 1, 1.E-5, &dum, &iter, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, periods[p], p);
      bm = maxi*step;
      am = bm - step;
      cm = bm + step;
      maxval = dbrent(am, bm, cm, -1, 1.E-5, &dum, &iter, Nsubharm, Nharm, subharmA, subharmB, harmA, harmB, fundA, fundB, periods[p], p);
      maxval = -maxval;
      amp[p] = maxval - minval;
    }
  /* Put the output terms into the right format*/
  if(outtype != KILLHARM_OUTTYPE_DEFAULT)
    {
      for(i=0;i<Nper;i++)
	{
	  d1 = sqrt(fundA[i]*fundA[i] + fundB[i]*fundB[i]);
	  switch(outtype)
	    {
	    case KILLHARM_OUTTYPE_AMPPHASE:
	      d2 = atan2(-fundA[i],fundB[i])/2.0/M_PI;
	      break;
	    case KILLHARM_OUTTYPE_AMPRADPHASE:
	      d2 = atan2(-fundA[i],fundB[i]);
	      break;
	    case KILLHARM_OUTTYPE_RPHI:
	      d2 = atan2(-fundA[i],fundB[i])/2.0/M_PI;
	      break;
	    case KILLHARM_OUTTYPE_RRADPHI:
	      d2 = atan2(-fundA[i],fundB[i]);
	      break;
	    }
	  fundA[i] = d1;
	  fundB[i] = d2;
	  for(j=0; j<Nsubharm;j++)
	    {
	      if(outtype == KILLHARM_OUTTYPE_AMPPHASE || outtype == KILLHARM_OUTTYPE_AMPRADPHASE)
		{
		  d1 = sqrt(subharmA[i][j]*subharmA[i][j] + subharmB[i][j]*subharmB[i][j]);
		  if(outtype == KILLHARM_OUTTYPE_AMPPHASE)
		    d2 = atan2(-subharmA[i][j],subharmB[i][j])/2.0/M_PI;
		  else
		    d2 = atan2(-subharmA[i][j],subharmB[i][j]);
		}
	      else
		{
		  d1 = sqrt(subharmA[i][j]*subharmA[i][j] + subharmB[i][j]*subharmB[i][j])/fundA[i];
		  if(outtype == KILLHARM_OUTTYPE_RPHI)
		    d2 = atan2(-subharmA[i][j],subharmB[i][j])/2.0/M_PI - fundB[i]/((double) (j+2));
		  else
		    d2 = atan2(-subharmA[i][j],subharmB[i][j]) - fundB[i]/((double) (j+2));
		}
	      subharmA[i][j] = d1;
	      subharmB[i][j] = d2;
	    }
	  for(j=0; j<Nharm;j++)
	    {
	      if(outtype == KILLHARM_OUTTYPE_AMPPHASE || outtype == KILLHARM_OUTTYPE_AMPRADPHASE)
		{
		  d1 = sqrt(harmA[i][j]*harmA[i][j] + harmB[i][j]*harmB[i][j]);
		  if(outtype == KILLHARM_OUTTYPE_AMPPHASE)
		    d2 = atan2(-harmA[i][j],harmB[i][j])/2.0/M_PI;
		  else
		    d2 = atan2(-harmA[i][j],harmB[i][j]);
		}
	      else
		{
		  d1 = sqrt(harmA[i][j]*harmA[i][j] + harmB[i][j]*harmB[i][j])/fundA[i];
		  if(outtype == KILLHARM_OUTTYPE_RPHI)
		    d2 = atan2(-harmA[i][j],harmB[i][j])/2.0/M_PI - fundB[i]*((double) (j+2));
		  else
		    d2 = atan2(-harmA[i][j],harmB[i][j]) - fundB[i]*((double) (j+2));
		}
	      harmA[i][j] = d1;
	      harmB[i][j] = d2;
	    }
	}
    }

  /*for(i=0;i<ncomp;i++)
    free(a[i]);
  free(a);
  free(b);
  free(indx);*/
  free(Avector);
  free(Bvector);
  free(w);
  for(i=0;i<ncomp;i++)
    free(v[i]);
  free(v);
  for(i=0;i<N;i++)
    {
      free(Amatrix[i]);
      free(DesignMatrix[i]);
    }
  free(DesignMatrix);
  free(Amatrix);
  free(isclip);
  free(magmodel);

}


