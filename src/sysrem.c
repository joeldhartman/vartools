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

/* This file contains functions to perform SYSREM on an ensemble of light curves, the sysrem algorithm is described in Tamuz, Mazeh & Zucker (2005, MNRAS, 356, 1466). */

#define TOLERANCE 1.0E-9

#define MOON_FREQ 0.212768702
#define MINMAGCHANGE 0.0001
#define MAXITS 1000
#define SIGMIN 0.003
#define INDEFVAL -9999

#ifdef isinf
#else
int isinf(double d)
{
  if(d == HUGE_VAL || d == -HUGE_VAL)
    return 1;
  else
    return 0;
}
#endif

double itersolve_chisqr_sysrem(double *star_coeff, double *image_coeff, int numlc, int numJD, double **mag, double **sig, double *mag_ave, double sigma_clip, double upperdiff, double saturation, int **clip)
{
  int j,k,n;
  long double S,dum;
  S = 0.0;
  n = 0;
  for(j=0;j<numlc;j++)
    for(k=0;k<numJD;k++)
      if(clip[j][k])
	{
	  dum = (long double) (mag[j][k] - mag_ave[j] - star_coeff[j]*image_coeff[k])/sig[j][k];
	  S += dum*dum;
	  n++;
	}
  return ((double) (S/ (double) n));
}

double chisqr_sysrem(double p[], int numlc, int numJD, double **mag, double **sig, double *mag_ave, double sigma_clip, double saturation)
{
  int j,k,n;
  long double S,dum;
  S = 0.0;
  n = 0;
  for(j=0;j<numlc;j++)
    for(k=0;k<numJD;k++)
      if(!isnan(mag[j][k]) && !isinf(mag[j][k]) && mag[j][k] > saturation && (ABS_(mag[j][k] - mag_ave[j])/sig[j][k] < sigma_clip || sigma_clip < 0))
	{
	  dum = (long double) (mag[j][k] - mag_ave[j] - p[j+1]*p[numlc+k+1])/sig[j][k];
	  S += dum*dum;
	}
  return((double) (S/(double) (numJD*numlc)));
}

void dchisqr_sysrem(double p[], double xi[], int numlc, int numJD, double **mag, double **sig, double *mag_ave, double sigma_clip, double saturation)
{
  int j,k;
  int n, num;
  long double S;
  num = 0;
  for(n=1;n<=numlc+numJD;n++)
    {
      S = 0;
      if(n <= numlc)
	{
	  j = n-1;
	  for(k=0;k<numJD;k++)
	    {
	      if(!isnan(mag[j][k]) && !isinf(mag[j][k]) && mag[j][k] > saturation && (ABS_(mag[j][k] - mag_ave[j])/sig[j][k] < sigma_clip || sigma_clip < 0))
		S += -2.0*(p[numlc+k+1]*(mag[j][k] - mag_ave[j] - p[n]*p[numlc+k+1])/(sig[j][k]*sig[j][k]));
	    }
	}
      else
	{
	  k = n - numlc - 1;
	  for(j=0;j<numlc;j++)
	    {
	      if(!isnan(mag[j][k]) && !isinf(mag[j][k]) && mag[j][k] > saturation && (ABS_(mag[j][k] - mag_ave[j])/sig[j][k] < sigma_clip || sigma_clip < 0))
		S += -2.0*(p[j+1]*(mag[j][k] - mag_ave[j] - p[j+1]*p[n])/(sig[j][k]*sig[j][k]));
	    }
	}
      xi[n] = (double) (S / ((double) numJD*numlc));
    }
}

int ncom=0;	/* defining declarations */
double *pcom=0,*xicom=0;

double f1dim_sysrem(x,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation)
double x, **mag, **sig, *mag_ave, sigma_clip,saturation;
int numlc, numJD;
{
	int j;
	double f,*xt;

	xt=(double *) malloc(ncom*sizeof(double))-1;
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	f=chisqr_sysrem(xt,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
	free(xt+1);
	return f;
}

double df1dim_sysrem(x,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation)
double x,**mag,**sig,*mag_ave,sigma_clip,saturation;
int numlc,numJD;
{
	int j;
	double df1=0.0;
	double *xt,*df;

	xt=(double *) malloc(ncom*sizeof(double))-1;
	df=(double *) malloc(ncom*sizeof(double))-1;
	for (j=1;j<=ncom;j++) xt[j]=pcom[j]+x*xicom[j];
	dchisqr_sysrem(xt,df,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
	for (j=1;j<=ncom;j++) df1 += df[j]*xicom[j];
	free(df+1);
	free(xt+1);
	return df1;
}

#define GOLD 1.618034
#define GLIMIT 100.0
#define TINY 1.0e-20
#define MAX(a,b) ((a) > (b) ? (a) : (b))
#define SHFT(a,b,c,d) (a)=(b);(b)=(c);(c)=(d);

void mnbrak_sysrem(ax,bx,cx,fa,fb,fc,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation)
double *ax,*bx,*cx,*fa,*fb,*fc,**mag,**sig,*mag_ave,sigma_clip,saturation;
int numlc,numJD;
{
	double ulim,u,r,q,fu,dum;

	*fa=f1dim_sysrem(*ax,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
	*fb=f1dim_sysrem(*bx,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
	if (*fb > *fa) {
		SHFT(dum,*ax,*bx,dum)
		SHFT(dum,*fb,*fa,dum)
	}
	*cx=(*bx)+GOLD*(*bx-*ax);
	*fc=f1dim_sysrem(*cx,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
	while (*fb > *fc) {
		r=(*bx-*ax)*(*fb-*fc);
		q=(*bx-*cx)*(*fb-*fa);
		u=(*bx)-((*bx-*cx)*q-(*bx-*ax)*r)/
			(2.0*SIGN(MAX(ABS_(q-r),TINY),q-r));
		ulim=(*bx)+GLIMIT*(*cx-*bx);
		if ((*bx-u)*(u-*cx) > 0.0) {
			fu=f1dim_sysrem(u,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
			if (fu < *fc) {
				*ax=(*bx);
				*bx=u;
				*fa=(*fb);
				*fb=fu;
				return;
			} else if (fu > *fb) {
				*cx=u;
				*fc=fu;
				return;
			}
			u=(*cx)+GOLD*(*cx-*bx);
			fu=f1dim_sysrem(u,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
		} else if ((*cx-u)*(u-ulim) > 0.0) {
			fu=f1dim_sysrem(u,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
			if (fu < *fc) {
				SHFT(*bx,*cx,u,*cx+GOLD*(*cx-*bx))
				SHFT(*fb,*fc,fu,f1dim_sysrem(u,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation))
			}
		} else if ((u-ulim)*(ulim-*cx) >= 0.0) {
			u=ulim;
			fu=f1dim_sysrem(u,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
		} else {
			u=(*cx)+GOLD*(*cx-*bx);
			fu=f1dim_sysrem(u,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
		}
		SHFT(*ax,*bx,*cx,u)
		SHFT(*fa,*fb,*fc,fu)
	}
}

#undef GOLD
#undef GLIMIT
#undef TINY
#undef MAX
#undef SIGN
#undef SHFT


#define ITMAX 100
#define ZEPS 1.0e-10
#define SIGN(a,b) ((b) > 0.0 ? ABS_(a) : -ABS_(a))
#define MOV3(a,b,c, d,e,f) (a)=(d);(b)=(e);(c)=(f);

double dbrent_sysrem(ax,bx,cx,tol,xmin,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation)
     double ax,bx,cx,tol,*xmin,**mag,**sig,*mag_ave,sigma_clip,saturation;
int numlc,numJD;
{
	int iter,ok1,ok2;
	double a,b,d,d1,d2,du,dv,dw,dx,e=0.0;
	double fu,fv,fw,fx,olde,tol1,tol2,u,u1,u2,v,w,x,xm;

	a=(ax < cx ? ax : cx);
	b=(ax > cx ? ax : cx);
	x=w=v=bx;
	fw=fv=fx=f1dim_sysrem(x,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
	dw=dv=dx=df1dim_sysrem(x,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
	for (iter=1;iter<=ITMAX;iter++) {
		xm=0.5*(a+b);
		tol1=tol*ABS_(x)+ZEPS;
		tol2=2.0*tol1;
		if (ABS_(x-xm) <= (tol2-0.5*(b-a))) {
			*xmin=x;
			return fx;
		}
		if (ABS_(e) > tol1) {
			d1=2.0*(b-a);
			d2=d1;
			if (dw != dx)  d1=(w-x)*dx/(dx-dw);
			if (dv != dx)  d2=(v-x)*dx/(dx-dv);
			u1=x+d1;
			u2=x+d2;
			ok1 = (a-u1)*(u1-b) > 0.0 && dx*d1 <= 0.0;
			ok2 = (a-u2)*(u2-b) > 0.0 && dx*d2 <= 0.0;
			olde=e;
			e=d;
			if (ok1 || ok2) {
				if (ok1 && ok2)
					d=(ABS_(d1) < ABS_(d2) ? d1 : d2);
				else if (ok1)
					d=d1;
				else
					d=d2;
				if (ABS_(d) <= ABS_(0.5*olde)) {
					u=x+d;
					if (u-a < tol2 || b-u < tol2)
						d=SIGN(tol1,xm-x);
				} else {
					d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
				}
			} else {
				d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
			}
		} else {
			d=0.5*(e=(dx >= 0.0 ? a-x : b-x));
		}
		if (ABS_(d) >= tol1) {
			u=x+d;
			fu=f1dim_sysrem(u,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
		} else {
			u=x+SIGN(tol1,d);
			fu=f1dim_sysrem(u,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
			if (fu > fx) {
				*xmin=x;
				return fx;
			}
		}
		du=df1dim_sysrem(u,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
		if (fu <= fx) {
			if (u >= x) a=x; else b=x;
			MOV3(v,fv,dv, w,fw,dw)
			MOV3(w,fw,dw, x,fx,dx)
			MOV3(x,fx,dx, u,fu,du)
		} else {
			if (u < x) a=u; else b=u;
			if (fu <= fw || w == x) {
				MOV3(v,fv,dv, w,fw,dw)
				MOV3(w,fw,dw, u,fu,du)
			} else if (fu < fv || v == x || v == w) {
				MOV3(v,fv,dv, u,fu,du)
			}
		}
	}
	fprintf(stderr,"Too many iterations in routine DBRENT\n");
	return -1.;
}

#undef ITMAX
#undef ZEPS
#undef SIGN
#undef MOV3


#define TOL 2.0e-4

void dlinmin_sysrem(p,xi,fret,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation)
double p[],xi[],*fret, **mag, **sig, *mag_ave, sigma_clip,saturation;
int numlc,numJD;
{
	int j,n;
	double xx,xmin,fx,fb,fa,bx,ax;
      	void mnbrak_sysrem();

	n = numlc+numJD;
	ncom=n;
	pcom=(double *) malloc(n*sizeof(double))-1;
	xicom=(double *) malloc(n*sizeof(double))-1;
	for (j=1;j<=n;j++) {
		pcom[j]=p[j];
		xicom[j]=xi[j];
	}
	ax=0.0;
	xx=1.0;
	bx=2.0;
	mnbrak_sysrem(&ax,&xx,&bx,&fa,&fx,&fb,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
	*fret=dbrent_sysrem(ax,xx,bx,TOL,&xmin,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
	for (j=1;j<=n;j++) {
		xi[j] *= xmin;
		p[j] += xi[j];
	}
	free(xicom+1);
	free(pcom+1);
}

#undef TOL

#define ITMAX 200
#define EPS 1.0e-10
#define FREEALL free(xi+1);free(h+1);free(g+1);

void frprmn_sysrem(p,ftol,iter,fret,numlc,numJD, mag,sig,mag_ave,sigma_clip,saturation)
     double p[],ftol,*fret, **mag, **sig, *mag_ave, sigma_clip,saturation;
     int *iter,numlc,numJD;
{
        int n;
	int j,its;
	double gg,gam,fp,dgg;
	double *g,*h,*xi;

	n = numlc+numJD;
	g= (double *) malloc(n*sizeof(double)) - 1;
	h= (double *) malloc(n*sizeof(double)) - 1;
	xi= (double *) malloc(n*sizeof(double)) - 1;
	fp=chisqr_sysrem(p,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
	dchisqr_sysrem(p,xi,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
	for (j=1;j<=n;j++) {
		g[j] = -xi[j];
		xi[j]=h[j]=g[j];
	}
	for (its=1;its<=ITMAX;its++) {
		*iter=its;
		dlinmin_sysrem(p,xi,fret,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
		if (2.0*ABS_(*fret-fp) <= ftol*(ABS_(*fret)+ABS_(fp)+EPS)) {
			FREEALL
			return;
		}
		fp=chisqr_sysrem(p,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
		dchisqr_sysrem(p,xi,numlc,numJD,mag,sig,mag_ave,sigma_clip,saturation);
		dgg=gg=0.0;
		for (j=1;j<=n;j++) {
			gg += g[j]*g[j];
/*		  dgg += xi[j]*xi[j];	*/
			dgg += (xi[j]+g[j])*xi[j];
		}
		if (gg == 0.0) {
			FREEALL
			return;
		}
		gam=dgg/gg;
		for (j=1;j<=n;j++) {
			g[j] = -xi[j];
			xi[j]=h[j]=g[j]+gam*h[j];
		}
	}
	fprintf(stderr,"Too many iterations in FRPRMN\n");
}

#undef ITMAX
#undef EPS
#undef FREEALL


void nrerror(error_text)
     char *error_text;
{
  void exit();
  fprintf(stderr," Run error....");
  fprintf(stderr,"%s\n",error_text);
  fprintf(stderr,"Goodbye ! \n");
}


#define TOL 1.0e-7

#define ALN2I 1.442695022
#define TINY 1.0e-5

void shell_sysrem(n,arr)
float arr[];
int n;
{
	int nn,m,j,i,lognb2;
	float t;

	lognb2=(log((double) n)*ALN2I+TINY);
	m=n;
	for (nn=1;nn<=lognb2;nn++) {
		m >>= 1;
		for (j=m+1;j<=n;j++) {
			i=j-m;
			t=arr[j];
			while (i >= 1 && arr[i] > t) {
				arr[i+m]=arr[i];
				i -= m;
			}
			arr[i+m]=t;
		}
	}
}

#undef ALN2I
#undef TINY


double stddev_sysrem(int n, double *data)
{
  int i,num;
  long double var1 = 0;
  long double var2 = 0;
  num = 0;
  for(i=0;i<n;i++){
    if(!isnan(data[i]) && !isinf(data[i])) num++;
  }
  for(i=0;i<n;i++){
    if(!isnan(data[i]) && !isinf(data[i]))
      {
	var1 += (long double) (data[i]/ (double) num);
	var2 += (long double) (data[i]*data[i]) / (double) num;
      }
  }
  return(sqrt((double)((var2) - ((var1*var1)))));
}

void mysort_sysrem(int N, double* data1, double* data2, double* data3)
{
  int i, j;
  double v, t1, t2, t3;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(data1[++i] < v && i < N) { }
    while(data1[--j] > v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
    t3 = data3[i]; data3[i] = data3[j]; data3[j] = t3;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  t3 = data3[i-1]; data3[i-1] = data3[0]; data3[0] = t3;
  mysort_sysrem(i-1,data1,data2, data3);
  mysort_sysrem(N-i,data1+i,data2+i, data3+i);
}

double median_sysrem(int n, double *mag_bin){
  int i;
  float *temp;
  double temp2;
  if(n == 1)
    return(mag_bin[0]);
  temp = (float *) malloc(2*(n+1)*sizeof(float));
  for(i=0;i<n;i++)
    temp[i+1] = (float) mag_bin[i];
    shell_sysrem(n, temp);
  if(n % 2 == 0){
    temp2 = (double) (temp[n/2]+temp[(n/2) + 1])/2.0;
    free(temp);
    return(temp2);
  }
  else {
    temp2 = (double) temp[(n/2) + 1];
    free(temp);
    return(temp2);
  }
}

void get_magave_std_sysrem(double **mag, double **sig, double *mag_ave, int numLC, int numJD, double sigma_clip1, int useweights)
{
  int j, k,num;
  double tempmean,std_val, stdval2;
  for(j=0;j<numLC;j++)
    {
      mag_ave[j] = 0.0;
      tempmean = 0.0;
      std_val = stddev_sysrem(numJD,mag[j]);
      num = 0;
      for(k=0;k<numJD;k++)
	if(!isnan(mag[j][k]) && !isinf(mag[j][k]))
	  num++;
      for(k=0;k<numJD;k++)
	if(!isnan(mag[j][k]) && !isinf(mag[j][k]))
	  tempmean += mag[j][k] / ((double) num);
      num = 0;
      for(k=0;k<numJD;k++)
	if(!isnan(mag[j][k]) && !isinf(mag[j][k]) && (ABS_(mag[j][k]-tempmean) < sigma_clip1*std_val || sigma_clip1 < 0))
	  num++;
      stdval2 = 0.;
      for(k=0;k<numJD;k++)
	if(!isnan(mag[j][k]) && !isinf(mag[j][k]) && (ABS_(mag[j][k]-tempmean) < sigma_clip1*std_val || sigma_clip1 < 0))
	  {
	    mag_ave[j] += mag[j][k]/((double) num);
	    stdval2 += mag[j][k]*mag[j][k]/((double) num);
	  }
      stdval2 = sqrt(stdval2 - (mag_ave[j]*mag_ave[j]));
      if(!useweights)
	{
	  for(k=0;k<numJD;k++)
	    sig[j][k] = 1.;
	}
      else if(useweights == 1)
	{
	  for(k=0;k<numJD;k++)
	    sig[j][k] = stdval2;
	}
    }
}

void decorr_star_sysrem(double *star_coeff,double **mag,double **sig,double *mag_ave,int numlc,int numJD,int **clip, double sigma_clip2,double *image_coeff,double saturation)
{
  int j, k, n;
  long double val1, val2;
  for(j=0;j<numlc;j++)
    {
      n = 0;
      val1 = 0.0;
      val2 = 0.0;
      for(k=0;k<numJD;k++)
	{
	  if(clip[j][k])
	    {
	      val1 += (long double) (mag[j][k]-mag_ave[j])*image_coeff[k]/(sig[j][k]*sig[j][k]);
	      val2 += (long double) (image_coeff[k]*image_coeff[k])/(sig[j][k]*sig[j][k]);
	      n++;
	    }
	}
      if(n > 0 && val2 > 0.0)
	star_coeff[j] = (double) val1/val2;
      else
	star_coeff[j] = 0.0;
    }
}

void decorr_image_sysrem(double *star_coeff,double **mag,double **sig,double *mag_ave,int numlc,int numJD,int **clip,double sigma_clip2,double upperdiff,double *image_coeff,double saturation)
{
  int j, k;
  long double av_val1, av_val2, val1, val2;

  int n;
  for(k=0; k<numJD; k++)
    {
      n = 0;
      val1 = 0.0;
      val2 = 0.0;
      av_val1 = 0.0;
      av_val2 = 0.0;
      n = 0;
      for(j=0;j<numlc;j++)
	{
	  if(clip[j][k])
	    {
	      val1 += (long double) (mag[j][k] - mag_ave[j])*star_coeff[j]/(sig[j][k]*sig[j][k]);
	      val2 += (long double) (star_coeff[j]*star_coeff[j])/(sig[j][k]*sig[j][k]);
	      n++;
	    }
	}
      if(n > 0 && val2 > 0.0)
	image_coeff[k] = (double) val1/val2;
      else
	image_coeff[k] = 0.0;
    }
}

void magcorr_sysrem(double *star_coeff,double *image_coeff,double **mag,double **sig,int numlc,int numJD)
{
  int j, k;
  for(j=0;j<numlc;j++)
    for(k=0;k<numJD;k++)
      if(!isnan(mag[j][k]) && !isinf(mag[j][k]))
	{
	  mag[j][k] -= star_coeff[j]*image_coeff[k];
	}
}

int iter_solve_sysrem(int starttype, int maxits, double convergeval, double *star_coeff, double* image_coeff,double **mag,double **sig,double *mag_ave,int numlc,int numJD,int **clip,double sigma_clip2, double upperdiff, double saturation)
{
  double converged = 0, image_coeff_stddev, image_coeff_ave;
  int iternum = 0, ngood;
  double oldchisqr,newchisqr;
  double *old_star_coeff, *old_image_coeff;
  long double var1, var2;
  int j,k;
  oldchisqr = itersolve_chisqr_sysrem(star_coeff, image_coeff, numlc, numJD, mag, sig, mag_ave, sigma_clip2, upperdiff, saturation,clip);
  old_star_coeff = (double *) malloc(numlc*sizeof(double));
  old_image_coeff = (double *) malloc(numJD*sizeof(double));
  while (iternum <= maxits && !converged)
    {
      if(starttype == 1)
	{
	  for(j = 0; j < numlc;j++)
	    old_star_coeff[j] = star_coeff[j];
	  decorr_star_sysrem(star_coeff,mag,sig,mag_ave,numlc,numJD,clip,sigma_clip2,image_coeff,saturation);
	  starttype = 0;
	}
      else
	{
	  for(k = 0; k < numJD; k++)
	    old_image_coeff[k] = image_coeff[k];
	  decorr_image_sysrem(star_coeff,mag,sig,mag_ave,numlc,numJD,clip,sigma_clip2,upperdiff,image_coeff,saturation);
	  /* adjust the image_coeffs to have stddev = 1. */
	  var1 = 0.;
	  var2 = 0.;
	  ngood = 0;
	  for(k = 0; k < numJD; k++)
	    {
	      if(!isnan(image_coeff[k]) && !isinf(image_coeff[k]))
		{
		  var1 += image_coeff[k];
		  var2 += image_coeff[k]*image_coeff[k];
		  ngood++;
		}
	    }
	  image_coeff_stddev = sqrt((double)(var2 / (long double) ngood - (var1*var1 / (long double) (ngood * ngood))));
	  image_coeff_ave = (double) (var1 / (long double) ngood);
	  for(k=0; k < numJD; k++)
	    {
	      if(!isnan(image_coeff[k]) && !isinf(image_coeff[k]))
		{
		  image_coeff[k] = (image_coeff[k] - image_coeff_ave)/image_coeff_stddev;
		}
	    }
	  starttype = 1;
	}
      iternum++;
      newchisqr = itersolve_chisqr_sysrem(star_coeff, image_coeff, numlc, numJD, mag, sig, mag_ave, sigma_clip2, upperdiff,saturation,clip);
      if(newchisqr > oldchisqr)
	{
	  if(starttype == 1)
	    for(k = 0; k < numJD; k++)
	      image_coeff[k] = old_image_coeff[k];
	  if(starttype == 0)
	    for(j = 0; j < numlc; j++)
	      star_coeff[j] = old_star_coeff[j];
	  converged = 1;
	}
      else if(oldchisqr - newchisqr < convergeval)
	converged = 1;
      oldchisqr = newchisqr;
    }
  decorr_star_sysrem(star_coeff,mag,sig,mag_ave,numlc,numJD,clip,sigma_clip2,image_coeff,saturation);
  free(old_star_coeff);
  free(old_image_coeff);
  magcorr_sysrem(star_coeff,image_coeff,mag,sig,numlc,numJD);
  return iternum;
}

void initialize_sysrem(_Sysrem *Sysrem, int numlcs, int matchstringid)
{
  FILE *dates;

  char *line;
  size_t line_size = MAXLEN;
  int i, j, k;
  void error2(int, char *);

  line = malloc(line_size);

  if((dates = fopen(Sysrem->dates_name,"r")) == NULL)
    error2(ERR_FILENOTFOUND,Sysrem->dates_name);
  Sysrem->Njd = 0;
  while(gnu_getline(&line,&line_size,dates) >= 0)
    Sysrem->Njd++;
  rewind(dates);

  if(matchstringid)
    {
      if((Sysrem->stringid = (char **) malloc(Sysrem->Njd * sizeof(char *))) == NULL ||
	 (Sysrem->stringid_idx = (int *) malloc(Sysrem->Njd * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0;i<Sysrem->Njd;i++)
	{
	  if((Sysrem->stringid[i] = (char *) malloc(MAXIDSTRINGLENGTH * sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	}
    }
  else
    {
      if((Sysrem->JD = (double *) malloc(Sysrem->Njd * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  if(Sysrem->Nsysrem_airmass > 0)
    {
      if((Sysrem->initial_X = (double **) malloc(Sysrem->Nsysrem_airmass * sizeof(double *))) == NULL)
	error(ERR_MEMALLOC);
    }
  if(Sysrem->Nsysrem_color > 0)
    {
      if((Sysrem->initial_colors = (double **) malloc(Sysrem->Nsysrem_color * sizeof(double *))) == NULL)
	error(ERR_MEMALLOC);
    }
  if((Sysrem->final_X = (double **) malloc(Sysrem->Nsysrem_total * sizeof(double *))) == NULL ||
     (Sysrem->mag_ave = (double *) malloc(numlcs * sizeof(double))) == NULL ||
     (Sysrem->rms_out = (double *) malloc(numlcs * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0;i<Sysrem->Nsysrem_airmass;i++)
    if((Sysrem->initial_X[i] = (double *) malloc(Sysrem->Njd * sizeof(double))) == NULL ||
       (Sysrem->final_X[i] = (double *) malloc(Sysrem->Njd * sizeof(double))) == NULL ||
       (Sysrem->final_colors[i] = (double *) malloc(numlcs * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  for(;i<Sysrem->Nsysrem_total;i++)
    if((Sysrem->final_X[i] = (double *) malloc(Sysrem->Njd * sizeof(double))) == NULL ||
       (Sysrem->final_colors[i] = (double *) malloc(numlcs * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  for(i=0;i<Sysrem->Nsysrem_color;i++)
    if((Sysrem->initial_colors[i] = (double *) malloc(numlcs * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);

  for(i=0;i<Sysrem->Nsysrem_color;i++)
    for(j=0;j<numlcs;j++)
      Sysrem->initial_colors[i][j] = Sysrem->initial_colors_readin[j][i];

  i = 0;
  while(gnu_getline(&line,&line_size,dates) >= 0)
    {
      j = 0;
      if(matchstringid)
	{
	  j += parseone(&line[j], (void *) (Sysrem->stringid[i]), VARTOOLS_TYPE_STRING);
	  Sysrem->stringid_idx[i] = i;
	}
      else
	j += parseone(&line[j], (void *) (&Sysrem->JD[i]), VARTOOLS_TYPE_DOUBLE);
      for(k=0;k<Sysrem->Nsysrem_airmass;k++)
	{
	  if(line[j] != '\0' && line[j] != '\n')
	    j += parseone(&line[j], (void *) (&Sysrem->initial_X[k][i]), VARTOOLS_TYPE_DOUBLE);
	  else
	    error2(ERR_INPUTMISSINGCOLUMN,Sysrem->dates_name);
	}
      i++;
    }
  fclose(dates);
  free(line);
  if(matchstringid)
    mysortstringint(Sysrem->Njd, MAXIDSTRINGLENGTH, Sysrem->stringid, Sysrem->stringid_idx);
}

void do_sysrem(ProgramData *p, _Sysrem *Sysrem, int numlc, int *Njd_in, double **t_in, double **mag_in, double **sig_in, char **lcnames, int matchstringid, char ***stringid, int **stringid_idx)
{
  FILE *lcout, *trendout;
  char outname[MAXLEN];
  double *JD, **initial_X, **final_X, **initial_colors, **final_colors, sigma_clip1, sigma_clip2, saturation, ftol, **decorr, *Avector, *A_errvector, delmag, initial_X_stddev, initial_X_ave;
  int Nsysrem_color, Nsysrem_airmass, Nsysrem_total, numjd, Nsysrem_done, Ngood;
  double **mag, **sig, *mag_corr1, *mag_corr2, *mag_ave, *old_starcoeff, *rms_out, *mag_out, outchi2val;
  int **clip, i, j, k, l, *numX, i1, i2, *order;
  long double var1, var2, var3, var4;

  JD = Sysrem->JD; initial_X = Sysrem->initial_X; final_X = Sysrem->final_X; initial_colors = Sysrem->initial_colors; final_colors = Sysrem->final_colors; sigma_clip1 = Sysrem->sigma_clip1; sigma_clip2 = Sysrem->sigma_clip2; saturation = Sysrem->saturation; Nsysrem_color = Sysrem->Nsysrem_color; Nsysrem_airmass = Sysrem->Nsysrem_airmass; Nsysrem_total = Sysrem->Nsysrem_total; numjd = Sysrem->Njd; mag_ave = Sysrem->mag_ave; rms_out = Sysrem->rms_out;

  /* Store memory for the mag/sig vectors to work on, this is necessary because the light curves may have different lengths */
  if((mag = (double **) malloc(numlc * sizeof(double *))) == NULL ||
     (mag_out = (double *) malloc(numjd * sizeof(double))) == NULL ||
     (sig = (double **) malloc(numlc * sizeof(double *))) == NULL ||
     (clip = (int **) malloc(numlc * sizeof(int *))) == NULL ||
     (mag_corr1 = (double *) malloc(numlc * sizeof(double))) == NULL ||
     (mag_corr2 = (double *) malloc(numjd * sizeof(double))) == NULL ||
     (old_starcoeff = (double *) malloc(numlc * sizeof(double))) == NULL ||
     (numX = (int *) malloc(numjd * sizeof(int))) == NULL ||
     (decorr = (double **) malloc(numjd * sizeof(double *))) == NULL ||
     (Avector = (double *) malloc((Nsysrem_total + 1) * sizeof(double))) == NULL ||
     (A_errvector = (double *) malloc((Nsysrem_total + 1)* sizeof(double))) == NULL ||
     (order = (int *) malloc(Nsysrem_total * sizeof(int))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0;i<numlc;i++)
    {
      if((mag[i] = (double *) malloc(numjd * sizeof(double))) == NULL ||
	 (sig[i] = (double *) malloc(numjd * sizeof(double))) == NULL ||
	 (clip[i] = (int *) malloc(numjd * sizeof(int))) == NULL)
	error(ERR_MEMALLOC);
    }
  for(i=0;i<numjd;i++)
    {
      if((decorr[i] = (double *) malloc(Nsysrem_total * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
  for(i=0;i<numjd;i++)
    numX[i] = 0;

  /* Set up the mag and sig vectors */
  for(j=0;j<numlc;j++)
    {
      k = 0;
      if(!matchstringid)
	{
	  for(i=0;i<Njd_in[j];i++)
	    {
	      while(t_in[j][i] > JD[k] + p->JDTOL && k < numjd){
		mag[j][k] = sqrt(-1);
		sig[j][k] = sqrt(-1);
		//X[j][k] = INDEFVAL - 1.0;
		k++;
	      }
	      if(k < numjd && t_in[j][i] > JD[k] - p->JDTOL && t_in[j][i] < JD[k] + p->JDTOL)
		{
		  if(!isnan(mag_in[j][i]) && !isinf(mag_in[j][i]))
		    {
		      mag[j][k] = mag_in[j][i];
		      sig[j][k] = sig_in[j][i];
		      if(sig[j][k] < SIGMIN)
			sig[j][k] = SIGMIN;
		      numX[k]++;
		    }
		  else
		    {
		      mag[j][k] = sqrt(-1);
		      sig[j][k] = sqrt(-1);
		    }
		  k++;
		}
	    }
	  for(; k<numjd; k++)
	    {
	      mag[j][k] = sqrt(-1);
	      sig[j][k] = sqrt(-1);
	    }
	}
      else
	{
	  for(i=0;i<Njd_in[j];i++)
	    {
	      while(k < numjd ? strncmp(stringid[j][stringid_idx[j][i]], Sysrem->stringid[Sysrem->stringid_idx[k]], MAXIDSTRINGLENGTH) > 0 : 0) {
		mag[j][Sysrem->stringid_idx[k]] = sqrt(-1);
		sig[j][Sysrem->stringid_idx[k]] = sqrt(-1);
		//X[j][k] = INDEFVAL - 1.0;
		k++;
	      }
	      if(k < numjd ? !strncmp(stringid[j][stringid_idx[j][i]], Sysrem->stringid[Sysrem->stringid_idx[k]], MAXIDSTRINGLENGTH) : 0)
		{
		  if(!isnan(mag_in[j][stringid_idx[j][i]]) && !isinf(mag_in[j][stringid_idx[j][i]]))
		    {
		      mag[j][Sysrem->stringid_idx[k]] = mag_in[j][stringid_idx[j][i]];
		      sig[j][Sysrem->stringid_idx[k]] = sig_in[j][stringid_idx[j][i]];
		      if(sig[j][Sysrem->stringid_idx[k]] < SIGMIN)
			sig[j][Sysrem->stringid_idx[k]] = SIGMIN;
		      numX[Sysrem->stringid_idx[k]]++;
		    }
		  else
		    {
		      mag[j][Sysrem->stringid_idx[k]] = sqrt(-1);
		      sig[j][Sysrem->stringid_idx[k]] = sqrt(-1);
		    }
		  k++;
		}
	    }
	  for(; k<numjd; k++)
	    {
	      mag[j][Sysrem->stringid_idx[k]] = sqrt(-1);
	      sig[j][Sysrem->stringid_idx[k]] = sqrt(-1);
	    }
	}
    }

  /* Get the light curve statistics */
  get_magave_std_sysrem(mag,sig,mag_ave,numlc,numjd,sigma_clip1,Sysrem->useweights);

  // Determine which points to clip in calculating chisqr for the iterative solution
  // Do this now so that these points don't end up wandering out of the calculation
  if(sigma_clip2 > 0)
    {
      for(j=0;j<numlc;j++)
	{
	  for(k=0;k<numjd;k++)
	    {
	      if(!isnan(mag[j][k]) && !isinf(mag[j][k]) && mag[j][k] > saturation && ABS_((mag[j][k] - mag_ave[j])/sig[j][k]) < sigma_clip2)
		clip[j][k] = 1;
	      else
		clip[j][k] = 0;
	    }
	}
    }
  else
    {
      for(j=0;j<numlc;j++)
	{
	  for(k=0;k<numjd;k++)
	    {
	      if(!isnan(mag[j][k]) && !isinf(mag[j][k]) && mag[j][k] > saturation)
		clip[j][k] = 1;
	      else
		clip[j][k] = 0;
	    }
	}
    }

  /* The sysrem algorithm begins here */
  Nsysrem_done = 0;
  l=0;
  // Now do the airmass corrections
  for(Nsysrem_done=0;Nsysrem_done < Nsysrem_airmass; Nsysrem_done++)
    {
      var1 = 0.;
      var2 = 0.;
      Ngood = 0;
      for(k=0;k<numjd;k++)
	{
	  if(!isnan(initial_X[Nsysrem_done][k]) && !isinf(initial_X[Nsysrem_done][k]))
	    {
	      var1 += (long double) initial_X[Nsysrem_done][k];
	      var2 += (long double) initial_X[Nsysrem_done][k]*initial_X[Nsysrem_done][k];
	      Ngood++;
	    }
	}
      initial_X_stddev = sqrt((double) (var2 / (long double) Ngood - (var1*var1 / (long double) (Ngood * Ngood))));
      initial_X_ave = (double) (var1 / (long double) Ngood);
      for(k=0;k<numjd;k++)
	mag_corr2[k] = (initial_X[Nsysrem_done][k] - initial_X_ave)/initial_X_stddev;
      decorr_star_sysrem(mag_corr1,mag,sig,mag_ave,numlc,numjd,clip,sigma_clip2,initial_X[Nsysrem_done],saturation);
      ftol = 1.0E-4;
      iter_solve_sysrem(1, MAXITS, ftol, mag_corr1, mag_corr2,mag,sig,mag_ave,numlc,numjd,clip,sigma_clip2,0.,saturation);

      /* Save the terms */
      var1 = 0.;
      var2 = 0.;
      Ngood = 0;
      for(k=0;k<numjd;k++)
	{
	  if(!isnan(mag_corr2[k]) && !isinf(mag_corr2[k]))
	    {
	      var1 += (long double) mag_corr2[k];
	      var2 += (long double) mag_corr2[k]*mag_corr2[k];
	      Ngood++;
	    }
	}
      initial_X_stddev = sqrt((double) (var2 / (long double) Ngood - (var1*var1 / (long double) (Ngood * Ngood))));
      initial_X_ave = (double) (var1 / (long double) Ngood);

      for(k=0;k<numjd;k++)
	{
	  final_X[l][k] = (mag_corr2[k] - initial_X_ave)/initial_X_stddev;
	  mag_corr2[k] = final_X[l][k];
	}

      for(k=0;k<numlc;k++)
	final_colors[l][k] = mag_corr1[k];
      l++;

      magcorr_sysrem(mag_corr1,mag_corr2,mag,sig,numlc,numjd);
    }

  //Next correct for color terms
  for(Nsysrem_done = 0; Nsysrem_done < Nsysrem_color; Nsysrem_done++)
    {
      var1 = 0;
      var2 = 0;
      for(j=0; j < numlc; j++)
	{
	  var1 += initial_colors[Nsysrem_done][j];
	  var2 += initial_colors[Nsysrem_done][j]*initial_colors[Nsysrem_done][j];
	}
      for(j=0;j<numlc;j++)
	mag_corr1[j] = (initial_colors[Nsysrem_done][j] - (var1/numlc))/sqrt((double)((var2/numlc)-(var1*var1/(numlc*numlc))));


      iter_solve_sysrem(0, MAXITS, ftol, mag_corr1, mag_corr2,mag,sig,mag_ave,numlc,numjd,clip,sigma_clip2,0.,saturation);

      /* Save the terms */
      var1 = 0.;
      var2 = 0.;
      Ngood = 0;
      for(k=0;k<numjd;k++)
	{
	  if(!isnan(mag_corr2[k]) && !isinf(mag_corr2[k]))
	    {
	      var1 += (long double) mag_corr2[k];
	      var2 += (long double) mag_corr2[k]*mag_corr2[k];
	      Ngood++;
	    }
	}
      initial_X_stddev = sqrt((double) (var2 / (long double) Ngood - (var1*var1 / (long double) (Ngood * Ngood))));
      initial_X_ave = (double) (var1 / (long double) Ngood);

      for(k=0;k<numjd;k++)
	{
	  final_X[l][k] = (mag_corr2[k] - initial_X_ave)/initial_X_stddev;
	  mag_corr2[k] = final_X[l][k];
	}

      for(k=0;k<numlc;k++)
	final_colors[l][k] = mag_corr1[k];
      l++;

      magcorr_sysrem(mag_corr1,mag_corr2,mag,sig,numlc,numjd);

    }

  /* Now decorrelate the light curves against the full trend set */
  for(k=0;k<numjd;k++)
    {
      for(j=0;j<Nsysrem_total;j++)
	decorr[k][j] = final_X[j][k];
    }
  for(j=0;j<Nsysrem_total;j++)
    order[j] = 1;
  for(k=0;k<numlc;k++)
    {
      docorr(mag[k], sig[k], numjd, Nsysrem_total, decorr, order, Avector, A_errvector, mag_ave[k], 1);
      for(j=0;j<numjd;j++)
	mag_out[j] = mag[k][j];
      mag_ave[k] = Avector[0];
      if(matchstringid)
	magcorr((void *) (Sysrem->stringid), VARTOOLS_TYPE_STRING, mag_out, sig[k], numjd, Nsysrem_total, decorr, order, Avector, &outchi2val, &mag_ave[k], mag_ave[k], 0, NULL, 1);
      else
	magcorr((void *) JD, VARTOOLS_TYPE_DOUBLE, mag_out, sig[k], numjd, Nsysrem_total, decorr, order, Avector, &outchi2val, &mag_ave[k], mag_ave[k], 0, NULL, 1);
      for(j=0;j<Nsysrem_total;j++)
	final_colors[j][k] = Avector[j+1];

      /* Correct the original light curve and print out the model if we're asked to */
      if(Sysrem->omodel)
	{
	  i1 = 0;
	  i2 = 0;
	  while(lcnames[k][i1] != '\0')
	    {
	      if(lcnames[k][i1] == '/')
		i2 = i1 + 1;
	      i1++;
	    }
	  sprintf(outname,"%s/%s%s",Sysrem->model_outdir,&lcnames[k][i2],Sysrem->model_suffix);
	  if((lcout = fopen(outname,"w")) == NULL)
	    error2(ERR_CANNOTWRITE,outname);
	}
      var1 = 0.;
      var2 = 0.;
      var3 = 0.;
      var4 = 0.;
      Ngood = 0;
      if(!matchstringid)
	{
	  for(i=0;i<Njd_in[k];i++)
	    {
	      l = 0;
	      while(t_in[k][i] > JD[l] + p->JDTOL && l < numjd){
		l++;
	      }
	      if(l < numjd && t_in[k][i] > JD[l] - p->JDTOL && t_in[k][i] < JD[l] + p->JDTOL)
		{
		  if(!isnan(mag_in[k][i]) && !isinf(mag_in[k][i]))
		    {
		      if(Sysrem->omodel)
			{
			  fprintf(lcout,"%f %f %f %f %d\n",JD[l],mag_in[k][i],mag_in[k][i] - mag_out[l] + mag_ave[k],sig_in[k][i], clip[k][l]);
			}
		      delmag = mag_out[l];
		      if(clip[k][l])
			{
			  var1 += (long double) delmag;
			  var2 += (long double) (delmag * delmag);
			  Ngood++;
			}
		      if(Sysrem->correctlc)
			{
			  mag_in[k][i] = delmag;
			}
		    }
		  l++;
		}
	    }
	}
      else
	{
	  for(i=0;i<Njd_in[k];i++)
	    {
	      l = 0;
	      while(l < numjd ? strncmp(stringid[k][stringid_idx[k][i]], Sysrem->stringid[Sysrem->stringid_idx[l]],MAXIDSTRINGLENGTH) > 0: 0){
		l++;
	      }
	      if(l < numjd ? !strncmp(stringid[k][stringid_idx[k][i]],Sysrem->stringid[Sysrem->stringid_idx[l]],MAXIDSTRINGLENGTH) : 0)
		{
		  if(!isnan(mag_in[k][stringid_idx[k][i]]) && !isinf(mag_in[k][stringid_idx[k][i]]))
		    {
		      if(Sysrem->omodel)
			{
			  fprintf(lcout,"%s %f %f %f %d\n",Sysrem->stringid[Sysrem->stringid_idx[l]],mag_in[k][stringid_idx[k][i]],mag_in[k][stringid_idx[k][i]] - mag_out[Sysrem->stringid_idx[l]] + mag_ave[k],sig_in[k][stringid_idx[k][i]], clip[k][Sysrem->stringid_idx[l]]);
			}
		      delmag = mag_out[Sysrem->stringid_idx[l]];
		      if(clip[k][Sysrem->stringid_idx[l]])
			{
			  var1 += (long double) delmag;
			  var2 += (long double) (delmag * delmag);
			  Ngood++;
			}
		      if(Sysrem->correctlc)
			{
			  mag_in[k][stringid_idx[k][i]] = delmag;
			}
		    }
		  l++;
		}
	    }
	}

      if(Sysrem->omodel)
	fclose(lcout);
      rms_out[k] = sqrt((double) (var2 / (long double) Ngood - (var1*var1 / (long double) (Ngood * Ngood))));

    }

  /* Print out the trends if we're asked to */
  if(Sysrem->otrend)
    {
      if((trendout = fopen(Sysrem->trends_outname,"w")) == NULL)
	error2(ERR_CANNOTWRITE,Sysrem->trends_outname);
      for(k=0;k<numjd;k++)
	{
	  if(matchstringid)
	    fprintf(trendout,"%s",Sysrem->stringid[k]);
	  else
	    fprintf(trendout,"%f",JD[k]);
	  for(j=0;j<Nsysrem_total;j++)
	    fprintf(trendout," %f",final_X[j][k]);
	  fprintf(trendout,"\n");
	}
      fclose(trendout);
    }

  /* Free all the non-static memory that we allocated */
  for(i=0;i<numlc;i++)
    {
      free(mag[i]);
      free(sig[i]);
      free(clip[i]);
    }
  for(i=0;i<numjd;i++)
    {
      free(decorr[i]);
    }
  free(mag);
  free(mag_out);
  free(sig);
  free(clip);
  free(mag_corr1);
  free(mag_corr2);
  free(old_starcoeff);
  free(numX);
  free(decorr);
  free(Avector);
  free(A_errvector);
  free(order);
}


