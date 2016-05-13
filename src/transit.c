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

/* This is a c implementation of the occultquad.f function for analytic transit light curves due to Mandel & Agol 2002 */

#define SQR_(A) ((A)*(A))

#ifndef PARALLEL
int NRVcurve;
double *RVcurveJD;
double *RVcurveRV;
double *RVcurvesig;
int dofitRV;
#else
typedef struct {
  int NRVcurve;
  double *RVcurveJD;
  double *RVcurveRV;
  double *RVcurvesig;
  int dofitRV;
} _TypeRVdata;
#endif

#define DEFAULTSIZEMODELRVCURVE 1000

void occultquad(double *z0, double u1, double u2, double p, double *muo1, double *mu0, int nz)
{
  /* C  This routine computes the lightcurve for occultation
C  of a quadratically limb-darkened source without microlensing.
C  Please cite Mandel & Agol (2002) if you make use of this routine
C  in your research.  Please report errors or bugs to agol@tapir.caltech.edu
  */
  int i;
  double *mu, *lambdad, *etad, *lambdae, lam, pi, x1, x2, x3, z, omega, kap0, kap1, q, Kk, Ek, Pk, n;

  double ellec(double);
  double ellk(double);
  double rc(double x, double y);
  double rj(double, double, double, double);
  double rf(double, double, double);

  if((mu = (double *) malloc(nz * sizeof(double))) == NULL ||
     (lambdad = (double *) malloc(nz * sizeof(double))) == NULL ||
     (etad = (double *) malloc(nz * sizeof(double))) == NULL ||
     (lambdae = (double *) malloc(nz * sizeof(double))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(2);
    }

  if(fabs(p - 0.5) < 1.E-3) p = 0.5;

  /*
z0 and muo1 and mu0 should be of length nz

C
C Input:
C
C rs   radius of the source (set to unity)
C z0   impact parameter in units of rs
C p    occulting star size in units of rs
C u1   linear    limb-darkening coefficient (gamma_1 in paper)
C u2   quadratic limb-darkening coefficient (gamma_2 in paper)
C
C Output:
C
C muo1 fraction of flux at each z0 for a limb-darkened source
C mu0  fraction of flux at each z0 for a uniform source
C
C Limb darkening has the form:
C  I(r)=[1-u1*(1-sqrt(1-(r/rs)^2))-u2*(1-sqrt(1-(r/rs)^2))^2]/(1-u1/3-u2/6)/pi
C
C To use this routine
C
C Now, compute pure occultation curve:
  */
  omega = 1.0 - u1/3.0 - u2/6.0;
  pi = acos(-1.);
  for(i=0;i<nz;i++)
    {
      z = z0[i];
      if(fabs(1.0 - z - p) < 1.E-5)
	z = z - 1.E-5;
      x1 = (p-z);
      x1 = x1*x1;
      x2 = (p+z);
      x2 = x2*x2;
      x3 = p*p - z*z;

      /* The source is unocculted:
	 Table 3, I.:
      */
      if(z >= 1.0+p) {
	lambdad[i] = 0.0;
	etad[i] = 0.0;
	lambdae[i] = 0.0;
      }

      /* The source is completely occulted:
	 Table 3, II.:
      */
      else if(p >= 1.0 && z <= p-1.0) {
	lambdad[i] = 1.0;
	etad[i] = 1.0;
	lambdae[i] = 1.0;
      }

      /* The source is partly occulted and the occulting object crosses the limb:
	 Equation (26):
      */
      else
	{
	  if(z >= fabs(1.0-p) && z <= 1.0+p) {
	    kap1 = acos(MIN_((1.0-p*p+z*z)/2.0/z,1.0));
	    kap0 = acos(MIN_((p*p+z*z-1.0)/2.0/p/z,1.0));
	    lambdae[i] = p*p*kap0+kap1;
	    lambdae[i] = (lambdae[i]-0.5*sqrt(MAX_(4.0*z*z-SQR_(1.0+z*z-p*p),0.0)))/pi;
	  }
	  /* The occulting object transit the source star (but doesn't completely cover it): */
	  if(z <= 1.0-p) {
	    lambdae[i] = p*p;
	  }
      /* The edge of the occulting star lies at the origin-special expressions in this case: */
	  if(fabs(z-p) < 0.0001*(z+p)) {
	    /* Table 3, Case V.: */
	    if(z >= 0.5) {
	      lam = 0.5*pi;
	      q = 0.5/p;
	      Kk = ellk(q);
	      Ek = ellec(q);
	      lambdad[i] = 1.0/3.0 + 16.0*p/9.0/pi*(2.0*p*p-1.0)*Ek - (32.0*p*p*p*p-20.0*p*p+3.0)/9.0/pi/p*Kk;
	      etad[i] = 1.0/2.0/pi*(kap1+p*p*(p*p+2.0*z*z)*kap0-(1.0+5.0*p*p+z*z)/4.0*sqrt((1.0-x1)*(x2-1.0)));
	      if(p == 0.5) {
		/* Case VIII: p=1/2, z=1/2 */
		lambdad[i] = 1.0/3.0 - 4.0/pi/9.0;
		etad[i] = 3.0/32.0;
	      }
	    }
	    else {
	      /* Table 3, Case VI.: */
	      lam=0.5*pi;
	      q = 2.0*p;
	      Kk=ellk(q);
	      Ek=ellec(q);
	      lambdad[i] = 1.0/3.0+2.0/9.0/pi*(4.0*(2.0*p*p-1.0)*Ek+(1.0-4.0*p*p)*Kk);
	      etad[i]=p*p/2.0*(p*p+2.0*z*z);
	    }
	  }
	  /* The occulting star partly occults the source and crosses the limb:
	     Table 3, Case III:
	  */
	  else if((z > 0.5 + fabs(p-0.5) && z < 1.0+p) || (p > 0.5 && z > fabs(1.0-p)*1.0001 && z < p)) {
	    lam=0.5*pi;
	    q = sqrt((1.0-SQR_(p-z))/4.0/z/p);
	    Kk=ellk(q);
	    Ek=ellec(q);
	    n=1.0/x1-1.0;
	    Pk=Kk-n/3.0*rj(0.0,1.0-q*q,1.0,1.0+n);
	    lambdad[i] = 1.0/9.0/pi/sqrt(p*z)*(((1.0-x2)*(2.0*x2+x1-3.0)-3.0*x3*(x2-2.0))*Kk+4.0*p*z*(z*z+7.0*p*p-4.0)*Ek-3.0*x3/x1*Pk);
	    if(z < p) lambdad[i]=lambdad[i]+2.0/3.0;
	    etad[i] = 1.0/2.0/pi*(kap1+p*p*(p*p+2.0*z*z)*kap0-(1.0+5.0*p*p+z*z)/4.0*sqrt((1.0-x1)*(x2-1.0)));
	  }
	  /* The occulting star transits the source:
	     Table 3, Case IV.:
	  */
	  else if(p <= 1.0 && z <= (1.0-p)*0.9999) {
	    lam=0.5*pi;
	    q = sqrt((x2-x1)/(1.0-x1));
	    Kk=ellk(q);
	    Ek=ellec(q);
	    n = x2/x1 - 1.0;
	    Pk = Kk - n/3.0*rj(0.0,1.0-q*q,1.0,1.0+n);
	    lambdad[i] = 2.0/9.0/pi/sqrt(1.0-x1)*((1.0-5.0*z*z+p*p+x3*x3)*Kk+(1.0-x1)*(z*z+7.0*p*p-4.0)*Ek-3.0*x3/x1*Pk);
	    if(z < p) lambdad[i] = lambdad[i] + 2.0/3.0;
	    if(fabs(p+z-1.0) <= 0.0001) {
	      lambdad[i] = 2./3.0/pi*acos(1.0-2.0*p)-4.0/9.0/pi*sqrt(p*(1.0-p))*(3.0+2.0*p-8.0*p*p);
	    }
	    etad[i]=p*p/2.0*(p*p+2.0*z*z);
	  }
	  /* The occulting star transits the source, with the limbs coinciding:
	     Table 3, Case IV.:
	  */
	  else if(p <= 1.0 && z > (1.0-p)*0.9999 && z <= (1.0-p)*1.0001) {
	    lam=0.5*pi;
	    q = sqrt(1.0 - 1.1*TINY);
	    Kk=ellk(q);
	    Ek=ellec(q);
	    n = x2/x1 - 1.0;
	    Pk = Kk - n/3.0*rj(0.0,1.0-1.1*TINY,1.0,1.0+n);
	    lambdad[i] = 2.0/9.0/pi/sqrt(1.0-x1)*((1.0-5.0*z*z+p*p+x3*x3)*Kk+(1.0-x1)*(z*z+7.0*p*p-4.0)*Ek-3.0*x3/x1*Pk);
	    if(z < p) lambdad[i] = lambdad[i] + 2.0/3.0;
	    if(fabs(p+z-1.0) <= 0.0001) {
	      lambdad[i] = 2./3.0/pi*acos(1.0-2.0*p)-4.0/9.0/pi*sqrt(p*(1.0-p))*(3.0+2.0*p-8.0*p*p);
	    }
	    etad[i]=p*p/2.0*(p*p+2.0*z*z);
	  }
	}
      /* Now, using equation (33): */
      muo1[i]=1.0-((1.0-u1-2.0*u2)*lambdae[i]+(u1+2.0*u2)*lambdad[i]+u2*etad[i])/omega;
      mu0[i]=1.0-lambdae[i];
    }
  free(mu);
  free(lambdad);
  free(etad);
  free(lambdae);
  return;
}

#ifdef ERRTOL
#undef ERRTOL
#endif
#define ERRTOL 0.08
#ifdef TINY
#undef TINY
#endif
#define TINY 1.5E-38
#ifdef BIG
#undef BIG
#endif
#define BIG 3.0E37

#define MAX3_(A,B,C) ((A) > (B) ? (MAX_(A,C)) : (MAX_(B,C)))
#define MIN3_(A,B,C) ((A) < (B) ? (MIN_(A,C)) : (MIN_(B,C)))

double rf(double x, double y, double z)
{
  double c1, c2, c3, c4, alamb, ave, delx, dely, delz, e2, e3, sqrtx, sqrty, sqrtz, xt, yt, zt, rf_val, THIRD;
  THIRD = 1.0/3.0;
  c1 = 1.0/24.0; c2 = 0.1; c3 = 3.0/44.0; c4 = 1.0/14.0;
  if(MIN3_(x,y,z) < 0.0 || MIN3_(x+y,x+z,y+z) < TINY || MAX3_(x,y,z) > BIG)
    {
      fprintf(stderr,"Invalid arguments in rf\nAborting!\n");
      exit(4);
    }
  xt = x;
  yt = y;
  zt = z;
  do
    {
      sqrtx = sqrt(xt);
      sqrty = sqrt(yt);
      sqrtz = sqrt(zt);
      alamb = sqrtx*(sqrty + sqrtz) + sqrty*sqrtz;
      xt = 0.25 * (xt + alamb);
      yt = 0.25 * (yt + alamb);
      zt = 0.25 * (zt + alamb);
      ave = THIRD * (xt + yt + zt);
      delx = (ave - xt)/ave;
      dely = (ave - yt)/ave;
      delz = (ave - zt)/ave;
    } while (MAX3_(fabs(delx),fabs(dely),fabs(delz)) > ERRTOL);
  e2 = delx*dely-delz*delz;
  e3 = delx*dely*delz;
  rf_val=(1.0 + (c1*e2 - c2 - c3*e3)*e2 + c4*e3)/sqrt(ave);
  return(rf_val);
}


#ifdef ERRTOL
#undef ERRTOL
#endif
#define ERRTOL 0.04
#ifdef TINY
#undef TINY
#endif
#define TINY 1.69e-38
#ifdef SQRTNY
#undef SQRTNY
#endif
#define SQRTNY 1.3e-19
#ifdef BIG
#undef BIG
#endif
#define BIG 3.0e37

double rc(double x, double y)
{
  double TNBG, COMP1, COMP2, THIRD, C1, C2, C3, C4, alamb, ave, s, w, xt, yt, rc_val;
  TNBG = TINY*BIG;
  COMP1 = 2.236/SQRTNY;
  COMP2 = TNBG*TNBG/25.0;
  THIRD = 1.0/3.0;
  C1 = 0.3;
  C2 = 1.0/7.0;
  C3 = 0.375;
  C4 = 9.0/22.0;

  if(x < 0 || y == 0 || (x + fabs(y)) < TINY || (x + fabs(y)) > BIG || (y < -COMP1 && x > 0 && x < COMP2))
    {
      fprintf(stderr,"Invalid arguments in rc\n");
      exit(10);
    }

  if(y > 0.0) {
    xt = x;
    yt = y;
    w = 1.;
  } else {
    xt = x - y;
    yt = -y;
    w = sqrt(x)/sqrt(xt);
  }
  do {
    alamb = 2.0*sqrt(xt)*sqrt(yt) + yt;
    xt = 0.25*(xt+alamb);
    yt = 0.25*(yt+alamb);
    ave = THIRD*(xt + yt + yt);
    s = (yt - ave) / ave;
  } while (fabs(s) > ERRTOL);
  rc_val = w*(1.0 + s*s*(C1 + s*(C2 + s*(C3 + s*C4))))/sqrt(ave);
  return rc_val;
} /* Copr 1986-92 Numerical Recipes Software */

#ifdef ERRTOL
#undef ERRTOL
#endif
#define ERRTOL 0.05
#ifdef TINY
#undef TINY
#endif
#define TINY 2.5e-13
#ifdef BIG
#undef BIG
#endif
#define BIG 9.0e11

//#define MIN3(A,B,C) ((A) <= (B) && (A) <= (C) ? (A) : ((B) <= (A) && (B) <= (C) ? (B) : (C)))
//#define MAX3(A,B,C) ((A) >= (B) && (A) >= (C) ? (A) : ((B) >= (A) && (B) >= (C) ? (B) : (C)))
#define MIN4(A,B,C,D) ((A) <= (B) && (A) <= (C) && (A) <= (D) ? (A) : ((B) <= (A) && (B) <= (C) && (B) <= (D) ? (B) : ((C) <= (A) && (C) <= (B) && (C) <= (D) ? (C) : (D))))
#define MAX4(A,B,C,D) ((A) >= (B) && (A) >= (C) && (A) >= (D) ? (A) : ((B) >= (A) && (B) >= (C) && (B) >= (D) ? (B) : ((C) >= (A) && (C) >= (B) && (C) >= (D) ? (C) : (D))))

double rj(double x, double y, double z, double p)
{
  double rj_val, C1, C2, C3, C4, C5, C6, C7, C8, a, alamb, alpha, ave, b, beta, delp, delx, dely, delz, ea, eb, ec, ed, ee, fac, pt, rcx, rho, sqrtx, sqrty, sqrtz, sum, tau, xt, yt, zt;

  C1 = 3.0/14.0; C2 = 1.0/3.0; C3 = 3.0/22.0; C4 = 3.0/26.0; C5 = 0.75*C3;
  C6 = 1.5*C4; C7 = 0.5*C2; C8 = C3 + C3;

  /* This program calls rc and rf */
  if(MIN3_(x,y,z) < 0. || MIN4(x+y,x+z,y+z,fabs(p)) < TINY || MAX4(x,y,z,fabs(p)) > BIG)
    {
      fprintf(stderr,"Invalid arguments in rj\n");
      exit(10);
    }
  sum = 0.0;
  fac = 1.0;
  if(p > 0.0) {
    xt = x;
    yt = y;
    zt = z;
    pt = p;
  }
  else {
    xt = MIN3_(x,y,z);
    zt = MAX3_(x,y,z);
    yt = x+y+z-xt-zt;
    a = 1.0/(yt-p);
    b = a*(zt-yt)*(yt-xt);
    pt = yt+b;
    rho = xt*zt/yt;
    tau = p*pt/yt;
    rcx = rc(rho,tau);
  }
  do {
    sqrtx = sqrt(xt);
    sqrty = sqrt(yt);
    sqrtz = sqrt(zt);
    alamb = sqrtx*(sqrty+sqrtz)+sqrty*sqrtz;
    alpha = SQR_(pt*(sqrtx+sqrty+sqrtz)+sqrtx*sqrty*sqrtz);
    beta = pt*SQR_(pt+alamb);
    sum = sum+fac*rc(alpha,beta);
    fac = 0.25*fac;
    xt = 0.25*(xt+alamb);
    yt = 0.25*(yt+alamb);
    zt = 0.25*(zt + alamb);
    pt = 0.25*(pt + alamb);
    ave = 0.2*(xt + yt + zt + pt + pt);
    delx = (ave - xt)/ave;
    dely = (ave - yt)/ave;
    delz = (ave - zt)/ave;
    delp = (ave - pt)/ave;
  } while (MAX4(fabs(delx),fabs(dely),fabs(delz),fabs(delp))>ERRTOL);
  ea = delx*(dely + delz) + dely*delz;
  eb = delx*dely*delz;
  ec = delp*delp;
  ed = ea - 3.0*ec;
  ee = eb + 2.0*delp*(ea - ec);
  rj_val = 3.0*sum+fac*(1.0+ed*(-C1+C5*ed-C6*ee)+eb*(C7+delp*(-C8+delp*C4))+delp*ea*(C2-delp*C3)-C2*delp*ec)/(ave*sqrt(ave));
  if(p <= 0.0)
    rj_val = a*(b*rj_val + 3.0*(rcx-rf(xt,yt,zt)));
  return rj_val;
} /* Copr. 1986 - 92 Numerical Recipes Software */

double ellec(double k)
{
  double m1, a1, a2, a3, a4, b1, b2, b3, b4, ee1, ee2, ellec_val;
  m1=1.0-k*k;
  a1=0.44325141463;
  a2=0.06260601220;
  a3=0.04757383546;
  a4=0.01736506451;
  b1=0.24998368310;
  b2=0.09200180037;
  b3=0.04069697526;
  b4=0.00526449639;
  ee1=1.0+m1*(a1+m1*(a2+m1*(a3+m1*a4)));
  ee2=m1*(b1+m1*(b2+m1*(b3+m1*b4)))*log(1.0/m1);
  ellec_val=ee1+ee2;
  return(ellec_val);
}

double ellk(double k)
{
  double a0,a1,a2,a3,a4,b0,b1,b2,b3,b4,ellk_val, ek1,ek2,m1;
  m1=1.0-k*k;
  a0=1.38629436112;
  a1=0.09666344259;
  a2=0.03590092383;
  a3=0.03742563713;
  a4=0.01451196212;
  b0=0.5;
  b1=0.12498593597;
  b2=0.06880248576;
  b3=0.03328355346;
  b4=0.00441787012;
  ek1=a0+m1*(a1+m1*(a2+m1*(a3+m1*a4)));
  ek2=(b0+m1*(b1+m1*(b2+m1*(b3+m1*b4))))*log(m1);
  ellk_val=ek1-ek2;
  return(ellk_val);
}

#define NMAX 513

#define MAX_(A,B) ((A) > (B) ? (A) : (B))
#define MIN_(A,B) ((A) < (B) ? (A) : (B))

void occultuniform(double *b0, double w, double *muo1, int nb)
{
  /*C
   This routine computes the lightcurve for occultation
C  of a uniform source without microlensing  (Mandel & Agol 2002).
C Input:
C
C  rs   radius of the source (set to unity)
C  b0   impact parameter in units of rs
C  w    occulting star size in units of rs
C
C Output:
C  muo1 fraction of flux at each b0 for a uniform source
C  */
  int i;
  double z, pi, lambdae, kap0, kap1;
  if(fabs(w-0.5) < 0.001) w = 0.5;
  pi = acos(-1.);

  /* Now, compute pure occulation curve: */
  for(i=0;i<nb;i++)
    {
      /* substitute z=b0[i] to shorten expressions */
      z = b0[i];
      /* The source is unocculted:
	 Table 3, I. */
      if(z >= 1.0+w) {
	muo1[i] = 1.;
	continue;
      }
      /* The source is completely occulted:
	 Table 3, II. */
      if(w >= 1.0 && z <= w-1.0) {
	muo1[i] = 0.0;
	continue;
      }
      /* The source is partly occulted and the occulting object crosses the limb:
	 Equation (26):
      */
      if(z >= fabs(1.0 - w) && z <= 1.0+w) {
	kap1 = acos(MIN_((1.0-w*w+z*z)/2.0/z,1.0));
	kap0 = acos(MIN_((w*w+z*z-1.0)/2.0/w/z,1.0));
	lambdae = w*w*kap0+kap1;
	lambdae = (lambdae-0.5*sqrt(MAX_(4.0*z*z-SQR_(1.0+z*z-w*w),0.0)))/pi;
	muo1[i] = 1.0 - lambdae;
	continue;
      }
      /* The occulting object transits the source star (but doesn't completely cover it): */
      if(z <= 1.0-w)
	muo1[i] = 1.0-w*w;
    }
  return;
}


void occultnl(double rl, double c1, double c2, double c3, double c4, double *b0, double *mulimb0, double **mulimbf, int nb)
{
     /*
C  This routine uses the results for a uniform source to
C  compute the lightcurve for a limb-darkened source
C  (5-1-02 notes)
C Input:
C   rl        radius of the lens   in units of the source radius
C   c1-c4     limb-darkening coefficients
C   b0        impact parameter normalized to source radius
C Output:
C  mulimb0 limb-darkened magnification
C  mulimbf lightcurves for each component
C
C  First, make grid in radius:
C  Call magnification of uniform source:

Memory for mulimbf and mulimb0 must be allocated before calling this function, they are mulimb0 is assumed to be of size [0....nb-1] while mulimbf is of size [0.....nb-1][0...4]
     */
     int i, j, nr, i1, i2;
     double pi, *bt0, *mulimb, *mulimbp, dt, sig, *mulimb1, *mulimbhalf, *mulimb3half, *mulimb2, sig1, sig2, omega, dmumax, fac, *mu, f1, f2, term1, term2;
     double *t, *th, *r;

     if((t = (double *) malloc(NMAX * sizeof(double))) == NULL ||
	(th = (double *) malloc(NMAX * sizeof(double))) == NULL ||
	(r = (double *) malloc(NMAX * sizeof(double))) == NULL)
       error(ERR_MEMALLOC);

     if((bt0 = (double *) malloc(nb * sizeof(double))) == NULL ||
	(mulimb = (double *) malloc(nb * sizeof(double))) == NULL ||
	(mulimbp = (double *) malloc(nb * sizeof(double))) == NULL ||
	(mulimb1 = (double *) malloc(nb * sizeof(double))) == NULL ||
	(mulimbhalf = (double *) malloc(nb * sizeof(double))) == NULL ||
	(mulimb3half = (double *) malloc(nb * sizeof(double))) == NULL ||
	(mulimb2 = (double *) malloc(nb * sizeof(double))) == NULL ||
	(mu = (double *) malloc(nb * sizeof(double))) == NULL)
       {
	 fprintf(stderr,"Memory Allocation Error\n");
	 exit(2);
       }

     pi=acos(-1.0);
     occultuniform(b0,rl,mulimb0,nb);
     i1 = nb-1;
     i2 = 0;
     fac = 0.0;
     for(i=0;i<nb;i++)
       {
	 bt0[i] = b0[i];
	 mulimbf[i][0] = 1.;
	 mulimbf[i][1] = 0.8;
	 mulimbf[i][2] = 2.0/3.0;
	 mulimbf[i][3] = 4.0/7.0;
	 mulimbf[i][4] = 0.5;
	 mulimb[i] = mulimb0[i];
	 if(mulimb0[i] != 1.0) {
	   i1 = MIN_(i1,i);
	   i2 = MAX_(i2,i);
	 }
	 fac = MAX_(fac,fabs(mulimb0[i]-1.));
       }
     omega = 4.0*((1.0-c1-c2-c3-c4)/4. + c1/5. + c2/6. + c3/7. + c4/8.);
     nr = 2;
     dmumax = 1.0;
     while(dmumax > fac*0.001 && nr < NMAX-1)
       {
	 for(i=i1;i<=i2;i++)
	   mulimbp[i]=mulimb[i];
	 nr = nr*2;
	 if(nr > NMAX) nr = (NMAX-1);
	 dt = 0.5*pi/(double)nr;
	 for(j=0;j<=nr;j++)
	   {
	     t[j] = dt*(double)(j);
	     th[j] = t[j]+0.5*dt;
	     r[j] = sin(t[j]);
	   }
	 sig = sqrt(cos(th[nr-1]));
	 for(i=i1;i<=i2;i++)
	   {
	     mulimbhalf[i] = sig*sig*sig*mulimb0[i]/(1.0-r[nr-1]);
	     mulimb1[i] = mulimbhalf[i]*sig;
	     mulimb3half[i] = mulimb1[i]*sig;
	     mulimb2[i] = mulimb3half[i]*sig;
	   }
	 for(j=1;j<nr;j++)
	   {
	     for(i=0;i<nb;i++)
	       {
		 b0[i]=bt0[i]/r[j];
	       }
	     occultuniform(b0,rl/r[j],mu,nb);
	     sig1 = sqrt(cos(th[j-1]));
	     sig2 = sqrt(cos(th[j]));
	     dmumax = 0.0;
	     for(i=i1;i<=i2;i++)
	       {
		 f1 = r[j]*r[j]*mu[i]/(r[j]-r[j-1]);
		 f2 = r[j]*r[j]*mu[i]/(r[j+1]-r[j]);
		 term1 = f1*sig1*sig1*sig1;
		 term2 = f2*sig2*sig2*sig2;
		 mulimbhalf[i] = mulimbhalf[i] + term1 - term2;
		 term1 *= sig1;
		 term2 *= sig2;
		 mulimb1[i] = mulimb1[i] + term1 - term2;
		 term1 *= sig1;
		 term2 *= sig2;
		 mulimb3half[i] = mulimb3half[i] + term1 - term2;
		 term1 *= sig1;
		 term2 *= sig2;
		 mulimb2[i] = mulimb2[i] + term1 - term2;
		 mulimb[i] = ((1.0 - c1-c2-c3-c4)*mulimb0[i] + c1*mulimbhalf[i]*dt + c2*mulimb1[i]*dt+c3*mulimb3half[i]*dt + c4*mulimb2[i]*dt) / omega;
		 if(mulimb[i] + mulimbp[i] != 0.0) {
		   dmumax = MAX_(dmumax,fabs(mulimb[i] - mulimbp[i])/(mulimb[i]+mulimbp[i]));
		 }
	       }
	   }
       }
     for(i=i1;i<=i2;i++)
       {
	 mulimbf[i][0] = mulimb0[i];
	 mulimbf[i][1] = mulimbhalf[i]*dt;
	 mulimbf[i][2] = mulimb1[i]*dt;
	 mulimbf[i][3] = mulimb3half[i]*dt;
	 mulimbf[i][4] = mulimb2[i]*dt;
	 mulimb0[i] = mulimb[i];
       }
     for(i=0;i<nb;i++)
       b0[i] = bt0[i];

     free(bt0);
     free(mulimb);
     free(mulimbp);
     free(mulimb1);
     free(mulimbhalf);
     free(mulimb3half);
     free(mulimb2);
     free(mu);
     free(t);
     free(th);
     free(r);
}

double getphaseofconjunction (double e, double omega)
{
  double E, nu, M, x, y, r, phi_c;
  nu = M_PI/2.0 - omega;
  while(nu < 0.)
    nu += 2.0*M_PI;
  while(nu > 2.0*M_PI)
    nu -= 2.0*M_PI;
  r = (1 - e*e)/(1 + e*cos(nu));
  x = e + r*cos(nu);
  y = sqrt(fabs(1 - x*x));
  if(nu > M_PI) y = -1*y;
  E = atan2(y,x);
  M = E - e*sin(E);
  phi_c = M/2.0/M_PI;
  return phi_c;
}

#define BISECTION_SAFETY_MAX_ITERATIONS 1000
#define ECC_ANOMALY_MAX_ERROR 1.0e-4  // Convergence threshold (average error is less than a third of this)


/* The code below is taken from J. Devor's Debil program */
double eccentricAnomaly (double M, double e)
{
  int counter = BISECTION_SAFETY_MAX_ITERATIONS ;
  double Emin = M - 1.0, Emax = M + 1.0, E = M ;
  double f, dfm, dE ;

  do
    {
      f = M + (e * sin(E)) - E ;  // May be optimized by setting cos=sqrt(1-sin^2)
      dfm = 1.0 - (e * cos(E)) ;  // Minus differential of f (always non-negative)
      dE = dfm * ECC_ANOMALY_MAX_ERROR ;

      if (f > dE)
	{
	  Emin = E ;
	  dE = Emax - Emin ;

	  if ((dfm * dE) > f)
	    E += (f / dfm) ;
	  else
	    E = 0.5 * (Emax + Emin) ;
	}
      else if (f < (-dE))
	{
	  Emax = E ;
	  dE = Emax - Emin ;

	  if ((dfm * dE) > -f)
	    E += (f / dfm) ;
	  else
	    E = 0.5 * (Emax + Emin) ;
	}
      else return (E) ;
    }
  while ((dE > ECC_ANOMALY_MAX_ERROR) && ((--counter) > 0)) ;

  return (E) ;  // should almost never exits here
}

double getimpactparameter(double phase, double sin_i, double a, double e, double omega, double phi_c, double p)
{
  const double meanAnomaly = 2.0 * M_PI * (phase + phi_c);
  double D2, D, E;

  if(e > 0)
    E = eccentricAnomaly(meanAnomaly, e);
  else
    E = meanAnomaly;
  D2 = SQR_(1.0 - (e * cos(E))) - SQR_((((cos(E) - e) * sin(omega)) + (sqrt(1.0 - (e*e)) * sin(E) * cos(omega))) * sin_i) ;
  D = sqrt(D2) * a;

  /* Check which side of the orbit the planet is, if the planet is behind the star set D to be something out of transit */
  if(sin(E + omega) < 0.)
    D = p + 2.0;
  return D;
}

void mandelagoltransitmodel(int Npoints, double *phase, double *outlc, int type, double *ldcoeffs, double sin_i, double a, double e, double p, double omega)
/* Returns a Mandel & Agol transit model evaluated at particular phases. 

   Npoints - the number of points in the light curve to evaluate the model at.
   phase - vector of phases to evaluate the model at (phase=0, and phase=1 is 
           the time of mid transit).
   outlc - vector to output the model fractional fluxes to (fractional flux
           is defined such that 1 = no transit).
   type - 0 for quadratic limb darkening (2 parameters in ldcoeffs), 1 for
          non-linear limb-darkening (4 parameters in ldcoeffs).
   ldcoeffs - array of limb-darkening coefficients.
   sin_i - sin of the inclination angle of the orbit. sin_i = 1 corresponds to
          an edge on orbit (with transits), sin_i = 0 corresponds to a face-on
          orbit (no transits)
   a - ratio of semi-major axis to star radius.
   e - eccentricity of the orbit.
   p - ratio of planet radius to star radius.
   omega - argument of periastron (in radians)
*/
{
  double *z, *mu0, *muo1, *mulimb0, **mulimbf, u1, u2, c1, c2, c3, c4, phi_c, cos_i;
  int i;

  if((z = (double *) malloc(Npoints * sizeof(double))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(3);
    }

  /* Get the impact parameters from the phase and orbital parameters */
  phi_c = getphaseofconjunction(e,omega);
  //printf("%f\n",phi_c);
  for(i=0;i<Npoints;i++)
    z[i] = getimpactparameter(phase[i],sin_i,a,e,omega,phi_c,p);

  /* Quadratic limb darkening */
  if(type == 0)
    {
      if((mu0 = (double *) malloc(Npoints * sizeof(double))) == NULL ||
	 (muo1 = (double *) malloc(Npoints * sizeof(double))) == NULL)
	{
	  fprintf(stderr,"Memory Allocation Error\n");
	  exit(3);
	}
      u1 = ldcoeffs[0];
      u2 = ldcoeffs[1];
      occultquad(z, u1, u2, p, muo1, mu0, Npoints);
      for(i=0;i<Npoints;i++)
	outlc[i] = muo1[i];
      free(mu0);
      free(muo1);
    }

  /* Non-linear limb darkening */
  else if(type == 1)
    {
      if((mulimb0 = (double *) malloc(Npoints * sizeof(double))) == NULL ||
	 (mulimbf = (double **) malloc(Npoints * sizeof(double *))) == NULL)
	{
	  fprintf(stderr,"Memory Allocation Error\n");
	  exit(3);
	}
      for(i=0;i<Npoints;i++)
	{
	  if((mulimbf[i] = (double *) malloc(5 * sizeof(double))) == NULL)
	    {
	      fprintf(stderr,"Memory Allocation Error\n");
	      exit(3);
	    }
	}
      c1 = ldcoeffs[0];
      c2 = ldcoeffs[1];
      c3 = ldcoeffs[2];
      c4 = ldcoeffs[3];
      occultnl(p, c1, c2, c3, c4, z, mulimb0, mulimbf, Npoints);
      for(i=0;i<Npoints;i++)
	{
	  outlc[i] = mulimb0[i];
	  free(mulimbf[i]);
	}
      free(mulimb0);
      free(mulimbf);
    }
  free(z);
}

void integratemandelagoltransitmodel(double exptime_phase, int Npoints, double *phase, double *outlc, int type, double *ldcoeffs, double sin_i, double a, double e, double p, double omega, int Nresamp)
/* Returns a Mandel & Agol transit model including integration over exposures

   exptime_phase - the exposure time to integrate over, in units of phase.
   Npoints - the number of points in the light curve to evaluate the model at.
   phase - vector of phases to evaluate the model at (phase=0, and phase=1 is 
           the time of mid transit).
   outlc - vector to output the model fractional fluxes to (fractional flux
           is defined such that 1 = no transit).
   type - 0 for quadratic limb darkening (2 parameters in ldcoeffs), 1 for
          non-linear limb-darkening (4 parameters in ldcoeffs).
   ldcoeffs - array of limb-darkening coefficients.
   sin_i - sin of the inclination angle of the orbit. sin_i = 1 corresponds to
          an edge on orbit (with transits), sin_i = 0 corresponds to a face-on
          orbit (no transits)
   a - ratio of semi-major axis to star radius.
   e - eccentricity of the orbit.
   p - ratio of planet radius to star radius.
   omega - argument of periastron (in radians)
   Nresamp - Number of points per integration at which to evaluate the M&A
             model
*/
{
  double *fluxtmp, *phasetmp;
  double dphasesamp, fluxsum;
  int i, j;

  if(Nresamp <= 1 || exptime_phase <= 0.) {
    mandelagoltransitmodel(Npoints, phase, outlc, type, ldcoeffs, sin_i, a, e, p, omega);
    return;
  }

  if((fluxtmp = (double *) malloc(Nresamp * sizeof(double))) == NULL ||
     (phasetmp = (double *) malloc(Nresamp * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  dphasesamp = exptime_phase/(Nresamp - 1);
  for(i=0; i < Npoints; i++) {
    for(j=0; j < Nresamp; j++) {
      phasetmp[j] = phase[i] - 0.5*exptime_phase + dphasesamp*j;
      if(phasetmp[j] < 0) phasetmp[j] += 1.;
      if(phasetmp[j] <= 1.) phasetmp[j] -= 1.;
    }
    mandelagoltransitmodel(Nresamp, phasetmp, fluxtmp, type, ldcoeffs, sin_i, a, e, p, omega);
    fluxsum = 0.;
    for(j=0; j < Nresamp; j++) {
      fluxsum += fluxtmp[j];
    }
    fluxsum = fluxsum / Nresamp;
    outlc[i] = fluxsum;
  }
  free(fluxtmp);
  free(phasetmp);
  return;
}

void getRVcurve(int N, double *phase, double *RV, double K, double e, double omega, double gamma)
{

  int i;
  double phi_c, meanAnomaly, E, trueAnomaly, cosE, sinE;

  phi_c = getphaseofconjunction(e,omega);

  for(i=0;i<N;i++)
    {

      meanAnomaly = 2.0 * M_PI * (phase[i] + phi_c);

      if(e > 0)
	E = eccentricAnomaly(meanAnomaly, e);
      else
	E = meanAnomaly;

      cosE = cos(E);
      sinE = sin(E);
      trueAnomaly = acos((cosE - e)/(1. - e*cosE));
      if(sinE < 0.)
	trueAnomaly = 2.0*M_PI - trueAnomaly;
      RV[i] = K*(cos(trueAnomaly + omega) + e*cos(omega)) + gamma;
    }
}

double chisqmandelagoltransit(double *p, int ma, int N, double *t, double *mag, double *sig, void *userparams)
{
#ifdef PARALLEL
  double *newlc = NULL, *phase = NULL, *phaseRV = NULL, *newRV = NULL;
  int size_newlc = 0, size_newRV = 0;
  _TypeRVdata *RVdata;
  int NRVcurve;
  double *RVcurveJD;
  double *RVcurveRV;
  double *RVcurvesig;
  int dofitRV;
#else
  static double *newlc, *phase, *phaseRV, *newRV;
  static int size_newlc = 0, size_newRV = 0;
#endif
  double period, T0, sin_i, cos_i, a, e, omega, r, ldcoeffs[4], mconst, chisqval, K, gamma, bimpact, N1, N2, Tc1, Tc2;
  int type, i;
  Tc1 = p[0]; Tc2 = p[1]; r = p[2]; a = p[3]; bimpact = p[4]; e = p[5]; omega = p[6]; mconst = p[7]; type = (int) p[8]; ldcoeffs[0] = p[9]; ldcoeffs[1] = p[10]; ldcoeffs[2] = p[11]; ldcoeffs[3] = p[12]; K = p[13]; gamma = p[14];

  N1 = p[15]; N2 = p[16];

  T0 = Tc1;
  period = (Tc2 - Tc1)/(N2 - N1);

  if(bimpact > 1.)
    bimpact = 1.;

  cos_i = bimpact*(1. + e*cos(omega))/(1. - e*e)/(a);
  if(cos_i > 1.)
    cos_i = 1.;
  else if(cos_i < -1.)
    cos_i = -1.;
  sin_i = sqrt(1. - cos_i*cos_i);

#ifdef PARALLEL
  RVdata = (_TypeRVdata *) userparams;
  NRVcurve = RVdata->NRVcurve;
  RVcurveJD = RVdata->RVcurveJD;
  RVcurveRV = RVdata->RVcurveRV;
  RVcurvesig = RVdata->RVcurvesig;
  dofitRV = RVdata->dofitRV;
#endif

  if(sin_i > 1.0)
    {
      p[4] = 0.0;
      sin_i = 1.0;
    }
  if(N > size_newlc)
    {
      if(!size_newlc)
	{
	  if((newlc = (double *) malloc(N * sizeof(double))) == NULL ||
	     (phase = (double *) malloc(N * sizeof(double))) == NULL)
	    {
	      fprintf(stderr,"Memory Allocation Error\n");
	      exit(3);
	    }
	}
      else
	{
	  if((newlc = (double *) realloc(newlc, N * sizeof(double))) == NULL ||
	     (phase = (double *) realloc(phase, N * sizeof(double))) == NULL)
	    {
	      fprintf(stderr,"Memory Allocation Error\n");
	      exit(3);
	    }
	}
      size_newlc = N;
    }
  /* Phase the light curve at the trial period/T0 */
  for(i=0;i<N;i++)
    {
      phase[i] = ((t[i] - T0)/period);
      phase[i] -= (double) floor(phase[i]);
      //phase[i] = ((t[i])/period) - (double) floor(((t[i])/period)) - T0;
      //while(phase[i] < 0.)
      //phase[i] += T0*100.*period;
      //phase[i] = phase[i] - (double) floor(phase[i]);
    }

  mandelagoltransitmodel(N, phase, newlc, type, ldcoeffs, sin_i, a, e, r, omega);

  chisqval = 0.;
  for(i=0;i<N;i++)
    {
      newlc[i] = -2.5*log(newlc[i])/log(10.0) + mconst;
      if(!isnan(mag[i]))
	chisqval += (mag[i] - newlc[i])*(mag[i] - newlc[i])/sig[i]/sig[i];
    }
  if(dofitRV)
    {
      if(NRVcurve > size_newRV)
	{
	  if(!size_newRV)
	    {
	      if((newRV = (double *) malloc(NRVcurve * sizeof(double))) == NULL ||
		 (phaseRV = (double *) malloc(NRVcurve * sizeof(double))) == NULL)
		{
		  fprintf(stderr,"Memory Allocation Error\n");
		  exit(3);
		}
	    }
	  else
	    {
	      if((newRV = (double *) realloc(newRV, NRVcurve * sizeof(double))) == NULL ||
		 (phaseRV = (double *) realloc(phaseRV, NRVcurve * sizeof(double))) == NULL)
		{
		  fprintf(stderr,"Memory Allocation Error\n");
		  exit(3);
		}
	    }
	  size_newRV = NRVcurve;
	}
      /* Phase the light curve at the trial period/T0 */
      for(i=0;i<NRVcurve;i++)
	{
	  phaseRV[i] = ((RVcurveJD[i] - T0)/period);
	  phaseRV[i] -= (double) floor(phaseRV[i]);
	  //phaseRV[i] = ((RVcurveJD[i])/period) - (double) floor(((RVcurveJD[i])/period)) - T0;
	  //while(phaseRV[i] < 0.)
	  //phaseRV[i] += T0*100.*period;
	  //phaseRV[i] = phaseRV[i] - (double) floor(phaseRV[i]);
	}

      getRVcurve(NRVcurve, phaseRV, newRV, K, e, omega, gamma);
      for(i=0;i<NRVcurve;i++)
	{
	  if(!isnan(RVcurveRV[i]))
	    chisqval += (RVcurveRV[i] - newRV[i])*(RVcurveRV[i] - newRV[i])/RVcurvesig[i]/RVcurvesig[i];
	}
    }

#ifdef PARALLEL
  if(newlc != NULL) free(newlc);
  if(phase != NULL) free(phase);
  if(phaseRV != NULL) free(phaseRV);
  if(newRV != NULL) free(newRV);
#endif
  return chisqval;
}

/* For a given set of orbital parameters, this function returns the minimum value of sin_i that yields transits */
double getminsini(double a, double e, double omega, double p)
{
  double phi_c, meanAnomaly;
  double D2, D, E;

  phi_c = getphaseofconjunction (e, omega);
  meanAnomaly = 2.0 * M_PI * (phi_c);

  if(e > 0)
    E = eccentricAnomaly(meanAnomaly, e);
  else
    E = meanAnomaly;
  D2 = (SQR_(1.0 - (e * cos(E))) - SQR_((1. + p)/a))/(SQR_((cos(E) - e)*sin(omega) + sqrt(1.0 - e*e)*sin(E)*cos(omega)));
  D = sqrt(D2);

  return D;
}

#define CONVERGENCELIMIT 0.00001
#define INITIALSTEP 0.05

void fitmandelagoltransit_amoeba(int N, double *t, double *mag, double *sig, double *P, double *T0, double *r, double *a, double *inc, double *bimpact, double *e, double *omega, double *mconst, int type, double *ldcoeffs, int fitephem, int fitr, int fita, int fitinclterm, int fite, int fitomega, int fitmconst, int *fitldcoeffs, double *chi2_, int correctlc, int omodel, char *modelname, int fitRV, char *RVfilename, char *omodelRVcurve, double *K, double *gamma, int fitK, int fitgamma, int refititer, int ophcurve, char *ophcurvename, double phmin, double phmax, double phstep, int ojdcurve, char *ojdcurvename, double jdstep, char *modelvarname, _Variable *modelvar, int threadid)
{
  /* This function fits a Mandel and Agol 2002 transit model to a light curve with t, mag and sig; The parameters are:
     P = period, T0 = time of transit minimum, r = planet to star radius ratio, a = ratio of semi-major axis to star radius, sin_i = sin inclination angle, e = eccentricity, omega = argument of periastron, mconst = out of transit magnitude, type = 0 for quadratic limb darkening, 1 for non-linear limb-darkening. ldcoeffs should be ldcoeffs[4]. */

  int l, j, k, ma, ngood, nfunk, nvar, sinivaryflag, Pvaryflag, fititer;
  double **p, *y, *delmag, *phase, ftol, meanval1, meanval2, phase0_init, phase0_final, minsini, qtran, Tbaseline, sin_i, T1, T2, N1, N2;
  int *ia, amoeba_val, i;
  double (*func)(double *, int, int, double *, double *, double *, void *);
  FILE *outfile, *outfile2, *RVfile;
  char *line;
  size_t line_size = MAXLEN;
  double *phaseRVmodel, *RVmodel, dphase;
  int sizemodelRVcurve;
  FILE *omodelRVcurvefile;

  double *tout, ph, jdtmp;
  int Nphase;

#ifdef PARALLEL
  _TypeRVdata RVdata;
  int sizeRVcurvevec = 0;
  int NRVcurve = 0;
  double *RVcurveJD = NULL;
  double *RVcurveRV = NULL;
  double *RVcurvesig = NULL;
  int dofitRV;

#else
  static int sizeRVcurvevec = 0;
#endif

  func = &chisqmandelagoltransit;

  dofitRV = fitRV;
  if(fitRV)
    {
      line = malloc(line_size);
      NRVcurve = 0;
      if((RVfile = fopen(RVfilename,"r")) == NULL)
	error2(ERR_FILENOTFOUND,RVfilename);
      while(gnu_getline(&line,&line_size,RVfile) >= 0)
	NRVcurve++;
      rewind(RVfile);
      if(NRVcurve > sizeRVcurvevec)
	{
	  if(!sizeRVcurvevec)
	    {
	      if((RVcurveJD = (double *) malloc(NRVcurve * sizeof(double))) == NULL ||
		 (RVcurveRV = (double *) malloc(NRVcurve * sizeof(double))) == NULL ||
		 (RVcurvesig = (double *) malloc(NRVcurve * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  else
	    {
	      if((RVcurveJD = (double *) realloc(RVcurveJD, NRVcurve * sizeof(double))) == NULL ||
		 (RVcurveRV = (double *) realloc(RVcurveRV, NRVcurve * sizeof(double))) == NULL ||
		 (RVcurvesig = (double *) realloc(RVcurvesig, NRVcurve * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	  sizeRVcurvevec = NRVcurve;
	}
      NRVcurve = 0;
      while(gnu_getline(&line,&line_size,RVfile) >= 0)
	{
	  sscanf(line,"%lf %lf %lf", &RVcurveJD[NRVcurve], &RVcurveRV[NRVcurve], &RVcurvesig[NRVcurve]);
	  NRVcurve++;
	}
      fclose(RVfile);
      free(line);
    }
  else
    {
      fitK = 0;
      fitgamma = 0;
    }

#ifdef PARALLEL
  RVdata.NRVcurve = NRVcurve;
  RVdata.RVcurveJD = RVcurveJD;
  RVdata.RVcurveRV = RVcurveRV;
  RVdata.RVcurvesig = RVcurvesig;
  RVdata.dofitRV = dofitRV;
#endif

  ma = 17;

  if((*P <= 0))
    {
      *chi2_ = 9999999.;
      return;
    }

  if((ia = (int *) malloc(ma * sizeof(int))) == NULL)
    error(ERR_MEMALLOC);

  if(fitephem) ia[0] = 1; else ia[0] = 0;
  if(fitephem) ia[1] = 1; else ia[1] = 0;
  if(fitr) ia[2] = 1; else ia[2] = 0;
  if(fita) ia[3] = 1; else ia[3] = 0;
  if(fitinclterm) ia[4] = 1; else ia[4] = 0;
  if(fite) ia[5] = 1; else ia[5] = 0;
  if(fitomega) ia[6] = 1; else ia[6] = 0;
  if(fitmconst) ia[7] = 1; else ia[7] = 0;
  ia[8] = 0;
  if(fitldcoeffs[0]) ia[9] = 1; else ia[9] = 0;
  if(fitldcoeffs[1]) ia[10] = 1; else ia[10] = 0;
  if(fitldcoeffs[2]) ia[11] = 1; else ia[11] = 0;
  if(fitldcoeffs[3]) ia[12] = 1; else ia[12] = 0;
  if(fitK) ia[13] = 1; else ia[13] = 0;
  if(fitgamma) ia[14] = 1; else ia[14] = 0;
  ia[15] = 0;
  ia[16] = 0;

  for(nvar=0,j=0;j<ma;j++)
    if(ia[j])
      nvar++;

  if((p = (double **) malloc((nvar + 1)* sizeof(double *))) == NULL ||
     (y = (double *) malloc((nvar + 1) * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  for(j=0;j<nvar+1;j++)
    if((p[j] = (double *) malloc(ma * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);

  /* If mconst is less than zero, set it equal to the average magnitude of the light curve */
  if(*mconst < 0)
    {
      for(j=0, l = 0; j < N; j++)
	{
	  if(!isnan(mag[j]))
	    {
	      l++;
	      (*mconst) += mag[j];
	    }
	}
      (*mconst) /= (double) l;
    }

  if(correctlc || omodel || ophcurve || ojdcurve || modelvarname != NULL)
    {
      for(j=0,meanval1=0.0,meanval2 = 0.0; j<N;j++)
	{
	  if(!isnan(mag[j]))
	    {
	      meanval1 += mag[j]/(sig[j]*sig[j]);
	      meanval2 += 1./(sig[j]*sig[j]);
	    }
	  meanval1 /= meanval2;
	}
    }

  /* Initialize the trial guesses */
  /* Note we'll treat the phase of point zero as the fitting variable for internal purposes */
  Tbaseline = t[N-1] - t[0];
  for(fititer = 0; fititer < 1; fititer++) {
    p[0][15] = 0.;
    p[0][16] = floor(Tbaseline / (*P));
    if(p[0][16] <= 0.)
      p[0][16] = 1.;
    p[0][0] = (*T0);
    p[0][1] = (*T0) + (*P)*p[0][16];
    phase0_init = (*T0)/(*P) - (double) floor((*T0)/(*P));
    //p[0][0] = *P; p[0][1] = phase0_init;
    p[0][2] = *r; p[0][3] = *a; p[0][4] = *bimpact; p[0][5] = *e; p[0][6] = *omega; p[0][7] = *mconst; p[0][8] = (double) type; p[0][9] = ldcoeffs[0]; p[0][10] = ldcoeffs[1]; p[0][11] = ldcoeffs[2]; p[0][12] = ldcoeffs[3];
    p[0][13] = *K; p[0][14] = *gamma;


    for(k=1;k<=nvar;k++)
      {
	sinivaryflag = 0;
	Pvaryflag = 0;
	for(l=0,j=0;j<ma;j++)
	  {
	    if(ia[j])
	      {
		if(l == k -1)
		  {
		    if(j != 1)
		      {
			if(p[0][j] != 0.)
			  p[k][j] = p[0][j]*(1 + INITIALSTEP);
			else
			  p[k][j] = INITIALSTEP;
			/* Make sure we don't try a sin_i value > 1.0, or a sin_i value that does not yield transits */
			if(j == 4)
			  {
			    /* We'll deal with sini later if we're varying it. */
			    sinivaryflag = 1;
			  }
			else if(j == 0 || j == 1)
			  {
			    /* We'll deal with period later as well if we're varying it */
			    qtran = 1./(*a)/M_PI;
			    p[k][j] = p[0][j] + qtran*(*P)*(*P)/Tbaseline;
			    Pvaryflag = 1;
			  }
		      }
		    else
		      p[k][j] = floor(p[0][j]) + (p[0][j] - floor(p[0][j]))*(1 + INITIALSTEP);

		  }
		else
		  p[k][j] = p[0][j];
		l++;
	      }
	    else
	      p[k][j] = p[0][j];
	  }
	if(sinivaryflag)
	  {
	    if(p[0][4] >= 1.)
	      p[0][4] = 1.;
	    else if(p[0][4] <= -1.)
	      p[0][4] = -1.;
	    if(p[0][4] >= 1. - INITIALSTEP) {
	      p[k][4] = p[0][4] - INITIALSTEP;
	    }
	    else {
	      p[k][4] = p[0][4] + INITIALSTEP;
	    }

	  }
      }
    for(k=0;k<nvar+1;k++)
#ifdef PARALLEL
      y[k] = (*func)(p[k],ma,N,t,mag,sig,(void *) (&RVdata));
#else
      y[k] = (*func)(p[k],ma,N,t,mag,sig, NULL);
#endif

      ftol = CONVERGENCELIMIT;

      /* Run amoeba if there are any variable terms */
      if(nvar > 0)
#ifdef PARALLEL
	amoeba_val = amoeba(p, y, ia, ma, ftol, func, &nfunk, 0, N, t, mag, sig, (void *) (&RVdata));
#else
      amoeba_val = amoeba(p, y, ia, ma, ftol, func, &nfunk, 0, N, t, mag, sig, NULL);
#endif
      else amoeba_val = 0;

  /* If amoeba didn't converge, then write out a garbage model */
      if(amoeba_val)
	{
	  *P = -1.; *T0 = -1.; *r = -1.; *a = -1.; *bimpact = -1.; *inc = -1.; *e = -1.; *omega = -1.; *mconst = -1.; ldcoeffs[0] = -1.; ldcoeffs[1] = -1.; ldcoeffs[2] = -1.; ldcoeffs[3] = -1.; *K = -1.; *gamma = -1.;
	  break;
	}
      else
	{
	  /* Find the minimum among the vertices */
	  k = 0;
	  for(j=0;j<=nvar;j++)
	    if(y[j] < y[k])
	      k = j;

	  T1 = p[k][0]; T2 = p[k][1];
	  N1 = p[k][15]; N2 = p[k][16];
	  *T0 = T1;
	  *P = (T2 - T1)/(N2 - N1);
	  //phase0_final = p[k][1];
	  //*T0 = (*T0 + (phase0_final - phase0_init)*p[k][0]);
	  //*P = p[k][0];
	  *r = p[k][2]; *a = p[k][3]; *bimpact = p[k][4]; *e = p[k][5]; *omega = p[k][6]; *mconst = p[k][7]; ldcoeffs[0] = p[k][9]; ldcoeffs[1] = p[k][10]; ldcoeffs[2] = p[k][11]; ldcoeffs[3] = p[k][12]; *K = p[k][13]; *gamma = p[k][14];
	  *inc = (*bimpact)*(1. + *e)*cos(*omega)/(1. - (*e)*(*e))/(*a);
	  *inc = 180.*acos((*inc))/M_PI;
	}
  }

  /* Remove the signal if we're doing that */
  if(correctlc || omodel || modelvarname != NULL)
    {
      if((delmag = (double *) malloc(N * sizeof(double))) == NULL ||
	 (phase = (double *) malloc(N * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      if(!amoeba_val)
	{
	  for(i=0;i<N;i++)
	    {
	      phase[i] = ((t[i]-(*T0))/(*P)) - (double) floor(((t[i]-(*T0))/(*P)));
	    }

	  sin_i = sin((*inc)*M_PI/180.);
	  mandelagoltransitmodel(N, phase, delmag, type, ldcoeffs, sin_i, *a, *e, *r, *omega);

	}

      if(omodel)
	{
	  if(!strncmp(modelname,"-",1) && strlen(modelname) == 1)
	    {
	      outfile = stdout;
	    }
	  else
	    {
	      if((outfile = fopen(modelname,"w")) == NULL)
		error2(ERR_CANNOTWRITE,modelname);
	    }
	  if(!amoeba_val)
	    {
	      for(j=0;j<N;j++)
		if(!isnan(mag[j]))
		  {
		    fprintf(outfile,"%f %f %f %f %f\n",t[j],mag[j],-2.5*log(delmag[j])/log(10.0) + *mconst,sig[j], phase[j]);
		  }
	    }
	  else
	    {
	      for(j=0;j<N;j++)
		if(!isnan(mag[j]))
		  fprintf(outfile,"%f %f 0. %f\n",t[j],mag[j],sig[j]);
	    }
	  if(outfile != stdout)
	    fclose(outfile);
	}
      if(modelvarname != NULL) {
	if(!amoeba_val)
	  {
	    for(j=0;j<N;j++)
	      (*((double ***) modelvar->dataptr))[threadid][j] = -2.5*log(delmag[j])/log(10.0) + *mconst;
	  }
	else
	  {
	    for(j=0;j<N;j++)
	      (*((double ***) modelvar->dataptr))[threadid][j] = 0.0;
	  }
      }
      if(correctlc && !amoeba_val)
	{
	  for(j=0;j<N;j++)
	    if(!isnan(mag[j]))
	      { mag[j] = mag[j] + 2.5*log(delmag[j])/log(10.0) - *mconst + meanval1; }
	}
      free(delmag);
      free(phase);
    }
  /* Output the phase curve if asked to */
    if(ophcurve)
    {
      if((outfile2 = fopen(ophcurvename,"w")) == NULL)
	error2(ERR_CANNOTWRITE,ophcurvename);

      fprintf(outfile2,"#Phase Mag_model\n");
      Nphase = ceil((phmax - phmin)/phstep)+1;

      if((delmag = (double *) malloc(Nphase * sizeof(double))) == NULL ||
	 (phase = (double *) malloc(Nphase * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0, ph=phmin;i<Nphase && ph <= phmax;i++, ph += phstep)
	{
	  if(ph < 0) {
	    phase[i] = ph + ceil(-ph);
	  } else if(ph > 1) {
	    phase[i] = ph - ceil(ph - 1.0);
	  } else {
	    phase[i] = ph;
	  }
	}
      Nphase = i;
      if(!amoeba_val)
	{
	  sin_i = sin((*inc)*M_PI/180.);
	  mandelagoltransitmodel(Nphase, phase, delmag, type, ldcoeffs, sin_i, *a, *e, *r, *omega);
	}
      else {
	for(i=0; i < Nphase; i++) delmag[i] = 1.;
      }
      for(i=0; i < Nphase; i++) {
	if((phase[i] < phmin && i > 0) || (phase[i] < phmin - JDTOL && i == 0)) {
	  phase[i] = phase[i] + ceil(phmin - phase[i]);
	} else if(i == 0 && phase[i] == phmin+1.0) {
	  phase[i] = phmin;
	} else if((phase[i] > phmax && i < (Nphase - 1)) || (phase[i] > phmax + JDTOL && i == Nphase - 1)) {
	  phase[i] = phase[i] - ceil(phase[i] - phmax);
	} else if(i == Nphase - 1 && phase[i] == phmax-1.0)
	  phase[i] = phmax;
	fprintf(outfile2,"%.17g %.17g\n", phase[i], -2.5*log(delmag[i])/log(10.0) + *mconst);
      }
      fclose(outfile2);
      free(delmag);
      free(phase);
    }

  /* Output the JD curve if asked to */
    if(ojdcurve)
    {
      if((outfile2 = fopen(ojdcurvename,"w")) == NULL)
	error2(ERR_CANNOTWRITE,ojdcurvename);

      fprintf(outfile2,"#Time Mag_model Phase\n");
      Nphase = ceil((t[N-1] - t[0])/jdstep)+1;
      if((delmag = (double *) malloc(Nphase * sizeof(double))) == NULL ||
	 (phase = (double *) malloc(Nphase * sizeof(double))) == NULL ||
	 (tout = (double *) malloc(Nphase * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0, jdtmp=t[0];i<Nphase && jdtmp <= t[N-1];i++, jdtmp += jdstep)
	{
	  tout[i] = jdtmp;
	  if(!amoeba_val)
	    {
	      ph = ((tout[i]-(*T0))/(*P)) - (double) floor(((tout[i]-(*T0))/(*P)));
	    }
	  else
	    ph = 0.0;
	  phase[i] = ph;
	}
      Nphase = i;
      if(!amoeba_val)
	{
	  sin_i = sin((*inc)*M_PI/180.);
	  mandelagoltransitmodel(Nphase, phase, delmag, type, ldcoeffs, sin_i, *a, *e, *r, *omega);
	}
      else {
	for(i=0; i < Nphase; i++) delmag[i] = 1.;
      }
      for(i=0; i < Nphase; i++) {
	fprintf(outfile2,"%.17g %.17g %.17g\n", tout[i], -2.5*log(delmag[i])/log(10.0) + *mconst, phase[i]);
      }
      fclose(outfile2);
      free(delmag);
      free(phase);
      free(tout);
    }

  if(fitRV && *P > 0.)
    {
      if((omodelRVcurvefile = fopen(omodelRVcurve,"w")) == NULL)
	error2(ERR_CANNOTWRITE,omodelRVcurve);

      sizemodelRVcurve = DEFAULTSIZEMODELRVCURVE;
      if((phaseRVmodel = (double *) malloc(sizemodelRVcurve * sizeof(double))) == NULL ||
	 (RVmodel = (double *) malloc(sizemodelRVcurve * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);

      dphase = 1./((double) sizemodelRVcurve);
      phaseRVmodel[0] = 0.;
      for(i=1;i<sizemodelRVcurve;i++)
	phaseRVmodel[i] = phaseRVmodel[i-1] + dphase;
      getRVcurve(sizemodelRVcurve, phaseRVmodel, RVmodel, *K, *e, *omega, *gamma);
      for(i=0; i < sizemodelRVcurve; i++)
	{
	  fprintf(omodelRVcurvefile, "%f %f\n",phaseRVmodel[i], RVmodel[i]);
	}
      fclose(omodelRVcurvefile);
      free(phaseRVmodel);
      free(RVmodel);
    }

  ngood = 0;
  for(j=0;j<N;j++)
    if(!isnan(mag[j]))
      ngood++;

  if(ngood > 9)
    *chi2_ = (y[k] / (ngood - nvar));
  else
    *chi2_ = -1.;
  free(ia);
  free(y);
  for(j=0;j<=nvar;j++)
    {
      free(p[j]);
    }
  free(p);

#ifdef PARALLEL
  if(RVcurveJD != NULL) free(RVcurveJD);
  if(RVcurveRV != NULL) free(RVcurveRV);
  if(RVcurvesig != NULL) free(RVcurvesig);
#endif

}




