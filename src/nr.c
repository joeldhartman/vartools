/* Various numerical recipes functions, many of these are scattered throughout
   the code, I will now try to collect these here. Modifications are made to
   index arrays and matrices from 0, and to use double precision. */
#include "commands.h"
#include "programdata.h"
#include "functions.h"

#ifndef EPS
#define EPS 1.0e-16
#endif
#ifndef FPMIN
#define FPMIN 1.0e-30
#endif
#ifndef MAXIT
#define MAXIT 10000
#endif
#ifndef XMIN
#define XMIN 2.0
#endif

#ifndef NUSE1
#define NUSE1 7
#endif
#ifndef NUSE2
#define NUSE2 8
#endif

/* Chebyshev evaluation: All arguments are input. c[0...m-1] is an array of
   Chebyshev coefficients, the first m elements of c output from chebft. */
double chebev(double a, double b, double c[], int m, double x)
{
  double d=0.0,dd=0.0,sv,y,y2;
  int j;

  if ((x-a)*(x-b) > 0.0) error2(ERR_NR,"x not in range in routine chebev");
  y2=2.0*(y=(2.0*x-a-b)/(b-a));
  for(j=m-1;j>=1;j--) {
    sv=d;
    d=y2*d-dd+c[j];
    dd=sv;
  }
  return y*d-dd+0.5*c[0];
}


/* Evaluates Gamma_1 and Gamma_2 by Chebyshev expanions for |x| <= 0.5; Also
   returns 1/Gamma(1+x) and 1/Gamma(1-x). */
void beschb(double x, double *gam1, double *gam2, double *gampl, double *gammi)
{
  double chebev(double a, double b, double c[], int m, double x);
  double xx;
  static double c1[] = {
    -1.142022680371168e0,6.5165112670737e-3,
    3.087090173086e-4,-3.4706269649e-6,6.9437664e-9,
    3.67795e-11,-1.356e-13};
  static double c2[] = {
    1.843740587300905e0,-7.68528408447867e-2,
    1.2719271366546e-3,-4.9717367042e-6,-3.31261198e-8,
    2.423096e-10,-1.702e-13,-1.49e-15};
  
  xx = 8.0*x*x-1.0;
  *gam1=chebev(-1.0,1.0,c1,NUSE1,xx);
  *gam2=chebev(-1.0,1.0,c2,NUSE2,xx);
  *gampl=*gam2-x*(*gam1);
  *gammi=*gam2+x*(*gam1);
}

/* Compute Modified Bessel function of fractional order xnu.
   ri = I_nu; rk = K_nu; Their derivatives are rip = I_prime_nu and
   rkp = K_prime_nu; We must have x positive, and xnu >= 0. The relative
   accuracy is within one or two sig. digits of EPS. FPMIN is a number close
   to the machine's smallest floating point number. */
void bessik(double x, double xnu, double *ri, double *rk, double *rip, double *rkp)
{
	void beschb(double x, double *gam1, double *gam2, double *gampl,
		double *gammi);
	int i,l,nl;
	double a,a1,b,c,d,del,del1,delh,dels,e,f,fact,fact2,ff,gam1,gam2,
		gammi,gampl,h,p,pimu,q,q1,q2,qnew,ril,ril1,rimu,rip1,ripl,
		ritemp,rk1,rkmu,rkmup,rktemp,s,sum,sum1,x2,xi,xi2,xmu,xmu2;

	if (x <= 0.0 || xnu < 0.0) error2(ERR_NR,"bad arguments in bessik");
	nl=(int)(xnu+0.5);
	xmu=xnu-nl;
	xmu2=xmu*xmu;
	xi=1.0/x;
	xi2=2.0*xi;
	h=xnu*xi;
	if (h < FPMIN) h=FPMIN;
	b=xi2*xnu;
	d=0.0;
	c=h;
	for (i=1;i<=MAXIT;i++) {
		b += xi2;
		d=1.0/(b+d);
		c=b+1.0/c;
		del=c*d;
		h=del*h;
		if (fabs(del-1.0) < EPS) break;
	}
	if (i > MAXIT) error2(ERR_NR,"x too large in bessik; try asymptotic expansion");
	ril=FPMIN;
	ripl=h*ril;
	ril1=ril;
	rip1=ripl;
	fact=xnu*xi;
	for (l=nl;l>=1;l--) {
		ritemp=fact*ril+ripl;
		fact -= xi;
		ripl=fact*ritemp+ril;
		ril=ritemp;
	}
	f=ripl/ril;
	if (x < XMIN) {
		x2=0.5*x;
		pimu=M_PI*xmu;
		fact = (fabs(pimu) < EPS ? 1.0 : pimu/sin(pimu));
		d = -log(x2);
		e=xmu*d;
		fact2 = (fabs(e) < EPS ? 1.0 : sinh(e)/e);
		beschb(xmu,&gam1,&gam2,&gampl,&gammi);
		ff=fact*(gam1*cosh(e)+gam2*fact2*d);
		sum=ff;
		e=exp(e);
		p=0.5*e/gampl;
		q=0.5/(e*gammi);
		c=1.0;
		d=x2*x2;
		sum1=p;
		for (i=1;i<=MAXIT;i++) {
			ff=(i*ff+p+q)/(i*i-xmu2);
			c *= (d/i);
			p /= (i-xmu);
			q /= (i+xmu);
			del=c*ff;
			sum += del;
			del1=c*(p-i*ff);
			sum1 += del1;
			if (fabs(del) < fabs(sum)*EPS) break;
		}
		if (i > MAXIT) error2(ERR_NR,"bessk series failed to converge");
		rkmu=sum;
		rk1=sum1*xi2;
	} else {
		b=2.0*(1.0+x);
		d=1.0/b;
		h=delh=d;
		q1=0.0;
		q2=1.0;
		a1=0.25-xmu2;
		q=c=a1;
		a = -a1;
		s=1.0+q*delh;
		for (i=2;i<=MAXIT;i++) {
			a -= 2*(i-1);
			c = -a*c/i;
			qnew=(q1-b*q2)/a;
			q1=q2;
			q2=qnew;
			q += c*qnew;
			b += 2.0;
			d=1.0/(b+a*d);
			delh=(b*d-1.0)*delh;
			h += delh;
			dels=q*delh;
			s += dels;
			if (fabs(dels/s) < EPS) break;
		}
		if (i > MAXIT) error2(ERR_NR,"bessik: failure to converge in cf2");
		h=a1*h;
		rkmu=sqrt(M_PI/(2.0*x))*exp(-x)/s;
		rk1=rkmu*(xmu+x+0.5-h)*xi;
	}
	rkmup=xmu*xi*rkmu-rk1;
	rimu=xi/(f*rkmu-rkmup);
	*ri=(rimu*ril1)/ril;
	*rip=(rimu*rip1)/ril;
	for (i=1;i<=nl;i++) {
		rktemp=(xmu+i)*xi2*rk1+rkmu;
		rkmu=rk1;
		rk1=rktemp;
	}
	*rk=rkmu;
	*rkp=xnu*xi*rkmu-rk1;
}
#undef EPS
#undef FPMIN
#undef MAXIT
#undef XMIN

double gammln(double xx)
{
  double x,y,tmp,ser;
  double cof[6]={76.18009172947146,-86.50532032941677,24.01409824083091,-1.231739572450155,0.1208650973866179e-2,-0.5395239384953e-5};
  int j;
  
  y=x=xx;
  tmp = x + 5.5;
  tmp -= (x+0.5)*log(tmp);
  ser=1.000000000190015;
  for(j=0;j<=5;j++) ser += cof[j]/++y;
  return -tmp+log(2.5066282746310005*ser/x);
}

#ifndef MAXIT
#define MAXIT 1000
#endif
#ifndef EPS
#define EPS 3.0e-7
#endif
#ifndef FPMIN
#define FPMIN 1.0e-30
#endif

double betacf(double a, double b, double x)
{
  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;
  qab = a+b;
  qap = a+1.0;
  qam = a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for(m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if(fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d = 1.0+aa*d;
    if(fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if(fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if(fabs(del-1.0) < EPS) break;
  }
  if(m > MAXIT) {
    error2(ERR_NR,"a or b too big, or MAXIT too small in betacf\n\n");
    exit(3);
  }
  return h;
}

double lbetacf(double a, double b, double x)
{
  int m, m2;
  double aa, c, d, del, h, qab, qam, qap;
  qab = a+b;
  qap = a+1.0;
  qam = a-1.0;
  c=1.0;
  d=1.0-qab*x/qap;
  if (fabs(d) < FPMIN) d=FPMIN;
  d=1.0/d;
  h=d;
  for(m=1;m<=MAXIT;m++) {
    m2=2*m;
    aa=m*(b-m)*x/((qam+m2)*(a+m2));
    d=1.0+aa*d;
    if (fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if(fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    h *= d*c;
    aa = -(a+m)*(qab+m)*x/((a+m2)*(qap+m2));
    d = 1.0+aa*d;
    if(fabs(d) < FPMIN) d=FPMIN;
    c=1.0+aa/c;
    if(fabs(c) < FPMIN) c=FPMIN;
    d=1.0/d;
    del=d*c;
    h *= del;
    if(fabs(del-1.0) < EPS) break;
  }
  if(m > MAXIT) {
    error2(ERR_NR,"a or b too big, or MAXIT too small in betacf\n\n");
    exit(3);
  }
  return log(h);
}

double betai(double a, double b, double x)
{
  double bt;
  if( x < 0.0 || x > 1.0) return 0.;
  if(x == 0.0 || x == 1.0) bt = 0.0;
  else
    bt = exp(gammln(a+b) - gammln(a) - gammln(b) + a*log(x) + b*log(1.0-x));
  if (x < (a + 1.0)/(a+b+2.0))
    return bt*betacf(a,b,x)/a;
  else
    return 1.0-bt*betacf(b,a,1.0-x)/b;
}

double log1minusbetai(double a, double b, double x)
{
  double bt;
  if( x < 0.0 || x > 1.0) return 0.;
  if(x == 0.0 || x == 1.0) bt = 0.0;
  else
    bt = gammln(a+b) - gammln(a) - gammln(b) + a*log(x) + b*log(1.0-x);
  if (x < (a + 1.0)/(a+b+2.0))
    return log(1. - exp(bt)*betacf(a,b,x)/a);
  else
    return bt + lbetacf(b,a,1.0-x) - log(b);
}

#undef MAXIT
#undef EPS
#undef FPMIN

#ifndef CHOLDC_TINY
#define CHOLDC_TINY 1.e-99
#endif

#ifndef ITMAX_GFUNC
#define ITMAX_GFUNC 100000
#endif

#ifndef EPS_GFUNC
#define EPS_GFUNC 3.0e-9
#endif
/* Incomplete Gamma Function P(a,x) evaluated by series representation, 
   taken from NR */
void gser(double *gamser, double a, double x, double *gln)
{
  int n;
  double sum, del, ap;
  *gln=gammln(a);
  if(x <= 0) {
    if( x < 0.0) {
      fprintf(stderr,"x less than 0 in routine gser\n");
    }
    *gamser = 0.0;
    return;
  } else {
    ap = a;
    del = sum = 1.0/a;
    for(n=1;n<=ITMAX_GFUNC;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if(fabs(del) < fabs(sum)*EPS_GFUNC) {
	*gamser = sum*exp(-x+a*log(x)-(*gln));
	return;
      }
    }
    fprintf(stderr,"a too large, ITMAX too small in routine gser\n");
    return;
  }
}

void lngser(double *lngamser, double a, double x, double *gln)
{
  int n;
  double sum, del, ap;
  *gln=gammln(a);
  if(x <= 0) {
    if( x < 0.0) {
      fprintf(stderr,"x less than 0 in routine gser\n");
    }
    *lngamser = -sqrt(-1.0);
    return;
  } else {
    ap = a;
    del = sum = 1.0/a;
    for(n=1;n<=ITMAX_GFUNC;n++) {
      ++ap;
      del *= x/ap;
      sum += del;
      if(fabs(del) < fabs(sum)*EPS_GFUNC) {
	*lngamser = log(sum) - x + a*log(x) - (*gln);
	return;
      }
    }
    fprintf(stderr,"a too large, ITMAX too small in routine gser\n");
    return;
  }
}

/* Incomplete Gamma FUnction Q(a,x) evaluated by its continued fraction
   representation, taken from NR */
#define FPMIN_GFUNC 1.0e-30
void gcf(double *gammcf, double a, double x, double *gln)
{
  int i;
  double an, b, c, d, del, h;
  *gln = gammln(a);
  b = x+1.0-a;
  c = 1.0/FPMIN_GFUNC;
  d = 1.0/b;
  h = d;
  for(i=1; i <= ITMAX_GFUNC; i++) {
    an = -i*(i-a);
    b += 2.0;
    d = an*d+b;
    if (fabs(d) < FPMIN_GFUNC) d = FPMIN_GFUNC;
    c = b+an/c;
    if (fabs(c) < FPMIN_GFUNC) c = FPMIN_GFUNC;
    d = 1.0/d;
    del = d*c;
    h *= del;
    if (fabs(del-1.0) < EPS_GFUNC) break;
  }
  if(i > ITMAX_GFUNC) {
    fprintf(stderr,"a too large, ITMAX too small in gcf\n");
  }
  *gammcf = exp(-x+a*log(x)-(*gln))*h;
}

void lngcf(double *lngammcf, double a, double x, double *gln)
{
  int i;
  double an, b, c, d, del, h;
  *gln = gammln(a);
  b = x+1.0-a;
  c = 1.0/FPMIN_GFUNC;
  d = 1.0/b;
  h = d;
  for(i=1; i <= ITMAX_GFUNC; i++) {
    an = -i*(i-a);
    b += 2.0;
    d = an*d+b;
    if (fabs(d) < FPMIN_GFUNC) d = FPMIN_GFUNC;
    c = b+an/c;
    if (fabs(c) < FPMIN_GFUNC) c = FPMIN_GFUNC;
    d = 1.0/d;
    del = d*c;
    h *= del;
    if (fabs(del-1.0) < EPS_GFUNC) break;
  }
  if(i > ITMAX_GFUNC) {
    fprintf(stderr,"a too large, ITMAX too small in gcf\n");
  }
  *lngammcf = -x+a*log(x)-(*gln) + log(h);
}

/* Incomplete Gamma Function, taken from NR */
double gammp(double a, double x)
{
  double gamser, gammcf, gln;
  if(x < 0.0 || a <= 0.0) {
    fprintf(stderr,"Invalid arguments in routine gammp\n");
    exit(4);
  }
  if ( x < (a+1.0)) {
    gser(&gamser,a,x,&gln);
    return gamser;
  } else {
    gcf(&gammcf, a, x, &gln);
    return 1.0 - gammcf;
  }
}

double lngammp(double a, double x)
{
  double lngamser, lngammcf, gln;
  if(x < 0.0 || a <= 0.0) {
    fprintf(stderr,"Invalid arguments in routine gammp\n");
    exit(4);
  }
  if ( x < (a+1.0)) {
    lngser(&lngamser,a,x,&gln);
    return lngamser;
  } else {
    lngcf(&lngammcf, a, x, &gln);
    return log(1.0 - exp(lngammcf));
  }
}

/* Gaussian Error Function, taken from NR */
double erff_(double x)
{
  return x < 0.0 ? -gammp(0.5,x*x) : gammp(0.5,x*x);
}

/* Return the inverse error function by bisection */
double inverff(double p)
{
  double xlow, xhigh, xtry, ptry;
  if(p == 0) {
    return 0.0;
  }
  if(p <= -1.0) {
    return -sqrt(-1.0);
  }
  else if(p >= 1.0) {
    return sqrt(1.0);
  }
  if(p < 0) {
    xhigh = 0.0;
    xlow = -0.5;
    while(erff_(xlow) > p) {
      xhigh = xlow;
      xlow = (xlow - 1.0)/2.;
    }
    xtry = (xlow + xhigh)/2.;
  }
  else {
    xlow = 0.0;
    xhigh = 0.5;
    while(erff_(xhigh) < p) {
      xlow = xhigh;
      xhigh = (xhigh + 1.0)/2.;
    }
    xtry = (xlow + xhigh)/2.;
  }
  do {
    ptry = erff_(xtry);
    if(fabs(ptry - p) < EPS_GFUNC)
      return xtry;
    if(ptry < p) {
      xlow = xtry;
      xtry = (xlow + xhigh)/2.;
    }
    else {
      xhigh = xtry;
      xtry = (xlow + xhigh)/2.;
    }
  } while(1);
    
}

/* Sparse Cholesky Decomposition for near center diagonal matrix */
void choldc_sparse_neardiag(double **a, int N, int *Nvec, double *p) {
  int i,j,k;
  double sum;
  for(i=0; i < N; i++) {
    for(j=i; j < N; j++) {
      if(a[i][j] == 0.) {
	sum = 0;
	break;
      }
      for(sum=a[i][j],k=i-1;k >= 0; k--) {
	if(Nvec[k] < i || Nvec[k] < j) {break;}
	else {sum -= a[i][k]*a[j][k];}
      }
      if (i == j) {
	if(sum <= 0.) {
	  sum = CHOLDC_TINY;
	}
	p[i] = sqrt(sum);
      } else a[j][i] = sum/p[i];
    }
    Nvec[i] = j;
    if(j < N) a[j][i] = 0.0;
  }
}

/* Multiply the cholesky decomposition stored in a by the vector stored in b,
   result is stored in vector x */
void cholmult_sparse_neardiag(double **a, int N, int *Nvec, double *p, double *b, double *x)
{
  int i, k;
  double sum;
  for(i=N-1; i >= 0; i--) {
    for(sum=b[i]*p[i],k=i-1;k>=0;k--) {if(Nvec[k] < i) break; sum += a[i][k]*b[k];}
    x[i]=sum;
  }
}

void cholsl_sparse_neardiag(double **a, int N, int *Nvec, double *p, double *b, double *x)
{
  int i, k;
  double sum;
  for(i=0; i < N; i++) {
    for(sum=b[i],k=i-1;k>=0;k--) {if(Nvec[k] < i) break; sum -= a[i][k]*x[k];}
    x[i]=sum/p[i];
  }
  for(i=N-1; i>=0; i--) {
    for(sum=x[i],k=i+1;k<Nvec[i];k++) sum -= a[k][i]*x[k];
    x[i]=sum/p[i];
  }
}

void cholsl(double **a, int N, double *p, double *b, double *x)
{
  int i, k;
  double sum;
  for(i=0; i < N; i++) {
    for(sum=b[i],k=i-1;k>=0;k--) sum -= a[i][k]*x[k];
    x[i]=sum/p[i];
  }
  for(i=N-1; i>=0; i--) {
    for(sum=x[i],k=i+1;k<N;k++) sum -= a[k][i]*x[k];
    x[i]=sum/p[i];
  }
}

/* Cholesky Decomposition, taken from NR */
void choldc(double **a, int N) {
  double *p;
  p = (double *) malloc(N * sizeof(double));
  int i,j,k;
  double sum;
  for(i=0; i < N; i++) {
    for(j=i; j < N; j++) {
      for(sum=a[i][j],k=i-1;k >= 0; k--) sum -= a[i][k]*a[j][k];
      if (i == j) {
	if(sum <= 0.) {
	  sum = CHOLDC_TINY;
	}
	p[i] = sqrt(sum);
      } else a[j][i] = sum/p[i];
    }
  }

  /* Fill out the rest of the matrix */
  for(i=0; i < N; i++) {
    a[i][i] = p[i];
    for(j=i + 1; j < N; j++) {
      a[i][j] = a[j][i];
    }
  }

  free(p);
}

#ifndef SWAP
#define SWAP(a,b) temp=(a);(a)=(b);(b)=temp
#endif

int gaussj_nob(double **a, int n)
{
  /* Linear equation solution by Gauss-Jordan elimination. a[0...n-1][0....n-1] is the input matrix. b[0....n-1][0...m-1] is input containing the m right-hand side vectors. On output, a is replaced by its matrix inverse, and b is replaced by the corresponding set of solution vectors. 

This algorithm is taken from Press et al. 1992, it has been modified slightly to put it into the program vartools by J. Hartman

*/
  int *indxc, *indxr, *ipiv;
  int i, icol, irow, j, k, l, ll;
  double big, dum, pivinv, temp;
  
  if((indxc = (int *) malloc(n * sizeof(int))) == NULL ||
     (indxr = (int *) malloc(n * sizeof(int))) == NULL ||
     (ipiv = (int *) malloc(n * sizeof(int))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(3);
    }

  for(j=0;j<n;j++) ipiv[j] = 0;
  for(i=0;i<n;i++) {
    big = 0.0;
    for(j=0;j<n;j++)
      if(ipiv[j] != 1)
	for(k=0;k<n;k++) {
	  if(ipiv[k] == 0) {
	    if (fabs(a[j][k]) >= big) {
	      big = fabs(a[j][k]);
	      irow = j;
	      icol = k;
	    }
	  }
	}
    ++(ipiv[icol]);
    if(irow != icol) {
      for(l=0;l<n;l++) { SWAP(a[irow][l],a[icol][l]); }
    }
    indxr[i] = irow;
    indxc[i] = icol;
    if(a[icol][icol] == 0.0)
      {
	free(ipiv);
	free(indxr);
	free(indxc);
	return(1);
      }
    pivinv=1.0/a[icol][icol];
    a[icol][icol] = 1.0;
    for(l=0;l<n;l++) a[icol][l] *= pivinv;
    for(ll=0;ll<n;ll++)
      if(ll != icol) {
	dum=a[ll][icol];
	a[ll][icol] = 0.0;
	for(l=0;l<n;l++) a[ll][l] -= a[icol][l]*dum;
      }
  }
  
  for(l=n-1;l>=0;l--) {
    if(indxr[l] != indxc[l])
      for(k=0;k<n;k++){ SWAP(a[k][indxr[l]],a[k][indxc[l]]); }
  }
  free(ipiv);
  free(indxr);
  free(indxc);
  return(0);
}
