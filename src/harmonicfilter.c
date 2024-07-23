/* This file contains functions which implement the -harmonicfilter
   command for VARTOOLS */
#include "commands.h"
#include "programdata.h"
#include "functions.h"

#ifndef TWOPI
#define TWOPI 6.28318530717958647692528676656
#endif

typedef struct {
  double r;
  double i;
} Complex;

void complex_product(Complex *a, Complex *b, Complex *ret) {
  /* a times b */
  ret->r = a->r*b->r - a->i*b->i;
  ret->i = a->r*b->i + a->i*b->r;
}

void complex_conjugate_product(Complex *a, Complex *b, Complex *ret) {
  /* a times b_conj */
  ret->r = a->r*b->r + a->i*b->i;
  ret->i = -a->r*b->i + a->i*b->r;
}

void complex_division(Complex *a, Complex *b, Complex *ret) {
  /* a / b */
  double amp;
  amp = b->r*b->r + a->i*b->i;
  ret->r = (a->r*b->r + a->i*b->i)/amp;
  ret->i = (-a->r*b->i + a->i*b->r)/amp;
}

void RAG_alg31(int m, int n, Complex *z, double *phase, double *w, 
	       Complex *g,
               Complex *c_prime, Complex *gamma, double *sigma) 
/* Executes algorithm 3.1 of Reichel, Ammar and Gragg, 1991,
   Mathematics of Computation, Volume 57, Number 195, Pages 273-289 

   Input: integers m and n;
          distinct nodes {z_k}_k=0^m-1 on the unit circle (i.e., z_k = exp(i*2.0*pi*t_k*freq))
          phases associated with z_k ( z_k = exp(i*phase[k]))
          sets of weights {w_k}_k=0^m-1 {w_k = 1./err_k} where sum (w_k^2) = 1.0;
          vector {g_k}_k=0^m-1 with g_k = (z_k^l)*f_k for observation f_k, and
              l = (n-1)/2;
   Returns: vector c' = {c^_j}_j=0^n-1 := Q*Dg; note that c' must be length n+1
            vector {gamma_j}_j=1^n-1  (assumed to be length n)
            vector {sigma_j}_j=0^n-1  (assumed to be length n)
                   
*/
{
  int i, j, k, j_prime;
  double *beta;
  Complex *alpha;
  Complex ztmp;
  double rho;
  Complex tau;

  Complex s1, s2, s3, s4;

  if((beta = (double *) malloc((n+1)*sizeof(double))) == NULL ||
     (alpha = (Complex *) malloc((n+1)*sizeof(Complex))) == NULL)
    error(ERR_MEMALLOC);

  sigma[0] = w[0];
  gamma[0].r = 1.; gamma[1].i = 0.;
  gamma[1].r = -z[0].r; gamma[1].i = -z[0].i;
  sigma[1] = 0.0;
  c_prime[0].r = w[0]*g[0].r; c_prime[0].i = w[0]*g[0].i;
  c_prime[1].r = 0.; c_prime[1].i = 0.;
  for(j=1; j < m; j++) {
    j_prime = j < n-1 ? j : n-1;
    for(k=j_prime+1; k >= 1; k--) {
      c_prime[k].r = c_prime[k-1].r;
      c_prime[k].i = c_prime[k-1].i;
    }
    c_prime[0].r = w[j]*g[j].r; c_prime[0].i = w[j]*g[j].i;
    beta[0] = sigma[0]; sigma[0] = sqrt(sigma[0]*sigma[0] + w[j]*w[j]);
    beta[0] = beta[0]/sigma[0]; alpha[0].r = -w[j]/sigma[0]; alpha[0].i = 0.;
    complex_conjugate_product(&(c_prime[0]),&(alpha[0]),&s1);
    s2.r = beta[0]*c_prime[1].r; s2.i = beta[0]*c_prime[1].i;
    s3.r = beta[0]*c_prime[0].r; s3.i = beta[0]*c_prime[0].i;
    complex_product(&(c_prime[1]),&(alpha[0]),&s4);
    c_prime[0].r = -s1.r + s2.r; c_prime[0].i = -s1.i + s2.i;
    c_prime[1].r = s3.r + s4.r; c_prime[1].i = s3.i + s4.i;
    if(j+1 < n) {
      complex_product(&(gamma[j]),&(z[j]),&(gamma[j+1]));
      gamma[j+1].r = -gamma[j+1].r;
      gamma[j+1].i = -gamma[j+1].i;
      sigma[j+1] = 0.0;
    }
    for(k=1; k <= j_prime; k++) {
      ztmp.r = cos((phase[j]*((double) (k-2))));
      ztmp.i = sin((phase[j]*((double) (k-2))));
      complex_conjugate_product(&ztmp,&(alpha[k-1]),&s1);
      complex_product(&(gamma[k]),&s1,&s2);
      tau.r = alpha[k-1].r + s2.r;
      tau.i = alpha[k-1].i + s2.i;
      rho = beta[k-1]*sqrt(sigma[k]*sigma[k] + tau.r*tau.r + tau.i*tau.i);
      complex_product(&tau,&(z[j]),&s1);
      alpha[k].r = beta[k-1]*s1.r/rho;
      alpha[k].i = beta[k-1]*s1.i/rho;
      beta[k] = beta[k-1]*sigma[k]/rho;
      ztmp.r = cos((-phase[j])*((double) (k-2)));
      ztmp.i = sin((-phase[j])*((double) (k-2)));
      complex_product(&alpha[k-1],&alpha[k-1],&s1);
      complex_product(&ztmp,&s1,&s2);
      gamma[k].r = beta[k-1]*beta[k-1]*gamma[k].r - s2.r;
      gamma[k].i = beta[k-1]*beta[k-1]*gamma[k].i - s2.i;
      sigma[k] = rho;

      complex_conjugate_product(&(c_prime[k]),&(alpha[k]),&s1);
      s2.r = beta[k]*c_prime[k+1].r; s2.i = beta[k]*c_prime[k+1].i;
      s3.r = beta[k]*c_prime[k].r; s3.i = beta[k]*c_prime[k].i;
      complex_product(&(c_prime[k+1]),&(alpha[k]),&s4);
      c_prime[k].r = -s1.r + s2.r; c_prime[k].i = -s1.i + s2.i;
      c_prime[k+1].r = s3.r + s4.r; c_prime[k+1].i = s3.i + s4.i;
    }
  }

  free(beta);
  free(alpha);
}

void complex_vector_conjugate_reversal(int n, Complex *r, Complex *r_out) {
  int i, j;
  for(i=n-1, j=0; i >= 0; i--, j++) {
    r_out[i].r = r[j].r;
    r_out[i].i = -r[j].i;
  }
}

void RAG_alg41(int n, Complex *gamma, double *sigma, Complex *b,
	       Complex *a)
/* Executes algorithm 4.1 of Reichel, Ammar and Gragg, 1991,
   Mathematics of Computation, Volume 57, Number 195, Pages 273-289 

   Input: integer n;
          vectors gamma, and sigma computed from RAG_alg31
          vector b is the same as vector c_prime computed from RAG_alg31
   Returns: vector a_j=0^n-1 = R^-1 b
*/
{
  int j, k;
  Complex *r0;
  Complex *r0_rev;
  Complex *r1;
  Complex s1;
  
  if((r0 = (Complex *) malloc(n * sizeof(Complex))) == NULL ||
     (r1 = (Complex *) malloc(n * sizeof(Complex))) == NULL ||
     (r0_rev = (Complex *) malloc(n * sizeof(Complex))) == NULL)
    error(ERR_MEMALLOC);

  for(j=0; j < n; j++) {
    a[j].r = 0.0; a[j].i = 0.0;
  }
  r0[0].r = 1./sigma[0];
  r0[0].i = 0.0;
  complex_product(&(r0[0]),&(b[0]),&(a[0]));
  for(j=1; j <= n-1; j++) {
    complex_vector_conjugate_reversal(j, r0, r0_rev);
    r1[j].r = 0.; r1[j].i = 0.;
    for(k=0; k < j; k++) {
      complex_product(&(r0_rev[k]),&(gamma[j]),&s1);
      r1[k].r = s1.r/sigma[j];
      r1[k].i = s1.i/sigma[j];
    }
    for(k=1; k <= j; k++) {
      r1[k].r += r0[k-1].r/sigma[j];
      r1[k].i += r0[k-1].i/sigma[j];
    }
    for(k=0; k <= j; k++) {
      complex_product(&(b[j]),&(r1[k]),&s1);
      a[k].r += s1.r;
      a[k].i += s1.i;
      fprintf(stderr,"%d %d %.17g %.17g\n", j, k, r1[k].r, r1[k].i);
      r0[k].r = r1[k].r;
      r0[k].i = r1[k].i;
    }
  }
  free(r0); free(r1); free(r0_rev);
}

void RAG_alg42(int m, int n, Complex *gamma, double *sigma, Complex *c_prime,
	       Complex *z, Complex *p)
/* Executes algorithm 4.2 of Reichel, Ammar and Gragg, 1991,
   Mathematics of Computation, Volume 57, Number 195, Pages 273-289 

   Input: integer m data points;
          integer n harmonics;
          vectors gamma, sigma and c_prime computed from RAG_alg31
          vector z containing the points at which to evaluate the polynomial
   Returns: vector p containing the polynomial evaluation
*/
{
  int i, j;
  Complex q, r;
  Complex s1_c;
  Complex s2_c;
  Complex s3_c;
  Complex s4_c;
  for(i=0; i < m; i++) {
    q.r = 1.0/sigma[0];
    q.i = 0.0;
    r.r = 1.0/sigma[0];
    r.i = 0.0;
    complex_product(&(c_prime[0]),&q,&(p[i]));
    for(j=1; j < n; j++) {
      complex_product(&(z[i]),&q,&s1_c);
      complex_product(&(gamma[j]),&r,&s2_c);
      complex_conjugate_product(&(z[i]),&(gamma[j]),&s3_c);
      complex_product(&s3_c,&q,&s4_c);
      q.r = (s1_c.r + s2_c.r)/sigma[j];
      q.i = (s1_c.i + s2_c.i)/sigma[j];
      r.r = (s4_c.r + r.r)/sigma[j];
      r.i = (s4_c.i + r.r)/sigma[j];
      complex_product(&(c_prime[j]),&q,&s1_c);
      p[i].r = p[i].r + s1_c.r;
      p[i].i = p[i].i + s1_c.i;
    }
  }
}

void fit_harmonic_series_RAG(int N, double *t, double *mag, double *err, double f0, int Nharm, double *avals, double *bvals)
/* Follows the procedure of 
   Reichel, Ammar and Gragg, 1991,
   Mathematics of Computation, Volume 57, Number 195, Pages 273-289 
   to fit a harmonic series to a light curve
*/
{
  double Nharm_orig_d;

  double Ntimesfreq, lcave, sumweight;

  int sizeNvecs = 0, sizeNharmvecs = 0, i;
  int Nharm_orig_i;

  int j, n;

  double var1, var2;

  double *weight = NULL, *phase = NULL, *sigma = NULL;

  Complex *z = NULL, *g = NULL, *gamma = NULL, *c_prime = NULL,
    *c_hat = NULL;

  if((z = (Complex *) malloc(N * sizeof(Complex))) == NULL ||
     (weight = (double *) malloc(N * sizeof(double))) == NULL ||
     (g = (Complex *) malloc(N * sizeof(Complex))) == NULL ||
     (phase = (double *) malloc(N * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  sumweight = 0.;
  for(i=0; i < N; i++) {
    weight[i] = 1./err[i]/err[i];
    sumweight += weight[i];
  }
  for(i=0; i < N; i++) {
    weight[i] = sqrt(weight[i]/sumweight);
  }
    
  var1 = 0.0; var2 = 0.0;
  for(i=0;i<N;i++)
    {
      var1 += (double) mag[i]*weight[i];
      var2 += (double) weight[i];
    }
  lcave = (double) (var1 / var2);


  sizeNharmvecs = 2*Nharm + 1;
  if((gamma = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL ||
     (sigma = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
     (c_prime = (Complex *) malloc((sizeNharmvecs+1)*sizeof(Complex))) == NULL ||
     (c_hat = (Complex *) malloc((sizeNharmvecs+1)*sizeof(Complex))) == NULL)
    error(ERR_MEMALLOC);

  /* change the frequency to angular frequency for the computation, also double Nharm to account for the complex polynomials being of double the order of the real space polynomials */
  
  f0 = TWOPI * f0;
  
  Nharm_orig_i = Nharm;
  Nharm_orig_d = (double) Nharm;

  for(i=0;i<N;i++)
    {
      z[i].r = cos(f0*t[i]);
      z[i].i = sin(f0*t[i]);
      phase[i] = f0*t[i];
      g[i].r = cos(Nharm_orig_d*f0*t[i])*(mag[i]-lcave);
      g[i].i = sin(Nharm_orig_d*f0*t[i])*(mag[i]-lcave);
    }

  RAG_alg31(N, sizeNharmvecs, z, phase, weight, g, c_prime, gamma, sigma);

  RAG_alg41(sizeNharmvecs, gamma, sigma, c_prime, c_hat);


  /* Now translate to the fourier coefficients */
  avals[0] = c_hat[Nharm_orig_i].r;
  for(i=1; i <= Nharm_orig_i; i++) {
    avals[i] = 2.0*c_hat[i+1].r;
    bvals[i] = -2.0*c_hat[i+1].i;
  }

  free(z); free(weight); free(g); free(phase);
  free(gamma); free(sigma); free(c_prime); free(c_hat);
}


void fit_harmonic_series_orthogonal_poly_complex2(int N, double *t, double *mag, double *err, double f0, int Nharm, double *avals, double *bvals) 
/* Uses the method of Schwarzenberg-Czerny 1996, ApJ, 460, L107 to fit a 
   harmonic series to the data via projection onto orthogonal polynomials

   This version uses the relations given in Jagels and Reichel 1993,
   Journal of Computational and Applied Mathematics, 46, 241 for the
   recurrence
*/
{
  double Nharm_orig_d;

  double var1, var2, var5;
  Complex var1_c, var2_c, var3_c, var4_c;
  Complex s1_c, s2_c, s3_c, s4_c, s5_c, s6_c;
  double delta;
  double sigma0, sigma1;
  Complex gamma0, gamma1;

  double Ntimesfreq, tmp1, tmp2, th_coeff1, th_coeff2, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, periodogram, periodogram_new, a, b, lcave, sumweight;

  int sizeNvecs = 0, sizeNharmvecs = 0, i;
  int Nharm_orig_i;

  int j, n;

  Complex *cn = NULL, *psi = NULL, *z = NULL, *phi = NULL, *phiconj = NULL,
    *zn = NULL, *alpha = NULL, *aN = NULL, *zcoeff_final = NULL,
    *aN2 = NULL, *aNconj = NULL, *aNconj2 = NULL;
 
  double *weight = NULL, *phi_norm = NULL;

  if((psi = (Complex *) malloc(N * sizeof(Complex))) == NULL ||
     (z = (Complex *) malloc(N * sizeof(Complex))) == NULL ||
     (zn = (Complex *) malloc(N * sizeof(Complex))) == NULL ||
     (phi = (Complex *) malloc(N * sizeof(Complex))) == NULL ||
     (phiconj = (Complex *) malloc(N * sizeof(Complex))) == NULL ||
     (weight = (double *) malloc(N * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  sumweight = 0.;
  for(i=0; i < N; i++) {
    weight[i] = 1./err[i]/err[i];
    sumweight += weight[i];
  }
  for(i=0; i < N; i++) {
    weight[i] = weight[i]/sumweight;
  }
    


  var1 = 0.0; var2 = 0.0;
  for(i=0;i<N;i++)
    {
      var1 += (double) mag[i]*weight[i];
      var2 += (double) weight[i];
    }
  lcave = (double) (var1 / var2);


  sizeNharmvecs = 2*Nharm + 1;
  if((cn = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL ||
     (alpha = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL ||
     (phi_norm = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
     (zcoeff_final = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL)
    error(ERR_MEMALLOC);
  
  if((aN = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL ||
     (aN2 = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL ||
     (aNconj = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL ||
     (aNconj2 = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL)
    error(ERR_MEMALLOC);
  
  for(i=0; i < sizeNharmvecs; i++) {
    aN[i].r = 0.0;
    aN[i].i = 0.0;
    aNconj[i].r = 0.0;
    aNconj[i].i = 0.0;
  }

  
  /* change the frequency to angular frequency for the computation, also double Nharm to account for the complex polynomials being of double the order of the real space polynomials */
  
  f0 = TWOPI * f0;

  Nharm_orig_i = Nharm;
  Nharm_orig_d = (double) Nharm;
  Nharm *= 2;

  th_coeff1 = (double) (N - Nharm - 1);
  th_coeff2 = (double) Nharm;

  Ntimesfreq = f0*Nharm_orig_d;
  for(i=0;i<N;i++)
    {
      z[i].r = cos(f0*t[i]);
      z[i].i = sin(f0*t[i]);
      zn[i].r = 1.;
      zn[i].i = 0.;
      psi[i].r = (mag[i]-lcave)*(cos(Ntimesfreq*t[i]));
      psi[i].i = (mag[i]-lcave)*(sin(Ntimesfreq*t[i]));
      phi[i].r = 1.;
      phi[i].i = 0.;
      phiconj[i].r = 1.;
      phiconj[i].i = 0.;
    }
  aN[0].r = 1.0;
  aN[0].i = 0.0;
  aNconj[0].r = 1.0;
  aNconj[0].i = 0.0;
  

  
  /* Now get the cn values using the recurrence algorithm */
  
  /* First get the n = 0 term */
  var1_c.r = 0.; var1_c.i = 0.;
  var2_c.r = 0.; var2_c.i = 0.;
  var3_c.r = 0.; var3_c.i = 0.;
  var4_c.r = 0.; var4_c.i = 0.;
  for(i=0;i<N;i++)
    {
      complex_conjugate_product(&(phi[i]),&(phi[i]),&s6_c);
      var4_c.r += weight[i]*s6_c.r;
      var4_c.i += weight[i]*s6_c.i;
    }
  var1 = (double) sqrt((double) var4_c.r);
  phi_norm[0] = var1;
  for(i=0;i<N;i++) {
    weight[i] = weight[i]/phi_norm[0];
  }

  cn[0].r = 0.; cn[0].i = 0.;
  for(i=0; i < N; i++) {
    complex_conjugate_product(&(psi[i]),&(phi[i]),&s1_c);
    cn[0].r += weight[i]*s1_c.r;
    cn[0].i += weight[i]*s1_c.i;
  }

  for(i=0; i <= Nharm; i++) {
    zcoeff_final[i].r = 0.;
    zcoeff_final[i].i = 0.;
  }
  complex_product(&(cn[0]),&(aN[0]),&s1_c);
  zcoeff_final[0].r = s1_c.r;
  zcoeff_final[0].i = s1_c.i;


  delta = 1.0;
  for(n=1; n <= Nharm; n++) {
    var1_c.r = 0.; var1_c.i = 0.;
    for(i=0; i < N; i++) {
      complex_product(&(z[i]),&(phi[i]),&s1_c);
      var1_c.r += weight[i]*s1_c.r;
      var1_c.i += weight[i]*s1_c.i;
    }
    gamma1.r = -var1_c.r/delta;
    gamma1.i = -var1_c.i/delta;
    sigma1 = sqrt(1.0 - (gamma1.r*gamma1.r + gamma1.i*gamma1.i));
    delta = delta*sigma1;
    cn[n].r = 0.;
    cn[n].i = 0.;
    for(i=0; i < N; i++) {
      complex_product(&gamma1, &(phiconj[i]), &s1_c);
      complex_product(&(z[i]), &(phi[i]), &(s2_c));
      complex_conjugate_product(&(s2_c),&gamma1,&s3_c);
      
      phi[i].r = (s2_c.r + s1_c.r)/sigma1;
      phi[i].i = (s2_c.i + s1_c.i)/sigma1;
      phiconj[i].r = (s3_c.r + phiconj[i].r)/sigma1;
      phiconj[i].i = (s3_c.i + phiconj[i].i)/sigma1;
      complex_conjugate_product(&(psi[i]),&(phi[i]),&s1_c);
      cn[n].r += s1_c.r;
      cn[n].i += s1_c.i;
    }
    for(i=0; i <= n; i++) {
      aN2[i].r = 0.; aN2[i].i = 0.;
      aNconj2[i].r = 0.; aNconj2[i].i = 0.;
    }
    for(i=0; i < n; i++) {
      complex_product(&gamma1,&(aNconj[i]),&s1_c);
      aN2[i].r += s1_c.r;
      aN2[i].i += s1_c.i;
      aN2[i+1].r += aN[i].r;
      aN2[i+1].i += aN[i].i;
      
      complex_conjugate_product(&(aN[i]),&gamma1,&s1_c);
      aNconj2[i+1].r += s1_c.r;
      aNconj2[i+1].i += s1_c.i;
      aNconj2[i].r += aNconj[i].r;
      aNconj2[i].i += aNconj[i].i;
    }
    for(i=0; i <= n; i++) {
      aN[i].r = aN2[i].r/sigma1;
      aN[i].i = aN2[i].i/sigma1;
      complex_product(&(cn[n]),&(aN[i]),&s1_c);
      zcoeff_final[i].r += s1_c.r;
      zcoeff_final[i].i += s1_c.i;
      aNconj[i].r = aNconj2[i].r/sigma1;
      aNconj[i].i = aNconj2[i].i/sigma1;
    }
  }
  
  /*  for(i=0; i <= Nharm; i++) {
    fprintf(stderr,"%d %.17g %.17g %.17g %.17g %.17g\n", i, zcoeff_final_r[i], zcoeff_final_i[i], phi_norm[i], c_r[i], c_i[i]);
    }*/
  
  avals[0] = zcoeff_final[Nharm_orig_i].r;
  bvals[0] = 0.0;

  for(i=1; i <= Nharm_orig_i; i++) {
    avals[i] = zcoeff_final[Nharm_orig_i+i].r+zcoeff_final[Nharm_orig_i-i].i;
    bvals[i] = zcoeff_final[Nharm_orig_i-i].i-zcoeff_final[Nharm_orig_i+i].i;
  }


  /*
  for(i=0; i <= Nharm; i++) {
    for(j=0; j <= (i < Nharm_orig_i ? i : Nharm_orig_i); j++) {
      avals[Nharm_orig_i-j] += 2.0*(aN_r[i][j]*c_r[i]/phi_norm[i]-aN_i[i][j]*c_i[i]/phi_norm[i]);
      bvals[Nharm_orig_i-j] += 2.0*(aN_r[i][j]*c_i[i]/phi_norm[i]+aN_i[i][j]*c_r[i]/phi_norm[i]);
    }
    }*/

  free(zcoeff_final);

  free(aN);
  free(aN2);
  free(aNconj);
  free(aNconj2);

  free(psi);
  free(z);
  free(zn);
  free(phiconj);
  free(phi);
  free(weight);
  free(cn);
  free(alpha);
  free(phi_norm);
}


void fit_harmonic_series_orthogonal_poly_complex(int N, double *t, double *mag, double *err, double f0, int Nharm, double *avals, double *bvals) 
/* Uses the method of Schwarzenberg-Czerny 1996, ApJ, 460, L107 to fit a 
   harmonic series to the data via projection onto orthogonal polynomials
*/
{
  double Nharm_orig_d;

  double var1, var2, var5;
  Complex var1_c, var2_c, var3_c, var4_c;
  Complex s1_c, s2_c, s3_c, s4_c, s5_c, s6_c;


  double Ntimesfreq, tmp1, tmp2, th_coeff1, th_coeff2, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, periodogram, periodogram_new, a, b, lcave, sumweight;

  int sizeNvecs = 0, sizeNharmvecs = 0, i;
  int Nharm_orig_i;

  int j, n;

  Complex *cn = NULL, *psi = NULL, *z = NULL, *phi = NULL, 
    *zn = NULL, *alpha = NULL, *aN = NULL, *zcoeff_final = NULL,
    *aN2 = NULL;
 
  double *weight = NULL, *phi_norm = NULL;

  if((psi = (Complex *) malloc(N * sizeof(Complex))) == NULL ||
     (z = (Complex *) malloc(N * sizeof(Complex))) == NULL ||
     (zn = (Complex *) malloc(N * sizeof(Complex))) == NULL ||
     (phi = (Complex *) malloc(N * sizeof(Complex))) == NULL ||
     (weight = (double *) malloc(N * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  sumweight = 0.;
  for(i=0; i < N; i++) {
    weight[i] = 1./err[i]/err[i];
    sumweight += weight[i];
  }
  for(i=0; i < N; i++) {
    weight[i] = weight[i]/sumweight;
  }
    


  var1 = 0.0; var2 = 0.0;
  for(i=0;i<N;i++)
    {
      var1 += (double) mag[i]*weight[i];
      var2 += (double) weight[i];
    }
  lcave = (double) (var1 / var2);


  sizeNharmvecs = 2*Nharm + 1;
  if((cn = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL ||
     (alpha = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL ||
     (phi_norm = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
     (zcoeff_final = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL)
    error(ERR_MEMALLOC);
  
  if((aN = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL ||
     (aN2 = (Complex *) malloc(sizeNharmvecs * sizeof(Complex))) == NULL)
    error(ERR_MEMALLOC);
  
  for(i=0; i < sizeNharmvecs; i++) {
    aN[i].r = 0.0;
    aN[i].i = 0.0;
  }

  
  /* change the frequency to angular frequency for the computation, also double Nharm to account for the complex polynomials being of double the order of the real space polynomials */
  
  f0 = TWOPI * f0;

  Nharm_orig_i = Nharm;
  Nharm_orig_d = (double) Nharm;
  Nharm *= 2;

  th_coeff1 = (double) (N - Nharm - 1);
  th_coeff2 = (double) Nharm;

  Ntimesfreq = f0*Nharm_orig_d;
  for(i=0;i<N;i++)
    {
      z[i].r = cos(f0*t[i]);
      z[i].i = sin(f0*t[i]);
      zn[i].r = 1.;
      zn[i].i = 0.;
      psi[i].r = (mag[i]-lcave)*(cos(Ntimesfreq*t[i]));
      psi[i].i = (mag[i]-lcave)*(sin(Ntimesfreq*t[i]));
      phi[i].r = 1.;
      phi[i].i = 0.;
    }
  
  /* Now get the cn values using the recurrence algorithm */
  
  /* First get the n = 0 term */
  var1_c.r = 0.; var1_c.i = 0.;
  var2_c.r = 0.; var2_c.i = 0.;
  var3_c.r = 0.; var3_c.i = 0.;
  var4_c.r = 0.; var4_c.i = 0.;
  for(i=0;i<N;i++)
    {
      complex_conjugate_product(&(phi[i]),&(phi[i]),&s6_c);
      var4_c.r += weight[i]*s6_c.r;
      var4_c.i += weight[i]*s6_c.i;
    }
  var1 = (double) sqrt((double) var4_c.r);
  phi_norm[0] = var1;
  for(i=0;i<N;i++) {
    phi[i].r = phi[i].r/phi_norm[0];
    phi[i].i = phi[i].i/phi_norm[0];
  }

  var1_c.r = 0.; var1_c.i = 0.;
  var2_c.r = 0.; var2_c.i = 0.;
  var3_c.r = 0.; var3_c.i = 0.;
  var4_c.r = 0.; var4_c.i = 0.;
  for(i=0;i<N;i++)
    {
      complex_product(&(z[i]),&(phi[i]),&s1_c);
      complex_conjugate_product(&s1_c,&(phi[i]),&s2_c);

      complex_conjugate_product(&(zn[i]),&(phi[i]),&s3_c);
      complex_conjugate_product(&s3_c,&(phi[i]),&s4_c);

      complex_conjugate_product(&(psi[i]),&(phi[i]),&s5_c);
      complex_conjugate_product(&(phi[i]),&(phi[i]),&s6_c);
      
      var1_c.r += weight[i]*s2_c.r;
      var1_c.i += weight[i]*s2_c.i;

      var2_c.r += weight[i]*s4_c.r;
      var2_c.i += weight[i]*s4_c.i;
      
      var3_c.r += weight[i]*s5_c.r;
      var3_c.i += weight[i]*s5_c.i;

      var4_c.r += weight[i]*s6_c.r;
      var4_c.i += weight[i]*s6_c.i;
    }
  
  complex_division(&var3_c,&var4_c,&(cn[0]));
  complex_division(&var1_c,&var2_c,&(alpha[0]));

  aN[0].r = 1.0/phi_norm[0];
  aN[0].i = 0.0;

  for(i=0; i <= Nharm; i++) {
    zcoeff_final[i].r = 0.;
    zcoeff_final[i].i = 0.;
  }
  complex_product(&(cn[0]),&(aN[0]),&(zcoeff_final[0]));
  
  /* Get the rest of the harmonics */
  for(n=1;n<=Nharm;n++)
    {
      complex_conjugate_product(&(alpha[n-1]),&(aN[n-1]),&s1_c);
      aN2[0].r = -s1_c.r;
      aN2[0].i = -s1_c.i;
      for(i=1; i <= n-1; i++) {
	complex_conjugate_product(&(alpha[n-1]),&(aN[n-1-i]),&s1_c);
	aN2[i].r = aN[i-1].r - s1_c.r;
	aN2[i].i = aN[i-1].i - s1_c.i;
      }
      aN2[n].r = aN[n-1].r;
      aN2[n].i = aN[n-1].i;

      for(i=0;i<N;i++)
	{
	  /* Get the new phi values */
	  complex_product(&(z[i]),&(phi[i]),&s1_c);
	  complex_conjugate_product(&(zn[i]),&(phi[i]),&s2_c);
	  complex_product(&(alpha[n-1]),&s2_c,&s3_c);
	  phi[i].r = s1_c.r - s3_c.r;
	  phi[i].i = s1_c.i - s3_c.i;
	  complex_product(&(z[i]),&(zn[i]),&s1_c);
	  zn[i].r = s1_c.r;
	  zn[i].i = s1_c.i;
	}
      var4_c.r = 0.; var4_c.i = 0.;
      for(i=0;i<N;i++)
	{
	  complex_conjugate_product(&(phi[i]),&(phi[i]),&s6_c);
	  var4_c.r += weight[i]*s6_c.r;
	  var4_c.i += weight[i]*s6_c.i;
	}
      var1 = (double) sqrt((double) var4_c.r);
      phi_norm[n] = var1;
      for(i=0;i<N;i++) {
	phi[i].r = phi[i].r/phi_norm[n];
	phi[i].i = phi[i].i/phi_norm[n];
      }

      for(i=0; i <= n; i++) {
	aN[i].r = aN2[i].r/phi_norm[n];
	aN[i].i = aN2[i].i/phi_norm[n];
	if(fabs(aN[i].r) < 1.0e-10) aN[i].r = 0.;
	if(fabs(aN[i].i) < 1.0e-10) aN[i].i = 0.;
	fprintf(stderr,"%d %d %.17g %.17g\n",n,i,aN[i].r,aN[i].i);
      }

      /* Calculate the new alpha_n and c values */
      var1_c.r = 0.; var1_c.i = 0.;
      var2_c.r = 0.; var2_c.i = 0.;
      var3_c.r = 0.; var3_c.i = 0.;
      var4_c.r = 0.; var4_c.i = 0.;
      for(i=0;i<N;i++)
	{
	  complex_product(&(z[i]),&(phi[i]),&s1_c);
	  complex_conjugate_product(&s1_c,&(phi[i]),&s2_c);

	  complex_conjugate_product(&(zn[i]),&(phi[i]),&s3_c);
	  complex_conjugate_product(&s3_c,&(phi[i]),&s4_c);
	  
	  complex_conjugate_product(&(psi[i]),&(phi[i]),&s5_c);
	  complex_conjugate_product(&(phi[i]),&(phi[i]),&s6_c);
	  
	  var1_c.r += weight[i]*s2_c.r;
	  var1_c.i += weight[i]*s2_c.i;
	  
	  var2_c.r += weight[i]*s4_c.r;
	  var2_c.i += weight[i]*s4_c.i;
	  
	  var3_c.r += weight[i]*s5_c.r;
	  var3_c.i += weight[i]*s5_c.i;
	  
	  var4_c.r += weight[i]*s6_c.r;
	  var4_c.i += weight[i]*s6_c.i;
	}
      var1 = (double) sqrt((double) var4_c.r);
      complex_division(&var3_c,&var4_c,&(cn[n]));
      complex_division(&var1_c,&var2_c,&(alpha[n]));
      for(i=0; i <= n; i++) {
	complex_product(&(cn[n]),&(aN[i]),&s1_c);
	zcoeff_final[i].r += s1_c.r;
	zcoeff_final[i].i += s1_c.i;
      }
    }
  
  /*  for(i=0; i <= Nharm; i++) {
    fprintf(stderr,"%d %.17g %.17g %.17g %.17g %.17g\n", i, zcoeff_final_r[i], zcoeff_final_i[i], phi_norm[i], c_r[i], c_i[i]);
    }*/
  
  avals[0] = zcoeff_final[Nharm_orig_i].r;
  bvals[0] = 0.0;

  for(i=1; i <= Nharm_orig_i; i++) {
    avals[i] = zcoeff_final[Nharm_orig_i+i].r+zcoeff_final[Nharm_orig_i-i].i;
    bvals[i] = zcoeff_final[Nharm_orig_i-i].i-zcoeff_final[Nharm_orig_i+i].i;
  }


  /*
  for(i=0; i <= Nharm; i++) {
    for(j=0; j <= (i < Nharm_orig_i ? i : Nharm_orig_i); j++) {
      avals[Nharm_orig_i-j] += 2.0*(aN_r[i][j]*c_r[i]/phi_norm[i]-aN_i[i][j]*c_i[i]/phi_norm[i]);
      bvals[Nharm_orig_i-j] += 2.0*(aN_r[i][j]*c_i[i]/phi_norm[i]+aN_i[i][j]*c_r[i]/phi_norm[i]);
    }
    }*/

  free(zcoeff_final);

  free(aN);
  free(aN2);

  free(psi);
  free(z);
  free(zn);
  free(phi);
  free(weight);
  free(cn);
  free(alpha);
  free(phi_norm);
}


void fit_harmonic_series_orthogonal_poly(int N, double *t, double *mag, double *err, double f0, int Nharm, double *avals, double *bvals) 
/* Uses the method of Schwarzenberg-Czerny 1996, ApJ, 460, L107 to fit a 
   harmonic series to the data via projection onto orthogonal polynomials
*/
{
  double Nharm_orig_d;

  double var1, var2, var1_r, var1_i, var2_r, var2_i, var3_r, var3_i, var4_r, var4_i, var5;
  double Ntimesfreq, tmp1, tmp2, th_coeff1, th_coeff2, s1, s2, s3, s4, s5, s6, s7, s8, s9, s10, s11, s12, periodogram, periodogram_new, a, b, lcave, sumweight;

  int sizeNvecs = 0, sizeNharmvecs = 0, i;
  int Nharm_orig_i;

  int j, n;

  double *c_r = NULL, *c_i = NULL, *psi_r = NULL, *psi_i = NULL, 
    *z_r = NULL, *z_i = NULL, *phi_r = NULL, *phi_i = NULL, *zn_r = NULL, 
    *zn_i = NULL, *alpha_r = NULL, *alpha_i = NULL, *weight = NULL,
    **aN_r = NULL, **aN_i = NULL, *phi_norm = NULL, *zcoeff_final_r = NULL,
    *zcoeff_final_i = NULL;

  if((psi_r = (double *) malloc(N * sizeof(double))) == NULL ||
     (psi_i = (double *) malloc(N * sizeof(double))) == NULL ||
     (z_r = (double *) malloc(N * sizeof(double))) == NULL ||
     (z_i = (double *) malloc(N * sizeof(double))) == NULL ||
     (zn_r = (double *) malloc(N * sizeof(double))) == NULL ||
     (zn_i = (double *) malloc(N * sizeof(double))) == NULL ||
     (phi_r = (double *) malloc(N * sizeof(double))) == NULL ||
     (phi_i = (double *) malloc(N * sizeof(double))) == NULL ||
     (weight = (double *) malloc(N * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  sumweight = 0.;
  for(i=0; i < N; i++) {
    weight[i] = 1./err[i]/err[i];
    sumweight += weight[i];
  }
  for(i=0; i < N; i++) {
    weight[i] = weight[i]/sumweight;
  }
    


  var1 = 0.0; var2 = 0.0;
  for(i=0;i<N;i++)
    {
      var1 += (double) mag[i]*weight[i];
      var2 += (double) weight[i];
    }
  lcave = (double) (var1 / var2);


  sizeNharmvecs = 2*Nharm + 1;
  if((c_r = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
     (c_i = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
     (alpha_r = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
     (alpha_i = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
     (phi_norm = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
     (zcoeff_final_r = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL ||
     (zcoeff_final_i = (double *) malloc(sizeNharmvecs * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);
  
  if((aN_r = (double **) malloc(sizeNharmvecs * sizeof(double *))) == NULL ||
     (aN_i = (double **) malloc(sizeNharmvecs * sizeof(double *))) == NULL)
    error(ERR_MEMALLOC);
  
  for(i=0; i < sizeNharmvecs; i++) {
    if((aN_r[i] = (double *) malloc((i+1) * sizeof(double))) == NULL ||
       (aN_i[i] = (double *) malloc((i+1) * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    for(j=0; j <= i; j++) {
      aN_r[i][j] = 0.0;
      aN_i[i][j] = 0.0;
    }
  }

  
  /* change the frequency to angular frequency for the computation, also double Nharm to account for the complex polynomials being of double the order of the real space polynomials */
  
  f0 = TWOPI * f0;

  Nharm_orig_i = Nharm;
  Nharm_orig_d = (double) Nharm;
  Nharm *= 2;

  th_coeff1 = (double) (N - Nharm - 1);
  th_coeff2 = (double) Nharm;

  Ntimesfreq = f0*Nharm_orig_d;
  for(i=0;i<N;i++)
    {
      z_r[i] = cos(f0*t[i]);
      z_i[i] = sin(f0*t[i]);
      zn_r[i] = 1.;
      zn_i[i] = 0.;
      psi_r[i] = (mag[i]-lcave)*(cos(Ntimesfreq*t[i]));
      psi_i[i] = (mag[i]-lcave)*(sin(Ntimesfreq*t[i]));
      phi_r[i] = 1.;
      phi_i[i] = 0.;
    }
  
  /* Now get the cn values using the recurrence algorithm */
  
  /* First get the n = 0 term */
  var1_r = 0.; var1_i = 0.;
  var2_r = 0.; var2_i = 0.;
  var3_r = 0.; var3_i = 0.;
  var4_r = 0.; var4_i = 0.;
  for(i=0;i<N;i++)
    {
      s1 = z_r[i]*phi_r[i];
      s2 = z_i[i]*phi_i[i];
      s3 = z_r[i]*phi_i[i];
      s4 = z_i[i]*phi_r[i];
      s5 = zn_r[i]*phi_r[i];
      s6 = zn_i[i]*phi_i[i];
      s7 = zn_r[i]*phi_i[i];
      s8 = zn_i[i]*phi_r[i];
      s9 = s1 - s2;
      s10 = s4 + s3;
      s11 = s5 + s6;
      s12 = s8 - s7;
      /*
	var1_r += (double) weight[i]*(s9*phi_r[i] + s10*phi_i[i]);
	var1_i += (double) weight[i]*(-s9*phi_i[i] + s10*phi_r[i]);
	var2_r += (double) weight[i]*(s11*phi_r[i] + s12*phi_i[i]);
	var2_i += (double) weight[i]*(-s11*phi_i[i] + s12*phi_r[i]);*/
      var1_r += (double) weight[i]*s9;
      var1_i += (double) weight[i]*s10;
      var2_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
      var3_r += (double) weight[i]*(psi_r[i]*phi_r[i] + psi_i[i]*phi_i[i]);
      var3_i += (double) weight[i]*(-psi_r[i]*phi_i[i] + psi_i[i]*phi_r[i]);
      var4_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
    }
  var1 = (double) sqrt((double) var4_r);
  phi_norm[0] = var1;
  c_r[0] = (double) var3_r / var1;
  c_i[0] = (double) var3_i / var1;
  alpha_r[0] = (double) (var1_r / (double) var2_r);
  alpha_i[0] = (double) (var1_i / (double) var2_r);
  
  aN_r[0][0] = 1.0;
  aN_i[0][0] = 0.0;
    
  /* Get the rest of the harmonics */
  for(n=1;n<=Nharm;n++)
    {
      aN_r[n][0] = -alpha_r[n-1]*aN_r[n-1][n-1]-alpha_i[n-1]*aN_i[n-1][n-1];
      aN_i[n][0] = -alpha_i[n-1]*aN_r[n-1][n-1]+alpha_r[n-1]*aN_i[n-1][n-1];
      for(i=1; i <= n-1; i++) {
	aN_r[n][i] = aN_r[n-1][i-1]-alpha_r[n-1]*aN_r[n-1][n-1-i]-alpha_i[n-1]*aN_i[n-1][n-1-i];
	aN_i[n][i] = aN_i[n-1][i-1]-alpha_i[n-1]*aN_r[n-1][n-1-i]+alpha_r[n-1]*aN_i[n-1][n-1-i];
      }
      aN_r[n][n] = aN_r[n-1][n-1];
      aN_i[n][n] = aN_i[n-1][n-1];

      for(i=0;i<N;i++)
	{
	  /* Get the new phi values */
	  s1 = zn_r[i]*phi_r[i];
	  s2 = zn_i[i]*phi_i[i];
	  s3 = zn_r[i]*phi_i[i];
	  s4 = zn_i[i]*phi_r[i];
	  s5 = s1 + s2;
	  s6 = s4 - s3;
	  tmp1 = z_r[i]*phi_r[i] - z_i[i]*phi_i[i] - (alpha_r[n-1]*s5 - alpha_i[n-1]*s6);
	  tmp2 = z_r[i]*phi_i[i] + z_i[i]*phi_r[i] - (alpha_r[n-1]*s6 + alpha_i[n-1]*s5);
	  phi_r[i] = tmp1;
	  phi_i[i] = tmp2;
	  /* Update zn */
	  tmp1 = zn_r[i]*z_r[i] - zn_i[i]*z_i[i];
	  tmp2 = zn_r[i]*z_i[i] + zn_i[i]*z_r[i];
	  zn_r[i] = tmp1;
	  zn_i[i] = tmp2;
	}
      /* Calculate the new alpha_n and c values */
      var1_r = 0.; var1_i = 0.;
      var2_r = 0.; var2_i = 0.;
      var3_r = 0.; var3_i = 0.;
      var4_r = 0.; var4_i = 0.;
      for(i=0;i<N;i++)
	{
	  s1 = z_r[i]*phi_r[i];
	  s2 = z_i[i]*phi_i[i];
	  s3 = z_r[i]*phi_i[i];
	  s4 = z_i[i]*phi_r[i];
	  s5 = zn_r[i]*phi_r[i];
	  s6 = zn_i[i]*phi_i[i];
	  s7 = zn_r[i]*phi_i[i];
	  s8 = zn_i[i]*phi_r[i];
	  s9 = s1 - s2;
	  s10 = s4 + s3;
	  s11 = s5 + s6;
	  s12 = s8 - s7;
	  /*var1_r += (double) weight[i]*(s9*phi_r[i] + s10*phi_i[i]);
	    var1_i += (double) weight[i]*(-s9*phi_i[i] + s10*phi_r[i]);
	    var2_r += (double) weight[i]*(s11*phi_r[i] + s12*phi_i[i]);
	    var2_i += (double) weight[i]*(-s11*phi_i[i] + s12*phi_r[i]);*/
	  var1_r += (double) weight[i]*s9;
	  var1_i += (double) weight[i]*s10;
	  var2_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
	  var3_r += (double) weight[i]*(psi_r[i]*phi_r[i] + psi_i[i]*phi_i[i]);
	  var3_i += (double) weight[i]*(-psi_r[i]*phi_i[i] + psi_i[i]*phi_r[i]);
	  var4_r += (double) weight[i]*(phi_r[i]*phi_r[i] + phi_i[i]*phi_i[i]);
	}
      var1 = (double) sqrt((double) var4_r);
      phi_norm[n] = var1;
      c_r[n] = (double) var3_r / var1;
      c_i[n] = (double) var3_i / var1;
      alpha_r[n] = (double) (var1_r / (double) var2_r);
      alpha_i[n] = (double) (var1_i / (double) var2_r);
    }
  for(i=0; i <= Nharm; i++) {
    zcoeff_final_r[i] = 0.;
    zcoeff_final_i[i] = 0.;
    for(j=i; j <= Nharm; j++) {
      /*
      zcoeff_final_r[i] += (c_r[j]*aN_r[j][i] - c_i[j]*aN_i[j][i])/phi_norm[j];
      zcoeff_final_i[i] += (c_r[j]*aN_i[j][i] + c_i[j]*aN_r[j][i])/phi_norm[j];
      */
      fprintf(stderr,"%d %d %.17g %.17g\n", i, j, aN_r[j][i]/phi_norm[j], aN_i[j][i]/phi_norm[j]);
      zcoeff_final_r[i] += (c_r[j]*aN_r[j][i] - c_i[j]*aN_i[j][i]);
      zcoeff_final_i[i] += (c_r[j]*aN_i[j][i] + c_i[j]*aN_r[j][i]);
    }
  }

  /*  for(i=0; i <= Nharm; i++) {
    fprintf(stderr,"%d %.17g %.17g %.17g %.17g %.17g\n", i, zcoeff_final_r[i], zcoeff_final_i[i], phi_norm[i], c_r[i], c_i[i]);
    }*/
  
  avals[0] = zcoeff_final_r[Nharm_orig_i];
  bvals[0] = 0.0;

  for(i=1; i <= Nharm_orig_i; i++) {
    avals[i] = zcoeff_final_r[Nharm_orig_i+i]+zcoeff_final_r[Nharm_orig_i-i];
    bvals[i] = zcoeff_final_i[Nharm_orig_i-i]-zcoeff_final_i[Nharm_orig_i+i];
  }


  /*
  for(i=0; i <= Nharm; i++) {
    for(j=0; j <= (i < Nharm_orig_i ? i : Nharm_orig_i); j++) {
      avals[Nharm_orig_i-j] += 2.0*(aN_r[i][j]*c_r[i]/phi_norm[i]-aN_i[i][j]*c_i[i]/phi_norm[i]);
      bvals[Nharm_orig_i-j] += 2.0*(aN_r[i][j]*c_i[i]/phi_norm[i]+aN_i[i][j]*c_r[i]/phi_norm[i]);
    }
    }*/

  free(zcoeff_final_r);
  free(zcoeff_final_i);

  for(i=0; i < sizeNharmvecs; i++) {
    free(aN_r[i]);
    free(aN_i[i]);
  }
  free(aN_r);
  free(aN_i);

  free(psi_r);
  free(psi_i);
  free(z_r);
  free(z_i);
  free(zn_r);
  free(zn_i);
  free(phi_r);
  free(phi_i);
  free(weight);
  free(c_r);
  free(c_i);
  free(alpha_r);
  free(alpha_i);
  free(phi_norm);
}

void doHarmonicFilter(ProgramData *p, _HarmonicFilter *c, int threadid, int lcid) {
  double *avals, *avals_orig = NULL;
  double *bvals, *bvals_orig = NULL;
  int Njd;
  double *t, *mag, *err, *mag_model;
  double maxfreq_calc, maxfreq_bandcut = -1, minfreq_bandcut = -1;
  double maxfreq = -1, minfreq = -1, minfreq_model = -1, maxfreq_model;
  double T, df, f, delmin, delt, twojdtol;
  double filter_val;
  int Nf, Nftot, Nfcalc, i, j, isuniform;
  FILE *outfile;
  char *lcoutname;

  double f0, xt, xdf, s0, c0, sdf, cdf, s, cval;

  Njd = p->NJD[threadid];
  t = p->t[threadid];
  mag = p->mag[threadid];
  err = p->sig[threadid];

  if(Njd < 3)
    return;

  T = t[Njd-1] - t[0];
  /* Determine the frequency spacing and number of frequencies */
  /* Also check if the sampling is uniform */
  isuniform = 1;
  twojdtol = 2.0*p->JDTOL;
  if(t[1] != t[0])
    delmin = t[1] - t[0];
  else
    delmin = T;
  for(i=2;i<Njd;i++)
    {
      delt = t[i] - t[i-1];
      if(isuniform) {
	if(delt < delmin - twojdtol || delt > delmin + twojdtol)
	  isuniform = 0;
      }
      if(delt < delmin && delt > 0.)
	delmin = t[i] - t[i-1];
    }
  maxfreq_calc = 1./(2.0*delmin);
  Nftot = floor(4.0 * T * maxfreq_calc);
  df = maxfreq_calc / (double) Nftot;
  if(2.0*Nftot+1 > Njd) {
    Nftot = (Njd-1)/2;
    maxfreq_calc = df*Nftot;
  }

  isuniform = 0;

  maxfreq_model = maxfreq_calc;

  /* Determine any upper or lower limits on the frequencies to compute */
  switch(c->filtertype) {
  case VARTOOLS_HARMONICFILTER_FULLSPEC:
    break;
  case VARTOOLS_HARMONICFILTER_LOWPASS:
    GetDoubleParameterValue(threadid, lcid, &maxfreq, 
			    c->maxfreq_source, c->maxfreq_fix,
			    c->maxfreq, c->maxfreq_linkedcolumn,
			    c->maxfreq_expr);
    if(!c->calc_full_spec) {
      if(maxfreq < maxfreq_calc)
	maxfreq_calc = maxfreq;
    }
    maxfreq_model = maxfreq;
    break;
  case VARTOOLS_HARMONICFILTER_HIGHPASS:
    GetDoubleParameterValue(threadid, lcid, &minfreq, 
			    c->minfreq_source, c->minfreq_fix,
			    c->minfreq, c->minfreq_linkedcolumn,
			    c->minfreq_expr);
    if(c->filter_exprstring == NULL && !c->calc_full_spec) {
      if(minfreq < maxfreq_calc)
	maxfreq_calc = minfreq;
    } else {
      minfreq_model = minfreq;
    }
    maxfreq_model = maxfreq_calc;
    break;
  case VARTOOLS_HARMONICFILTER_BANDPASS:
    GetDoubleParameterValue(threadid, lcid, &minfreq, 
			    c->minfreq_source, c->minfreq_fix,
			    c->minfreq, c->minfreq_linkedcolumn,
			    c->minfreq_expr);
    GetDoubleParameterValue(threadid, lcid, &maxfreq, 
			    c->maxfreq_source, c->maxfreq_fix,
			    c->maxfreq, c->maxfreq_linkedcolumn,
			    c->maxfreq_expr);
    if(!c->calc_full_spec) {
      if(maxfreq_calc < maxfreq)
	maxfreq_calc = maxfreq;
    }
    minfreq_model = minfreq;
    maxfreq_model = maxfreq;
    break;
  case VARTOOLS_HARMONICFILTER_BANDCUT:
    GetDoubleParameterValue(threadid, lcid, &minfreq, 
			    c->minfreq_source, c->minfreq_fix,
			    c->minfreq, c->minfreq_linkedcolumn,
			    c->minfreq_expr);
    GetDoubleParameterValue(threadid, lcid, &maxfreq, 
			    c->maxfreq_source, c->maxfreq_fix,
			    c->maxfreq, c->maxfreq_linkedcolumn,
			    c->maxfreq_expr);
    if(c->filter_exprstring == NULL && !c->calc_full_spec) {
      if(maxfreq_calc < maxfreq)
	maxfreq_calc = maxfreq;
      minfreq_model = minfreq;
      maxfreq_model = maxfreq;
      maxfreq_bandcut = -1;
      minfreq_bandcut = -1;
    } else {
      maxfreq_bandcut = maxfreq;
      minfreq_bandcut = minfreq;
    }
    break;
  }
  
  Nfcalc = floor(maxfreq_calc/df);

  if(!isuniform && !c->forcefft) {
    /* For non-uniform sampling and not using forcefft, calculate the
       fourier series using the projection onto orthogonal
       polynomials */
    if((avals = (double *) malloc((Nfcalc + 1)*sizeof(double))) == NULL ||
       (bvals = (double *) malloc((Nfcalc + 1)*sizeof(double))) == NULL ||
       (mag_model = (double *) malloc(Njd * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);

    if(c->filter_exprstring != NULL && c->ofourier) {
      if((avals_orig = (double *) malloc((Nfcalc + 1)*sizeof(double))) == NULL ||
	 (bvals_orig = (double *) malloc((Nfcalc + 1)*sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
    
    fit_harmonic_series_RAG(Njd, t, mag, err, df, Nfcalc, avals, bvals);

    if(avals_orig != NULL) {
      for(i=0; i < Nfcalc+1; i++) {
	avals_orig[i] = avals[i];
	bvals_orig[i] = bvals[i];
      }
    }
    
    if(c->filter_exprstring != NULL) {
      /* Apply the analytic filter to the fourier series */
      for(i=1,f=df;i < Nfcalc+1;i++, f += df) {
	SetVariable_Value_Double(lcid, threadid, 0, c->freq_var, f);
	filter_val = EvaluateExpression(lcid, threadid, 0, c->filter_expr);
	avals[i] *= filter_val;
	bvals[i] *= filter_val;
      }
    }

    /* Output the Fourier Coefficients if requested */
    if(c->ofourier) {
      if((lcoutname = (char *) malloc(MAXLEN)) == NULL)
	error(ERR_MEMALLOC);
      GetOutputFilename(lcoutname, p->lcnames[lcid], c->ofourier_dir,
			"fouriercoeffs", c->ofourier_format, lcid);
      if((outfile = fopen(lcoutname,"w")) == NULL) {
	error2(ERR_CANNOTWRITE,lcoutname);
      }
      fprintf(outfile,"#Frequency");
      if(avals_orig != NULL) {
	fprintf(outfile," CosCoeff_orig SinCoeff_orig CosCoeff_filter SinCoeff_filter\n");
      } else {
	fprintf(outfile," CosCoeff SinCoeff\n");
      }
      for(i=0,f=0; i <= Nfcalc; i++, f += df) {
	fprintf(outfile, "%.17g", f);
	if(avals_orig != NULL) {
	  fprintf(outfile," %.17g %.17g %.17g %.17g\n", avals_orig[i], bvals_orig[i], avals[i], bvals[i]);
	} else {
	  fprintf(outfile," %.17g %.17g\n", avals[i], bvals[i]);
	}
      }
      fclose(outfile);
      free(lcoutname);
      if(avals_orig != NULL) free(avals_orig);
      if(bvals_orig != NULL) free(bvals_orig);
    }

    /* Calculate the time series from these fourier coefficients */
    if(c->filtertype != VARTOOLS_HARMONICFILTER_BANDCUT ||
       (c->filtertype == VARTOOLS_HARMONICFILTER_BANDCUT && minfreq_bandcut < 0 && maxfreq_bandcut < 0)) {
      for(i=0; i < Njd; i++) {
	mag_model[i] = avals[0];
	if(Nfcalc >= 1) {
	  f0 = df;
	  xt = TWOPI * t[i];
	  xdf = xt*df;
	  s0 = sin(xdf);
	  c0 = cos(xdf);
	  sdf = s0;
	  cdf = c0;
	  if(f0 >= minfreq_model && f0 <= maxfreq_model) {
	    mag_model[i] += avals[1]*c0 + bvals[1]*s0;
	  }
	  for(j=2; j <= Nfcalc; j++) {
	    f0 += df;
	    if(f0 > maxfreq_model)
	      break;
	    cval = c0*cdf - s0*sdf;
	    s = s0*cdf + c0*sdf;
	    if(f0 >= minfreq_model) {
	      mag_model[i] += avals[j]*cval + bvals[j]*s;
	    }
	    c0 = cval;
	    s0 = s;
	  }
	}
      }
    } else {
      for(i=0; i < Njd; i++) {
	mag_model[i] = avals[0];
	if(Nfcalc >= 1) {
	  f0 = df;
	  xt = TWOPI * t[i];
	  xdf = xt*df;
	  s0 = sin(xdf);
	  c0 = cos(xdf);
	  sdf = s0;
	  cdf = c0;
	  if(f0 <= minfreq_bandcut || f0 >= maxfreq_bandcut) {
	    mag_model[i] += avals[1]*c0 + bvals[1]*s0;
	  }
	  for(j=2; j <= Nfcalc; j++) {
	    f0 += df;
	    cval = c0*cdf - s0*sdf;
	    s = s0*cdf + c0*sdf;
	    if(f0 <= minfreq_bandcut || f0 >= maxfreq_bandcut) {
	      mag_model[i] += avals[j]*cval + bvals[j]*s;
	    }
	    c0 = cval;
	    s0 = s;
	  }
	}
      }
    }

    /* Either replace the light curve with the model or subtract the
       model, depending on how the filter is to be applied */
    switch(c->filtertype) {
    case VARTOOLS_HARMONICFILTER_FULLSPEC:
      /* Replace the light curve with the model */
	for(i=0; i < Njd; i++) {
	  mag[i] = mag_model[i];
	}
	break;
    case VARTOOLS_HARMONICFILTER_HIGHPASS:
      if(c->filter_exprstring != NULL) {
	/* Replace the light curve with the model */
	for(i=0; i < Njd; i++) {
	  mag[i] = mag_model[i];
	}
      } else {
	/* Subtract the model from the light curve, but save the 0 term */
	for(i=0; i < Njd; i++) {
	  mag[i] = mag[i] - mag_model[i] + avals[0];
	}
      }
      break;
    case VARTOOLS_HARMONICFILTER_LOWPASS:
      /* Replace the light curve with the model */
      for(i=0; i < Njd; i++) {
	mag[i] = mag_model[i];
      }
      break;
    case VARTOOLS_HARMONICFILTER_BANDPASS:
      /* Replace the light curve with the model */
      for(i=0; i < Njd; i++) {
	mag[i] = mag_model[i];
      }
      break;
    case VARTOOLS_HARMONICFILTER_BANDCUT:
      if(minfreq_bandcut < 0 && maxfreq_bandcut < 0) {
	/* Subtract the model from the light curve, but save the 0 term */
	for(i=0; i < Njd; i++) {
	  mag[i] = mag[i] - mag_model[i] + avals[0];
	}
      }
      else {
	/* Replace the light curve with the model */
	for(i=0; i < Njd; i++) {
	  mag[i] = mag_model[i];
	}
      }
      break;
    }

    free(avals);
    free(bvals);
    free(mag_model);
  }
  /******* TBD: implement the FFT-based filter ********/
}

int ParseHarmonicFilterCommand(int *iret, int argc, char **argv, ProgramData *p,
			       _HarmonicFilter *c, int cnum)
/* Parse the command line for the "-harmonicfilter" command */
{
  int i, j, k;
  FILE *infile;
  char *line = NULL;
  size_t line_size = MAXLEN;
  int colnum;
  int size_tout = 0;
  
  i = *iret;
  if(i >= argc)
    return(1);

  if(!strcmp(argv[i],"highpass")) {
    c->filtertype = VARTOOLS_HARMONICFILTER_HIGHPASS;
    i++;
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "minfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->minfreq_source),
				    0, (void *) (&(c->minfreq_fix)),
				    (void *) (&(c->minfreq)),
				    (void *) (&(c->minfreq_linkedcolumn)),
				    (void *) (&(c->minfreq_exprstring)),
				    "HARMONICFILTER_MINFREQ",
				    0, NULL)) {
      *iret = i-1; return 1;}
  } else if(!strcmp(argv[i],"lowpass")) {
    c->filtertype = VARTOOLS_HARMONICFILTER_LOWPASS;
    i++;
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "maxfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->maxfreq_source),
				    0, (void *) (&(c->maxfreq_fix)),
				    (void *) (&(c->maxfreq)),
				    (void *) (&(c->maxfreq_linkedcolumn)),
				    (void *) (&(c->maxfreq_exprstring)),
				    "HARMONICFILTER_MAXFREQ",
				    0, NULL)) {
      *iret = i-1; return 1;}
  } else if(!strcmp(argv[i],"bandpass")) {
    c->filtertype = VARTOOLS_HARMONICFILTER_BANDPASS;
    i++;
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "minfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->minfreq_source),
				    0, (void *) (&(c->minfreq_fix)),
				    (void *) (&(c->minfreq)),
				    (void *) (&(c->minfreq_linkedcolumn)),
				    (void *) (&(c->minfreq_exprstring)),
				    "HARMONICFILTER_MINFREQ",
				    0, NULL)) {
      *iret = i; return 1;}
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "maxfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->maxfreq_source),
				    0, (void *) (&(c->maxfreq_fix)),
				    (void *) (&(c->maxfreq)),
				    (void *) (&(c->maxfreq_linkedcolumn)),
				    (void *) (&(c->maxfreq_exprstring)),
				    "HARMONICFILTER_MAXFREQ",
				    0, NULL)) {
      *iret = i-1; return 1;}
  } else if(!strcmp(argv[i],"bandcut")) {
    c->filtertype = VARTOOLS_HARMONICFILTER_BANDCUT;
    i++;
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "minfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->minfreq_source),
				    0, (void *) (&(c->minfreq_fix)),
				    (void *) (&(c->minfreq)),
				    (void *) (&(c->minfreq_linkedcolumn)),
				    (void *) (&(c->minfreq_exprstring)),
				    "HARMONICFILTER_MINFREQ",
				    0, NULL)) {
      *iret = i-1; return 1;}
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "maxfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->maxfreq_source),
				    0, (void *) (&(c->maxfreq_fix)),
				    (void *) (&(c->maxfreq)),
				    (void *) (&(c->maxfreq_linkedcolumn)),
				    (void *) (&(c->maxfreq_exprstring)),
				    "HARMONICFILTER_MAXFREQ",
				    0, NULL)) {
      *iret = i-1; return 1;}
  } else {*iret = i-1; return 1;}

  *iret = i-1;
  return 0;
}
