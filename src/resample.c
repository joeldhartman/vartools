#include "commands.h"
#include "programdata.h"
#include "functions.h"

#ifdef _HAVE_GSL
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#endif

typedef struct {
  double yp1;
  double ypn;
  int nbreaks;
  int spline_order;
} _InterpParams;

void SetupResampleExpression(ProgramData *p, _Resample *c)
{
  if(c->use_near_far) {
    if(c->minsep_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
      c->minsep_expr = ParseExpression(c->minsep_exprstring, p);
    }
  }

  if(c->use_file) return;

  if(c->tstart_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    c->tstart_expr = ParseExpression(c->tstart_exprstring, p);
  }

  if(c->tstop_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    c->tstop_expr = ParseExpression(c->tstop_exprstring, p);
  }
  
  if(c->Nresamp_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    c->Nresamp_expr = ParseExpression(c->Nresamp_exprstring, p);
  }
  
  if(c->delt_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
    c->delt_expr = ParseExpression(c->delt_exprstring, p);
  }
  
    
}

void InterpolateVec_Linear(int N_in, double *t_in, double *vec_in,
			   int N_out, double *t_out, double *vec_out)
/* Performs linear interpolation. Assumes both t_in and t_out are sorted. */
{
  int i, j;
  double h;
  
  for(i=0, j=0; i < N_out; i++) {
    while(j < N_in - 1) {
      if(t_out[i] < t_in[j+1]) {
	break;
      }
      j++;
    }
    if(j == 0 && t_out[i] < t_in[j]) {
	/* Point needs to be extrapolated to the left */
	h = t_in[1] - t_in[0];
	vec_out[i] = vec_in[0] + (t_out[i] - t_in[0])*((vec_in[1]-vec_in[0])/h);
    }
    else if(j == N_in - 1) {
      /* Point needs to be extrapolated to the right */
	h = t_in[j] - t_in[j-1];
	vec_out[i] = vec_in[j] + (t_out[i] - t_in[j])*((vec_in[j]-vec_in[j-1])/h);
    }
    else {
      h = t_in[j+1] - t_in[j];
      vec_out[i] = vec_in[j] + (t_out[i] - t_in[j])*((vec_in[j+1]-vec_in[j])/h);
    }
  }
  return;
}

void InterpolateErr_Linear(int N_in, double *t_in, double *vec_in,
			   int N_out, double *t_out, double *vec_out)
/* Returns Errors from linear interpolation. Assumes both t_in and t_out are sorted, and that input uncertainties are uncorrelated. */
{
  int i, j;
  double h, c;
  
  for(i=0, j=0; i < N_out; i++) {
    while(j < N_in - 1) {
      if(t_out[i] < t_in[j+1]) {
	break;
      }
      j++;
    }
    if(j == 0 && t_out[i] < t_in[j]) {
      /* Point needs to be extrapolated to the left */
      h = t_in[1] - t_in[0];
      c = (t_out[i] - t_in[0])/h;
      vec_out[i] = sqrt((1.0-c)*(1.0-c)*vec_in[0]*vec_in[0] + c*c*vec_in[1]*vec_in[1]);
    }
    else if(j == N_in - 1) {
      /* Point needs to be extrapolated to the right */
      h = t_in[j] - t_in[j-1];
      c = (t_out[i] - t_in[j])/h;
      vec_out[i] = sqrt((1.0+c)*(1.0+c)*vec_in[j]*vec_in[j] + c*c*vec_in[j-1]*vec_in[j-1]);
    }
    else {
      h = t_in[j+1] - t_in[j];
      c = (t_out[i] - t_in[j])/h;
      vec_out[i] = sqrt((1.0-c)*(1.0-c)*vec_in[j]*vec_in[j] + c*c*vec_in[j+1]*vec_in[j+1]);
    }
  }
  return;
}

#ifdef _HAVE_GSL

/* Perform BSpline interpolation using GSL functions. The code is based on the
   B-Spline example code supplied with GSL */
void InterpolateVec_BSpline(int N_in, double *t_in, double *vec_in, int N_out,
			    double *t_out, double *vec_out, int nbreaks, 
			    int spline_order, int use_err, double *err_in,
			    double *err_out)
{
  int ncoeffs, i, j, Nshift;
  double Bj, yi, yerr;
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  gsl_multifit_linear_workspace *mw;
  gsl_vector *x, *y, *c, *w, *breakpts;
  gsl_matrix *X, *cov;
  double chisq;
  double t0, t1, span, nbreaks2;
  
  if(nbreaks - 2 > N_in) {nbreaks = N_in + 2;}

  ncoeffs = nbreaks - 2 + spline_order;
  bw = gsl_bspline_alloc(spline_order, nbreaks);
  B = gsl_vector_alloc(ncoeffs);
  x = gsl_vector_alloc(N_in);
  y = gsl_vector_alloc(N_in);
  X = gsl_matrix_alloc(N_in, ncoeffs);
  c = gsl_vector_alloc(ncoeffs);
  w = gsl_vector_alloc(N_in);
  breakpts = gsl_vector_alloc(nbreaks);
  cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
  mw = gsl_multifit_linear_alloc(N_in, ncoeffs);
  
  for(i=0; i < N_in; i++) {
    gsl_vector_set(x, i, t_in[i]);
    gsl_vector_set(y, i, vec_in[i]);
    if(use_err)
      gsl_vector_set(w, i, 1.0/(err_in[i]*err_in[i]));
    else
      gsl_vector_set(w, i, 1.0);
  }

  t0 = (t_out[0] < t_in[0] ? t_out[0] : t_in[0]);
  t1 = (t_out[N_out-1] > t_in[N_in-1] ? t_out[N_out-1] : t_in[N_in-1]);
  span = t1-t0;
  t0 = t0 - 0.1*(span/nbreaks);
  t1 = t1 + 0.1*(span/nbreaks);

  nbreaks2 = nbreaks - 2;
  if(nbreaks2 > 0) {
    Nshift = N_in/nbreaks2;
  }
  gsl_vector_set(breakpts, 0, t0);
  gsl_vector_set(breakpts, (nbreaks-1), t1);
  j = Nshift/2;
  for(i=1; i < (nbreaks - 1); i++) {
    gsl_vector_set(breakpts, i, t_in[j]);
    j += Nshift;
  }
    
  gsl_bspline_knots(breakpts, bw);

  /* Setup the fit matrix X */
  for(i = 0; i < N_in; i++)
    {
      /* compute B_j(t_in[i]) for all j */
      gsl_bspline_eval(t_in[i], B, bw);
      
      /* fill in row i of X */
      for(j=0; j < ncoeffs; j++) {
	Bj = gsl_vector_get(B, j);
	gsl_matrix_set(X, i, j, Bj);
      }
    }

  /* Do the fit */
  gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);

  /* Output the smoothed curve */
  for(i=0; i < N_out; i++) {
    gsl_bspline_eval(t_out[i], B, bw);
    gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
    if(use_err) {
      err_out[i] = yerr;
    }
    vec_out[i] = yi;
  }
  
  /* Clean-up */
  gsl_bspline_free(bw);
  gsl_vector_free(B);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(w);
  gsl_matrix_free(cov);
  gsl_vector_free(breakpts);
  gsl_multifit_linear_free(mw);
  return;
}

void InterpolateVec_BSpline_AutoBreaks(int N_in, double *t_in, double *vec_in, 
				       int N_out, double *t_out, 
				       double *vec_out, int *nbreaks_out, 
				       int spline_order, int use_err, 
				       double *err_in, double *err_out)
/* This version varies the number of knots to use in the B-spline until
   chi2/dof < 1, at which point it stops. The number of knots used is returned
   in nbreaks_out. */
{
  int ncoeffs, i, j;
  double Bj, yi, yerr;
  gsl_bspline_workspace *bw;
  gsl_vector *B;
  gsl_multifit_linear_workspace *mw;
  gsl_vector *x, *y, *c, *w, *breakpts;
  gsl_matrix *X, *cov;
  double chisq, chisqdof;
  int nbreaks, dof, nbreaks2, Nshift;
  double t0, t1, span;

  x = gsl_vector_alloc(N_in);
  y = gsl_vector_alloc(N_in);
  w = gsl_vector_alloc(N_in);
  
  for(i=0; i < N_in; i++) {
    gsl_vector_set(x, i, t_in[i]);
    gsl_vector_set(y, i, vec_in[i]);
    if(use_err)
      gsl_vector_set(w, i, 1.0/(err_in[i]*err_in[i]));
    else
      gsl_vector_set(w, i, 1.0);
  }
  
  nbreaks = 1;

  do {

    nbreaks += 10;

    ncoeffs = nbreaks - 2 + spline_order;
    bw = gsl_bspline_alloc(spline_order, nbreaks);
    B = gsl_vector_alloc(ncoeffs);
    X = gsl_matrix_alloc(N_in, ncoeffs);
    c = gsl_vector_alloc(ncoeffs);
    cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
    breakpts = gsl_vector_alloc(nbreaks);
    mw = gsl_multifit_linear_alloc(N_in, ncoeffs);

    t0 = (t_out[0] < t_in[0] ? t_out[0] : t_in[0]);
    t1 = (t_out[N_out-1] > t_in[N_in-1] ? t_out[N_out-1] : t_in[N_in-1]);
    span = t1-t0;
    t0 = t0 - 0.1*(span/nbreaks);
    t1 = t1 + 0.1*(span/nbreaks);

    nbreaks2 = nbreaks - 2;
    if(nbreaks2 > 0) {
      Nshift = N_in/nbreaks2;
    }
    gsl_vector_set(breakpts, 0, t0);
    gsl_vector_set(breakpts, (nbreaks-1), t1);
    j = Nshift/2;
    for(i=1; i < (nbreaks - 1); i++) {
      gsl_vector_set(breakpts, i, t_in[j]);
      j += Nshift;
    }
    
    gsl_bspline_knots(breakpts, bw);

    
    /* Setup the fit matrix X */
    for(i = 0; i < N_in; i++)
      {
	/* compute B_j(t_in[i]) for all j */
	gsl_bspline_eval(t_in[i], B, bw);
	
	/* fill in row i of X */
	for(j=0; j < ncoeffs; j++) {
	  Bj = gsl_vector_get(B, j);
	  gsl_matrix_set(X, i, j, Bj);
	}
      }
    
    /* Do the fit */
    gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);
   
    dof = N_in - ncoeffs;
    
    chisqdof = chisq/dof;

    /* Clean-up */
    gsl_bspline_free(bw);
    gsl_vector_free(B);
    gsl_matrix_free(X);
    gsl_vector_free(c);
    gsl_vector_free(breakpts);
    gsl_matrix_free(cov);
    gsl_multifit_linear_free(mw);
  } while(chisqdof > 1 && nbreaks - 2 < N_in);
  
  /* Adopt the next smallest number of breaks, and use that
     B-spline solution */
  if(nbreaks > 12) nbreaks -= 10;

  ncoeffs = nbreaks - 2 + spline_order;
  bw = gsl_bspline_alloc(spline_order, nbreaks);
  B = gsl_vector_alloc(ncoeffs);
  X = gsl_matrix_alloc(N_in, ncoeffs);
  c = gsl_vector_alloc(ncoeffs);
  cov = gsl_matrix_alloc(ncoeffs, ncoeffs);
  mw = gsl_multifit_linear_alloc(N_in, ncoeffs);
  breakpts = gsl_vector_alloc(nbreaks);
  
  t0 = (t_out[0] < t_in[0] ? t_out[0] : t_in[0]);
  t1 = (t_out[N_out-1] > t_in[N_in-1] ? t_out[N_out-1] : t_in[N_in-1]);
  span = t1-t0;
  t0 = t0 - 0.1*(span/nbreaks);
  t1 = t1 + 0.1*(span/nbreaks);

  nbreaks2 = nbreaks - 2;
  if(nbreaks2 > 0) {
    Nshift = N_in/nbreaks2;
  }
  gsl_vector_set(breakpts, 0, t0);
  gsl_vector_set(breakpts, (nbreaks-1), t1);
  j = Nshift/2;
  for(i=1; i < (nbreaks - 1); i++) {
    gsl_vector_set(breakpts, i, t_in[j]);
    j += Nshift;
  }
  
  gsl_bspline_knots(breakpts, bw);
  
  /* Setup the fit matrix X */
  for(i = 0; i < N_in; i++)
    {
      /* compute B_j(t_in[i]) for all j */
      gsl_bspline_eval(t_in[i], B, bw);
      
      /* fill in row i of X */
      for(j=0; j < ncoeffs; j++) {
	Bj = gsl_vector_get(B, j);
	gsl_matrix_set(X, i, j, Bj);
      }
    }
    
  /* Do the fit */
  gsl_multifit_wlinear(X, w, y, c, cov, &chisq, mw);  

  /* Output the smoothed curve */
  for(i=0; i < N_out; i++) {
    gsl_bspline_eval(t_out[i], B, bw);
    gsl_multifit_linear_est(B, c, cov, &yi, &yerr);
    if(use_err) {
      err_out[i] = yerr;
    }
    vec_out[i] = yi;
  }

  *nbreaks_out = nbreaks;
    
  /* Clean-up */
  gsl_bspline_free(bw);
  gsl_vector_free(x);
  gsl_vector_free(y);
  gsl_vector_free(w);
  gsl_vector_free(B);
  gsl_matrix_free(X);
  gsl_vector_free(c);
  gsl_vector_free(breakpts);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free(mw);

  return;
}

#endif


void InterpolateVec_Spline(int N_in, double *t_in, double *vec_in,
			   int N_out, double *t_out, double *vec_out,
			   double yp1, double ypn,
			   double **y2_in, double **u_in, int *size_y2_in)
{
  int j;
  double *y2;
  double *u;
  if(N_in > (*size_y2_in)) {
    if(!(*size_y2_in)) {
      if((*y2_in = (double *) malloc(N_in * sizeof(double))) == NULL ||
	 (*u_in = (double *) malloc(N_in * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    } else {
      if((*y2_in = (double *) realloc(*y2_in, N_in*sizeof(double))) == NULL ||
	 (*u_in = (double *) realloc(*u_in, N_in*sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
    *size_y2_in = N_in;
  }
  y2 = *y2_in;
  u = *u_in;
  spline(t_in, vec_in, N_in, yp1, ypn, y2, u);
  for(j=0; j < N_out; j++) {
    splint(t_in, vec_in, y2, N_in, t_out[j], &(vec_out[j]));
  }
  return;
}

void InterpolateVec_SplineMonotonic(int N_in, double *t_in, double *vec_in,
				    int N_out, double *t_out, double *vec_out,
				    double **y2_in, int *size_y2_in)
{
  int j;
  double *y2;
  if(N_in > (*size_y2_in)) {
    if(!(*size_y2_in)) {
      if((*y2_in = (double *) malloc(N_in * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    } else {
      if((*y2_in = (double *) realloc(*y2_in, N_in*sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }
    *size_y2_in = N_in;
  }
  y2 = *y2_in;
  spline_monotonic(N_in, t_in, vec_in, y2);
  for(j=0; j < N_out; j++) {
    vec_out[j] = splint_monotonic(N_in, t_in, vec_in, y2, t_out[j]);
  }
  return;
}

void InterpolateVec_Nearest(int N_in, double *t_in, double *vec_in,
			    int N_out, double *t_out, double *vec_out)
/* Assumes both t_in and t_out are sorted. Copies into vec_out the elements
   of vec_in for with t_in is closest to t_out. In cases of equal distance, the
   value from the lower t_in is copied */
{
  int i, j;
  double d;
  for(i=0, j=0; i < N_out; i++) {
    d = fabs(t_out[i] - t_in[j]);
    while(j < N_in - 1) {
      if(fabs(t_out[i] - t_in[j+1]) >= d) {
	break;
      }
      d = fabs(t_out[i] - t_in[j+1]);
      j++;
    }
    vec_out[i] = vec_in[j];
  }
}

void InterpolateVec_Nearest_String(int N_in, double *t_in, char **vec_in,
				   int N_out, double *t_out, char **vec_out)
/* Assumes both t_in and t_out are sorted. Copies into vec_out the elements
   of vec_in for with t_in is closest to t_out. In cases of equal distance, the
   value from the lower t_in is copied */
{
  int i, j;
  double d;
  for(i=0, j=0; i < N_out; i++) {
    d = fabs(t_out[i] - t_in[j]);
    while(j < N_in - 1) {
      if(fabs(t_out[i] - t_in[j+1]) >= d) {
	break;
      }
      d = fabs(t_out[i] - t_in[j+1]);
      j++;
    }
    sprintf(vec_out[i],"%s",vec_in[j]);
  }
}

void InterpolateVec_Nearest_Char(int N_in, double *t_in, char *vec_in,
				   int N_out, double *t_out, char *vec_out)
/* Assumes both t_in and t_out are sorted. Copies into vec_out the elements
   of vec_in for with t_in is closest to t_out. In cases of equal distance, the
   value from the lower t_in is copied */
{
  int i, j;
  double d;
  for(i=0, j=0; i < N_out; i++) {
    d = fabs(t_out[i] - t_in[j]);
    while(j < N_in - 1) {
      if(fabs(t_out[i] - t_in[j+1]) >= d) {
	break;
      }
      d = fabs(t_out[i] - t_in[j+1]);
      j++;
    }
    vec_out[i] = vec_in[j];
  }
}

/* Performs separate modes of interpolation for points within a
   pre-determined distance from an observation, for points outside
   the distance, and for points that need to be extrapolated. */
void InterpolateVec_Gap_Extrap(int N_in, double *t_in, double *vec_in, 
			int N_out, double *t_out, double *vec_out, 
                        double minsep, int method_close, _InterpParams *p_close,
			int method_far, _InterpParams *p_far,
			int method_extrap, _InterpParams *p_extrap,
			int use_err, double *err_in, double *err_out)
{
  double *vec_out_close, *vec_out_far, *vec_out_extrap;
  double *t_out_close, *t_out_far, *t_out_extrap;
  double *err_out_close = NULL, *err_out_far = NULL, *err_out_extrap = NULL;
  int N_out_close, N_out_far, N_out_extrap;
  char *which_method;
  int i, j, j1, j2, j3;
  double d;
  double *y2 = NULL, *u_spline = NULL;
  int size_y2 = 0;
  double *v_tmp, *t_tmp, *e_tmp;
  int method_tmp, N_tmp;
  _InterpParams *p_tmp;

  if((vec_out_close = (double *) malloc(N_out*sizeof(double))) == NULL ||
     (vec_out_far = (double *) malloc(N_out*sizeof(double))) == NULL ||
     (vec_out_extrap = (double *) malloc(N_out*sizeof(double))) == NULL ||
     (t_out_close = (double *) malloc(N_out*sizeof(double))) == NULL ||
     (t_out_far = (double *) malloc(N_out*sizeof(double))) == NULL ||
     (t_out_extrap = (double *) malloc(N_out*sizeof(double))) == NULL ||
     (which_method = (char *) malloc(N_out*sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  
  if(use_err) {
    if((err_out_close = (double *) malloc(N_out*sizeof(double))) == NULL ||
       (err_out_far = (double *) malloc(N_out*sizeof(double))) == NULL ||
       (err_out_extrap = (double *) malloc(N_out*sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
  }
  
  N_out_close = 0;
  N_out_far = 0;
  N_out_extrap = 0;
  for(i=0, j=0; i < N_out; i++) {
    while(j < N_in - 1) {
      if(t_out[i] < t_in[j+1]) {
	break;
      }
      j++;
    }
    if((j == 0 && t_out[i] < t_in[j]) || j == N_in - 1) {
      /* Point needs to be extrapolated */
      if(method_extrap >= 0) {
	which_method[i] = 2;
	t_out_extrap[N_out_extrap] = t_out[i];
	N_out_extrap++;
      } else {
	if(fabs(t_out[i] - t_in[j]) < minsep) {
	  which_method[i] = 0;
	  t_out_close[N_out_close] = t_out[i];
	  N_out_close++;
	}
	else {
	  which_method[i] = 1;
	  t_out_far[N_out_far] = t_out[i];
	  N_out_far++;
	}
      }
    }
    else {
      if(fabs(t_in[j+1]-t_out[i]) < minsep && fabs(t_out[i]-t_in[j]) < minsep) {
	/* Use the close interpolation */
	which_method[i] = 0;
	t_out_close[N_out_close] = t_out[i];
	N_out_close++;
      } else {
	/* Use the far interpolation */
	which_method[i] = 1;
	t_out_far[N_out_far] = t_out[i];
	N_out_far++;
      }
    }
  }

  /* Run the close, far, and extrapolation interpolations */
  for(i=0; i < 3; i++) {
    if(i == 0) {
      v_tmp = vec_out_close; t_tmp = t_out_close; e_tmp = err_out_close;
      method_tmp = method_close; N_tmp = N_out_close; p_tmp = p_close;
    }
    else if(i == 1) {
      v_tmp = vec_out_far; t_tmp = t_out_far; e_tmp = err_out_far;
      method_tmp = method_far; N_tmp = N_out_far; p_tmp = p_far;
    }
    else if(i == 2) {
      v_tmp = vec_out_extrap; t_tmp = t_out_extrap; e_tmp = err_out_extrap;
      method_tmp = method_extrap; N_tmp = N_out_extrap; p_tmp = p_extrap;
    }
    if(N_tmp > 0) {
      switch(method_tmp) {
      case VARTOOLS_RESAMPLE_NEAREST:
	InterpolateVec_Nearest(N_in, t_in, vec_in,
			       N_tmp, t_tmp, v_tmp);
	if(use_err) {
	  InterpolateVec_Nearest(N_in, t_in, err_in,
				 N_tmp, t_tmp, e_tmp);
	}
	break;
      case VARTOOLS_RESAMPLE_LINEAR:
	InterpolateVec_Linear(N_in, t_in, vec_in,
			      N_tmp, t_tmp, v_tmp);
	if(use_err) {
	  InterpolateErr_Linear(N_in, t_in, err_in,
				N_tmp, t_tmp, e_tmp);
	}
	break;
      case VARTOOLS_RESAMPLE_SPLINE:
	InterpolateVec_Spline(N_in, t_in, vec_in,
			      N_tmp, t_tmp, v_tmp,
			      p_tmp->yp1, p_tmp->ypn, &y2, &u_spline,
			      &size_y2);
	if(use_err) {
	  InterpolateErr_Linear(N_in, t_in, err_in,
				N_tmp, t_tmp, e_tmp);
	}
	break;
      case VARTOOLS_RESAMPLE_SPLINEMONOTONIC:
	InterpolateVec_SplineMonotonic(N_in, t_in, vec_in,
				       N_tmp, t_tmp, v_tmp,
				       &y2, &size_y2);
	if(use_err) {
	  InterpolateErr_Linear(N_in, t_in, err_in,
				N_tmp, t_tmp, e_tmp);
	}
	break;
#ifdef _HAVE_GSL
      case VARTOOLS_RESAMPLE_BSPLINE:
	if(p_tmp->nbreaks >= 2) {
	  InterpolateVec_BSpline(N_in, t_in, vec_in,
				 N_tmp, t_tmp, v_tmp,
				 p_tmp->nbreaks, p_tmp->spline_order,
				 use_err, err_in, e_tmp);
	} else {
	  InterpolateVec_BSpline_AutoBreaks(N_in, t_in, vec_in,
					    N_tmp, t_tmp,
					    v_tmp, &(p_tmp->nbreaks),
					    p_tmp->spline_order, use_err,
					    err_in, e_tmp);
	}
	break;
#endif
      default:
	error(ERR_CODEERROR);
      }
    }
  }

  /* Copy the interpolated data into the output vectors */
  j1 = 0; j2 = 0; j3 = 0;
  for(i=0; i < N_out; i++) {
    if(which_method[i] == 0) {
      vec_out[i] = vec_out_close[j1];
      j1++;
    }
    else if(which_method[i] == 1) {
      vec_out[i] = vec_out_far[j2];
      j2++;
    }
    else if(which_method[i] == 2) {
      vec_out[i] = vec_out_extrap[j3];
      j3++;
    }
  }
  if(use_err) {
    j1 = 0; j2 = 0; j3 = 0;
    for(i=0; i < N_out; i++) {
      if(which_method[i] == 0) {
	err_out[i] = err_out_close[j1];
	j1++;
      }
      else if(which_method[i] == 1) {
	err_out[i] = err_out_far[j2];
	j2++;
      }
      else if(which_method[i] == 2) {
	err_out[i] = err_out_extrap[j3];
	j3++;
      }
    }
  }

  /* Clean up */
  if(vec_out_close != NULL) free(vec_out_close);
  if(vec_out_far != NULL) free(vec_out_far);
  if(vec_out_extrap != NULL) free(vec_out_extrap);
  if(t_out_close != NULL) free(t_out_close);
  if(t_out_far != NULL) free(t_out_far);
  if(t_out_extrap != NULL) free(t_out_extrap);
  if(which_method != NULL) free(which_method);
  if(err_out_close != NULL) free(err_out_close);
  if(err_out_far != NULL) free(err_out_far);
  if(err_out_extrap != NULL) free(err_out_extrap);
  if(y2 != NULL) free(y2);
  if(u_spline != NULL) free(u_spline);
  
  return;
}


void InterpolateLC(ProgramData *p, int threadid, int lcid, int interp_mode, _InterpParams *params, int Ninterp, double *t_interp, int use_near_far, double minsep, int interp_mode_far, _InterpParams *params_far, int use_extrap, int interp_mode_extrap, _InterpParams *params_extrap)
/* Interpolate the light curve and associated data, string and char
   data will be subject to nearest neighbor interpolation, all other
   data types will be cast to double, interpolated, and recast to
   their original type.

   interp_mode = 
      VARTOOLS_RESAMPLE_NEAREST - take the value of the nearest neighbor.
                                  Rounds down for ties.
      VARTOOLS_RESAMPLE_LINEAR - linear interpolation
      VARTOOLS_RESAMPLE_SPLINE - cubic spline
      VARTOOLS_RESAMPLE_SPLINEMONOTONIC - cubic spline, forced to be monotonic
                                          between sample points.
      VARTOOLS_RESAMPLE_BSPLINE - Basis Spline interpolation.

   params - structure with method-specific parameters for interpolation.

   Ninterp - number of points in the vector of times to resample at.
   t_interp - vector of new times to resample the light curve at.

   use_near_far - 1 if we are using different methods for points close to
                    observations and points far from observations.
                  0 if not.
   minsep - if use_near_far is 1, this is the time difference to distinguish
            between near and far points.

   interp_mode_far - if use_near_far is 1, this is the method to use for the
            far points (the near points will be interpolated using the
            interp_mode method).  Options are the same as for interp_mode.

   params_far - structure with method-specific parameters for the far point
                interpolation method.

   use_extrap - 1 if we are using a different method for extrapolations.
                0 if not.

   interp_mode_extrap - if use_extrap is 1, this is the method to use for
                        extrapolating points. OPtions are the same as for
                        interp_mode.

   params_extrap - structure with method-specific parameters for the 
                   extrapolation.
*/
{
  int i, j, k, Nc;
  double ***dblptr;
  double ****dblptr2;
  short ***shortptr;
  short ****shortptr2;
  int ***intptr;
  int ****intptr2;
  char *charvec_in;
  char *charvec_out = NULL;
  int sizecharvec_out = 0;
  char ***charptr;
  char ****charptr2;
  char **stringvec_in;
  char **stringvec_out = NULL;
  int sizestringvec_out = 0;
  char ****stringptr;
  char *****stringptr2;
  float ***floatptr;
  float ****floatptr2;
  long ***longptr;
  long ****longptr2;
  double *vec_out, *vec_in, *vec_out_err = NULL, *y2 = NULL, *u_spline = NULL;
  double yp1, ypn;
  int size_y2 = 0;
  int is_time_vec, is_err_vec, is_id_vec, is_mag_vec;
  int N_in, u;
  int nbreaks, spline_order;
  _DataFromLightCurve *d;
  va_list varlist;
  /*_Variable *tvar;
  _Variable *magvar;
  _Variable *sigvar;
  _Variable *stringidvar;*/
  
  int interp_mode_near;

  yp1 = params->yp1; ypn = params->ypn;
  nbreaks = params->nbreaks; spline_order = params->spline_order;

  if(use_near_far || use_extrap) {
    if(!use_near_far) {
      /* No distinction between near and far modes. We can make this happen 
	 by setting the minimum separation to a negative number in which
         case all points will be "far" points. */
      minsep = -1.0;
      interp_mode_far = interp_mode;
      params_far = params;
    } 
    interp_mode_near = interp_mode;
    if(!use_extrap) {
      interp_mode_extrap = -1;
    }
    interp_mode = VARTOOLS_RESAMPLE_MULTIPLE;
  }

  N_in = p->NJD[threadid];
  if(N_in < 2) {
    error(ERR_INTERP_NOTENOUGHPOINTS);
  }
  if(Ninterp < 1) {
    /* The light curve is being deleted. Just set NJD to 0 and return */
    p->NJD[threadid] = 0;
    return;
  }
  if((vec_out = (double *) malloc(Ninterp * sizeof(double))) == NULL ||
     (vec_in = (double *) malloc(N_in * sizeof(double))) == NULL) {
    error(ERR_MEMALLOC);
  }

  /* Grow the light curve vectors to store the interpolated data if
     needed */
  /*tvar = NULL; magvar = NULL; sigvar = NULL; stringidvar = NULL;*/
  if(Ninterp > p->NJD[threadid]) {
    MemAllocDataFromLightCurveMidProcess(p, threadid, Ninterp);
    /*for(i=0; i < p->NDefinedVariables; i++) {
      if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_LC && (
	 p->DefinedVariables[i]->datatype == VARTOOLS_TYPE_DOUBLE ||
	 p->DefinedVariables[i]->datatype == VARTOOLS_TYPE_CONVERTJD)) {
	if(p->t[threadid] == (*((double ***) p->DefinedVariables[i]->dataptr))[threadid]) {
	  tvar = p->DefinedVariables[i];
	}
	if(p->mag[threadid] == (*((double ***) p->DefinedVariables[i]->dataptr))[threadid]) {
	  magvar = p->DefinedVariables[i];
	}
	if(p->sig[threadid] == (*((double ***) p->DefinedVariables[i]->dataptr))[threadid]) {
	  sigvar = p->DefinedVariables[i];
	}
      } 
      else if(p->DefinedVariables[i]->vectortype == VARTOOLS_VECTORTYPE_LC &&
	      p->DefinedVariables[i]->datatype == VARTOOLS_TYPE_STRING) {
	if(p->stringid[threadid] == (*((char ****) p->DefinedVariables[i]->dataptr))[threadid])
	  stringidvar = p->DefinedVariables[i];
      }
    }
    MemAllocDataFromLightCurve(p, threadid, Ninterp);
    if(tvar != NULL)
      p->t[threadid] = (*((double ***) tvar->dataptr))[threadid];
    if(magvar != NULL)
      p->mag[threadid] = (*((double ***) magvar->dataptr))[threadid];
    if(sigvar != NULL)
      p->sig[threadid] = (*((double ***) sigvar->dataptr))[threadid];
    if(stringidvar != NULL)
    p->stringid[threadid] = (*((char ****) stringidvar->dataptr))[threadid];*/
  }
    
  /* If we're doing or Gap-interpolation, B-splines, then take care of
     the magnitude interpolation first */
  if(interp_mode == VARTOOLS_RESAMPLE_MULTIPLE) {
    if((vec_out_err = (double *) malloc(Ninterp * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < N_in; i++) {
      vec_in[i] = p->mag[threadid][i];
    }
    InterpolateVec_Gap_Extrap(N_in, p->t[threadid], vec_in,
			      Ninterp, t_interp, vec_out, minsep,
			      interp_mode_near, params,
			      interp_mode_far, params_far,
			      interp_mode_extrap, params_extrap,
			      1, p->sig[threadid], vec_out_err);
    for(j=0; j < Ninterp; j++) {
      p->sig[threadid][j] = vec_out_err[j];
      p->mag[threadid][j] = vec_out[j];
    }
    free(vec_out_err);
  }    
#ifdef _HAVE_GSL
  else if(interp_mode == VARTOOLS_RESAMPLE_BSPLINE) {
    if((vec_out_err = (double *) malloc(Ninterp * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < N_in; i++) {
      vec_in[i] = p->mag[threadid][i];
    }
    if(nbreaks >= 2) {
      InterpolateVec_BSpline(N_in, p->t[threadid], vec_in, Ninterp,
			     t_interp, vec_out, nbreaks, spline_order,
			     1, p->sig[threadid], vec_out_err);
    } else {
      InterpolateVec_BSpline_AutoBreaks(N_in, p->t[threadid], vec_in, Ninterp,
					t_interp, vec_out, &nbreaks, 
					spline_order, 1, p->sig[threadid], 
					vec_out_err);
    }
    for(j=0; j < Ninterp; j++) {
      p->sig[threadid][j] = vec_out_err[j];
      p->mag[threadid][j] = vec_out[j];
    }
    free(vec_out_err);
  }
#endif

  /* Cycle through the light curve data copying each input vector into
     vec_in. Cast the data to double, call the appropriate
     interpolation function for that vector, and copy back to the
     light curve array.  Special handling needs to be done for the
     time, error, and id vectors, as well as for string or char
     data */
  for(k=0;k<p->NDataFromLightCurve;k++)
    {
      d = &(p->DataFromLightCurve[k]);
      Nc = d->Ncolumns;
      for(u=0; u < (Nc == 0 ? 1 : Nc); u++) {
	is_time_vec = 0;
	is_err_vec = 0;
	is_id_vec = 0;
	is_mag_vec = 0;
	if(Nc == 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr = (double ***) d->dataptr;
	    if((*dblptr)[threadid] == p->t[threadid]) {
	      is_time_vec = 1;
	      break;
	    }
	    else if((*dblptr)[threadid] == p->sig[threadid]) {
	      is_err_vec = 1;
	      break;
	    }
	    else if((*dblptr)[threadid] == p->mag[threadid]) {
	      is_mag_vec = 1;
	    }
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (*dblptr)[threadid][j];
	    }
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    dblptr = (double ***) d->dataptr;
	    if((*dblptr)[threadid] == p->t[threadid]) {
	      is_time_vec = 1;
	      break;
	    }
	    else if((*dblptr)[threadid] == p->sig[threadid]) {
	      is_err_vec = 1;
	      break;
	    }
	    else if((*dblptr)[threadid] == p->mag[threadid]) {
	      is_mag_vec = 1;
	    }
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (*dblptr)[threadid][j];
	    }
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr = (char ****) d->dataptr;
	    stringvec_in = ((*stringptr)[threadid]);
	    if(p->stringid != NULL) {
	      if((*stringptr)[threadid] == p->stringid[threadid]) {
		is_id_vec = 1;
		break;
	      }
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr = (int ***) d->dataptr;
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (double) ((*intptr)[threadid][j]);
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr = (short ***) d->dataptr;
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (double) ((*shortptr)[threadid][j]);
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr = (long ***) d->dataptr;
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (double) ((*longptr)[threadid][j]);
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr = (float ***) d->dataptr;
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (double) ((*floatptr)[threadid][j]);
	    }
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr = (char ***) d->dataptr;
	    charvec_in = ((*charptr)[threadid]);
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	} else if(Nc > 0) {
	  switch(d->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    dblptr2 = (double ****) d->dataptr;
	    if((*dblptr2)[threadid][u] == p->t[threadid]) {
	      is_time_vec = 1;
	      break;
	    }
	    else if((*dblptr2)[threadid][u] == p->sig[threadid]) {
	      is_err_vec = 1;
	      break;
	    }
	    else if((*dblptr2)[threadid][u] == p->mag[threadid]) {
	      is_mag_vec = 1;
	    }
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (*dblptr2)[threadid][u][j];
	    }
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    dblptr2 = (double ****) d->dataptr;
	    if((*dblptr2)[threadid][u] == p->t[threadid]) {
	      is_time_vec = 1;
	      break;
	    }
	    else if((*dblptr2)[threadid][u] == p->sig[threadid]) {
	      is_err_vec = 1;
	      break;
	    }
	    else if((*dblptr2)[threadid][u] == p->mag[threadid]) {
	      is_mag_vec = 1;
	    }
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (*dblptr2)[threadid][u][j];
	    }
	    break;
	  case VARTOOLS_TYPE_STRING:
	    stringptr2 = (char *****) d->dataptr;
	    stringvec_in = ((*stringptr2)[threadid][u]);
	    if(p->stringid != NULL) {
	      if((*stringptr2)[threadid][u] == p->stringid[threadid]) {
		is_id_vec = 1;
		break;
	      }
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    intptr2 = (int ****) d->dataptr;
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (double) ((*intptr2)[threadid][u][j]);
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    shortptr2 = (short ****) d->dataptr;
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (double) ((*shortptr2)[threadid][u][j]);
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    longptr2 = (long ****) d->dataptr;
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (double) ((*longptr2)[threadid][u][j]);
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    floatptr2 = (float ****) d->dataptr;
	    for(j=0; j < N_in; j++) {
	      vec_in[j] = (double) ((*floatptr2)[threadid][u][j]);
	    }
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    charptr2 = (char ****) d->dataptr;
	    charvec_in = ((*charptr2)[threadid][u]);
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	}
	/* Do the interpolation for any data that is not a special case */
	if(!is_time_vec && !is_err_vec && d->datatype != VARTOOLS_TYPE_STRING
	   && d->datatype != VARTOOLS_TYPE_CHAR) {
	  switch(interp_mode) {
	  case VARTOOLS_RESAMPLE_NEAREST:
	    InterpolateVec_Nearest(N_in, p->t[threadid], vec_in,
				   Ninterp, t_interp, vec_out);
	    break;
	  case VARTOOLS_RESAMPLE_LINEAR:
	    InterpolateVec_Linear(N_in, p->t[threadid], vec_in,
				  Ninterp, t_interp, vec_out);
	    break;
	  case VARTOOLS_RESAMPLE_SPLINE:
	    InterpolateVec_Spline(N_in, p->t[threadid], vec_in,
				  Ninterp, t_interp, vec_out,
				  yp1, ypn,
				  &y2, &u_spline, &size_y2);
	    break;
	  case VARTOOLS_RESAMPLE_SPLINEMONOTONIC:
	    InterpolateVec_SplineMonotonic(N_in, p->t[threadid], vec_in,
					   Ninterp, t_interp, vec_out,
					   &y2, &size_y2);
	    break;
#ifdef _HAVE_GSL
	  case VARTOOLS_RESAMPLE_BSPLINE:
	    if(!is_mag_vec) {
	      InterpolateVec_BSpline(N_in, p->t[threadid], vec_in, Ninterp,
				     t_interp, vec_out, nbreaks, spline_order,
				     0, NULL, NULL);
	    }
	    break;
#endif
	  case VARTOOLS_RESAMPLE_MULTIPLE:
	    if(!is_mag_vec) {
	      InterpolateVec_Gap_Extrap(N_in, p->t[threadid], vec_in,
					Ninterp, t_interp, vec_out, minsep,
					interp_mode_near, params,
					interp_mode_far, params_far,
					interp_mode_extrap, params_extrap,
					0, NULL, NULL);
	    }
	    break;
	  default:
	    error(ERR_CODEERROR);
	  }
	  /* Copy the data from vec_out back to the light curve array */
	  if(Nc == 0) {
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      dblptr = (double ***) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		(*dblptr)[threadid][j] = vec_out[j];
	      }
	      break;
	    case VARTOOLS_TYPE_CONVERTJD:
	      dblptr = (double ***) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		(*dblptr)[threadid][j] = vec_out[j];
	      }
	      break;
	    case VARTOOLS_TYPE_INT:
	      intptr = (int ***) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		((*intptr)[threadid][j]) = (int) vec_out[j];
	      }
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      shortptr = (short ***) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		((*shortptr)[threadid][j]) = (short) vec_out[j];
	      }
	      break;
	    case VARTOOLS_TYPE_LONG:
	      longptr = (long ***) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		((*longptr)[threadid][j]) = (long) vec_out[j];
	      }
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      floatptr = (float ***) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		((*floatptr)[threadid][j]) = (float) vec_out[j];
	      }
	      break;
	    default:
	      error(ERR_BADTYPE);
	    }
	  } else if(Nc > 0) {
	    switch(d->datatype) {
	    case VARTOOLS_TYPE_DOUBLE:
	      dblptr2 = (double ****) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		(*dblptr2)[threadid][u][j] = vec_out[j];
	      }
	      break;
	    case VARTOOLS_TYPE_CONVERTJD:
	      dblptr2 = (double ****) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		(*dblptr2)[threadid][u][j] = vec_out[j];
	      }
	      break;
	    case VARTOOLS_TYPE_INT:
	      intptr2 = (int ****) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		((*intptr2)[threadid][u][j]) = (int) vec_out[j];
	      }
	      break;
	    case VARTOOLS_TYPE_SHORT:
	      shortptr2 = (short ****) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		((*shortptr2)[threadid][u][j]) = (short) vec_out[j];
	      }
	      break;
	    case VARTOOLS_TYPE_LONG:
	      longptr2 = (long ****) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		((*longptr2)[threadid][u][j]) = (long) vec_out[j];
	      }
	      break;
	    case VARTOOLS_TYPE_FLOAT:
	      floatptr2 = (float ****) d->dataptr;
	      for(j=0; j < Ninterp; j++) {
		((*floatptr2)[threadid][u][j]) = (float) vec_out[j];
	      }
	      break;
	    default:
	      error(ERR_BADTYPE);
	    }
	  }
	}
	else if(d->datatype == VARTOOLS_TYPE_STRING) {
	  if(!sizestringvec_out) {
	    if((stringvec_out = (char **) malloc(Ninterp * sizeof(char *))) == NULL)
	      error(ERR_MEMALLOC);
	    for(j=0; j < Ninterp; j++) {
	      if((stringvec_out[j] = (char *) malloc(MAXLEN*sizeof(char))) == NULL)
		error(ERR_MEMALLOC);
	    }
	    sizestringvec_out = Ninterp;
	  }
	  InterpolateVec_Nearest_String(N_in, p->t[threadid], stringvec_in,
					Ninterp, t_interp, stringvec_out);
	  for(j=0; j < Ninterp; j++) {
	    sprintf(stringvec_in[j],"%s",stringvec_out[j]);
	  }
	}
	else if(d->datatype == VARTOOLS_TYPE_CHAR) {
	  if(!sizecharvec_out) {
	    if((charvec_out = (char *) malloc(Ninterp * sizeof(char))) == NULL)
	      error(ERR_MEMALLOC);
	    sizestringvec_out = Ninterp;
	  }
	  InterpolateVec_Nearest_Char(N_in, p->t[threadid], charvec_in,
					Ninterp, t_interp, charvec_out);
	  for(j=0; j < Ninterp; j++) {
	    charvec_in[j] = charvec_out[j];
	  }
	}
      }
    }
  /* Resample the errors */
  switch(interp_mode) {
  case VARTOOLS_RESAMPLE_NEAREST:
    InterpolateVec_Nearest(N_in, p->t[threadid], p->sig[threadid],
			   Ninterp, t_interp, vec_out);
    break;
  case VARTOOLS_RESAMPLE_SPLINE:
  case VARTOOLS_RESAMPLE_SPLINEMONOTONIC:
  case VARTOOLS_RESAMPLE_LINEAR:
    /* Fudge the errors assuming linear interpolation when doing pure
       cubic spline interpolation. This is wrong, but one wouldn't be
       doing spline interpolation if one cared about having correct
       uncertainties anyway. */
    InterpolateErr_Linear(N_in, p->t[threadid], p->sig[threadid],
			  Ninterp, t_interp, vec_out);
    break;
#ifdef _HAVE_GSL
  case VARTOOLS_RESAMPLE_BSPLINE:
    break;
#endif
  case VARTOOLS_RESAMPLE_MULTIPLE:
    break;
  default:
    error(ERR_CODEERROR);
  }

  /* Resample the times */
  for(j=0; j < Ninterp; j++) {
    p->t[threadid][j] = t_interp[j];
  }
  p->NJD[threadid] = Ninterp;

  /* Reset the string id index if needed */
  if(p->readimagestring) {
    for(j=0; j < Ninterp; j++) {
      p->stringid_idx[threadid][j] = j;
    }
    mysortstringint(p->NJD[threadid], MAXIDSTRINGLENGTH,
		    p->stringid[threadid], p->stringid_idx[threadid]);
  }
    
  
  if(y2 != NULL)
    free(y2);
  if(u_spline != NULL)
    free(u_spline);
  free(vec_in);
  free(vec_out);
  if(charvec_out != NULL) free(charvec_out);
  if(stringvec_out != NULL) {
    for(j=0; j < Ninterp; j++) {
      free(stringvec_out[j]);
    }
    free(stringvec_out);
  }
  return;
}

void DoResample(ProgramData *p, _Resample *c, int threadid, int lcid) 
/* This function executes the -resample command for the light curve #lcid 
   being processed in thread #threadid. The function generates the vector of
   times to resample the light curve at, and then calls the InterpolateLC
   function using an additional arguments as needed, depending on the
   interpolation method being used. */ 
{
  double *t_out = NULL, tstart, tstop, delt;
  double *delt_in_vals = NULL;
  int size_tout, i, k, j;
  int Ninterp, Nresamp;
  FILE *infile;
  char *line = NULL;
  size_t line_size = MAXLEN;
  double minsep, min_true_delt, med_true_delt;

  _InterpParams params;
  _InterpParams params_far;
  _InterpParams params_extrap;

  /* Sort the input light curve in time, and merge any equal points */
  if(sortlcbytime(p->NJD[threadid], p->t[threadid], threadid, p)) {
    mergeequallctimes(p, threadid);
  }

  /* Check if the light curve is too short for interpolation */
  if(p->NJD[threadid] < 2) {
    /* Cannot interpolate */
    error(ERR_INTERP_NOTENOUGHPOINTS);
  }


  /* Setup the vector of times to resample the light curve at */
  
  if(c->use_file) {
    /* The times were specified in a file */
    if(c->resample_filename_source != VARTOOLS_SOURCE_FIXED) {
      /* Each light curve has its own file to use */
      if((infile = fopen(c->resample_filenames[lcid],"r")) == NULL) {
	error2(ERR_CANNOTOPEN,c->resample_filenames[lcid]);
      }
      line = malloc(line_size);
      size_tout = p->NJD[threadid];
      if((t_out = (double *) malloc(size_tout * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      Ninterp = 0;
      while(gnu_getline(&line,&line_size,infile) >= 0) {
	if(line[0] != '#') {
	  if(Ninterp >= size_tout) {
	    size_tout *= 2;
	    if((t_out = (double *) realloc(t_out, size_tout * sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  k = 1;
	  j = 0;
	  while( k < c->t_column) {
	    j += skipone(&line[j]);
	    k++;
	  }
	  if(line[j] != '\0' && line[j] != '\n')
	    {
	      j += parseone(&line[j], (void *) (&(t_out[Ninterp])), VARTOOLS_TYPE_DOUBLE);
	    }
	  else {
	    error2(ERR_INPUTMISSINGCOLUMN, c->resample_filenames[lcid]);
	  }
	  Ninterp++;
	}
      }
      fclose(infile);
    } else {
      /* All light curves use the same file, get the number of times, and
	 the vector from the command structure */
      Ninterp = c->N_resamp;
      t_out = c->t_resamp;
    }
  } else {
    /* Get the start time */
    if(c->tstart_source == VARTOOLS_SOURCE_FIXED) {
      tstart = c->tstart_fix;
    } else if(c->tstart_source == VARTOOLS_SOURCE_COMPUTED) {
      tstart = p->t[threadid][0];
    } else if(c->tstart_source == VARTOOLS_SOURCE_INLIST) {
      tstart = c->tstart[lcid];
    } else if(c->tstart_source == VARTOOLS_SOURCE_PRIORCOLUMN) {
      getoutcolumnvalue(c->tstart_linkedcolumn, threadid, lcid,
			VARTOOLS_TYPE_DOUBLE, &tstart);
    } else if(c->tstart_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
      tstart = EvaluateExpression(lcid, threadid, 0, c->tstart_expr);
    }

    /* Get the stop time */
    if(c->tstop_source == VARTOOLS_SOURCE_FIXED) {
      tstop = c->tstop_fix;
    } else if(c->tstop_source == VARTOOLS_SOURCE_COMPUTED) {
      tstop = p->t[threadid][p->NJD[threadid]-1];
    } else if(c->tstop_source == VARTOOLS_SOURCE_INLIST) {
      tstop = c->tstop[lcid];
    } else if(c->tstop_source == VARTOOLS_SOURCE_PRIORCOLUMN) {
      getoutcolumnvalue(c->tstop_linkedcolumn, threadid, lcid,
			VARTOOLS_TYPE_DOUBLE, &tstop);
    } else if(c->tstop_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
      tstop = EvaluateExpression(lcid, threadid, 0, c->tstop_expr);
    }

    if(tstop <= tstart) {
      /* light curve is being deleted */
      p->NJD[threadid] = 0;
      return;
    }

    if(c->Nresamp_source > -1) {
      /* The number of points to use has been specified */
      if(c->Nresamp_source == VARTOOLS_SOURCE_FIXED) {
	Nresamp = c->Nresamp_fix;
      }
      else if (c->Nresamp_source == VARTOOLS_SOURCE_INLIST) {
	Nresamp = c->Nresamp[lcid];
      } else if(c->Nresamp_source == VARTOOLS_SOURCE_PRIORCOLUMN) {
	getoutcolumnvalue(c->Nresamp_linkedcolumn, threadid, lcid,
			  VARTOOLS_TYPE_INT, &Nresamp);
      } else if(c->Nresamp_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	Nresamp = (int) EvaluateExpression(lcid, threadid, 0, c->Nresamp_expr);
      }

      if(Nresamp > 1) {
	delt = (tstop - tstart)/(Nresamp - 1);
      } else {
	delt = 2*tstop - tstart;
      }
    }
    else {
      /* The separation has been specified */
      if(c->delt_source == VARTOOLS_SOURCE_FIXED) {
	delt = c->delt_fix;
      }
      else if(c->delt_source = VARTOOLS_SOURCE_COMPUTED) {
	/* Find the minimum non-zero separation in the lc */
	delt = p->t[threadid][1]-p->t[threadid][0];
	for(j=2; j < p->NJD[threadid]; j++) {
	  if(p->t[threadid][j]-p->t[threadid][j-1] < delt) {
	    delt = p->t[threadid][j]-p->t[threadid][j-1];
	  }
	}
      }
      else if(c->delt_source == VARTOOLS_SOURCE_INLIST) {
	delt = c->delt[lcid];
      } else if(c->delt_source == VARTOOLS_SOURCE_PRIORCOLUMN) {
	getoutcolumnvalue(c->delt_linkedcolumn, threadid, lcid,
			  VARTOOLS_TYPE_DOUBLE, &delt);
      } else if(c->delt_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
	delt = EvaluateExpression(lcid, threadid, 0, c->delt_expr);
      }

    }
    Ninterp = floor((tstop - tstart)/delt) + 1;
    if((t_out = (double *) malloc(Ninterp * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    for(i=1, t_out[0] = tstart; i < Ninterp; i++) {
      t_out[i] = t_out[i-1] + delt;
    }
  }

  /* Depending on the interpolation scheme, find optional arguments */
  switch(c->resample_method) {
  case VARTOOLS_RESAMPLE_NEAREST:
  case VARTOOLS_RESAMPLE_LINEAR:
  case VARTOOLS_RESAMPLE_SPLINEMONOTONIC:
    /* No additional parameters */
    break;
  case VARTOOLS_RESAMPLE_SPLINE:
    params.yp1 = c->yp1; params.ypn = c->ypn;
    break;
  case VARTOOLS_RESAMPLE_BSPLINE:
    params.nbreaks = c->bspline_nbreaks;
    params.spline_order = c->bspline_order;
    break;
  default:
    error(ERR_CODEERROR);
    break;
  }
  
  /* Check if we're doing different near/far interpolation */
  if(c->use_near_far) {
    switch(c->resample_method_far) {
    case VARTOOLS_RESAMPLE_NEAREST:
    case VARTOOLS_RESAMPLE_LINEAR:
    case VARTOOLS_RESAMPLE_SPLINEMONOTONIC:
      /* No additional parameters */
      break;
    case VARTOOLS_RESAMPLE_SPLINE:
      params_far.yp1 = c->yp1_far; params_far.ypn = c->ypn_far;
      break;
    case VARTOOLS_RESAMPLE_BSPLINE:
      params_far.nbreaks = c->bspline_nbreaks_far;
      params_far.spline_order = c->bspline_order_far;
      break;
    default:
      error(ERR_CODEERROR);
      break;
    }

    /* Determine the separation between near and far points */
    /* Get the start time */
    if(c->minsep_source == VARTOOLS_SOURCE_FIXED) {
      minsep = c->minsep_fix;
    } else if(c->minsep_source == VARTOOLS_SOURCE_INLIST) {
      minsep = c->minsep[lcid];
    } else if(c->minsep_source == VARTOOLS_SOURCE_PRIORCOLUMN) {
      getoutcolumnvalue(c->minsep_linkedcolumn, threadid, lcid,
			VARTOOLS_TYPE_DOUBLE, &minsep);
    } else if(c->minsep_source == VARTOOLS_SOURCE_EVALEXPRESSION) {
      minsep = EvaluateExpression(lcid, threadid, 0, c->minsep_expr);
    } else if(c->minsep_source < 0) {
      /*** user gave the "frac_min_sep", "frac_med_sep" or "percentile_sep"
	   option ****/
      if((delt_in_vals = (double *) malloc((p->NJD[threadid]-1)*sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
      for(i=1; i < p->NJD[threadid]; i++) {
	delt_in_vals[i-1] = p->t[threadid][i] - p->t[threadid][i-1];
      }
      if(c->minsep_source == -1) {
	/* frac_min_sep */
	minsep = c->frac_min_sep_val*(getminimum((p->NJD[threadid]-1),delt_in_vals));
      }
      else if(c->minsep_source == -2) {
	/* frac_med_sep */
	minsep = c->frac_med_sep_val*(median(p->NJD[threadid]-1,delt_in_vals));
      }
      else if(c->minsep_source == -3) {
	/* percentile_sep */
	minsep = percentile(p->NJD[threadid]-1,delt_in_vals,c->percentile_sep);
      }
    }  
  }

  /* Check if we're using a different method for extrapolation */
  if(c->use_extrap) {
    switch(c->resample_method_extrap) {
    case VARTOOLS_RESAMPLE_NEAREST:
    case VARTOOLS_RESAMPLE_LINEAR:
    case VARTOOLS_RESAMPLE_SPLINEMONOTONIC:
      /* No additional parameters */
      break;
    case VARTOOLS_RESAMPLE_SPLINE:
      params_extrap.yp1 = c->yp1_extrap; params_extrap.ypn = c->ypn_extrap;
      break;
    case VARTOOLS_RESAMPLE_BSPLINE:
      params_extrap.nbreaks = c->bspline_nbreaks_extrap;
      params_extrap.spline_order = c->bspline_order_extrap;
      break;
    default:
      error(ERR_CODEERROR);
      break;
    }
  }

  
  InterpolateLC(p, threadid, lcid, c->resample_method, &params, Ninterp, t_out,
		c->use_near_far, minsep, c->resample_method_far, &params_far, 
		c->use_extrap, c->resample_method_extrap, &params_extrap);
    

  if(line != NULL)
    free(line);
  if(t_out != NULL)
    free(t_out);
  if(delt_in_vals != NULL)
    free(delt_in_vals);
  return;
}

int ParseResampleMethod(int *iret, int argc, char **argv, ProgramData *p,
			_Resample *c, int cnum, int *resample_method,
			double *yp1, double *ypn, int *bspline_order,
			int *bspline_nbreaks)
{

  int i;
  i = *iret;

  if(!strcmp(argv[i],"nearest")) {
    *(resample_method) = VARTOOLS_RESAMPLE_NEAREST;
  }
  else if(!strcmp(argv[i],"linear")) {
    *(resample_method) = VARTOOLS_RESAMPLE_LINEAR;
  }
  else if(!strcmp(argv[i],"spline")) {
    *(resample_method) = VARTOOLS_RESAMPLE_SPLINE;
    *(yp1) = 1.0e30;
    *(ypn) = 1.0e30;
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"left")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	*yp1 = atof(argv[i]);
      } else 
	i--;
    } else
      i--;
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"right")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	*ypn = atof(argv[i]);
      } else
	i--;
    } else
      i--;
  }
  else if(!strcmp(argv[i],"splinemonotonic")) {
    *resample_method = VARTOOLS_RESAMPLE_SPLINEMONOTONIC;
  }
#ifdef _HAVE_GSL
  else if(!strcmp(argv[i],"bspline")) {
    *resample_method = VARTOOLS_RESAMPLE_BSPLINE;
    *bspline_order = 4;
    *bspline_nbreaks = 15;
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"nbreaks")){
	  i++;
	  if(i >= argc) {*iret = i; return 1;}
	  *bspline_nbreaks = atoi(argv[i]);
      } else
	i--;
    }
    else
      i--;
    i++;
    if(i < argc) {
      if(!strcmp(argv[i],"order")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	*bspline_order = atoi(argv[i]) + 1;
	if(*bspline_order < 2) {
	  error2(ERR_INVALID_PARAMETERVALUE,"-resample, the order for the bspline must be >= 1");
	}
      } else
	i--;
    } else
      i--;
  }
#endif
  else {
    *iret = i; return 1;
  }

  *iret = i; return 0;

}

int ParseResampleCommand(int *iret, int argc, char **argv, ProgramData *p,
			 _Resample *c, int cnum)
/* Parse the command line for the "-resample" command */
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

  if(ParseResampleMethod(&i, argc, argv, p, c, cnum, &(c->resample_method),
			 &(c->yp1), &(c->ypn), &(c->bspline_order),
			 &(c->bspline_nbreaks))) {
    *iret = i;
    return 1;
  }


  /* Set the default resample options */
  c->tstart_source = VARTOOLS_SOURCE_COMPUTED;
  c->tstop_source = VARTOOLS_SOURCE_COMPUTED;
  c->Nresamp_source = -1;
  c->use_file = 0;
  c->delt_source = VARTOOLS_SOURCE_COMPUTED;
  
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"file")) {
      c->use_file = 1;
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if(!strcmp(argv[i],"fix")) {
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->resample_filename_source = VARTOOLS_SOURCE_FIXED;
	c->resample_filename_fix = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
	sprintf(c->resample_filename_fix,"%s",argv[i]);
	colnum = 1;
	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"column")) {
	    i++;
	    if(i >= argc) {*iret = i; return 1;}
	    colnum = atoi(argv[i]);
	  }
	  else i--;
	} else i--;
	/* Read in the input times from the file */
	if((infile = fopen(c->resample_filename_fix,"r")) == NULL) {
	  error2(ERR_CANNOTOPEN,c->resample_filename_fix);
	}
	line = malloc(line_size);
	size_tout = 1024;
	if((c->t_resamp = (double *) malloc(size_tout * sizeof(double))) == NULL)
	  error(ERR_MEMALLOC);
	c->N_resamp = 0;
	while(gnu_getline(&line,&line_size,infile) >= 0) {
	  if(line[0] != '#') {
	    if(c->N_resamp >= size_tout) {
	      size_tout *= 2;
	      if((c->t_resamp = (double *) realloc(c->t_resamp, size_tout * sizeof(double))) == NULL)
		error(ERR_MEMALLOC);
	    }
	    k = 1;
	    j = 0;
	    while( k < colnum) {
	      j += skipone(&line[j]);
	      k++;
	    }
	    if(line[j] != '\0' && line[j] != '\n')
	      {
		j += parseone(&line[j], (void *) (&(c->t_resamp[c->N_resamp])), VARTOOLS_TYPE_DOUBLE);
	      }
	    else {
	      error2(ERR_INPUTMISSINGCOLUMN, c->resample_filename_fix);
	    }
	    c->N_resamp++;
	  }
	}
	fclose(infile);
      } 
      else if(!strcmp(argv[i],"list")) {
	c->resample_filename_source = VARTOOLS_SOURCE_INLIST;
	c->t_column = 1;
	colnum = 0;
	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"listcolumn")) {
	    i++;
	    if(i >= argc) {*iret = i; return 1;}
	    colnum = atoi(argv[i]);
	  } else
	    i--;
	} else
	  i--;
	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"tcolumn")) {
	    i++;
	    if(i >= argc) {*iret = i; return 1;}
	    c->t_column = atoi(argv[i]);
	  } else
	    i--;
	} else
	  i--;
	RegisterDataFromInputList(p,
				  (void *) (&(c->resample_filenames)),
				  VARTOOLS_TYPE_STRING,
				  0, cnum, 0, 0, NULL, colnum,
				  "RESAMPLE_TIMEFILE");
      } else {
	*iret = i; return 1;
      }
    } else
      i--;
  } else
    i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"tstart")) {
      if(c->use_file) {*iret = i; return 1;}
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if(!strcmp(argv[i],"fix")) {
	c->tstart_source = VARTOOLS_SOURCE_FIXED;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->tstart_fix = atof(argv[i]);
      }
      else if(!strcmp(argv[i],"fixcolumn")) {
	c->tstart_source = VARTOOLS_SOURCE_PRIORCOLUMN;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	increaselinkedcols(p, &(c->tstart_linkedcolumn), argv[i], cnum);
      }
      else if(!strcmp(argv[i],"list")) {
	c->tstart_source = VARTOOLS_SOURCE_INLIST;
	colnum = 0;
	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"column")) {
	    i++;
	    if(i >= argc) {*iret = i; return 1;}
	    colnum = atoi(argv[i]);
	  } else i--;
	} else i--;
	RegisterDataFromInputList(p,
				  (void *) (&(c->tstart)),
				  VARTOOLS_TYPE_DOUBLE,
				  0, cnum, 0, 0, NULL, colnum,
				  "RESAMPLE_TSTART");
      }
      else if(!strcmp(argv[i],"expr")) {
	c->tstart_source = VARTOOLS_SOURCE_EVALEXPRESSION;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->tstart_exprstring = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
	sprintf(c->tstart_exprstring,"%s",argv[i]);
      }
      else { *iret = i; return 1; }
    } else i--;
  } else i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"tstop")) {
      if(c->use_file) {*iret = i; return 1;}
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if(!strcmp(argv[i],"fix")) {
	c->tstop_source = VARTOOLS_SOURCE_FIXED;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->tstop_fix = atof(argv[i]);
      }
      else if(!strcmp(argv[i],"fixcolumn")) {
	c->tstop_source = VARTOOLS_SOURCE_PRIORCOLUMN;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	increaselinkedcols(p, &(c->tstop_linkedcolumn), argv[i], cnum);
      }
      else if(!strcmp(argv[i],"list")) {
	c->tstop_source = VARTOOLS_SOURCE_INLIST;
	colnum = 0;
	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"column")) {
	    i++;
	    if(i >= argc) {*iret = i; return 1;}
	    colnum = atoi(argv[i]);
	  } else i--;
	} else i--;
	RegisterDataFromInputList(p,
				  (void *) (&(c->tstop)),
				  VARTOOLS_TYPE_DOUBLE,
				  0, cnum, 0, 0, NULL, colnum,
				  "RESAMPLE_TSTOP");
      }
      else if(!strcmp(argv[i],"expr")) {
	c->tstop_source = VARTOOLS_SOURCE_EVALEXPRESSION;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->tstop_exprstring = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
	sprintf(c->tstop_exprstring,"%s",argv[i]);
      }
      else { *iret = i; return 1; }
    } else i--;
  } else i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"delt")) {
      if(c->use_file) {*iret = i; return 1;}
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if(!strcmp(argv[i],"fix")) {
	c->delt_source = VARTOOLS_SOURCE_FIXED;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->delt_fix = atof(argv[i]);
      }
      else if(!strcmp(argv[i],"fixcolumn")) {
	c->delt_source = VARTOOLS_SOURCE_PRIORCOLUMN;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	increaselinkedcols(p, &(c->delt_linkedcolumn), argv[i], cnum);
      }
      else if(!strcmp(argv[i],"list")) {
	c->delt_source = VARTOOLS_SOURCE_INLIST;
	colnum = 0;
	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"column")) {
	    i++;
	    if(i >= argc) {*iret = i; return 1;}
	    colnum = atoi(argv[i]);
	  } else i--;
	} else i--;
	RegisterDataFromInputList(p,
				  (void *) (&(c->delt)),
				  VARTOOLS_TYPE_DOUBLE,
				  0, cnum, 0, 0, NULL, colnum,
				  "RESAMPLE_DELTA_T");
      }
      else if(!strcmp(argv[i],"expr")) {
	c->delt_source = VARTOOLS_SOURCE_EVALEXPRESSION;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->delt_exprstring = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
	sprintf(c->delt_exprstring,"%s",argv[i]);
      }
      else { *iret = i; return 1; }
    } else i--;
  } else i--;

  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"Npoints")) {
      if(c->use_file) {*iret = i; return 1;}
      if(c->delt_source != VARTOOLS_SOURCE_COMPUTED) {*iret = i; return 1;}
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if(!strcmp(argv[i],"fix")) {
	c->Nresamp_source = VARTOOLS_SOURCE_FIXED;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->Nresamp_fix = atoi(argv[i]);
      }
      else if(!strcmp(argv[i],"fixcolumn")) {
	c->Nresamp_source = VARTOOLS_SOURCE_PRIORCOLUMN;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	increaselinkedcols(p, &(c->Nresamp_linkedcolumn), argv[i], cnum);
      }
      else if(!strcmp(argv[i],"list")) {
	c->Nresamp_source = VARTOOLS_SOURCE_INLIST;
	colnum = 0;
	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"column")) {
	    i++;
	    if(i >= argc) {*iret = i; return 1;}
	    colnum = atoi(argv[i]);
	  } else i--;
	} else i--;
	RegisterDataFromInputList(p,
				  (void *) (&(c->Nresamp)),
				  VARTOOLS_TYPE_INT,
				  0, cnum, 0, 0, NULL, colnum,
				  "RESAMPLE_NPOINTS");
      }
      else if(!strcmp(argv[i],"expr")) {
	c->Nresamp_source = VARTOOLS_SOURCE_EVALEXPRESSION;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->Nresamp_exprstring = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
	sprintf(c->Nresamp_exprstring,"%s",argv[i]);
      }
      else { *iret = i; return 1; }
    } else i--;
  } else i--;

  c->use_near_far = 0;
  c->resample_method_far = -1;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"gaps")) {
      c->use_near_far = 1;
      
      i++;
      if(i >= argc) {*iret = i; return 1;}
      if(!strcmp(argv[i],"fix")) {
	c->minsep_source = VARTOOLS_SOURCE_FIXED;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->minsep_fix = atof(argv[i]);
      }
      else if(!strcmp(argv[i],"fixcolumn")) {
	c->minsep_source = VARTOOLS_SOURCE_PRIORCOLUMN;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	increaselinkedcols(p, &(c->minsep_linkedcolumn), argv[i], cnum);
      }
      else if(!strcmp(argv[i],"list")) {
	c->minsep_source = VARTOOLS_SOURCE_INLIST;
	colnum = 0;
	i++;
	if(i < argc) {
	  if(!strcmp(argv[i],"column")) {
	    i++;
	    if(i >= argc) {*iret = i; return 1;}
	    colnum = atoi(argv[i]);
	  } else i--;
	} else i--;
	RegisterDataFromInputList(p,
				  (void *) (&(c->minsep)),
				  VARTOOLS_TYPE_DOUBLE,
				  0, cnum, 0, 0, NULL, colnum,
				  "RESAMPLE_MINSEP");
      }
      else if(!strcmp(argv[i],"expr")) {
	c->minsep_source = VARTOOLS_SOURCE_EVALEXPRESSION;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->minsep_exprstring = (char *) malloc((strlen(argv[i])+1)*sizeof(char));
	sprintf(c->minsep_exprstring,"%s",argv[i]);
      }
      else if(!strcmp(argv[i],"frac_min_sep")) {
	c->minsep_source = -1;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->frac_min_sep_val = atof(argv[i]);
      }
      else if(!strcmp(argv[i],"frac_med_sep")) {
	c->minsep_source = -2;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->frac_med_sep_val = atof(argv[i]);
      }
      else if(!strcmp(argv[i],"percentile_sep")) {
	c->minsep_source = -3;
	i++;
	if(i >= argc) {*iret = i; return 1;}
	c->percentile_sep = atof(argv[i]);
      }
      else { *iret = i; return 1; }
      
      i++;
      if(ParseResampleMethod(&i, argc, argv, p, c, cnum, 
			     &(c->resample_method_far),
			     &(c->yp1_far), &(c->ypn_far), 
			     &(c->bspline_order_far),
			     &(c->bspline_nbreaks_far))) {
	*iret = i;
	return 1;
      }
	 
    } else {
      i--;
    }
  } else
    i--;

  c->use_extrap = 0;
  c->resample_method_extrap = -1;
  i++;
  if(i < argc) {
    if(!strcmp(argv[i],"extrap")) {
      i++;
      c->use_extrap = 1;
      if(ParseResampleMethod(&i, argc, argv, p, c, cnum, 
			     &(c->resample_method_extrap),
			     &(c->yp1_extrap), &(c->ypn_extrap), 
			     &(c->bspline_order_extrap),
			     &(c->bspline_nbreaks_extrap))) {
	*iret = i;
	return 1;
      }
    } else
      i--;
  } else
    i--;
      
  *iret = i; return 0;
}
