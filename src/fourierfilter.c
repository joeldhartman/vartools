/* This file contains functions which implement the -fourierfilter
   command for VARTOOLS.

   The entire implementation depends on GSL's mixed-radix complex FFT
   (`gsl_fft_complex_forward` / `gsl_fft_complex_inverse`).  Everything
   below the `#ifdef _HAVE_GSL` guard is built only when GSL was
   detected at configure time.  If GSL is absent, stubs at the bottom
   of this file supply the ParseFourierFilterCommand / doFourierFilter /
   SetupFourierFilterExpression symbols that the rest of vartools
   expects; they emit a clear error message and do nothing else so the
   command is effectively unavailable without breaking the link. */
#include "commands.h"
#include "programdata.h"
#include "functions.h"
#include <ctype.h>
#include <string.h>

#ifdef _HAVE_GSL
#include <gsl/gsl_fft_complex.h>

#ifndef TWOPI
#define TWOPI 6.28318530717958647692528676656
#endif

#ifndef PI
#define PI 3.14159265358979323846
#endif

/* Taper types for -fourierfilter taper <name>.  NONE means no taper is
   applied (brick-wall cut).  All listed tapers produce a transition
   that goes 0 -> 1 monotonically across the taper half-window, so they
   can be used interchangeably at a band edge; they differ in the
   trade-off between transition steepness and sidelobe suppression. */
#define VARTOOLS_FOURIERFILTER_TAPER_NONE     0
#define VARTOOLS_FOURIERFILTER_TAPER_LINEAR   1
#define VARTOOLS_FOURIERFILTER_TAPER_COSINE   2
#define VARTOOLS_FOURIERFILTER_TAPER_BLACKMAN 3
#define VARTOOLS_FOURIERFILTER_TAPER_KAISER   4

/* Sentinel values for _FourierFilter::resample_source and gapbreak_source
   beyond the standard VARTOOLS_SOURCE_* enum.  The _DELMIN value for
   resample_source means "use the minimum dt in the light curve".  The
   _FRAC_* and _PERCENTILE values for gapbreak_source follow the
   convention used by -resample's "gaps" keyword. */
#define VARTOOLS_FOURIERFILTER_RESAMPLE_DELMIN        -1
#define VARTOOLS_FOURIERFILTER_GAPBREAK_FRAC_MIN_SEP  -1
#define VARTOOLS_FOURIERFILTER_GAPBREAK_FRAC_MED_SEP  -2
#define VARTOOLS_FOURIERFILTER_GAPBREAK_PERCENTILE    -3

/* padmode values */
#define VARTOOLS_FOURIERFILTER_PAD_WRAP     0
#define VARTOOLS_FOURIERFILTER_PAD_REFLECT  1
#define VARTOOLS_FOURIERFILTER_PAD_ZERO     2

/* Modified Bessel function I_0(x) via Abramowitz & Stegun 9.8 — a
   7th-order polynomial approximation accurate to ~1.6e-7 for all x.
   Used only by the Kaiser taper, so pulling in gsl_sf_bessel_I0 is
   not worth the link dependency. */
static double _bessel_I0(double x) {
  double ax, y, ans;
  ax = fabs(x);
  if (ax < 3.75) {
    y = x / 3.75;
    y = y * y;
    ans = 1.0 + y*(3.5156229 + y*(3.0899424 + y*(1.2067492 +
          y*(0.2659732 + y*(0.0360768 + y*0.0045813)))));
  } else {
    y = 3.75 / ax;
    ans = (exp(ax) / sqrt(ax)) *
          (0.39894228 + y*(0.01328592 + y*(0.00225319 + y*(-0.00157565 +
          y*(0.00916281 + y*(-0.02057706 + y*(0.02635537 + y*(-0.01647633 +
          y*0.00392377))))))));
  }
  return ans;
}

/* Edge taper: returns a value in [0,1] that rises from 0 to 1 as f
   moves from (edge - delta) to (edge + delta), using the requested
   shape.  direction=+1 means the response goes 0 -> 1 with increasing
   f (appropriate for a highpass edge or the low-frequency edge of a
   bandpass); direction=-1 reverses it (lowpass or the high-frequency
   edge of a bandpass).  Outside the transition window the value is
   clamped to 0 (rejection side) or 1 (pass side). */
static double _edge_taper(double f, double edge, int direction,
                          int taper_type, double delta, double beta) {
  double x, u, ans;
  if (delta <= 0.0 || taper_type == VARTOOLS_FOURIERFILTER_TAPER_NONE) {
    /* Degenerate: brick wall */
    return (direction > 0 ? (f >= edge ? 1.0 : 0.0)
                          : (f <= edge ? 1.0 : 0.0));
  }
  x = (double)direction * (f - edge);   /* ranges -delta..0..+delta in window */
  if (x <= -delta) return 0.0;
  if (x >=  delta) return 1.0;
  u = (x + delta) / (2.0 * delta);      /* 0..1 across transition */
  switch (taper_type) {
  case VARTOOLS_FOURIERFILTER_TAPER_LINEAR:
    return u;
  case VARTOOLS_FOURIERFILTER_TAPER_COSINE:
    return 0.5 * (1.0 - cos(PI * u));
  case VARTOOLS_FOURIERFILTER_TAPER_BLACKMAN:
    return 0.42 - 0.5*cos(PI*u) + 0.08*cos(TWOPI*u);
  case VARTOOLS_FOURIERFILTER_TAPER_KAISER:
    /* Kaiser half-window: I_0(beta*sqrt(1-(1-u)^2)) / I_0(beta).
       Reaches exactly 1 at u=1 and approaches 0 at u=0 for large beta. */
    ans = _bessel_I0(beta * sqrt(1.0 - (1.0 - u) * (1.0 - u))) / _bessel_I0(beta);
    return ans;
  default:
    return 1.0;
  }
}

/* Composite frequency-domain mask for the whole -fourierfilter response.
   Returns 1 inside the pass band (as defined by filtertype, minfreq,
   maxfreq), 0 in the stop band, and a taper-weighted value in the
   transition region.  For BANDCUT the logic is inverted.  */
static double _fourierfilter_mask(double f, int filtertype,
                                  double minfreq, double maxfreq,
                                  int taper_type, double delta, double beta) {
  double lo, hi;
  switch (filtertype) {
  case VARTOOLS_FOURIERFILTER_FULLSPEC:
    return 1.0;   /* taper is a no-op without cut edges */
  case VARTOOLS_FOURIERFILTER_HIGHPASS:
    return _edge_taper(f, minfreq, +1, taper_type, delta, beta);
  case VARTOOLS_FOURIERFILTER_LOWPASS:
    return _edge_taper(f, maxfreq, -1, taper_type, delta, beta);
  case VARTOOLS_FOURIERFILTER_BANDPASS:
    lo = _edge_taper(f, minfreq, +1, taper_type, delta, beta);
    hi = _edge_taper(f, maxfreq, -1, taper_type, delta, beta);
    return lo * hi;
  case VARTOOLS_FOURIERFILTER_BANDCUT:
    lo = _edge_taper(f, minfreq, -1, taper_type, delta, beta);
    hi = _edge_taper(f, maxfreq, +1, taper_type, delta, beta);
    /* At f below minfreq, lo=1, hi=0 -> return 1.  At f above maxfreq,
       lo=0, hi=1 -> return 1.  Inside the cut band both go through
       their transitions; 1 - min(lo, hi) would also work but the
       disjoint-edge formulation generalizes cleanly. */
    return 1.0 - (1.0 - lo) * (1.0 - hi);
  default:
    return 1.0;
  }
}


/* Linear interpolation of v_src(t_src) onto the points t_dst.  Beyond
   the edges the nearest value is extrapolated flat.  t_src must be
   strictly monotonic increasing.  O(N_src + N_dst) - the two time
   arrays share a common forward walk. */
static void _ff_linear_interp(int N_src, const double *t_src,
                              const double *v_src,
                              int N_dst, const double *t_dst, double *v_dst)
{
  int i, j = 0;
  double x, u;
  for(i = 0; i < N_dst; i++) {
    x = t_dst[i];
    while(j + 1 < N_src && t_src[j + 1] < x) j++;
    if(x <= t_src[0]) {
      v_dst[i] = v_src[0];
    } else if(j + 1 >= N_src) {
      v_dst[i] = v_src[N_src - 1];
    } else {
      u = (x - t_src[j]) / (t_src[j + 1] - t_src[j]);
      v_dst[i] = v_src[j] * (1.0 - u) + v_src[j + 1] * u;
    }
  }
}

/* Core FFT-based filter kernel.

   Operates on a uniformly-sampled signal of length N at spacing
   `delta`.  The signal's time origin is implicit - samples are
   treated as being at (0, delta, 2*delta, ..., (N-1)*delta), so the
   FFT's natural bin grid 1/(N*delta) aligns with df used by the rest
   of the code.  Coefficients are computed at that grid regardless of
   any absolute-time offset in the caller's data.

   Steps: mean-subtract, optional edge-pad (padmode), complex FFT,
   multiply each bin by the band+taper+filterexpr mask (or `1-mask`
   when invert_mask=1, used by subtract-mode highpass/bandcut), complex
   IFFT, discard padding, add mean back.  Writes the resulting length-N
   signal into signal_out.

   If coeffs_a_out / coeffs_b_out are non-NULL, they are filled with
   the *pre-mask* cos/sin coefficients at bins 0..Ntot/2 inclusive,
   where Ntot is the padded FFT length.  Caller uses these for the
   "ofourier" coefficient dump and labels bin k at frequency k/(Ntot*delta).

   Returns Ntot, the FFT length used including any padding. */
static int _ff_fft_apply_mask(
    int N, const double *signal_in, double dc_level, double delta,
    int filtertype, double minfreq, double maxfreq,
    int taper_type, double taper_delta, double taper_beta,
    int has_filterexpr, _Variable *freq_var, _Expression *filter_expr,
    int padmode, double padfrac,
    int invert_mask,
    int lcid, int threadid,
    double *signal_out,
    double *coeffs_a_out, double *coeffs_b_out)
{
  int Npad, Ntot, k, i;
  double *tmpdata = NULL;
  gsl_fft_complex_wavetable *wt;
  gsl_fft_complex_workspace *ws;
  double df_fft, f, mask, eff_mask;

  if(N < 2 || delta <= 0.0) {
    for(i = 0; i < N; i++) signal_out[i] = dc_level;
    return 0;
  }

  /* Clamp pad length for reflect (can't mirror more than we have);
     zero-padding can be arbitrary but there's little value beyond ~N. */
  Npad = (padmode == VARTOOLS_FOURIERFILTER_PAD_WRAP)
         ? 0
         : (int)(padfrac * (double) N);
  if(padmode == VARTOOLS_FOURIERFILTER_PAD_REFLECT && Npad > N - 1)
    Npad = N - 1;
  if(Npad < 0) Npad = 0;
  Ntot = N + 2 * Npad;

  if((tmpdata = (double *) malloc(2 * Ntot * sizeof(double))) == NULL)
    vt_error(ERR_MEMALLOC);

  /* Pack: left-pad, signal (mean-subtracted), right-pad.  For reflect
     we mirror about indices 0 and N-1; for zero we use 0.0 which
     corresponds to the signal's mean level after subtraction. */
  for(i = 0; i < Ntot; i++) {
    double val;
    if(i < Npad) {
      if(padmode == VARTOOLS_FOURIERFILTER_PAD_REFLECT) {
        int src = Npad - i;   /* mirror about 0: i=Npad-1 -> src=1 */
        if(src < 0) src = 0;
        if(src >= N) src = N - 1;
        val = signal_in[src] - dc_level;
      } else {
        val = 0.0;
      }
    } else if(i < Npad + N) {
      val = signal_in[i - Npad] - dc_level;
    } else {
      if(padmode == VARTOOLS_FOURIERFILTER_PAD_REFLECT) {
        int src = 2 * N - 2 - (i - Npad);   /* mirror about N-1 */
        if(src < 0) src = 0;
        if(src >= N) src = N - 1;
        val = signal_in[src] - dc_level;
      } else {
        val = 0.0;
      }
    }
    tmpdata[2*i]     = val;
    tmpdata[2*i + 1] = 0.0;
  }

  wt = gsl_fft_complex_wavetable_alloc(Ntot);
  ws = gsl_fft_complex_workspace_alloc(Ntot);
  gsl_fft_complex_forward(tmpdata, 1, Ntot, wt, ws);

  /* Optionally extract pre-mask cos/sin coefficients at bins 0..Ntot/2
     for the ofourier output.  Convention: signal_in[n] ~=
     sum_k [a_k cos + b_k sin] with arguments 2*pi*k*n/Ntot (the FFT's
     implicit origin-at-zero time convention).  The k=0 DC coefficient
     is set to dc_level (the mean we subtracted), not Re(X[0])/Ntot. */
  if(coeffs_a_out != NULL && coeffs_b_out != NULL) {
    int Nhalf = Ntot / 2;
    coeffs_a_out[0] = dc_level;
    coeffs_b_out[0] = 0.0;
    for(k = 1; k <= Nhalf; k++) {
      /* Factor 2/Ntot for conjugate-paired bins, 1/Ntot for the
         self-conjugate Nyquist bin (even Ntot only). */
      double factor = 2.0;
      if(Ntot % 2 == 0 && k == Ntot / 2) factor = 1.0;
      coeffs_a_out[k] =  factor * tmpdata[2*k]     / (double) Ntot;
      coeffs_b_out[k] = -factor * tmpdata[2*k + 1] / (double) Ntot;
    }
  }

  /* Apply mask symmetrically (at |f|) so conjugate symmetry is
     preserved and IFFT is strictly real. */
  df_fft = 1.0 / ((double) Ntot * delta);
  for(k = 0; k < Ntot; k++) {
    f = (k <= Ntot / 2) ? ((double) k * df_fft)
                        : ((double) (Ntot - k) * df_fft);
    mask = _fourierfilter_mask(f, filtertype, minfreq, maxfreq,
                               taper_type, taper_delta, taper_beta);
    if(has_filterexpr) {
      SetVariable_Value_Double(lcid, threadid, 0, freq_var, f);
      mask *= EvaluateExpression(lcid, threadid, 0, filter_expr);
    }
    eff_mask = invert_mask ? (1.0 - mask) : mask;
    tmpdata[2*k]     *= eff_mask;
    tmpdata[2*k + 1] *= eff_mask;
  }

  gsl_fft_complex_inverse(tmpdata, 1, Ntot, wt, ws);
  gsl_fft_complex_wavetable_free(wt);
  gsl_fft_complex_workspace_free(ws);

  /* Extract the central N samples and add dc_level back. */
  for(i = 0; i < N; i++) signal_out[i] = tmpdata[2*(i + Npad)] + dc_level;

  free(tmpdata);
  return Ntot;
}

/* Resample + filter one LC segment.  Interpolates the (possibly non-
   uniformly-sampled) signal onto a uniform grid at spacing `delta`,
   hands it to _ff_fft_apply_mask, then interpolates the result back to
   the original sample times.  Writes the filtered model (pass band
   plus dc_level, or reject band plus dc_level when invert_mask=1) into
   mag_model_out[0..N_orig-1].  Caller decides replace vs subtract.
   Returns the padded FFT length Ntot (Nu + 2*Npad).  The caller uses
   that to derive the segment's FFT bin spacing df_fft = 1/(Ntot*delta)
   and bin-count statistics. */
static int _ff_resample_filter_segment(
    int N_orig, const double *t_orig, const double *mag_orig,
    double delta, double dc_level,
    int filtertype, double minfreq, double maxfreq,
    int taper_type, double taper_delta, double taper_beta,
    int has_filterexpr, _Variable *freq_var, _Expression *filter_expr,
    int padmode, double padfrac,
    int invert_mask,
    int lcid, int threadid,
    double *mag_model_out,
    double *coeffs_a_out, double *coeffs_b_out)
{
  int Nu, Ntot, i;
  double t0, t1;
  double *t_uniform = NULL, *mag_uniform = NULL, *filt_uniform = NULL;

  if(N_orig < 2 || delta <= 0.0) {
    for(i = 0; i < N_orig; i++) mag_model_out[i] = dc_level;
    return 0;
  }
  t0 = t_orig[0];
  t1 = t_orig[N_orig - 1];
  Nu = (int) floor((t1 - t0) / delta) + 1;
  if(Nu < 4) Nu = 4;

  if((t_uniform    = (double *) malloc(Nu * sizeof(double))) == NULL ||
     (mag_uniform  = (double *) malloc(Nu * sizeof(double))) == NULL ||
     (filt_uniform = (double *) malloc(Nu * sizeof(double))) == NULL)
    vt_error(ERR_MEMALLOC);

  for(i = 0; i < Nu; i++) t_uniform[i] = t0 + (double) i * delta;
  _ff_linear_interp(N_orig, t_orig, mag_orig, Nu, t_uniform, mag_uniform);

  Ntot = _ff_fft_apply_mask(
      Nu, mag_uniform, dc_level, delta,
      filtertype, minfreq, maxfreq,
      taper_type, taper_delta, taper_beta,
      has_filterexpr, freq_var, filter_expr,
      padmode, padfrac,
      invert_mask,
      lcid, threadid,
      filt_uniform,
      coeffs_a_out, coeffs_b_out);

  _ff_linear_interp(Nu, t_uniform, filt_uniform, N_orig, t_orig, mag_model_out);

  free(t_uniform); free(mag_uniform); free(filt_uniform);
  return Ntot;
}

/* Resolve the resample-delta value for this LC, honoring the
   fix/var/expr forms and the special "delmin" keyword.  Returns a
   positive value on success or 0.0 on failure (caller should warn). */
static double _ff_get_resample_delta(ProgramData *p, _FourierFilter *c,
                                     int threadid, int lcid,
                                     const double *t, int Njd)
{
  int i;
  double delta = 0.0, dt, dtmin;
  switch(c->resample_source) {
  case VARTOOLS_FOURIERFILTER_RESAMPLE_DELMIN:
    /* "delmin" — use the minimum dt in the LC, skipping duplicates
       (dt == 0) and any separation smaller than 1e-12*T to guard
       against timestamp precision issues. */
    dtmin = 0.0;
    if(Njd >= 2) {
      double guard = 1e-12 * (t[Njd - 1] - t[0]);
      if(guard <= 0.0) guard = 0.0;
      for(i = 1; i < Njd; i++) {
        dt = t[i] - t[i - 1];
        if(dt > guard && (dtmin == 0.0 || dt < dtmin)) dtmin = dt;
      }
    }
    delta = dtmin;
    break;
  case VARTOOLS_SOURCE_FIXED:
    delta = c->resample_delta_fix;
    break;
  case VARTOOLS_SOURCE_PRIORCOLUMN:
    getoutcolumnvalue(c->resample_delta_linkedcolumn, threadid, lcid,
                      VARTOOLS_TYPE_DOUBLE, &delta);
    break;
  case VARTOOLS_SOURCE_EVALEXPRESSION:
    delta = EvaluateExpression(lcid, threadid, 0, c->resample_delta_expr);
    break;
  }
  return delta;
}

/* Resolve the gap-break threshold for this LC.  Returns infinity when
   gapbreak is disabled (no splitting).  Fractional and percentile
   sources scan the current LC's dt values. */
static double _ff_get_gap_threshold(ProgramData *p, _FourierFilter *c,
                                    int threadid, int lcid,
                                    const double *t, int Njd)
{
  int i;
  double thresh, *dts, mean_dt;
  if(!c->gapbreak_enabled) return HUGE_VAL;
  switch(c->gapbreak_source) {
  case VARTOOLS_SOURCE_FIXED:
    return c->gapbreak_fix;
  case VARTOOLS_SOURCE_PRIORCOLUMN:
    thresh = HUGE_VAL;
    getoutcolumnvalue(c->gapbreak_linkedcolumn, threadid, lcid,
                      VARTOOLS_TYPE_DOUBLE, &thresh);
    return thresh;
  case VARTOOLS_SOURCE_EVALEXPRESSION:
    return EvaluateExpression(lcid, threadid, 0, c->gapbreak_expr);
  case VARTOOLS_FOURIERFILTER_GAPBREAK_FRAC_MIN_SEP:
  case VARTOOLS_FOURIERFILTER_GAPBREAK_FRAC_MED_SEP:
  case VARTOOLS_FOURIERFILTER_GAPBREAK_PERCENTILE:
    if(Njd < 2) return HUGE_VAL;
    dts = (double *) malloc((Njd - 1) * sizeof(double));
    if(dts == NULL) vt_error(ERR_MEMALLOC);
    for(i = 1; i < Njd; i++) dts[i - 1] = t[i] - t[i - 1];
    if(c->gapbreak_source == VARTOOLS_FOURIERFILTER_GAPBREAK_FRAC_MIN_SEP) {
      mean_dt = dts[0];
      for(i = 1; i < Njd - 1; i++) if(dts[i] < mean_dt) mean_dt = dts[i];
      thresh = c->gapbreak_frac_min * mean_dt;
    } else if(c->gapbreak_source == VARTOOLS_FOURIERFILTER_GAPBREAK_FRAC_MED_SEP) {
      thresh = c->gapbreak_frac_med * median(Njd - 1, dts);
    } else {
      thresh = percentile(Njd - 1, dts, c->gapbreak_percentile);
    }
    free(dts);
    return thresh;
  }
  return HUGE_VAL;
}

void doFourierFilter(ProgramData *p, _FourierFilter *c, int threadid, int lcid) {
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

  /* Initialize output columns to safe defaults so a short-LC early return
     doesn't leave garbage in the per-LC output row. */
  c->mean_mag[lcid] = 0.0;
  c->rms_in[lcid] = 0.0;
  c->rms_out[lcid] = 0.0;
  c->nfreqcalc[lcid] = 0;
  c->nfreqfilt[lcid] = 0;

  if(Njd < 3)
    return;

  /* Record input RMS before any filtering (simple unweighted RMS of
     non-NaN mag values, matching the convention of -rms). */
  {
    double _sm = 0, _sm2 = 0;
    int _n = 0;
    for(i = 0; i < Njd; i++) {
      if(!isnan(mag[i])) {
        _sm += mag[i]; _sm2 += mag[i]*mag[i]; _n++;
      }
    }
    if(_n > 1) {
      double _mean = _sm / _n;
      double _var = _sm2 / _n - _mean*_mean;
      c->rms_in[lcid] = (_var > 0) ? sqrt(_var) : 0.0;
    }
  }

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
  /* Use the FFT's natural bin grid: Nftot = Njd/2 positive-frequency
     bins up to Nyquist at spacing df = 1/(Njd*delmin).  This is the
     grid the non-resample FFT path runs on when padmode=wrap; the
     filter-type switch below then narrows maxfreq_calc based on the
     user's minfreq/maxfreq so Nfcalc is bounded by the requested band.
     The resample path recomputes its own df_fft inside
     _ff_fft_apply_mask from the caller-supplied delta. */
  Nftot = Njd / 2;
  df = 1.0 / ((double) Njd * delmin);
  if(2.0*Nftot+1 > Njd) {
    Nftot = (Njd-1)/2;
    maxfreq_calc = df*Nftot;
  }

  maxfreq_model = maxfreq_calc;

  /* Determine any upper or lower limits on the frequencies to compute */
  switch(c->filtertype) {
  case VARTOOLS_FOURIERFILTER_FULLSPEC:
    break;
  case VARTOOLS_FOURIERFILTER_LOWPASS:
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
  case VARTOOLS_FOURIERFILTER_HIGHPASS:
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
  case VARTOOLS_FOURIERFILTER_BANDPASS:
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
  case VARTOOLS_FOURIERFILTER_BANDCUT:
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

  /* Resample path (Phase R2).  When the user supplies `resample <delta>`
     we interpolate the LC onto a uniform grid, FFT, apply the filter
     mask, IFFT, and interpolate back — the canonical way to run a
     Fourier filter on non-uniformly-sampled data.  Optional
     `gapbreak <spec>` splits at large gaps and filters each segment
     independently; for segments whose DC bin is masked out
     (highpass/bandpass), all segments share the overall-LC weighted
     mean so inter-segment jumps don't appear in the output. */
  if(c->resample_enabled) {
    double delta, gap_thresh, lcave_overall, sum, sumw, w;
    int seg_start, seg_end, seg_n, total_nfreqcalc, total_nfreqfilt;

    delta = _ff_get_resample_delta(p, c, threadid, lcid, t, Njd);
    if(delta <= 0.0) {
      if(!c->nowarn) fprintf(stderr,
              "Warning: -fourierfilter resample delta is <= 0 (got %g) "
              "for LC %d; skipping filter for this LC.\n", delta, lcid);
      c->mean_mag[lcid] = 0.0;
      c->nfreqcalc[lcid] = 0;
      c->nfreqfilt[lcid] = 0;
      c->rms_out[lcid] = c->rms_in[lcid];
      return;
    }
    gap_thresh = _ff_get_gap_threshold(p, c, threadid, lcid, t, Njd);

    /* Warn once per LC if the largest *within-segment* gap exceeds
       1/minfreq — the filter's passband edge period is shorter than a
       gap that will be bridged by interpolation, so the reconstruction
       loses fidelity in/around that gap.  Gaps larger than gap_thresh
       are excluded from this check because gapbreak splits the LC at
       those and each segment is filtered independently — nothing is
       interpolated across them. */
    if(minfreq > 0.0 &&
       (c->filtertype == VARTOOLS_FOURIERFILTER_HIGHPASS ||
        c->filtertype == VARTOOLS_FOURIERFILTER_BANDPASS ||
        c->filtertype == VARTOOLS_FOURIERFILTER_BANDCUT)) {
      double within_seg_maxgap = 0.0;
      for(i = 1; i < Njd; i++) {
        double dt = t[i] - t[i - 1];
        if(dt > gap_thresh) continue;   /* gapbreak splits here */
        if(dt > within_seg_maxgap) within_seg_maxgap = dt;
      }
      if(within_seg_maxgap > 1.0 / minfreq && !c->nowarn) {
        fprintf(stderr,
                "Warning: -fourierfilter: within-segment sampling gap "
                "(%g) exceeds 1/minfreq (%g); the filter spans a "
                "frequency whose period is shorter than an "
                "interpolated gap — the result in/around that gap is "
                "inherently poorly defined.  Proceeding anyway.\n",
                within_seg_maxgap, 1.0 / minfreq);
      }
    }

    /* Weighted mean of the whole LC, used as the shared DC offset
       across all segments so highpass/bandpass don't show jumps. */
    sum = 0.0; sumw = 0.0;
    for(i = 0; i < Njd; i++) {
      w = 1.0 / (err[i] * err[i]);
      sum  += w * mag[i];
      sumw += w;
    }
    lcave_overall = (sumw > 0.0) ? (sum / sumw) : 0.0;

    if((mag_model = (double *) malloc(Njd * sizeof(double))) == NULL)
      vt_error(ERR_MEMALLOC);

    /* For ofourier output we capture the kernel's pre-mask cos/sin
       coefficients from the first segment.  Multi-segment LCs (when
       gapbreak fires) have a different bin grid per segment so we
       can't combine them into a single coefficient file; in that
       case we emit a warning and skip the file. */
    int    ofourier_seg_idx = -1;
    int    ofourier_Ntot = 0;
    double ofourier_delta = 0.0;
    double *ofourier_a = NULL, *ofourier_b = NULL;

    /* Segment loop.  A gap is detected when the inter-sample dt
       exceeds gap_thresh (HUGE_VAL when gapbreak is off, so the whole
       LC is one segment). */
    seg_start = 0;
    total_nfreqcalc = 0;
    total_nfreqfilt = 0;
    int seg_count = 0;
    while(seg_start < Njd) {
      int Ntot_seg;
      double *coeffs_a_arg = NULL, *coeffs_b_arg = NULL;
      seg_end = seg_start + 1;
      while(seg_end < Njd && (t[seg_end] - t[seg_end - 1]) <= gap_thresh) {
        seg_end++;
      }
      seg_n = seg_end - seg_start;

      /* Allocate ofourier coefficient buffers for the first segment
         only; estimate the upper bound on Ntot_seg from the segment
         length plus padding. */
      if(c->ofourier && seg_count == 0) {
        int Nu_est = (int) floor((t[seg_end - 1] - t[seg_start]) / delta) + 1;
        if(Nu_est < 4) Nu_est = 4;
        int Ntot_max = Nu_est + 2 * (int)(c->padfrac * (double) Nu_est) + 1;
        ofourier_a = (double *) malloc(Ntot_max * sizeof(double));
        ofourier_b = (double *) malloc(Ntot_max * sizeof(double));
        if(ofourier_a == NULL || ofourier_b == NULL) vt_error(ERR_MEMALLOC);
        coeffs_a_arg = ofourier_a;
        coeffs_b_arg = ofourier_b;
      }

      Ntot_seg = _ff_resample_filter_segment(
          seg_n, &t[seg_start], &mag[seg_start],
          delta, lcave_overall,
          c->filtertype, minfreq, maxfreq,
          c->taper_type, c->taper_deltafreq, c->taper_beta,
          c->filter_exprstring != NULL, c->freq_var, c->filter_expr,
          c->padmode, c->padfrac,
          0,                                    /* invert_mask: always replace-mode for resample path */
          lcid, threadid,
          &mag_model[seg_start],
          coeffs_a_arg, coeffs_b_arg);

      if(c->ofourier && seg_count == 0 && Ntot_seg > 0) {
        ofourier_seg_idx = 0;
        ofourier_Ntot    = Ntot_seg;
        ofourier_delta   = delta;
      }
      /* Segment FFT bin count = Ntot/2 positive-frequency bins up to
         Nyquist (excluding DC).  Filter-band bin count = bins with
         f < maxfreq_model in that segment's df_fft = 1/(Ntot*delta). */
      if(Ntot_seg > 0) {
        total_nfreqcalc += Ntot_seg / 2;
        total_nfreqfilt += (int) floor(maxfreq_model *
                                       (double) Ntot_seg * delta);
      }
      seg_start = seg_end;
      seg_count++;
    }

    /* Replace the LC with the pass-band reconstruction. */
    for(i = 0; i < Njd; i++) mag[i] = mag_model[i];
    free(mag_model);

    /* Write the ofourier file if requested.  Single-segment runs use
       the captured first-segment coefficients directly; multi-segment
       runs (gapbreak fired) print a warning and skip the file because
       there's no single coefficient grid that covers all segments. */
    if(c->ofourier) {
      if(seg_count > 1) {
        if(!c->nowarn) {
          fprintf(stderr,
                  "Warning: -fourierfilter ofourier with `gapbreak` "
                  "split the LC into %d segments, each with its own "
                  "FFT bin grid; no single coefficient file can "
                  "represent them all and the file is skipped for "
                  "this LC.\n", seg_count);
        }
      } else if(ofourier_seg_idx == 0 && ofourier_a != NULL && ofourier_Ntot > 0) {
        int Nhalf = ofourier_Ntot / 2;
        double df_fft_seg = 1.0 / ((double) ofourier_Ntot * ofourier_delta);
        char *lcoutname;
        FILE *outfile;
        if((lcoutname = (char *) malloc(MAXLEN)) == NULL) vt_error(ERR_MEMALLOC);
        GetOutputFilename(lcoutname, p->lcnames[lcid], c->ofourier_dir,
                          "fouriercoeffs", c->ofourier_format, lcid);
        if((outfile = fopen(lcoutname,"w")) == NULL) {
          vt_error2(ERR_CANNOTWRITE, lcoutname);
        }
        fprintf(outfile, "#Frequency");
        if(c->filter_exprstring != NULL) {
          fprintf(outfile, " CosCoeff_orig SinCoeff_orig CosCoeff_filter SinCoeff_filter\n");
        } else {
          fprintf(outfile, " CosCoeff SinCoeff\n");
        }
        for(i = 0; i <= Nhalf; i++) {
          double f_i = (double) i * df_fft_seg;
          if(c->filter_exprstring != NULL) {
            double fv;
            SetVariable_Value_Double(lcid, threadid, 0, c->freq_var, f_i);
            fv = EvaluateExpression(lcid, threadid, 0, c->filter_expr);
            fprintf(outfile, "%.17g %.17g %.17g %.17g %.17g\n",
                    f_i, ofourier_a[i], ofourier_b[i],
                    ofourier_a[i] * fv, ofourier_b[i] * fv);
          } else {
            fprintf(outfile, "%.17g %.17g %.17g\n",
                    f_i, ofourier_a[i], ofourier_b[i]);
          }
        }
        fclose(outfile);
        free(lcoutname);
      }
    }
    if(ofourier_a != NULL) free(ofourier_a);
    if(ofourier_b != NULL) free(ofourier_b);

    c->mean_mag[lcid]  = lcave_overall;
    c->nfreqcalc[lcid] = total_nfreqcalc;
    c->nfreqfilt[lcid] = total_nfreqfilt;

    /* Post-filter RMS on the filtered mag array. */
    {
      double _sm = 0, _sm2 = 0;
      int _n = 0;
      for(i = 0; i < Njd; i++) {
        if(!isnan(mag[i])) {
          _sm += mag[i]; _sm2 += mag[i]*mag[i]; _n++;
        }
      }
      if(_n > 1) {
        double _mean = _sm / _n;
        double _var = _sm2 / _n - _mean*_mean;
        c->rms_out[lcid] = (_var > 0) ? sqrt(_var) : 0.0;
      }
    }
    return;
  }

  /* Non-resample path: the data is (or is being treated as) uniformly
     sampled at step delmin.  Hand it to the FFT kernel directly.  For
     genuinely non-uniform data without `resample`, warn and skip
     this LC — the direct-DFT fallback was removed because it was
     mathematically invalid on non-uniform grids (the cos/sin basis is
     not orthogonal, so reconstructions systematically over-count).
     Skipping (rather than aborting the command) keeps batch pipelines
     running when a single bad LC turns up. */
  if(!isuniform && !c->forcefft) {
    if(!c->nowarn) {
      fprintf(stderr,
              "Warning: -fourierfilter on non-uniformly-sampled data "
              "requires `resample <delta>`.  Without resampling there "
              "is no mathematically valid Fourier decomposition on a "
              "non-uniform grid — skipping the filter for this light "
              "curve.  Add `resample delmin` (or another delta) to "
              "interpolate onto a uniform grid first.\n");
    }
    c->mean_mag[lcid]  = 0.0;
    c->nfreqcalc[lcid] = 0;
    c->nfreqfilt[lcid] = 0;
    c->rms_out[lcid]   = c->rms_in[lcid];
    return;
  }

  {
    int invert_mask, Ntot;
    double lcave, sum, sumw, w;
    double *coeffs_a_raw = NULL, *coeffs_b_raw = NULL;

    if(c->forcefft && !isuniform && !c->nowarn) {
      fprintf(stderr,
              "Warning: -fourierfilter forcefft was requested on "
              "non-uniformly-sampled data; the FFT treats samples as "
              "evenly spaced so the coefficients will not be "
              "meaningful.  For non-uniform data prefer the `resample "
              "<delta>` path instead.\n");
    }

    /* Warn on bandpass/bandcut if the taper windows from the two
       edges overlap — one warning per LC. */
    if(c->taper_type != VARTOOLS_FOURIERFILTER_TAPER_NONE &&
       (c->filtertype == VARTOOLS_FOURIERFILTER_BANDPASS ||
        c->filtertype == VARTOOLS_FOURIERFILTER_BANDCUT) &&
       maxfreq > minfreq &&
       c->taper_deltafreq > 0.5 * (maxfreq - minfreq) && !c->nowarn) {
      fprintf(stderr,
              "Warning: -fourierfilter taper deltafreq (%g) exceeds "
              "half the band width (%g); the two edge tapers overlap "
              "and the pass/reject plateau is reduced to a curved peak.\n",
              c->taper_deltafreq, 0.5 * (maxfreq - minfreq));
    }

    /* Weighted mean: the DC level that gets subtracted before the FFT
       and restored after the IFFT. */
    sum = 0.0; sumw = 0.0;
    for(i = 0; i < Njd; i++) {
      w = 1.0 / (err[i] * err[i]);
      sum  += w * mag[i];
      sumw += w;
    }
    lcave = (sumw > 0.0) ? (sum / sumw) : 0.0;

    /* Subtract-mode (highpass and bandcut without filterexpr/fullspec):
       the kernel computes the reject-band model and the caller
       subtracts.  This preserves the original high-frequency content
       exactly except for FP round-off in the FFT round-trip.  In all
       other cases we replace with the pass-band model. */
    invert_mask =
        ((c->filtertype == VARTOOLS_FOURIERFILTER_HIGHPASS ||
          c->filtertype == VARTOOLS_FOURIERFILTER_BANDCUT) &&
         c->filter_exprstring == NULL &&
         !c->calc_full_spec);

    if((mag_model = (double *) malloc(Njd * sizeof(double))) == NULL)
      vt_error(ERR_MEMALLOC);

    /* For ofourier we also want the pre-mask cos/sin coefficients.
       The kernel writes Ntot/2+1 of them when asked; we allocate the
       worst case (padded FFT length). */
    if(c->ofourier) {
      int Ntot_max = Njd + 2 * (int)(c->padfrac * (double) Njd) + 1;
      coeffs_a_raw = (double *) malloc(Ntot_max * sizeof(double));
      coeffs_b_raw = (double *) malloc(Ntot_max * sizeof(double));
      if(coeffs_a_raw == NULL || coeffs_b_raw == NULL) vt_error(ERR_MEMALLOC);
    }

    /* Run FFT -> mask -> IFFT. */
    Ntot = _ff_fft_apply_mask(
        Njd, mag, lcave, delmin,
        c->filtertype, minfreq, maxfreq,
        c->taper_type, c->taper_deltafreq, c->taper_beta,
        c->filter_exprstring != NULL, c->freq_var, c->filter_expr,
        c->padmode, c->padfrac,
        invert_mask,
        lcid, threadid,
        mag_model,
        coeffs_a_raw, coeffs_b_raw);

    /* Apply the model: replace for the pass-band case; subtract
       (keeping the DC term) for the reject-band case. */
    if(invert_mask) {
      for(i = 0; i < Njd; i++) {
        mag[i] = mag[i] - (mag_model[i] - lcave);
      }
    } else {
      for(i = 0; i < Njd; i++) mag[i] = mag_model[i];
    }

    /* ofourier file: write the pre-mask cos/sin coefficients that the
       kernel returned, plus the post-filterexpr coefficients if
       filterexpr is set.  Frequencies are at the padded-FFT grid
       df_fft = 1/(Ntot*delmin), which equals df only when padmode =
       WRAP (the default) and becomes finer with reflect/zero. */
    if(c->ofourier && Ntot > 0) {
      int Nhalf = Ntot / 2;
      double df_fft = 1.0 / ((double) Ntot * delmin);
      char *lcoutname;
      FILE *outfile;
      if((lcoutname = (char *) malloc(MAXLEN)) == NULL)
        vt_error(ERR_MEMALLOC);
      GetOutputFilename(lcoutname, p->lcnames[lcid], c->ofourier_dir,
                        "fouriercoeffs", c->ofourier_format, lcid);
      if((outfile = fopen(lcoutname,"w")) == NULL) {
        vt_error2(ERR_CANNOTWRITE, lcoutname);
      }
      fprintf(outfile, "#Frequency");
      if(c->filter_exprstring != NULL) {
        fprintf(outfile, " CosCoeff_orig SinCoeff_orig CosCoeff_filter SinCoeff_filter\n");
      } else {
        fprintf(outfile, " CosCoeff SinCoeff\n");
      }
      for(i = 0; i <= Nhalf; i++) {
        f = (double) i * df_fft;
        if(c->filter_exprstring != NULL) {
          /* Post-filter value is pre-filter multiplied by the
             user's filterexpr evaluated at this bin's frequency. */
          double fv;
          SetVariable_Value_Double(lcid, threadid, 0, c->freq_var, f);
          fv = EvaluateExpression(lcid, threadid, 0, c->filter_expr);
          fprintf(outfile, "%.17g %.17g %.17g %.17g %.17g\n",
                  f, coeffs_a_raw[i], coeffs_b_raw[i],
                  coeffs_a_raw[i] * fv, coeffs_b_raw[i] * fv);
        } else {
          fprintf(outfile, "%.17g %.17g %.17g\n",
                  f, coeffs_a_raw[i], coeffs_b_raw[i]);
        }
      }
      fclose(outfile);
      free(lcoutname);
    }
    if(coeffs_a_raw != NULL) free(coeffs_a_raw);
    if(coeffs_b_raw != NULL) free(coeffs_b_raw);

    c->mean_mag[lcid]  = lcave;
    /* Nfreqcalc = positive-frequency bins up to Nyquist (Ntot/2).
       Nfreqfilt = bins inside the filter pass band (below maxfreq_model). */
    c->nfreqcalc[lcid] = (Ntot > 0) ? (Ntot / 2) : 0;
    c->nfreqfilt[lcid] = (Ntot > 0) ?
        (int) floor(maxfreq_model * (double) Ntot * delmin) : 0;

    /* Post-filter RMS — recompute on the (now filtered) mag array. */
    {
      double _sm = 0, _sm2 = 0;
      int _n = 0;
      for(i = 0; i < Njd; i++) {
        if(!isnan(mag[i])) {
          _sm += mag[i]; _sm2 += mag[i]*mag[i]; _n++;
        }
      }
      if(_n > 1) {
        double _mean = _sm / _n;
        double _var = _sm2 / _n - _mean*_mean;
        c->rms_out[lcid] = (_var > 0) ? sqrt(_var) : 0.0;
      }
    }

    free(mag_model);
  }
}

/* Recursive walk over a parsed expression tree that rejects any
   reference to a variable of vectortype VARTOOLS_VECTORTYPE_LC.  The
   filter expression is evaluated in frequency space with obs_index=0,
   so LC-vector variables make no physical sense — silently using the
   value at observation 0 would give surprising results. */
static int _fourierfilter_expr_has_lc_var(_Expression *e, const char **badname)
{
  if(e == NULL) return 0;
  if(e->op1type == VARTOOLS_OPERANDTYPE_VARIABLE && e->op1_variable != NULL) {
    if(e->op1_variable->vectortype == VARTOOLS_VECTORTYPE_LC) {
      if(badname) *badname = e->op1_variable->varname;
      return 1;
    }
  }
  if(e->op2type == VARTOOLS_OPERANDTYPE_VARIABLE && e->op2_variable != NULL) {
    if(e->op2_variable->vectortype == VARTOOLS_VECTORTYPE_LC) {
      if(badname) *badname = e->op2_variable->varname;
      return 1;
    }
  }
  if(e->op1type == VARTOOLS_OPERANDTYPE_EXPRESSION
     && _fourierfilter_expr_has_lc_var((_Expression *) e->op1_expression, badname))
    return 1;
  if(e->op2type == VARTOOLS_OPERANDTYPE_EXPRESSION
     && _fourierfilter_expr_has_lc_var((_Expression *) e->op2_expression, badname))
    return 1;
  /* Function-call arguments also need walking; exp(-mag*f), sin(mag), etc.
     should all reject. */
  if(e->op1type == VARTOOLS_OPERANDTYPE_FUNCTION
     && e->op1_functioncall != NULL) {
    _FunctionCall *fc = (_FunctionCall *) e->op1_functioncall;
    int a;
    for(a = 0; a < fc->Nexpr; a++) {
      if(_fourierfilter_expr_has_lc_var(fc->arguments[a], badname))
        return 1;
    }
  }
  if(e->op2type == VARTOOLS_OPERANDTYPE_FUNCTION
     && e->op2_functioncall != NULL) {
    _FunctionCall *fc = (_FunctionCall *) e->op2_functioncall;
    int a;
    for(a = 0; a < fc->Nexpr; a++) {
      if(_fourierfilter_expr_has_lc_var(fc->arguments[a], badname))
        return 1;
    }
  }
  return 0;
}

/* Whole-word in-place substring replacement on src, writing into dst
   (length-dst_sz).  Only identifier characters (alnum or '_') on the
   boundaries count as part of a word.  Returns 0 on success, -1 if
   dst would overflow. */
static int _fourierfilter_subst_whole_word(const char *src,
                                           const char *from,
                                           const char *to,
                                           char *dst, size_t dst_sz)
{
  size_t flen = strlen(from), tlen = strlen(to);
  size_t wr = 0;
  const char *p = src;
  while(*p) {
    if(strncmp(p, from, flen) == 0) {
      int before_ok = (p == src)
                      || !(isalnum((unsigned char) *(p-1)) || *(p-1) == '_');
      int after_ok = (p[flen] == '\0')
                     || !(isalnum((unsigned char) p[flen]) || p[flen] == '_');
      if(before_ok && after_ok) {
        if(wr + tlen >= dst_sz) return -1;
        memcpy(dst + wr, to, tlen);
        wr += tlen;
        p += flen;
        continue;
      }
    }
    if(wr + 1 >= dst_sz) return -1;
    dst[wr++] = *p++;
  }
  dst[wr] = '\0';
  return 0;
}

/* Called from analytic.c::CompileAllExpressions once per
   -fourierfilter command that supplied a "filterexpr" string.  Builds
   a unique stump variable, substitutes the user-visible frequency name
   with the stump in the expression source, parses, and validates that
   only per-star / scalar / constant inputs are referenced. */
void SetupFourierFilterExpression(ProgramData *p, _FourierFilter *c, int cnum)
{
  char stump[64];
  char substituted[4*MAXLEN];
  const char *badname = NULL;

  /* Unique stump name, keyed by command index.  Keep it short, all-
     lowercase, starting with '__ff_' so it's recognizable in error
     messages. */
  snprintf(stump, sizeof(stump), "__ff_freq_%d", cnum);

  /* Allocate the INTERNALSCALAR variable and its per-thread storage. */
  c->freq_var = CreateVariable(p, stump, VARTOOLS_TYPE_DOUBLE,
                               VARTOOLS_VECTORTYPE_INTERNALSCALAR, NULL);
  RegisterDataFromLightCurve(p, c->freq_var->dataptr, VARTOOLS_TYPE_DOUBLE,
                             0, 0, 0, 0, 0, NULL, c->freq_var, -1, stump);

  /* Substitute whole-word occurrences of the user's freq variable name
     with the stump name so the expression parser resolves to the
     INTERNALSCALAR we just created — even if the user's name (default
     "f") collides with another variable in scope. */
  if(_fourierfilter_subst_whole_word(c->filter_exprstring,
                                     c->freq_varname, stump,
                                     substituted, sizeof(substituted))) {
    vt_error2(ERR_CODEERROR,
              "-fourierfilter filterexpr: substituted expression exceeds "
              "internal buffer size");
  }

  c->filter_expr = ParseExpression(substituted, p);

  if(_fourierfilter_expr_has_lc_var(c->filter_expr, &badname)) {
    fprintf(stderr,
            "\nError in -fourierfilter (command #%d): filterexpr may not "
            "reference light-curve-vector variables.\n"
            "'%s' has vectortype VECTORTYPE_LC.  Only per-star scalars, "
            "constants, and the frequency variable (%s) are allowed.\n\n",
            cnum, badname ? badname : "<unknown>", c->freq_varname);
    vt_error(ERR_CODEERROR);
  }
}


int ParseFourierFilterCommand(int *iret, int argc, char **argv, ProgramData *p,
			       _FourierFilter *c, int cnum)
/* Parse the command line for the "-fourierfilter" command */
{
  int i, j, k;
  FILE *infile;
  char *line = NULL;
  size_t line_size = MAXLEN;
  int colnum;
  int size_tout = 0;

  /* Zero-initialize all optional struct fields before any parser branch
     runs.  Without this they would hold arbitrary stack bytes when the
     caller malloc()'s a fresh _FourierFilter, which is what used to
     leave calc_full_spec / forcefft / ofourier / filter_expr /
     freq_var garbage-valued when the user did not explicitly set
     them.  filtertype is assigned by the filter-keyword branches below
     and is always overwritten before return. */
  c->filtertype = VARTOOLS_FOURIERFILTER_FULLSPEC;
  c->maxfreq = NULL;
  c->maxfreq_fix = 0.0;
  c->maxfreq_source = 0;
  c->maxfreq_linkedcolumn = NULL;
  c->maxfreq_expr = NULL;
  c->maxfreq_exprstring = NULL;
  c->minfreq = NULL;
  c->minfreq_fix = 0.0;
  c->minfreq_source = 0;
  c->minfreq_linkedcolumn = NULL;
  c->minfreq_expr = NULL;
  c->minfreq_exprstring = NULL;
  c->filter_expr = NULL;
  c->filter_exprstring = NULL;
  c->calc_full_spec = 0;
  c->forcefft = 0;
  c->taper_type = VARTOOLS_FOURIERFILTER_TAPER_NONE;
  c->taper_deltafreq = 0.0;
  c->taper_beta = 0.0;
  c->resample_enabled = 0;
  c->resample_source = VARTOOLS_SOURCE_FIXED;
  c->resample_delta_fix = 0.0;
  c->resample_delta = NULL;
  c->resample_delta_linkedcolumn = NULL;
  c->resample_delta_expr = NULL;
  c->resample_delta_exprstring = NULL;
  c->gapbreak_enabled = 0;
  c->gapbreak_source = VARTOOLS_SOURCE_FIXED;
  c->gapbreak_fix = 0.0;
  c->gapbreak_frac_min = 0.0;
  c->gapbreak_frac_med = 0.0;
  c->gapbreak_percentile = 0.0;
  c->gapbreak_threshold = NULL;
  c->gapbreak_linkedcolumn = NULL;
  c->gapbreak_expr = NULL;
  c->gapbreak_exprstring = NULL;
  c->ofourier = 0;
  c->ofourier_dir = NULL;
  c->ofourier_formatflag = 0;
  c->ofourier_format = NULL;
  c->nowarn = 0;
  c->padmode = VARTOOLS_FOURIERFILTER_PAD_WRAP;
  c->padfrac = 0.0;
  c->freq_var = NULL;

  /* Default frequency-variable name for filterexpr; overridden by
     "filterexpr <expr> freqvar <name>". */
  snprintf(c->freq_varname, sizeof(c->freq_varname), "f");

  i = *iret;
  if(i >= argc)
    return(1);

  if(!strcmp(argv[i],"full") || !strcmp(argv[i],"fullspec")) {
    /* "full" / "fullspec" filter type — reconstruct the full Fourier
       series with no band limits.  Mostly useful paired with
       "filterexpr" for an analytic filter applied across the whole
       spectrum. */
    c->filtertype = VARTOOLS_FOURIERFILTER_FULLSPEC;
    i++;
  } else if(!strcmp(argv[i],"highpass")) {
    c->filtertype = VARTOOLS_FOURIERFILTER_HIGHPASS;
    i++;
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "minfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->minfreq_source),
				    0, (void *) (&(c->minfreq_fix)),
				    (void *) (&(c->minfreq)),
				    (void *) (&(c->minfreq_linkedcolumn)),
				    (void *) (&(c->minfreq_exprstring)),
				    "FOURIERFILTER_MINFREQ",
				    0, NULL)) {
      *iret = i-1; return 1;}
  } else if(!strcmp(argv[i],"lowpass")) {
    c->filtertype = VARTOOLS_FOURIERFILTER_LOWPASS;
    i++;
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "maxfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->maxfreq_source),
				    0, (void *) (&(c->maxfreq_fix)),
				    (void *) (&(c->maxfreq)),
				    (void *) (&(c->maxfreq_linkedcolumn)),
				    (void *) (&(c->maxfreq_exprstring)),
				    "FOURIERFILTER_MAXFREQ",
				    0, NULL)) {
      *iret = i-1; return 1;}
  } else if(!strcmp(argv[i],"bandpass")) {
    c->filtertype = VARTOOLS_FOURIERFILTER_BANDPASS;
    i++;
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "minfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->minfreq_source),
				    0, (void *) (&(c->minfreq_fix)),
				    (void *) (&(c->minfreq)),
				    (void *) (&(c->minfreq_linkedcolumn)),
				    (void *) (&(c->minfreq_exprstring)),
				    "FOURIERFILTER_MINFREQ",
				    0, NULL)) {
      *iret = i; return 1;}
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "maxfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->maxfreq_source),
				    0, (void *) (&(c->maxfreq_fix)),
				    (void *) (&(c->maxfreq)),
				    (void *) (&(c->maxfreq_linkedcolumn)),
				    (void *) (&(c->maxfreq_exprstring)),
				    "FOURIERFILTER_MAXFREQ",
				    0, NULL)) {
      *iret = i-1; return 1;}
  } else if(!strcmp(argv[i],"bandcut")) {
    c->filtertype = VARTOOLS_FOURIERFILTER_BANDCUT;
    i++;
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "minfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->minfreq_source),
				    0, (void *) (&(c->minfreq_fix)),
				    (void *) (&(c->minfreq)),
				    (void *) (&(c->minfreq_linkedcolumn)),
				    (void *) (&(c->minfreq_exprstring)),
				    "FOURIERFILTER_MINFREQ",
				    0, NULL)) {
      *iret = i-1; return 1;}
    if(ParseParameterBuiltInCommand(p, cnum, &i, argv, argc, "maxfreq", 1,
				    VARTOOLS_TYPE_DOUBLE, &(c->maxfreq_source),
				    0, (void *) (&(c->maxfreq_fix)),
				    (void *) (&(c->maxfreq)),
				    (void *) (&(c->maxfreq_linkedcolumn)),
				    (void *) (&(c->maxfreq_exprstring)),
				    "FOURIERFILTER_MAXFREQ",
				    0, NULL)) {
      *iret = i-1; return 1;}
  } else {*iret = i-1; return 1;}

  /* Optional trailing keywords (any subset, any order):
       "filterexpr" <expr> ["freqvar" <name>]
       "fullspec"
       "forcefft"
       "ofourier" <outdir> ["nameformat" <fmt>]
     Stop as soon as we see an argument that doesn't match any of these
     keywords — the caller-level parser picks up from there. */
  while(i < argc) {
    if(!strcmp(argv[i],"filterexpr")) {
      i++;
      if(i >= argc) { *iret = i-1; return 1; }
      if((c->filter_exprstring = (char *) malloc(strlen(argv[i])+1)) == NULL)
        vt_error(ERR_MEMALLOC);
      strcpy(c->filter_exprstring, argv[i]);
      i++;
      /* Optional "freqvar <name>" right after the expression. */
      if(i < argc && !strcmp(argv[i],"freqvar")) {
        i++;
        if(i >= argc) { *iret = i-1; return 1; }
        snprintf(c->freq_varname, sizeof(c->freq_varname), "%s", argv[i]);
        i++;
      }
    } else if(!strcmp(argv[i],"fullspec")) {
      c->calc_full_spec = 1;
      i++;
    } else if(!strcmp(argv[i],"forcefft")) {
      c->forcefft = 1;
      i++;
    } else if(!strcmp(argv[i],"taper")) {
      /* taper <linear|cosine|tukey|hann|blackman|kaiser>
         deltafreq <value> [beta <value>]
         The function-name keyword comes first so the parser can report
         a useful error if the user mistypes.  deltafreq is required;
         beta is only consumed after "kaiser". */
      i++;
      if(i >= argc) { *iret = i-1; return 1; }
      if(!strcmp(argv[i],"linear")) {
        c->taper_type = VARTOOLS_FOURIERFILTER_TAPER_LINEAR;
      } else if(!strcmp(argv[i],"cosine") || !strcmp(argv[i],"tukey") ||
                !strcmp(argv[i],"hann")) {
        c->taper_type = VARTOOLS_FOURIERFILTER_TAPER_COSINE;
      } else if(!strcmp(argv[i],"blackman")) {
        c->taper_type = VARTOOLS_FOURIERFILTER_TAPER_BLACKMAN;
      } else if(!strcmp(argv[i],"kaiser")) {
        c->taper_type = VARTOOLS_FOURIERFILTER_TAPER_KAISER;
      } else {
        fprintf(stderr,
                "Error parsing -fourierfilter: unknown taper function '%s'.  "
                "Supported: linear, cosine (tukey, hann), blackman, kaiser.\n",
                argv[i]);
        *iret = i; return 1;
      }
      i++;
      if(i >= argc || strcmp(argv[i],"deltafreq")) {
        fprintf(stderr,
                "Error parsing -fourierfilter: taper requires \"deltafreq "
                "<value>\" after the taper function name.\n");
        *iret = i-1; return 1;
      }
      i++;
      if(i >= argc) { *iret = i-1; return 1; }
      c->taper_deltafreq = atof(argv[i]);
      if(c->taper_deltafreq <= 0.0) {
        fprintf(stderr,
                "Error parsing -fourierfilter: taper deltafreq must be "
                "positive; got %g.\n", c->taper_deltafreq);
        *iret = i; return 1;
      }
      i++;
      /* "beta" parameter: required after kaiser, silently ignored
         otherwise (belt-and-suspenders: accept and warn). */
      if(i < argc && !strcmp(argv[i],"beta")) {
        i++;
        if(i >= argc) { *iret = i-1; return 1; }
        c->taper_beta = atof(argv[i]);
        if(c->taper_type != VARTOOLS_FOURIERFILTER_TAPER_KAISER) {
          fprintf(stderr,
                  "Warning: -fourierfilter taper beta parameter is only "
                  "meaningful with kaiser; ignoring for this taper.\n");
        }
        if(c->taper_type == VARTOOLS_FOURIERFILTER_TAPER_KAISER &&
           c->taper_beta <= 0.0) {
          fprintf(stderr,
                  "Error parsing -fourierfilter: kaiser beta must be "
                  "positive; got %g.\n", c->taper_beta);
          *iret = i; return 1;
        }
        i++;
      } else if(c->taper_type == VARTOOLS_FOURIERFILTER_TAPER_KAISER) {
        fprintf(stderr,
                "Error parsing -fourierfilter: kaiser taper requires "
                "\"beta <value>\" after deltafreq.\n");
        *iret = i-1; return 1;
      }
      /* Taper on mode=full is meaningless: there are no cut edges to
         soften.  Warn at parse time so the user hears about it once
         rather than per-LC. */
      if(c->filtertype == VARTOOLS_FOURIERFILTER_FULLSPEC) {
        fprintf(stderr,
                "Warning: -fourierfilter taper is a no-op with mode=full "
                "(no cut edges to soften); ignoring.\n");
      }
    } else if(!strcmp(argv[i],"resample")) {
      /* resample <"delmin" | "fix" val | "var" name | "expr" e>
         The resampled-FFT path interpolates the LC onto a uniform grid
         at the requested step size, filters via FFT, and interpolates
         back to the original sample times.  Required before -fourierfilter
         can run on non-uniformly-sampled data without a warning. */
      c->resample_enabled = 1;
      i++;
      if(i >= argc) { *iret = i-1; return 1; }
      if(!strcmp(argv[i],"delmin")) {
        c->resample_source = VARTOOLS_FOURIERFILTER_RESAMPLE_DELMIN;
        i++;
      } else if(!strcmp(argv[i],"fix")) {
        c->resample_source = VARTOOLS_SOURCE_FIXED;
        i++;
        if(i >= argc) { *iret = i-1; return 1; }
        c->resample_delta_fix = atof(argv[i]);
        if(c->resample_delta_fix <= 0.0) {
          fprintf(stderr,
                  "Error parsing -fourierfilter: resample fix value must "
                  "be positive; got %g.\n", c->resample_delta_fix);
          *iret = i; return 1;
        }
        i++;
      } else if(!strcmp(argv[i],"var")) {
        c->resample_source = VARTOOLS_SOURCE_PRIORCOLUMN;
        i++;
        if(i >= argc) { *iret = i-1; return 1; }
        increaselinkedcols(p, &(c->resample_delta_linkedcolumn), argv[i], cnum);
        i++;
      } else if(!strcmp(argv[i],"expr")) {
        c->resample_source = VARTOOLS_SOURCE_EVALEXPRESSION;
        i++;
        if(i >= argc) { *iret = i-1; return 1; }
        c->resample_delta_exprstring =
            (char *) malloc((strlen(argv[i]) + 1) * sizeof(char));
        sprintf(c->resample_delta_exprstring, "%s", argv[i]);
        i++;
      } else {
        fprintf(stderr,
                "Error parsing -fourierfilter: resample requires one of "
                "\"delmin\", \"fix\", \"var\", or \"expr\"; got '%s'.\n",
                argv[i]);
        *iret = i; return 1;
      }
      /* Optional trailing "gapbreak <spec>" clause. */
      if(i < argc && !strcmp(argv[i],"gapbreak")) {
        c->gapbreak_enabled = 1;
        i++;
        if(i >= argc) { *iret = i-1; return 1; }
        if(!strcmp(argv[i],"fix")) {
          c->gapbreak_source = VARTOOLS_SOURCE_FIXED;
          i++;
          if(i >= argc) { *iret = i-1; return 1; }
          c->gapbreak_fix = atof(argv[i]);
          i++;
        } else if(!strcmp(argv[i],"expr")) {
          c->gapbreak_source = VARTOOLS_SOURCE_EVALEXPRESSION;
          i++;
          if(i >= argc) { *iret = i-1; return 1; }
          c->gapbreak_exprstring =
              (char *) malloc((strlen(argv[i]) + 1) * sizeof(char));
          sprintf(c->gapbreak_exprstring, "%s", argv[i]);
          i++;
        } else if(!strcmp(argv[i],"frac_min_sep")) {
          c->gapbreak_source = VARTOOLS_FOURIERFILTER_GAPBREAK_FRAC_MIN_SEP;
          i++;
          if(i >= argc) { *iret = i-1; return 1; }
          c->gapbreak_frac_min = atof(argv[i]);
          i++;
        } else if(!strcmp(argv[i],"frac_med_sep")) {
          c->gapbreak_source = VARTOOLS_FOURIERFILTER_GAPBREAK_FRAC_MED_SEP;
          i++;
          if(i >= argc) { *iret = i-1; return 1; }
          c->gapbreak_frac_med = atof(argv[i]);
          i++;
        } else if(!strcmp(argv[i],"percentile_sep")) {
          c->gapbreak_source = VARTOOLS_FOURIERFILTER_GAPBREAK_PERCENTILE;
          i++;
          if(i >= argc) { *iret = i-1; return 1; }
          c->gapbreak_percentile = atof(argv[i]);
          i++;
        } else {
          fprintf(stderr,
                  "Error parsing -fourierfilter: gapbreak requires one of "
                  "\"fix\", \"expr\", \"frac_min_sep\", \"frac_med_sep\", "
                  "or \"percentile_sep\"; got '%s'.\n", argv[i]);
          *iret = i; return 1;
        }
      }
    } else if(!strcmp(argv[i],"ofourier")) {
      i++;
      if(i >= argc) { *iret = i-1; return 1; }
      c->ofourier = 1;
      if((c->ofourier_dir = (char *) malloc(strlen(argv[i])+1)) == NULL)
        vt_error(ERR_MEMALLOC);
      strcpy(c->ofourier_dir, argv[i]);
      i++;
      if(i < argc && !strcmp(argv[i],"nameformat")) {
        i++;
        if(i >= argc) { *iret = i-1; return 1; }
        c->ofourier_formatflag = 1;
        if((c->ofourier_format = (char *) malloc(strlen(argv[i])+1)) == NULL)
          vt_error(ERR_MEMALLOC);
        strcpy(c->ofourier_format, argv[i]);
        i++;
      }
    } else if(!strcmp(argv[i],"padmode")) {
      /* padmode <wrap|reflect|zero> [padfrac <value>]
         Edge handling applied before the FFT in both the direct-FFT
         and the resample paths. */
      i++;
      if(i >= argc) { *iret = i-1; return 1; }
      if(!strcmp(argv[i],"wrap")) {
        c->padmode = VARTOOLS_FOURIERFILTER_PAD_WRAP;
      } else if(!strcmp(argv[i],"reflect")) {
        c->padmode = VARTOOLS_FOURIERFILTER_PAD_REFLECT;
      } else if(!strcmp(argv[i],"zero")) {
        c->padmode = VARTOOLS_FOURIERFILTER_PAD_ZERO;
      } else {
        fprintf(stderr,
                "Error parsing -fourierfilter: unknown padmode '%s'; "
                "expected \"wrap\", \"reflect\", or \"zero\".\n", argv[i]);
        *iret = i; return 1;
      }
      i++;
      /* Default padding fraction: 0.5 for reflect/zero, 0 for wrap. */
      c->padfrac = (c->padmode == VARTOOLS_FOURIERFILTER_PAD_WRAP) ? 0.0 : 0.5;
      if(i < argc && !strcmp(argv[i],"padfrac")) {
        i++;
        if(i >= argc) { *iret = i-1; return 1; }
        c->padfrac = atof(argv[i]);
        if(c->padfrac < 0.0) {
          fprintf(stderr,
                  "Error parsing -fourierfilter: padfrac must be >= 0; "
                  "got %g.\n", c->padfrac);
          *iret = i; return 1;
        }
        i++;
      }
    } else if(!strcmp(argv[i],"nowarn")) {
      c->nowarn = 1;
      i++;
    } else {
      break;
    }
  }

  *iret = i-1;
  return 0;
}

#else  /* ! _HAVE_GSL: stubs so vartools still links without GSL. */

int ParseFourierFilterCommand(int *iret, int argc, char **argv, ProgramData *p,
                              _FourierFilter *c, int cnum)
{
  fprintf(stderr,
          "Error: -fourierfilter requires vartools to be built with GSL "
          "(the FFT-based Fourier filter uses gsl_fft_complex_forward / "
          "gsl_fft_complex_inverse).  This vartools build does not have "
          "GSL; the -fourierfilter command is unavailable.\n");
  return 1;
}

void doFourierFilter(ProgramData *p, _FourierFilter *c, int threadid, int lcid)
{
  /* Should not be reachable: ParseFourierFilterCommand fails first.
     Provide a safe stub anyway that sets output columns to zero and
     emits a one-time error. */
  static int warned = 0;
  if(!warned) {
    fprintf(stderr,
            "Error: -fourierfilter is unavailable without GSL — skipping.\n");
    warned = 1;
  }
  if(c->mean_mag)  c->mean_mag[lcid]  = 0.0;
  if(c->rms_in)    c->rms_in[lcid]    = 0.0;
  if(c->rms_out)   c->rms_out[lcid]   = 0.0;
  if(c->nfreqcalc) c->nfreqcalc[lcid] = 0;
  if(c->nfreqfilt) c->nfreqfilt[lcid] = 0;
}

void SetupFourierFilterExpression(ProgramData *p, _FourierFilter *c, int cnum)
{
  /* Parse already failed without GSL; this is never reached.  No-op. */
}

#endif /* _HAVE_GSL */
