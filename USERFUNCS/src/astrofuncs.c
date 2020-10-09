#include "../../src/vartools.h"

/* This is the source code for a user-defined function to be used with
   vartools.

   This library defines a number of functions that can be included in
   vartools analytic expressions.
*/

void astrofuncs_Initialize(ProgramData *p) 
/* This function is used to initialize the library. Every library loadable
   with the name $LIBNAME.so that is loaded with the VARTOOLS -F option must 
   contain a function with the name $LIBNAME_Initialize which is used to
   register the new analytic functions with VARTOOLS.
*/
{
  double astrofuncs_eccentricAnomaly(double *);
  double astrofuncs_meanAnomaly(double *);
  double astrofuncs_meanAnomalyConjunction(double *);
  double astrofuncs_transitquadLD(double *);
  double astrofuncs_transitnonlinLD(double *);
  double astrofuncs_BroadeningProfile(double *);
  double astrofuncs_TransitProjectedX(double *);
  double astrofuncs_TransitProjectedY(double *);
  double astrofuncs_RV(double *);
  double astrofuncs_RVM(double *);
  double astrofuncs_RVdt(double *);
  double astrofuncs_RVdtp(double *);

  /* Use the VARTOOLS_RegisterUserFunction procedure to register each
     function that this library
     provides. VARTOOLS_RegisterUserFunction takes as input a pointer
     to the ProgramData structure p, the name of the function to use
     in the analytic expression evaluator, the number of arguments
     required by the function, a pointer to the function, and then either 0
     if no additional help will be provided for this function, or 1 if you
     will provide some text to describe the purpose of the function and its
     arguments. If you give one, then you will need to provide 1 + 2*Narg
     additional arguments. The first is a string giving a brief statement
     of the purpose of the function. Following this you should provide strings
     giving a variable name for each argument, and a brief description of it.*/

  VARTOOLS_RegisterUserFunction(p, "EccentricAnomaly", 2, 
				&astrofuncs_eccentricAnomaly, 1, 
				"returns the eccentric anomaly in radians", 
				"M", "mean anomaly in radians", 
				"e", "eccentricity");

  VARTOOLS_RegisterUserFunction(p, "MeanAnomaly", 2, 
				&astrofuncs_meanAnomaly, 1, 
				"returns the mean anomaly in radians", 
				"dt", "the time since periastron", 
				"P", "the orbital period");

  VARTOOLS_RegisterUserFunction(p, "MeanAnomalyConjunction", 4, 
				&astrofuncs_meanAnomalyConjunction, 1, 
				"returns the mean anomaly in radians", 
				"dt", "time since conjunction (or transit)", 
				"P", "the orbital period", 
				"e", "eccentricity", 
				"omega", "argument of periastron in degrees");

  VARTOOLS_RegisterUserFunction(p, "TransitQuadLD", 9, 
				&astrofuncs_transitquadLD, 1,
				"returns the relative flux of a source "
				"in transit using the Mandel & Agol 2002 "
				"semi-analytic transit model for quadratic "
				"limb darkening (=1 means out of transit, "
				"<1 means in transit).",
				"dt", "time since transit center", 
				"P", "orbital period", 
				"b", "impact parameter at conjunction "
				     "normalized to the stellar radius", 
				"Rp/R*", "planet to stellar radius ratio", 
				"a/R*", "semi-major axis in units of the "
				        "stellar radius", 
				"e", "eccentricity", 
				"omega", "argument of periastron", 
				"u1", "first limb darkening coefficient "
				      "for a quadratic law", 
				"u2", "second limb darkening coefficient "
				      "for a quadratic law");

  VARTOOLS_RegisterUserFunction(p, "TransitNonlinLD", 11, 
				&astrofuncs_transitnonlinLD, 1,
				"returns the relative flux of a source "
				"in transit using the Mandel & Agol 2002 "
				"semi-analytic transit model for a 4 parameter "
				"non-linear limb darkening law "
				"(=1 means out of transit, "
				"<1 means in transit).",
				"dt", "time since transit center", 
				"P", "orbital period", 
				"b", "impact parameter at conjunction "
				     "normalized to the stellar radius", 
				"Rp/R*", "planet to stellar radius ratio", 
				"a/R*", "semi-major axis in units of the "
				        "stellar radius", 
				"e", "eccentricity", 
				"omega", "argument of periastron", 
				"a1", "first limb darkening coefficient "
				      "for a non-linear law", 
				"a2", "second limb darkening coefficient "
				      "for a non-linear law",
				"a3", "third limb darkening coefficient "
				      "for a non-linear law",
				"a4", "fourth limb darkening coefficient "
				      "for a non-linear law");

  VARTOOLS_RegisterUserFunction(p, "BroadeningProfile", 12, 
				&astrofuncs_BroadeningProfile, 1,
				"returns the (distorted) line broadening "
				"function at a "
				"given wavelength for a star with a "
				"transiting planet",
				"delv", "(wl - wl0)*c_light/wl0, where wl is "
				        "the wavelength to return the "
                                        "broadening function at, wl0 is the "
                                        "central wavelength of the line, and "
                                        "c_light is the speed of light in "
				        "km/s.",
				"dt", "time since transit center",
				"P", "orbital period",
				"lambda", "projected obliquity angle in "
				          "degrees",
				"vsini", "projected equatorial rotation "
				         "velocity of the star, in km/s",
				"b", "impact parameter of the planet at "
                                     "conjunction, normalized to the stellar "
 				     "radius.",
                                "Rp/R*", "planet to stellar radius ratio",
                                "a/R*", "semi-major axis in units of the "
				        "stellar radius",
                                "e", "eccentricity",
                                "omega", "argument of periastron in degrees",
                                "u1", "first limb darkening coefficient "
				      "for a quadratic law", 
				"u2", "second limb darkening coefficient "
				      "for a quadratic law");

  VARTOOLS_RegisterUserFunction(p, "TransitProjectedX", 8, 
				&astrofuncs_TransitProjectedX, 1,
				"returns the sky-projected X position of the "
				"center of a planet in its orbit in front of "
				"the star, in units of the stellar radius. The "
				"coordinate system used has the rotation axis "
				"of the star along the Y direction.",
				"dt", "time since transit center",
				"P", "orbital period",
				"lambda", "sky-projected obliquity angle in "
				          "degrees",
				"b", "impact parameter of the planet at "
                                     "conjunction, normalized to the stellar "
 				     "radius.",
                                "Rp/R*", "planet to stellar radius ratio",
                                "a/R*", "semi-major axis in units of the "
				        "stellar radius",
                                "e", "eccentricity",
                                "omega", "argument of periastron in degrees");

  VARTOOLS_RegisterUserFunction(p, "TransitProjectedY", 8, 
				&astrofuncs_TransitProjectedY, 1,
				"returns the sky-projected Y position of the "
				"center of a planet in its orbit in front of "
				"the star, in units of the stellar radius. The "
				"coordinate system used has the rotation axis "
				"of the star along the Y direction.",
				"dt", "time since transit center",
				"P", "orbital period",
				"lambda", "sky-projected obliquity angle in "
				          "degrees",
				"b", "impact parameter of the planet at "
                                     "conjunction, normalized to the stellar "
 				     "radius.",
                                "Rp/R*", "planet to stellar radius ratio",
                                "a/R*", "semi-major axis in units of the "
				        "stellar radius",
                                "e", "eccentricity",
                                "omega", "argument of periastron in degrees");

  VARTOOLS_RegisterUserFunction(p, "RV_E", 4, 
				&astrofuncs_RV, 1,
				"returns the radial velocity of an object "
				"on a Keplerian orbit given the eccentric "
				"anomaly as the input time unit. The radial "
				"velocity will be in the same units as K.",
				"E", "eccentric anomaly in radians",
				"e", "eccentricity",
				"omega", "argument of periastron, in degrees.",
				"K", "RV semi-amplitude");

  VARTOOLS_RegisterUserFunction(p, "RV_M", 4, 
				&astrofuncs_RVM, 1,
				"returns the radial velocity of an object "
				"on a Keplerian orbit given the mean "
				"anomaly as the input time unit. The radial "
				"velocity will be in the same units as K.",
				"M", "mean anomaly in radians",
				"e", "eccentricity",
				"omega", "argument of periastron, in degrees.",
				"K", "RV semi-amplitude");

  VARTOOLS_RegisterUserFunction(p, "RV_dt", 5, 
				&astrofuncs_RVdt, 1,
				"returns the radial velocity of an object "
				"on a Keplerian orbit given the time from "
				"conjunction of its companion as the input "
				"time unit. E.g., this is the "
				"radial velocity of a transiting planet host "
				"star given the time since transit. The radial "
				"velocity will be in the same units as K.",
				"dt", "time from conjunction",
				"P", "orbital period",
				"e", "eccentricity",
				"omega", "argument of periastron, in degrees.",
				"K", "RV semi-amplitude");

  VARTOOLS_RegisterUserFunction(p, "RV_dtp", 5, 
				&astrofuncs_RVdtp, 1,
				"returns the radial velocity of an object "
				"on a Keplerian orbit given the time from "
				"periastron as the "
				"input time unit. The radial "
				"velocity will be in the same units as K.",
				"dtp", "time from periastron",
				"P", "orbital period",
				"e", "eccentricity",
				"omega", "argument of periastron, in degrees.",
				"K", "RV semi-amplitude");
}

#define BISECTION_SAFETY_MAX_ITERATIONS 1000
#define ECC_ANOMALY_MAX_ERROR 1.0e-4  // Convergence threshold (average error is less than a third of this)
double astrofuncs_eccentricAnomaly(double *param) 
/* This function returns the eccentric anomaly E given a mean anomaly M 
   and eccentricity e by solving Kepler's equation.

   param[0] = mean anomaly in radians
   param[1] = eccentricity

   return value = eccentric anomaly in radians
*/
{
  double M, e;
  int counter = BISECTION_SAFETY_MAX_ITERATIONS ;
  double Emin, Emax, E;
  double f, dfm, dE ;

  M = param[0];
  e = param[1];
  Emin = M - 1.0; Emax = M + 1.0; E = M;
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

double astrofuncs_meanAnomaly(double *param) 
/* This function returns the mean anomaly M given the time since periastron dt
   and the orbital period P.

   param[0] = time since periastron dt = (t-T)
   param[1] = orbital period P

   return value = mean anomaly in radians
*/
{
  double phase;
  phase = param[0]/param[1] - floor(param[0]/param[1]);
  return (2.0*M_PI*phase);
}

double astrofuncs_phaseofconjunction (double e, double omega){
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


double astrofuncs_meanAnomalyConjunction(double *param) 
/* This function returns the mean anomaly M given the time since
   conjunction (or transit) d, the orbital period P, the eccentricity
   e and the argument of periastron omega.

   param[0] = time since conjunction dt = (t-Tc)
   param[1] = orbital period P
   param[2] = eccentricity e
   param[3] = argument of periastron, omega, in degrees.

   return value = mean anomaly in radians
*/
{
  double phase, phaseC, M;
  phase = param[0]/param[1] - floor(param[0]/param[1]);

  phaseC = astrofuncs_phaseofconjunction(param[2],(param[3]*M_PI/180.0));
  M = (2.0*M_PI*(phase + phaseC));
  if(M > 2.0*M_PI) {
    M -= 2.0*M_PI;
  }
  return M;
}

#define SQR_(A) ((A)*(A))

double astrofuncs_getinstantimpactparameter_transit(double phase, double sin_i, double a, double e, double omega, double phi_c, double p)
{
  const double meanAnomaly = 2.0 * M_PI * (phase + phi_c);
  double D2, D, E;
  double params[2];

  if(e > 0) {
    params[0] = meanAnomaly;
    params[1] = e;
    E = astrofuncs_eccentricAnomaly(params);
  }
  else
    E = meanAnomaly;
  D2 = SQR_(1.0 - (e * cos(E))) - SQR_((((cos(E) - e) * sin(omega)) + (sqrt(1.0 - (e*e)) * sin(E) * cos(omega))) * sin_i) ;
  D = sqrt(D2) * a;

  /* Check which side of the orbit the planet is, if the planet is behind the star set D to be something out of transit */
  if(sin(E + omega) < 0.)
    D = p + 2.0;
  return D;
}


double astrofuncs_transitquadLD(double *param)
/* This function returns the relative flux of a source in transit using
   the Mandel & Agol 2002 semi-analytic transit model.  
   
   param[0] = time since transit center dt = (t-T0)
   param[1] = orbital period P
   param[2] = normalized impact parameter (normalized to R*);
   param[3] = Rp/R*;
   param[4] = a/R*;
   param[5] = eccentricity;
   param[6] = omega (argument of periastron, in degrees).
   param[7] = u1 (first limb darkening coefficient)
   param[8] = u2 (second limb darkening coefficient)

   return value = relative flux (=1 means out of transit, <1 in transit)
*/
{
  double dt, P, bimp, rprstar, arstar, eccen, omega, u1, u2;
  double phase_c, omega_rad, phase, sin_i, cos_i, z, flux1, flux2;

  dt = param[0]; P = param[1]; bimp = param[2]; rprstar = param[3];
  arstar = param[4]; eccen = param[5]; omega = param[6]; u1 = param[7];
  u2 = param[8];

  omega_rad = omega*M_PI/180.0;

  phase_c = astrofuncs_phaseofconjunction(eccen, omega_rad);

  phase = (dt/P) - floor(dt/P);

  cos_i = bimp*(1. + eccen)*cos(omega_rad)/(1. - eccen*eccen)/arstar;

  sin_i = sqrt(1. - cos_i*cos_i);

  z = astrofuncs_getinstantimpactparameter_transit(phase, sin_i, arstar, eccen, omega_rad, phase_c, rprstar);

  VARTOOLS_occultquad(&z, u1, u2, rprstar, &flux1, &flux2, 1);
  
  return flux1;
}

double astrofuncs_transitnonlinLD(double *param)
/* This function returns the relative flux of a source in transit using
   the Mandel & Agol 2002 semi-analytic transit model.  
   
   param[0] = time since transit center dt = (t-T0)
   param[1] = orbital period P
   param[2] = normalized impact parameter (normalized to R*);
   param[3] = Rp/R*;
   param[4] = a/R*;
   param[5] = eccentricity;
   param[6] = omega (argument of periastron, in degrees).
   param[7] = a1 (first limb darkening coefficient)
   param[8] = a2 (second limb darkening coefficient)
   param[9] = a3 (third limb darkening coefficient)
   param[10] = a4 (fourth limb darkening coefficient)

   return value = relative flux (=1 means out of transit, <1 in transit)
*/
{
  double dt, P, bimp, rprstar, arstar, eccen, omega, a1, a2, a3, a4;
  double phase_c, omega_rad, phase, sin_i, cos_i, z, flux1;
  double **flux2;

  flux2 = (double **) malloc(2*sizeof(double *));
  flux2[0] = (double *) malloc(5*sizeof(double));
  flux2[1] = (double *) malloc(5*sizeof(double));

  dt = param[0]; P = param[1]; bimp = param[2]; rprstar = param[3];
  arstar = param[4]; eccen = param[5]; omega = param[6]; a1 = param[7];
  a2 = param[8]; a3 = param[9]; a4 = param[10];

  omega_rad = omega*M_PI/180.0;

  phase_c = astrofuncs_phaseofconjunction(eccen, omega_rad);

  phase = (dt/P) - floor(dt/P);

  cos_i = bimp*(1. + eccen)*cos(omega_rad)/(1. - eccen*eccen)/arstar;

  sin_i = sqrt(1. - cos_i*cos_i);

  z = astrofuncs_getinstantimpactparameter_transit(phase, sin_i, arstar, eccen, omega_rad, phase_c, rprstar);

  VARTOOLS_occultnl(rprstar, a1, a2, a3, a4, &z, &flux1, flux2, 1);

  free(flux2[0]);
  free(flux2[1]);
  
  free(flux2);

  return flux1;
}

static double astrofuncs_y2lim(double rprs, double delv, double vp, double yp) {
  if(delv >= 1.0 || delv <= -1.0) return 0.0;
  if(delv-vp >= rprs || delv-vp <= -rprs) return 0.0;
  if(yp + sqrt(rprs*rprs - (delv - vp)*(delv - vp)) >= sqrt(1.0-delv*delv)) return sqrt(1.0-delv*delv);
  if(yp + sqrt(rprs*rprs - (delv - vp)*(delv - vp)) <= -sqrt(1.0-delv*delv)) return -sqrt(1.0-delv*delv);
  return yp + sqrt(rprs*rprs - (delv-vp)*(delv - vp));
}
static double astrofuncs_y1lim(double rprs, double delv, double vp, double yp) {
  if(delv >= 1.0 || delv <= -1.0) return 0.0;
  if(delv-vp >= rprs || delv-vp <= -rprs) return 0.0;
  if(yp - sqrt(rprs*rprs - (delv - vp)*(delv - vp)) >= sqrt(1.0-delv*delv)) return sqrt(1.0-delv*delv);
  if(yp - sqrt(rprs*rprs - (delv - vp)*(delv - vp)) <= -sqrt(1.0-delv*delv)) return -sqrt(1.0-delv*delv);
  return yp - sqrt(rprs*rprs - (delv-vp)*(delv - vp));
}
static double astrofuncs_Gv_broad(double delv, double u1, double u2) {
  if(delv >= 1.0 || delv <= -1.0) return 0.0;
  return 2.0*((1.0-u1-u2)*sqrt(1.0-delv*delv)+(u1+2.0*u2)*0.5*((1.0-delv*delv)*M_PI/2.0)-(2.0/3.0)*u2*(pow((1.0-delv*delv),1.5)));
}
static double astrofuncs_Gv_broad_norm(double u1, double u2) {
  return M_PI*(1.0 - u1/3. - u2/6.);
}
static double astrofuncs_Kv_broad(double delv, double rprs, double vp, double yp, double u1, double u2) {
  double y1_, y2_;
  if(delv >= 1.0 || delv <= -1.0) return 0.0;
  if(delv-vp <= -rprs || delv-vp >= rprs) return 0.0;
  y1_ = astrofuncs_y1lim(rprs,delv,vp,yp);
  y2_ = astrofuncs_y2lim(rprs,delv,vp,yp);
  return ((1.0-u1-u2)*(y2_-y1_)+(u1+2.0*u2)*0.5*((y2_*(1.0-delv*delv-y2_*y2_)+(1.0-delv*delv)*asin((y2_)/sqrt(1.0-delv*delv)))-(y1_*(1.0-delv*delv-y1_*y1_)+(1.0-delv*delv)*asin((y1_)/sqrt(1.0-delv*delv))))-u2*(1.0-delv*delv)*(y2_-y1_)+(u2/3.0)*(pow(y2_,3.0)-pow(y1_,3.0)));
}

int astrofuncs_solve_kepler_ME_arome(double M, double ecc, double *E)
/* return the eccentric anomaly computed from the mean anomaly
 
   @param M   (in)  mean anomaly (rad)
   @param ecc (in)  eccentricity
   @param E   (out) eccentric anomaly (rad)

   This function is taken from the AROME package to ensure consistent
   notation with other RM routines.
*/
{
   double ERRTOL = 1.E-10;
   double diff;
   double res;
   int niter;
   int niter_max = 10;
   
   res = M + ecc*sin(M)/(1.E0-ecc*cos(M));
   niter = 0;
   do 
   {
     diff = ( res - M - ecc*sin(res) ) / ( 1.E0 - ecc*cos(res) );
     res -= diff;
     niter++;
   } while (fabs(diff) > ERRTOL && niter < niter_max);
   
   if (fabs(diff) > ERRTOL)
     fprintf(stderr,"could not solve the Kepler Eq. for M = %.15E and ecc = %.15E"
	     " (error = %.15E) in solve_kepler_ME. You may want to increase niter_max.", 
	     M, ecc, diff);
   
   *E = res;
   return 0;
}


double astrofuncs_BroadeningProfile(double *param)
/* This function computes the Line Broadening Function at delv = (wl - wl0)*c_light/wl0. Some of the orbit calculations are taken from the AROME package by ***.

   param[0] = delv
   param[1] = time since transit center dt = (t-T0)
   param[2] = orbital period P
   param[3] = lambda in degrees (projected obliquity angle)
   param[4] = vsini (km/s)
   param[5] = normalized impact parameter (normalized to R*);
   param[6] = Rp/R*;
   param[7] = a/R*;
   param[8] = eccentricity;
   param[9] = omega (argument of periastron, in degrees).
   param[10] = u1 (first limb darkening coefficient)
   param[11] = u2 (second limb darkening coefficient)


   return value = relative flux (=1 means out of transit, <1 in transit)
*/
{
  double dt, P, bimp, rprstar, arstar, eccen, omega, u1, u2;
  double phase_c, omega_rad, phase, sin_i, cos_i, z, x, y, flux1, flux2;
  double delv, lambda, vsini, mean_anom, eccen_anom, true_anom, r_vec;
  double node, true_lat;
  double normval, broadout;

  delv = param[0]; dt = param[1]; P = param[2]; lambda = param[3];
  vsini = param[4]; bimp = param[5]; rprstar = param[6];
  arstar = param[7]; eccen = param[8]; omega = param[9]; u1 = param[10];
  u2 = param[11];

  omega_rad = omega*M_PI/180.0;

  phase_c = astrofuncs_phaseofconjunction(eccen, omega_rad);

  phase = (dt/P) - floor(dt/P);

  cos_i = bimp*(1. + eccen)*cos(omega_rad)/(1. - eccen*eccen)/arstar;

  sin_i = sqrt(1. - cos_i*cos_i);

  z = astrofuncs_getinstantimpactparameter_transit(phase, sin_i, arstar, eccen, omega_rad, phase_c, rprstar);

  VARTOOLS_occultquad(&z, u1, u2, rprstar, &flux1, &flux2, 1);


  normval = astrofuncs_Gv_broad_norm(u1,u2)*flux1;

  mean_anom = (2.0 * M_PI * (phase + phase_c));
  astrofuncs_solve_kepler_ME_arome(mean_anom, eccen, &eccen_anom);
  true_anom = 2.E0*atan(sqrt((1.E0+eccen)/(1.E0-eccen))*tan(eccen_anom/2.E0));

  r_vec = arstar*(1.E0-eccen*eccen)/(1.E0+eccen*cos(true_anom));
  node = M_PI + lambda*M_PI/180.;
  true_lat = omega_rad + true_anom;
  x = r_vec*(cos(node)*cos(true_lat)-sin(node)*sin(true_lat)*cos_i);
  y = r_vec*(sin(node)*cos(true_lat)+cos(node)*sin(true_lat)*cos_i);

  broadout = (astrofuncs_Gv_broad((delv/vsini), u1, u2) - astrofuncs_Kv_broad((delv/vsini), rprstar, x, y, u1, u2))/normval;

  return broadout;
}

double astrofuncs_TransitProjectedX(double *param)
  /* This function computes the sky-projected X position of the center of a planet in units of the stellar radius and in a coordinate system where the rotation axis of the star is vertical 

   param[0] = time since transit center dt = (t-T0)
   param[1] = orbital period P
   param[2] = lambda in degrees (projected obliquity angle)
   param[3] = normalized impact parameter (normalized to R*);
   param[4] = Rp/R*;
   param[5] = a/R*;
   param[6] = eccentricity;
   param[7] = omega (argument of periastron, in degrees).

   return value = Projected x position of the center of the planet in its orbit in front of the star, in units of the stellar radius.
*/
{
  double dt, P, bimp, rprstar, arstar, eccen, omega;
  double phase_c, omega_rad, phase, sin_i, cos_i, z, x, y;
  double lambda, mean_anom, eccen_anom, true_anom, r_vec;
  double node, true_lat;

  dt = param[0]; P = param[1]; lambda = param[2];
  bimp = param[3]; rprstar = param[4];
  arstar = param[5]; eccen = param[6]; omega = param[7];

  omega_rad = omega*M_PI/180.0;

  phase_c = astrofuncs_phaseofconjunction(eccen, omega_rad);

  phase = (dt/P) - floor(dt/P);

  cos_i = bimp*(1. + eccen)*cos(omega_rad)/(1. - eccen*eccen)/arstar;

  sin_i = sqrt(1. - cos_i*cos_i);

  z = astrofuncs_getinstantimpactparameter_transit(phase, sin_i, arstar, eccen, omega_rad, phase_c, rprstar);

  mean_anom = (2.0 * M_PI * (phase + phase_c));
  astrofuncs_solve_kepler_ME_arome(mean_anom, eccen, &eccen_anom);
  true_anom = 2.E0*atan(sqrt((1.E0+eccen)/(1.E0-eccen))*tan(eccen_anom/2.E0));

  r_vec = arstar*(1.E0-eccen*eccen)/(1.E0+eccen*cos(true_anom));
  node = M_PI + lambda*M_PI/180.;
  true_lat = omega_rad + true_anom;
  x = r_vec*(cos(node)*cos(true_lat)-sin(node)*sin(true_lat)*cos_i);

  return x;
}

double astrofuncs_TransitProjectedY(double *param)
/* This function computes the sky-projected Y position of the center of a planet in units of the stellar radius and in a coordinate system where the rotation axis of the star is vertical 

   param[0] = time since transit center dt = (t-T0)
   param[1] = orbital period P
   param[2] = lambda in degrees (projected obliquity angle)
   param[3] = normalized impact parameter (normalized to R*);
   param[4] = Rp/R*;
   param[5] = a/R*;
   param[6] = eccentricity;
   param[7] = omega (argument of periastron, in degrees).

   return value = Projected Y position of the center of the planet in its orbit in front of the star, in units of the stellar radius.
*/
{
  double dt, P, bimp, rprstar, arstar, eccen, omega;
  double phase_c, omega_rad, phase, sin_i, cos_i, z, x, y;
  double lambda, mean_anom, eccen_anom, true_anom, r_vec;
  double node, true_lat;

  dt = param[0]; P = param[1]; lambda = param[2];
  bimp = param[3]; rprstar = param[4];
  arstar = param[5]; eccen = param[6]; omega = param[7];

  omega_rad = omega*M_PI/180.0;

  phase_c = astrofuncs_phaseofconjunction(eccen, omega_rad);

  phase = (dt/P) - floor(dt/P);

  cos_i = bimp*(1. + eccen)*cos(omega_rad)/(1. - eccen*eccen)/arstar;

  sin_i = sqrt(1. - cos_i*cos_i);

  z = astrofuncs_getinstantimpactparameter_transit(phase, sin_i, arstar, eccen, omega_rad, phase_c, rprstar);

  mean_anom = (2.0 * M_PI * (phase + phase_c));
  astrofuncs_solve_kepler_ME_arome(mean_anom, eccen, &eccen_anom);
  true_anom = 2.E0*atan(sqrt((1.E0+eccen)/(1.E0-eccen))*tan(eccen_anom/2.E0));

  r_vec = arstar*(1.E0-eccen*eccen)/(1.E0+eccen*cos(true_anom));
  node = M_PI + lambda*M_PI/180.;
  true_lat = omega_rad + true_anom;
  x = r_vec*(cos(node)*cos(true_lat)-sin(node)*sin(true_lat)*cos_i);
  y = r_vec*(sin(node)*cos(true_lat)+cos(node)*sin(true_lat)*cos_i);

  return y;
}

double astrofuncs_RV(double *param)
  /* This function computes the radial velocity of an object on a 
     Keplerian orbit given eccentric anomaly as the input time unit

     param[0] = eccentric anomaly in radians
     param[1] = eccentricity
     param[2] = omega (argument of periastron, in degrees);
     param[3] = K - RV semi-amplitude.

     return value = Radial velocity in the same units as K.
*/
{
  double E, e, omega, K, cosE, sinE, nu, RVout;

  E = param[0]; e = param[1]; omega = param[2]; K = param[3];

  cosE = cos(E); sinE = sin(E);
  nu = acos((cosE - e)/(1. - e*cosE));
  if(sinE < 0.)
    nu = 2.0*M_PI - nu;
  
  omega = M_PI*omega/180.0;

  RVout = K*(cos(nu + omega) + e*cos(omega));

  return RVout;
}


double astrofuncs_RVM(double *param)
  /* This function computes the radial velocity of an object on a 
     Keplerian orbit given mean anomaly as the input time unit

     param[0] = mean anomaly in radians
     param[1] = eccentricity
     param[2] = omega (argument of periastron, in degrees);
     param[3] = K - RV semi-amplitude.

     return value = Radial velocity in the same units as K.
*/
{
  double M, E, e, omega, K, RVout;

  double param2[4];

  M = param[0]; e = param[1]; omega = param[2]; K = param[3];

  E = astrofuncs_eccentricAnomaly(param);

  param2[0] = E; param2[1] = e; param2[2] = omega; param2[3] = K;

  RVout = astrofuncs_RV(&(param2[0]));
  return RVout;
}

double astrofuncs_RVdt(double *param)
  /* This function computes the radial velocity of an object on a 
     Keplerian orbit given time from conjunction as the input time unit

     param[0] = dt = (T - T0), where T0 is the time of conjunction
     param[1] = P = orbital period
     param[2] = eccentricity
     param[3] = omega (argument of periastron, in degrees);
     param[4] = K - RV semi-amplitude.

     return value = Radial velocity in the same units as K.
*/
{
  double M, E, e, omega, K, RVout, P, dt;

  double param2[4];

  dt = param[0]; P = param[1]; e = param[2]; omega = param[3]; K = param[4];

  param2[0] = dt; param2[1] = P; param2[2] = e; param2[3] = omega;

  M = astrofuncs_meanAnomalyConjunction(&(param2[0]));

  param2[0] = M; param2[1] = e; param2[2] = omega; param2[3] = K;

  RVout = astrofuncs_RVM(&(param2[0]));
  return RVout;
}

double astrofuncs_RVdtp(double *param)
  /* This function computes the radial velocity of an object on a 
     Keplerian orbit given time from periastron as the input time unit

     param[0] = dtp = (T - TP), where TP is the time of periastron
     param[1] = P = orbital period
     param[2] = eccentricity
     param[3] = omega (argument of periastron, in degrees);
     param[4] = K - RV semi-amplitude.

     return value = Radial velocity in the same units as K.
*/
{
  double M, E, e, omega, K, RVout, P, dt;

  double param2[4];

  dt = param[0]; P = param[1]; e = param[2]; omega = param[3]; K = param[4];

  param2[0] = dt; param2[1] = P;

  M = astrofuncs_meanAnomaly(&(param2[0]));

  param2[0] = M; param2[1] = e; param2[2] = omega; param2[3] = K;

  RVout = astrofuncs_RVM(&(param2[0]));
  return RVout;
}

