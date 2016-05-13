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
#include "commands.h"
#include "programdata.h"
#include "functions.h"

#ifdef _HAVE_CSPICE
#include <SpiceUsr.h>
#endif

/* Speed of light in km/s */
#define C_LIGHT 299792.458
#define ARCSECTORADIANS 4.84813681109536e-6
#define MASPERYRTORADPERDAY 1.3273475e-11
#define KM_PER_AU 1.49597870691e8
#define SEC_PER_DAY 86400

/* Julian Date corresponding to 2000 JAN 01 12:00:00 */
#define JD_2000 2451545.0
#define JD_2000_INT 2451545

double mjd2jd(double deltat, long *tint)
{
  *tint = *tint + 2400000;
  deltat += 0.5;
  if(deltat >= 1.0) {
    *tint = *tint + 1;
    deltat -= 1.0;
  }
  return deltat;
}

double jd2mjd(double deltat, long *tint)
{
  *tint = *tint - 2400000;
  deltat -= 0.5;
  if(deltat < 0.0) {
    *tint = *tint - 1;
    deltat += 1.0;
  }
  return deltat;
}

/* Precess coordinates from epochin to epochout and apply the
   proper motion correction */
/* Precession is based on the precess.pro and premat.pro IDL
   routines by W. Landsman */
void coordprecess(double *rainout, double *decinout,
		    double epochin, double epochout,
		    double jd,
		    double mu_ra, double mu_dec)
/* rainout and decinout in radians
   mu_ra and mu_dec in mas per year */
{
  double ra, dec;
  double cd, sd, ca, sa;
  double M[3][3], t1, t2, a, b, c;
  double x[3], x2[3];
  double jd_epoch_in;
  double sina, sinb, sinc, cosa, cosb, cosc;
  int i, j;

  ra = *rainout; dec = *decinout;

  /* Apply the proper motion */
  if(mu_ra != 0. || mu_dec != 0.) {
    jd_epoch_in = JD_2000 + (epochin - epochout)*365.25;
    ra = ra + (jd - jd_epoch_in)*mu_ra*MASPERYRTORADPERDAY/cos(dec);
    dec = dec + (jd - jd_epoch_in)*mu_dec*MASPERYRTORADPERDAY;
  }

  ca = cos(ra); cd = cos(dec);
  sa = sin(ra); sd = sin(dec);

  if(epochin != epochout) {
    /* Get the precession matrix, assume FK5 */
    t1 = 0.001*(epochout - epochin);
    t2 = 0.001*(epochin - 2000.0);
    a = t1*ARCSECTORADIANS * (23062.181 + t2*(139.656 + 0.0139*t2) + t1*(30.188 - 0.344*t2+17.998*t1));
    b = t1*t1*ARCSECTORADIANS*(79.280 + 0.410*t2 + 0.205*t1) + a;
    c = ARCSECTORADIANS*t1*(20043.109 - t2*(85.33 + 0.217*t2) + t1*(-42.665 - 0.217*t2 - 41.833*t2));
    sina = sin(a); sinb = sin(b); sinc = sin(c);
    cosa = cos(a); cosb = cos(b); cosc = cos(c);
    M[0][0] = cosa*cosb*cosc-sina*sinb;
    M[0][1] = sina*cosb+cosa*sinb*cosc;
    M[0][2] = cosa*sinc;
    M[1][0] = -cosa*sinb-sina*cosb*cosc;
    M[1][1] = cosa*cosb-sina*sinb*cosc;
    M[1][2] = -sina*sinc;
    M[2][0] = -cosb*sinc;
    M[2][1] = -sinb*sinc;
    M[2][2] = cosc;

    /* Apply the precession matrix */
    x[0] = cd*ca; x[1] = cd*sa; x[2] = sd;
    for(i=0; i < 3; i++) {
      x2[i] = M[0][i]*x[0] + M[1][i]*x[1] + M[2][i]*x[2];
    }
    ra = atan2(x2[1],x2[0]);
    dec = asin(x2[2]);
  }
  *rainout = ra;
  *decinout = dec;
}

#define JD_1800 2378495.0
#define JD_2050 2469807.5

void get_EM_orbital_elements(double jd,
		     double *a, double *adot,
		     double *e, double *edot,
		     double *i, double *idot,
		     double *L, double *Ldot,
		     double *om, double *omdot,
		     double *Om, double *Omdot)
/* Returns the approximate Kepler Elements and rates for the Earth-Moon
   Barycenter from http://ssd.jpl.nasa.gov/txt/p_elem_t2.txt.
*/
{
  if(jd > JD_1800 && jd < JD_2050) {
    *a = 1.00000261; *adot = 0.00000562;
    *e = 0.01671123; *edot = -0.00004392;
    *i = -0.00001531; *idot = -0.01294668;
    *L = 100.46457166; *Ldot = 35999.37244981;
    *om = 102.93768193; *omdot = 0.32327364;
    *Om = 0.; *Omdot = 0.;
  }
  else {
    *a = 1.00000018; *adot = -0.00000003;
    *e = 0.01673163; *edot = -0.00003661;
    *i = -0.00054346; *idot = -0.01337178;
    *L = 100.46691572; *Ldot = 35999.37306329;
    *om = 102.93005885; *omdot = 0.31795260;
    *Om = -5.11260389; *Omdot = -0.24123856;
  }
  return;
}

#ifdef _HAVE_CSPICE
double jdtdb2jdutc(double jdtdbdiff, long tint) {
  SpiceChar jdstr[MAXLEN];
  SpiceDouble et;
  double jdutc, jdutcdiff, delta;
  long jdutc_int1, jdutc_int2;
  et = jdtdbdiff*SEC_PER_DAY + (double) ((tint - JD_2000_INT)*SEC_PER_DAY);
  deltet_c(et, "ET", &delta);
  jdutcdiff = jdtdbdiff - delta/SEC_PER_DAY;
  return jdutcdiff;
}
double jdutc2jdtdb(double jdutcdiff, long tint) {
  SpiceChar jdstr[MAXLEN];
  SpiceDouble et;
  double jdtdb, jdtdbdiff, jdutc_sec, delta;
  jdutc_sec = jdutcdiff*SEC_PER_DAY + (double)((tint - JD_2000_INT)*SEC_PER_DAY);
  deltet_c(jdutc_sec, "UTC", &delta);
  jdtdbdiff = jdutcdiff + delta/SEC_PER_DAY;
  return jdtdbdiff;
}
#endif

/* Gives an approximate correction to the helio-center accurate to
   ####TBD#### seconds. Assume the observations were made at the
   Earth-Moon Barycenter (this can lead to an error of ~30 milliseconds
   which is much less than the 8 second error made by not correcting
   to the Solar-System Barycenter) */
double jd2hjd(double jddiff, long tint, double ra, double dec, double epoch,
	      double mu_ra, double mu_dec, int timesys)
{
  double aEM, adotEM, eEM, edotEM,
    iEM, idotEM, LEM, LdotEM,
    omEM, omdotEM, OmEM, OmdotEM;

  double omperi, M, E;

  double cd, sd, ca, sa, Tcent, xprime[3], xecl[3], xeq[3];
  double co, so, cO, sO, ci, si;
  double delta_Romer, hjddiffout;

  double jd, ephemtime;

  jd = (double) (jddiff + tint);

  ephemtime = jd;

#ifdef _HAVE_CSPICE
  if(timesys == TIMESYSTEM_UTC) {
    ephemtime = (double) ( jdutc2jdtdb(jddiff, tint) + tint);
  }
#endif


  get_EM_orbital_elements(ephemtime, &aEM, &adotEM, &eEM, &edotEM,
			  &iEM, &idotEM, &LEM, &LdotEM,
			  &omEM, &omdotEM, &OmEM, &OmdotEM);

  ra = ra*M_PI/180.;
  dec = dec*M_PI/180.;

  /* Precess the coordinates and apply the proper motion
     if necessary */
  if(epoch != 2000.0 || mu_ra != 0. || mu_dec != 0.) {
    coordprecess(&ra, &dec, epoch, 2000.0, ephemtime, mu_ra, mu_dec);
  }
  cd = cos(dec); sd = sin(dec); ca = cos(ra); sa = sin(ra);

  /* Get the orbital elements for the Earth-Moon barycenter
     at the time of observation */
  Tcent = (ephemtime - JD_2000)/36525.;
  aEM += Tcent*adotEM; eEM += Tcent*edotEM; iEM += Tcent*idotEM;
  LEM += Tcent*LdotEM; omEM += Tcent*omdotEM; OmEM += Tcent*OmdotEM;

  /* Get the eccentric anomaly */
  omperi = omEM - OmEM;
  M = LEM - omEM;
  M = M - 360.*floor(M/360.);
  if(M > 180.)
    M -= 360.;
  M = M_PI*M/180.;
  omperi = omperi * M_PI/180.;
  OmEM = OmEM * M_PI/180.;
  iEM = iEM * M_PI/180.;
  E = eccentricAnomaly(M, eEM);

  /* Get the heliocentric coordinates of the Earth-Moon Barycenter in the
     EM orbital plane */
  xprime[0] = aEM*(cos(E)-eEM);
  xprime[1] = aEM*sqrt(1. - eEM*eEM)*sin(E);
  xprime[2] = 0.;

  /* Get the coordinates in the J2000 ecliptic plane */
  co = cos(omperi); so = sin(omperi);
  cO = cos(OmEM); sO = sin(OmEM);
  ci = cos(iEM); si = sin(iEM);
  xecl[0] = (co*cO - so*sO*ci)*xprime[0]
    + (-so*cO-co*sO*ci)*xprime[1];
  xecl[1] = (co*sO+so*cO*ci)*xprime[0]
    + (-so*sO+co*cO*ci)*xprime[1];
  xecl[2] = so*si*xprime[0] + co*si*xprime[1];

  /* Coordinates in the J2000 frame */
  xeq[0] = xecl[0];
  xeq[1] = 0.917482139208287*xecl[1] - 0.397776978008764*xecl[2];
  xeq[2] = 0.397776978008764*xecl[1] + 0.917482139208287*xecl[2];

  /* Calculate the Romer Delay (eq. 2 of Eastman, Siverd and Gaudi, 2010) */
  delta_Romer = (xeq[0]*cd*ca + xeq[1]*cd*sa + xeq[2]*sd)*KM_PER_AU/C_LIGHT/SEC_PER_DAY;

  hjddiffout = jddiff + delta_Romer;

  return(hjddiffout);
}

/* The inverse of jd2hjd */
double hjd2jd(double hjddiff, long tint, double ra, double dec, double epoch,
	      double mu_ra, double mu_dec, int timesys)
{
  double aEM, adotEM, eEM, edotEM,
    iEM, idotEM, LEM, LdotEM,
    omEM, omdotEM, OmEM, OmdotEM;
  double omperi, M, E;

  double cd, sd, ca, sa, Tcent, xprime[3], xecl[3], xeq[3];
  double co, so, cO, sO, ci, si;
  double delta_Romer, jddiffout, hjd, jdoutdiffold;
  double ephemtime, raorig, decorig;
  int numiter;

  hjd = (double) (hjddiff + tint);

  numiter = 0;
  ephemtime = hjd;

#ifdef _HAVE_CSPICE
  if(timesys == TIMESYSTEM_UTC) {
    ephemtime = (double) ( jdutc2jdtdb(hjddiff, tint) + tint);
  }
#endif

  jdoutdiffold = 0.;

  ra = ra*M_PI/180.;
  dec = dec*M_PI/180.;
  raorig = ra;
  decorig = dec;
  do {
    if(numiter > 0) {
      jdoutdiffold = jddiffout;
    }
    numiter++;

    get_EM_orbital_elements(ephemtime, &aEM, &adotEM, &eEM, &edotEM,
			    &iEM, &idotEM, &LEM, &LdotEM,
			    &omEM, &omdotEM, &OmEM, &OmdotEM);


    /* Precess the coordinates and apply the proper motion
       if necessary */
    ra = raorig;
    dec = decorig;
    if(epoch != 2000.0 || mu_ra != 0. || mu_dec != 0.) {
      coordprecess(&ra, &dec, epoch, 2000.0, ephemtime, mu_ra, mu_dec);
    }
    cd = cos(dec); sd = sin(dec); ca = cos(ra); sa = sin(ra);

    /* Get the orbital elements for the Earth-Moon barycenter
       at the time of observation */
    Tcent = (ephemtime - JD_2000)/36525.;
    aEM += Tcent*adotEM; eEM += Tcent*edotEM; iEM += Tcent*idotEM;
    LEM += Tcent*LdotEM; omEM += Tcent*omdotEM; OmEM += Tcent*OmdotEM;

    /* Get the eccentric anomaly */
    omperi = omEM - OmEM;
    M = LEM - omEM;
    M = M - 360.*floor(M/360.);
    if(M > 180.)
      M -= 360.;
    M = M_PI*M/180.;
    omperi = omperi * M_PI/180.;
    OmEM = OmEM * M_PI/180.;
    iEM = iEM * M_PI/180.;
    E = eccentricAnomaly(M, eEM);

    /* Get the heliocentric coordinates of the Earth-Moon Barycenter in the
       EM orbital plane */
    xprime[0] = aEM*(cos(E)-eEM);
    xprime[1] = aEM*sqrt(1. - eEM*eEM)*sin(E);
    xprime[2] = 0.;

    /* Get the coordinates in the J2000 ecliptic plane */
    co = cos(omperi); so = sin(omperi);
    cO = cos(OmEM); sO = sin(OmEM);
    ci = cos(iEM); si = sin(iEM);
    xecl[0] = (co*cO - so*sO*ci)*xprime[0]
      + (-so*cO-co*sO*ci)*xprime[1];
    xecl[1] = (co*sO+so*cO*ci)*xprime[0]
      + (-so*sO+co*cO*ci)*xprime[1];
    xecl[2] = so*si*xprime[0] + co*si*xprime[1];

    /* Coordinates in the J2000 frame */
    xeq[0] = xecl[0];
    xeq[1] = 0.917482139208287*xecl[1] - 0.397776978008764*xecl[2];
    xeq[2] = 0.397776978008764*xecl[1] + 0.917482139208287*xecl[2];

    /* Calculate the Romer Delay (eq. 2 of Eastman, Siverd and Gaudi, 2010) */
    delta_Romer = (xeq[0]*cd*ca + xeq[1]*cd*sa + xeq[2]*sd)*KM_PER_AU/C_LIGHT/SEC_PER_DAY;

    jddiffout = hjddiff - delta_Romer;

    ephemtime = (double) (jddiffout + tint);
#ifdef _HAVE_CSPICE
    if(timesys == TIMESYSTEM_UTC) {
      ephemtime = (double) ( jdutc2jdtdb(jddiffout, tint) + tint);
    }
#endif
  } while(numiter == 1 || fabs(86400*1000*(jddiffout - jdoutdiffold)) > 0.1);

  return(jddiffout);
}

#ifdef _HAVE_CSPICE

/* Use the JPL NAIF cspice library to convert from Julian Date
   measured at an observatory on the earth to Barycentric Julian Date */
double jd2bjd(double jddiff, long tint, double ra, double dec, double epoch,
	      double obslat, double obslong, double obsalt,
	      double mu_ra, double mu_dec, int timesys) {
  SpiceDouble et;
  SpiceDouble lt;
  SpiceDouble state[6];
  SpiceDouble abc[3];
  SpiceDouble jpos[3], epos[3], rotate[3][3];
  SpiceInt n;
  SpiceDouble equatr, polar, f;

  double cd, sd, ca, sa, delta_Romer, bjdoutdiff;
  double jd, jddiffephem;

  if(timesys == TIMESYSTEM_UTC) {
    jddiffephem = jdutc2jdtdb(jddiff, tint);
  }
  else {
    jddiffephem = jddiff;
  }



  jd = jddiffephem + tint;

  ra = ra*M_PI/180.;
  dec = dec*M_PI/180.;

  /* Precess the coordinates and apply the proper motion
     if necessary */
  if(epoch != 2000.0 || mu_ra != 0. || mu_dec != 0.) {
    coordprecess(&ra, &dec, epoch, 2000.0, jd, mu_ra, mu_dec);
  }
  cd = cos(dec); sd = sin(dec); ca = cos(ra); sa = sin(ra);


  /* Seconds since JD 2000.0 */
  et = jddiffephem*SEC_PER_DAY + ((double) ((tint - JD_2000_INT)*SEC_PER_DAY));

  /* Get the J2000 frame x, y, z coordinates of the Earth relative
     to the Solar-System Barycenter */
  spkezr_c( "Earth", et, "J2000", "None", "SSB", state, &lt);

  /* Get the J2000 frame x, y, z coordinates of the observatory
     relative to the center of the Earth */
  if(obsalt > -9998.) {
    bodvcd_c ( 399, "RADII", 3, &n, abc); /* Radii of the earth */
    equatr = abc[0]; polar = abc[2]; f = (equatr - polar)/equatr;
    georec_c(obslong, obslat, (obsalt/1000.), equatr, f, epos);
    pxform_c("IAU_EARTH", "J2000", et, rotate);
    mxv_c(rotate, epos, jpos);

    /* Get the vector from the SSB to the position of the observer */
    state[0] += jpos[0]; state[1] += jpos[1]; state[2] += jpos[2];
  }

  /* Calculate the Romer Delay (eq. 2 of Eastman, Siverd and Gaudi, 2010) */
  delta_Romer = (state[0]*cd*ca + state[1]*cd*sa + state[2]*sd)/C_LIGHT/SEC_PER_DAY;

  bjdoutdiff = jddiff + delta_Romer;

  /* We do not apply the Shapiro or Einstein corrections */

  return(bjdoutdiff);

}

/* JD at an observatory given a BJD */
double bjd2jd(double bjddiff, long tint, double ra, double dec, double epoch,
	      double obslat, double obslong, double obsalt,
	      double mu_ra, double mu_dec, int timesys) {
  SpiceDouble et;
  SpiceDouble lt;
  SpiceDouble state[6];
  SpiceDouble abc[3];
  SpiceDouble jpos[3], epos[3], rotate[3][3];
  SpiceInt n;
  SpiceDouble equatr, polar, f;

  double cd, sd, ca, sa, delta_Romer, jdoutdiff, bjd, checkdelt, jddiffephem;
  double jdoutdiffold, raorig, decorig, jdephem;
  int numiter;

  bjd = bjddiff + tint;

  raorig = ra*M_PI/180.;
  decorig = dec*M_PI/180.;


  /* The input time for the JPL ephemeris is JD, not BJD, so iterate the
     correction until the change in the output JD is less than 0.1 ms. */
  numiter = 0;
  if(timesys == TIMESYSTEM_UTC) {
    jddiffephem = jdutc2jdtdb(bjddiff, tint);
  }
  else {
    jddiffephem = bjddiff;
  }
  jdoutdiffold = 0.;
  do {
    if(numiter > 0) {
      jdoutdiffold = jdoutdiff;
    }
    numiter++;


    /* Precess the coordinates and apply the proper motion
       if necessary */
    jdephem = jddiffephem + tint;
    ra = raorig; dec = decorig;
    if(epoch != 2000.0 || mu_ra != 0. || mu_dec != 0.) {
      coordprecess(&ra, &dec, epoch, 2000.0, jdephem, mu_ra, mu_dec);
    }
    cd = cos(dec); sd = sin(dec); ca = cos(ra); sa = sin(ra);

    /* Seconds since JD 2000.0 */
    et = jddiffephem*SEC_PER_DAY + ((double) ((tint - JD_2000_INT)*SEC_PER_DAY));

    /* Get the J2000 frame x, y, z coordinates of the Solar-System
       Barycenter relative to the Earth */
    spkezr_c( "SSB", et, "J2000", "None", "Earth", state, &lt);

    /* Get the J2000 frame x, y, z coordinates of the observatory
       relative to the center of the Earth */
    if(obsalt > -9998.) {
      bodvcd_c ( 399, "RADII", 3, &n, abc); /* Radii of the earth */
      equatr = abc[0]; polar = abc[2]; f = (equatr - polar)/equatr;
      georec_c(obslong, obslat, (obsalt/1000.), equatr, f, epos);
      pxform_c("IAU_EARTH", "J2000", et, rotate);
      mxv_c(rotate, epos, jpos);

      /* Get the vector from the observer to the SSB */
      state[0] -= jpos[0]; state[1] -= jpos[1]; state[2] -= jpos[2];
    }

    /* Calculate the Romer Delay (eq. 2 of Eastman, Siverd and Gaudi, 2010) */
    delta_Romer = (state[0]*cd*ca + state[1]*cd*sa + state[2]*sd)/C_LIGHT/SEC_PER_DAY;

    jdoutdiff = bjddiff + delta_Romer;

    if(timesys == TIMESYSTEM_UTC) {
      jddiffephem = jdutc2jdtdb(jdoutdiff, tint);
    }
    else {
      jddiffephem = jdoutdiff;
    }
  } while(numiter == 1 || fabs(86400*1000*(jdoutdiff - jdoutdiffold)) > 0.1);

  /* We do not apply the Shapiro or Einstein corrections */

  return(jdoutdiff);

}


/* Load the ephemeris, leap-second, and planetary physical data kernels */
void load_cspice_kernels(char *ephemfile, char *leapsecfile, char *planetdatafile, int inputtimetype, int outputtimetype, int inputsys, int outputsys) {
  char *ephemfile_env, *leapsecfile_env, *planetdatafile_env;
  if(inputtimetype == TIMETYPE_BJD || outputtimetype == TIMETYPE_BJD) {
    if(ephemfile[0] == '\0') {
      ephemfile_env = getenv("CSPICE_EPHEM_FILE");
      if(ephemfile_env == NULL)
	error(ERR_NO_EPHEM_FILE);
      else {
	sprintf(ephemfile,"%s",ephemfile_env);
      }
    }
    furnsh_c(ephemfile);
    if(planetdatafile[0] == '\0') {
      planetdatafile_env = getenv("CSPICE_PLANETDATA_FILE");
      if(planetdatafile_env == NULL)
	error(ERR_NO_PLANETDATA_FILE);
      else {
	sprintf(planetdatafile,"%s",planetdatafile_env);
      }
    }
    furnsh_c(planetdatafile);
  }
  if(inputsys != outputsys || (inputsys == TIMESYSTEM_UTC &&
			       (inputtimetype == TIMETYPE_BJD ||
				inputtimetype == TIMETYPE_HJD ||
				outputtimetype == TIMETYPE_BJD ||
				outputtimetype == TIMETYPE_HJD))) {
    if(leapsecfile[0] == '\0') {
      leapsecfile_env = getenv("CSPICE_LEAPSEC_FILE");
      if(leapsecfile_env == NULL)
	error(ERR_NO_LEAPSEC_FILE);
      else {
	sprintf(leapsecfile,"%s",leapsecfile_env);
      }
    }
    furnsh_c(leapsecfile);
  }
}

#endif

typedef struct {
  char code[256];
  char fullname[MAXLEN];
  double obslat;
  double obslong;
  double obsalt;
} _ConvertTime_Observatory;

#define NUMOBSERVATORIES 21

/* Get the Geocentric coordinates for a given observatory obslat and
   obslong are in degrees (longitude is measured to the East from 0 to
   360) and obsalt is in meters.
 */
void ParseObservatoryCode(char *code, double *obslat, double *obslong, double *obsalt) {

  int i;

  _ConvertTime_Observatory obslist[NUMOBSERVATORIES] = {
    {"ctio", "Cerro Tololo Inter-American Observatory, Chile", -30.165, 289.185, 2215.},
    {"dao", "Dominion Astronomical Observatory, British Columbia", 48.52, 236.5833, 74.},
    {"flwo", "Fred Lawrence Whipple Observatory, Arizona", 31.6883, 249.11505, 2608.},
    {"haleakala", "Haleakala Observatories, Hawaii", 20.7075, 203.7442, 3055},
    {"hess", "High Energy Stereoscopic System Site, Namibia", -23.2716667, 16.50, 1800},
    {"kpno", "Kitt Peak National Observatory, Arizona", 31.9533, 248.38335, 1925.},
    {"lapalma", "Roque de los Muchachose, La Palma", 28.75833, 342.12, 2326.},
    {"lasilla", "European Southern Observatory, La Silla", -29.257, 289.2705, 2347.},
    {"lco", "Las Campanas Observatory, Chile", -29.00833, 289.30005, 2282.},
    {"lick", "Lick Observatory, California", 37.3433, 238.36335, 1290.},
    {"lowell", "Lowell Observatory, Arizona", 35.096667, 248.465, 2198.},
    {"maunakea", "Mauna Kea Observatory, Hawaii", 19.8267, 204.5283, 4215.},
    {"mcdonald", "McDonald Observatory, Texas", 30.6717, 255.9783, 2075.},
    {"mgio", "Mt. Graham International Observatory, Arizona", 32.7013028, 250.10799, 3191},
    {"ohp", "Observatoire de Haute-Provence, France", 43.9308333, 5.71333, 650},
    {"palomar", "Palomar Observatory, California", 33.35667, 243.13665, 1706.},
    {"piszkes", "Piszkesteto Mountain Station, Hungary", 47.918, 19.895, 940.},
    {"saao", "South African Astronomical Observatory, South Africa", 20.8107, -31.6205556, 1798.},
    {"sso", "Siding Springs Observatory, Australia", -31.277039, 149.066085, 1149.},
    {"vlt", "Very Large Telescope, Cerro Paranal, Chile", -24.625, 289.5966, 2635.},
    {"wise", "Wise Observatory, Israel", 30.5958333, 34.7633333, 875}};

  if(!strcmp(code,"show-codes")) {
    printf("# List of available Observatories\n");
    printf("#Code   ----    Name\n");
    printf("#\t\tLatitude[deg]   Longitude[deg_east]  Altitude[m]\n");
    for(i=0; i < NUMOBSERVATORIES; i++) {
      printf("%12s - %s\n", obslist[i].code, obslist[i].fullname);
      printf("\t\t%8.4f\t%8.4f\t%4f\n\n", obslist[i].obslat, obslist[i].obslong, obslist[i].obsalt);
    }
    exit(1);
  } else {
    for(i=0; i < NUMOBSERVATORIES; i++) {
      if(!strcmp(code,obslist[i].code)) {
	*obslat = obslist[i].obslat;
	*obslong = obslist[i].obslong;
	*obsalt = obslist[i].obsalt;
	return;
      }
    }
  }
  error2(ERR_UNKNOWNOBSERVATORY,code);
}

/* Main Function for converting time from one system to another */
void converttime(int N, double *t, int lc, int lcreal, _ConvertTime *c)
{
  int i;
  double rain, raout, decin, decout, ppm_mu_ra_in, ppm_mu_ra_out,
    ppm_mu_dec_in, ppm_mu_dec_out, epochin, epochout, obslat, obslong,
    obsalt, tdouble_add, delt;
  long tint;
  long tint_in, tint_out, tint_sub1, tint_sub2;

  /* Fill in the ra/dec, PPM, and observatory coordinates */
  if((c->inputtimetype == TIMETYPE_HJD ||
      c->inputtimetype == TIMETYPE_BJD ||
      c->outputtimetype == TIMETYPE_HJD ||
      c->outputtimetype == TIMETYPE_BJD))
    {
      if(c->radec_source == VARTOOLS_SOURCE_FIXED) {
	raout = c->raval_fix;
	decout = c->decval_fix;
      }
      else if(c->radec_source == VARTOOLS_SOURCE_INLIST) {
	raout = c->ravals[lcreal][0];
	decout = c->decvals[lcreal][0];
      }
      else
	error(ERR_CODEERROR);
      epochout = c->radecepoch;
      if(c->useinput_radec) {
	if(c->inputradec_source == VARTOOLS_SOURCE_FIXED) {
	  rain = c->inputraval_fix;
	  decin = c->inputdecval_fix;
	}
	else if(c->inputradec_source == VARTOOLS_SOURCE_INLIST) {
	  rain = c->inputravals[lcreal][0];
	  decin = c->inputdecvals[lcreal][0];
	}
	else
	  error(ERR_CODEERROR);
	epochin = c->inputradecepoch;
      }
      else {
	rain = raout; decin = decout; epochin = epochout;
      }
      if(c->useppm) {
	if(c->ppm_source == VARTOOLS_SOURCE_FIXED) {
	  ppm_mu_ra_out = c->ppm_mu_ra_fix;
	  ppm_mu_dec_out = c->ppm_mu_dec_fix;
	}
	else if(c->ppm_source == VARTOOLS_SOURCE_INLIST) {
	  ppm_mu_ra_out = c->ppm_mu_ra_vals[lcreal][0];
	  ppm_mu_dec_out = c->ppm_mu_dec_vals[lcreal][0];
	}
	else
	  error(ERR_CODEERROR);
      } else {
	ppm_mu_ra_out = 0.; ppm_mu_dec_out = 0.;
      }
      if(c->useinputppm) {
	if(c->inputppm_source == VARTOOLS_SOURCE_FIXED) {
	  ppm_mu_ra_in = c->inputppm_mu_ra_fix;
	  ppm_mu_dec_in = c->inputppm_mu_dec_fix;
	}
	else if(c->inputppm_source == VARTOOLS_SOURCE_INLIST) {
	  ppm_mu_ra_in = c->inputppm_mu_ra_vals[lcreal][0];
	  ppm_mu_dec_in = c->inputppm_mu_dec_vals[lcreal][0];
	}
	else
	  error(ERR_CODEERROR);
      } else {
	ppm_mu_ra_in = ppm_mu_ra_out; ppm_mu_dec_in = ppm_mu_dec_out;
      }
#ifdef _HAVE_CSPICE
      if(c->source_obs_coords == VARTOOLS_SOURCE_FIXED) {
	obslat = c->obslat_fixval;
	obslong = c->obslong_fixval;
	obsalt = c->obsalt_fixval;
      }
      else if(c->source_obs_coords == VARTOOLS_SOURCE_INLIST) {
	obslat = c->obslat_listvals[lcreal];
	obslong = c->obslong_listvals[lcreal];
	obsalt = c->obsalt_listvals[lcreal];
      }
      else if(c->source_obs_coords == VARTOOLS_SOURCE_NONE) {
	obslat = -999.; obslong = -999.; obsalt = -99999.;
      }
#endif
    }


  /* Do the time conversion */
  switch(c->inputtimetype) {
  case TIMETYPE_MJD:
    switch(c->outputtimetype) {
    case TIMETYPE_MJD:

      /* MJD -> MJD, no system change */
      if(c->inputsys == c->outputsys) {
	for(i=0; i < N; i++) {
	  t[i] = t[i] + (c->inputsubtractval - c->outputsubtractval);
	}
      }
#ifdef _HAVE_CSPICE

      /* MJD -> MJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = mjd2jd(delt,&tint);
	  delt = jdutc2jdtdb(delt,tint);
	  delt = jd2mjd(delt,&tint);
	  t[i] = delt + tint - c->outputsubtractval;
	}
      }

      /* MJD -> MJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = mjd2jd(delt,&tint);
	  delt = jdtdb2jdutc(delt,tint);
	  delt = jd2mjd(delt,&tint);
	  t[i] = (delt + tint)
	    - c->outputsubtractval;
	}
      }
#endif
      else
	error(ERR_CODEERROR);
      break;

    case TIMETYPE_JD:

      /* MJD -> JD, no system change */
      if(c->inputsys == c->outputsys) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = mjd2jd(delt,&tint);
	  t[i] = delt+tint - c->outputsubtractval;
	}
      }
#ifdef _HAVE_CSPICE

      /* MJD -> JD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = mjd2jd(delt,&tint);
	  t[i] = jdutc2jdtdb(delt,tint)+tint
	    - c->outputsubtractval;
	}
      }

      /* MJD -> JD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = mjd2jd(delt,&tint);
	  t[i] = jdtdb2jdutc(delt,tint)+tint
	    - c->outputsubtractval;
	}
      }
#endif
      else
	error(ERR_CODEERROR);
      break;

    case TIMETYPE_HJD:

      /* MJD -> HJD, no system change */
      if(c->inputsys == c->outputsys) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = mjd2jd(delt,&tint);
	  t[i] = jd2hjd(delt, tint, raout,
			decout, epochout, ppm_mu_ra_out, ppm_mu_dec_out,
			c->inputsys) + tint
	    - c->outputsubtractval;
	}
      }
#ifdef _HAVE_CSPICE

      /* MJD -> HJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = mjd2jd(delt,&tint);
	  t[i] = jd2hjd(jdutc2jdtdb(delt,tint),tint,
			raout, decout, epochout, ppm_mu_ra_out,
			ppm_mu_dec_out,
			TIMESYSTEM_TDB) + tint
	    - c->outputsubtractval;
	}
      }

      /* MJD -> HJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = mjd2jd(delt,&tint);
	  t[i] = jd2hjd(jdtdb2jdutc(delt,tint), tint,
			raout, decout, epochout, ppm_mu_ra_out,
			ppm_mu_dec_out,
			TIMESYSTEM_UTC) + tint
	    - c->outputsubtractval;
	}
      }
#endif
      else
	error(ERR_CODEERROR);
      break;
#ifdef _HAVE_CSPICE

    case TIMETYPE_BJD:

      /* MJD -> BJD, no system change */
      if(c->inputsys == c->outputsys) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = mjd2jd(delt,&tint);
	    t[i] = jd2bjd(delt, tint, raout,
			  decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  c->inputsys) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = mjd2jd(delt,&tint);
	    t[i] = jd2bjd(delt, tint, raout,
			  decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  c->inputsys) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* MJD -> BJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = mjd2jd(delt,&tint);
	    t[i] = jd2bjd(jdutc2jdtdb(delt,tint), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_TDB) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = mjd2jd(delt,&tint);
	    t[i] = jd2bjd(jdutc2jdtdb(delt,tint), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_TDB) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* MJD -> BJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = mjd2jd(delt,&tint);
	    t[i] = jd2bjd(jdtdb2jdutc(delt,tint), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_UTC) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = mjd2jd(delt,&tint);
	    t[i] = jd2bjd(jdtdb2jdutc(delt,tint), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_UTC)
	      - c->outputsubtractval;
	  }
	}
      }
      else
	error(ERR_CODEERROR);
      break;
#endif
    default:
      error(ERR_CODEERROR);
    }
    break;

  case TIMETYPE_JD:
    switch(c->outputtimetype) {

    case TIMETYPE_MJD:

      /* JD -> MJD, no system change */
      if(c->inputsys == c->outputsys) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = jd2mjd(delt,&tint);
	  t[i] = delt + tint - c->outputsubtractval;
	}
      }
#ifdef _HAVE_CSPICE

      /* JD -> MJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = jdutc2jdtdb(delt,tint);
	  delt = jd2mjd(delt,&tint);
	  t[i] = delt + tint
	    - c->outputsubtractval;
	}
      }

      /* JD -> MJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = jdtdb2jdutc(delt,tint);
	  delt = jd2mjd(delt,&tint);
	  t[i] = delt + tint
	    - c->outputsubtractval;
	}
      }
#endif
      else
	error(ERR_CODEERROR);
      break;

    case TIMETYPE_JD:

      /* JD -> JD, no system change */
      if(c->inputsys == c->outputsys) {
	for(i=0; i < N; i++)
	  t[i] = t[i] + c->inputsubtractval - c->outputsubtractval;
      }
#ifdef _HAVE_CSPICE

      /* JD -> JD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  t[i] = jdutc2jdtdb(delt,tint) + tint
	    - c->outputsubtractval;
	}
      }

      /* JD -> JD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  t[i] = jdtdb2jdutc(delt,tint) + tint
	    - c->outputsubtractval;
	}
      }
#endif
      else
	error(ERR_CODEERROR);
      break;

    case TIMETYPE_HJD:

      /* JD -> HJD, no system change */
      if(c->inputsys == c->outputsys) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  t[i] = jd2hjd(delt, tint, raout,
			decout, epochout, ppm_mu_ra_out, ppm_mu_dec_out,
			c->inputsys) + tint
	    - c->outputsubtractval;
	}
      }
#ifdef _HAVE_CSPICE

      /* JD -> HJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  t[i] = jd2hjd(jdutc2jdtdb(delt,tint), tint,
			raout, decout, epochout, ppm_mu_ra_out,
			ppm_mu_dec_out,
			TIMESYSTEM_TDB) + tint
	    - c->outputsubtractval;
	}
      }

      /* JD -> HJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  t[i] = jd2hjd(jdtdb2jdutc(delt,tint), tint,
			raout, decout, epochout, ppm_mu_ra_out,
			ppm_mu_dec_out,
			TIMESYSTEM_UTC) + tint
	    - c->outputsubtractval;
	}
      }
#endif
      else
	error(ERR_CODEERROR);
      break;

#ifdef _HAVE_CSPICE

    case TIMETYPE_BJD:

      /* JD -> BJD, no system change */
      if(c->inputsys == c->outputsys) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd(delt, tint, raout,
			  decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  c->inputsys) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd(delt, tint, raout,
			  decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  c->inputsys) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* JD -> BJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd(jdutc2jdtdb(delt,tint), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_TDB) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd(jdutc2jdtdb(delt,tint), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_TDB) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* JD -> BJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd(jdtdb2jdutc(delt,tint), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_UTC) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd(jdtdb2jdutc(delt,tint), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_UTC) + tint
	      - c->outputsubtractval;
	  }
	}
      }
      else
	error(ERR_CODEERROR);
      break;
#endif
    default:
      error(ERR_CODEERROR);
    }
    break;

  case TIMETYPE_HJD:
    switch(c->outputtimetype) {

    case TIMETYPE_MJD:

      /* HJD -> MJD, no system change */
      if(c->inputsys == c->outputsys) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = hjd2jd(delt, tint, rain,
			decin, epochin, ppm_mu_ra_in,
			ppm_mu_dec_in,
			c->inputsys);
	  delt = jd2mjd(delt,&tint);
	  t[i] = delt + tint
	    - c->outputsubtractval;
	}
      }
#ifdef _HAVE_CSPICE

      /* HJD -> MJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = hjd2jd(jdutc2jdtdb(delt,tint), tint,
			rain, decin, epochin, ppm_mu_ra_in,
			ppm_mu_dec_in,
			TIMESYSTEM_TDB);
	  delt = jd2mjd(delt,&tint);
	  t[i] = delt + tint
	    - c->outputsubtractval;
	}
      }

      /* HJD -> MJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  delt = hjd2jd(jdtdb2jdutc(delt,tint), tint,
			rain, decin, epochin, ppm_mu_ra_in,
			ppm_mu_dec_in,
			TIMESYSTEM_UTC);
	  delt = jd2mjd(delt,&tint);
	  t[i] = delt + tint
	    - c->outputsubtractval;
	}
      }
#endif
      else
	error(ERR_CODEERROR);
      break;

    case TIMETYPE_JD:

      /* HJD -> JD, no system change */
      if(c->inputsys == c->outputsys) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  t[i] = hjd2jd(delt, tint, rain,
			decin, epochin, ppm_mu_ra_in, ppm_mu_dec_in,
			c->inputsys) + tint
	    - c->outputsubtractval;
	}
      }
#ifdef _HAVE_CSPICE

      /* HJD -> JD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  t[i] = hjd2jd(jdutc2jdtdb(delt, tint),tint,
			rain, decin, epochin, ppm_mu_ra_in,
			ppm_mu_dec_in,
			TIMESYSTEM_TDB) + tint
	    - c->outputsubtractval;
	}
      }

      /* HJD -> JD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	for(i=0; i < N; i++) {
	  tint = floor(t[i] + c->inputsubtractval);
	  delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	  t[i] = hjd2jd(jdtdb2jdutc(delt,tint), tint,
			rain, decin, epochin, ppm_mu_ra_in,
			ppm_mu_dec_in, TIMESYSTEM_UTC) + tint
	    - c->outputsubtractval;
	}
      }
#endif
      else
	error(ERR_CODEERROR);
      break;

    case TIMETYPE_HJD:

      /* HJD -> HJD, no system change */
      if(c->inputsys == c->outputsys) {
	if(rain != raout || decin != decout || epochin != epochout ||
	   ppm_mu_ra_in != ppm_mu_ra_out || ppm_mu_dec_in != ppm_mu_dec_out) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2hjd(hjd2jd(delt,tint ,
				 rain, decin, epochin, ppm_mu_ra_in,
				 ppm_mu_dec_in,c->inputsys), tint,
			  raout, decout, epochout, ppm_mu_ra_out,
			  ppm_mu_dec_out,
			  c->inputsys) + tint - c->outputsubtractval;
	  }
	}
	else {
	  for(i=0; i < N; i++)
	    t[i] = t[i] + c->inputsubtractval - c->outputsubtractval;
	}
      }
#ifdef _HAVE_CSPICE

      /* HJD -> HJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	if(rain != raout || decin != decout || epochin != epochout ||
	   ppm_mu_ra_in != ppm_mu_ra_out || ppm_mu_dec_in != ppm_mu_dec_out) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2hjd(hjd2jd(jdutc2jdtdb(delt,tint),tint,
				 rain, decin, epochin, ppm_mu_ra_in,
				 ppm_mu_dec_in,TIMESYSTEM_TDB),tint,
			  raout, decout, epochout, ppm_mu_ra_out,
			  ppm_mu_dec_out,
			  TIMESYSTEM_TDB) + tint - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jdutc2jdtdb(delt,tint) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* HJD -> HJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	if(rain != raout || decin != decout || epochin != epochout ||
	   ppm_mu_ra_in != ppm_mu_ra_out || ppm_mu_dec_in != ppm_mu_dec_out) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2hjd((hjd2jd(jdtdb2jdutc(delt,tint),tint,
				  rain, decin, epochin, ppm_mu_ra_in,
				  ppm_mu_dec_in,TIMESYSTEM_UTC)),tint,
			  raout, decout, epochout, ppm_mu_ra_out,
			  ppm_mu_dec_out,
			  TIMESYSTEM_UTC) + tint - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jdtdb2jdutc(delt,tint) + tint
	      - c->outputsubtractval;
	  }
	}
      }
#endif
      else
	error(ERR_CODEERROR);
      break;

#ifdef _HAVE_CSPICE

    case TIMETYPE_BJD:

      /* HJD -> BJD, no system change */
      if(c->inputsys == c->outputsys) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd((hjd2jd(delt, tint, rain,
				  decin, epochin, ppm_mu_ra_in,
				  ppm_mu_dec_in,c->inputsys)),tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,c->inputsys) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd((hjd2jd(delt, tint, rain,
				  decin, epochin, ppm_mu_ra_in,
				  ppm_mu_dec_in,c->inputsys)), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,c->inputsys) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* HJD -> BJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd((hjd2jd((jdutc2jdtdb(delt,tint)), tint,
				  rain,
				  decin, epochin, ppm_mu_ra_in,
				  ppm_mu_dec_in,TIMESYSTEM_TDB)), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,TIMESYSTEM_TDB) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd((hjd2jd((jdutc2jdtdb(delt,tint)), tint,
				  rain,
				  decin, epochin, ppm_mu_ra_in,
				  ppm_mu_dec_in,TIMESYSTEM_TDB)), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,TIMESYSTEM_TDB) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* HJD -> BJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd((hjd2jd((jdtdb2jdutc(delt,tint)), tint,
				  rain,
				  decin, epochin, ppm_mu_ra_in,
				  ppm_mu_dec_in,TIMESYSTEM_UTC)), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,TIMESYSTEM_UTC) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd((hjd2jd((jdtdb2jdutc(delt,tint)), tint,
				  rain,
				  decin, epochin, ppm_mu_ra_in,
				  ppm_mu_dec_in,TIMESYSTEM_UTC)), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out, ppm_mu_dec_out,TIMESYSTEM_UTC) + tint
	      - c->outputsubtractval;
	  }
	}
      }
      else
	error(ERR_CODEERROR);
      break;
#endif
    default:
      error(ERR_CODEERROR);
    }
    break;
#ifdef _HAVE_CSPICE
  case TIMETYPE_BJD:
    switch(c->outputtimetype) {

    case TIMETYPE_MJD:

      /* BJD -> MJD, no system change */
      if(c->inputsys == c->outputsys) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = bjd2jd(delt, tint, rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, c->inputsys);
	    delt = jd2mjd(delt,&tint);
	    t[i] = delt + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = bjd2jd(delt, tint, rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, c->inputsys);
	    delt = jd2mjd(delt,&tint);
	    t[i] = delt + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* BJD -> MJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = bjd2jd((jdutc2jdtdb(delt,tint)), tint,
			  rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, TIMESYSTEM_TDB);
	    delt = jd2mjd(delt,&tint);
	    t[i] = delt + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = bjd2jd(jdutc2jdtdb(delt,tint), tint,
			  rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, TIMESYSTEM_TDB);
	    delt = jd2mjd(delt,&tint);
	    t[i] = delt + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* BJD -> MJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = bjd2jd(jdtdb2jdutc(delt,tint), tint,
			  rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, TIMESYSTEM_UTC);
	    delt = jd2mjd(delt,&tint);
	    t[i] = delt + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    delt = bjd2jd(jdtdb2jdutc(delt,tint), tint,
			  rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, TIMESYSTEM_UTC);
	    delt = jd2mjd(delt,&tint);
	    t[i] = delt + tint
	      - c->outputsubtractval;
	  }
	}
      }
      else
	error(ERR_CODEERROR);
      break;

    case TIMETYPE_JD:

      /* BJD -> JD, no system change */
      if(c->inputsys == c->outputsys) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = bjd2jd(delt, tint, rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, c->inputsys) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = bjd2jd(delt, tint, rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, c->inputsys) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* BJD -> JD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = bjd2jd(jdutc2jdtdb(delt,tint), tint, rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, TIMESYSTEM_TDB) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = bjd2jd(jdutc2jdtdb(delt,tint), tint, rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, TIMESYSTEM_TDB) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* BJD -> JD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = bjd2jd(jdtdb2jdutc(delt,tint), tint, rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, TIMESYSTEM_UTC) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = bjd2jd(jdtdb2jdutc(delt,tint), tint, rain,
			  decin, epochin,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_in, ppm_mu_dec_in, TIMESYSTEM_UTC) + tint
	      - c->outputsubtractval;
	  }
	}
      }
      else
	error(ERR_CODEERROR);
      break;

    case TIMETYPE_HJD:

      /* BJD -> HJD, no system change */
      if(c->inputsys == c->outputsys) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2hjd((bjd2jd(delt, tint, rain,
				  decin, epochin,
				  obslat, obslong, obsalt,
				  ppm_mu_ra_in, ppm_mu_dec_in,
				  c->inputsys)), tint,
			  raout, decout, epochout,
			  ppm_mu_ra_out, ppm_mu_dec_out, c->inputsys) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2hjd((bjd2jd(delt, tint, rain,
				  decin, epochin,
				  obslat, obslong, obsalt,
				  ppm_mu_ra_in, ppm_mu_dec_in,
				  c->inputsys)), tint,
			  raout, decout, epochout,
			  ppm_mu_ra_out, ppm_mu_dec_out, c->inputsys) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* BJD -> HJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2hjd((bjd2jd((jdutc2jdtdb(delt,tint)), tint,
				  rain,
				  decin, epochin,
				  obslat, obslong, obsalt,
				  ppm_mu_ra_in, ppm_mu_dec_in,
				  TIMESYSTEM_TDB)), tint,
			  raout, decout, epochout,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_TDB) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2hjd((bjd2jd((jdutc2jdtdb(delt,tint)), tint,
				  rain,
				  decin, epochin,
				  obslat, obslong, obsalt,
				  ppm_mu_ra_in, ppm_mu_dec_in,
				  TIMESYSTEM_TDB)), tint,
			  raout, decout, epochout,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_TDB) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* HJD -> BJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	if(c->source_obs_coords != VARTOOLS_SOURCE_LC) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2hjd((bjd2jd((jdtdb2jdutc(delt,tint)), tint,
				  rain,
				  decin, epochin,
				  obslat, obslong, obsalt,
				  ppm_mu_ra_in, ppm_mu_dec_in,
				  TIMESYSTEM_UTC)), tint,
			  raout, decout, epochout,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_UTC) + tint
	      - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    obslat = c->obslat_lcvals[lc][i];
	    obslong = c->obslong_lcvals[lc][i];
	    obsalt = c->obsalt_lcvals[lc][i];
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2hjd((bjd2jd((jdtdb2jdutc(delt,tint)), tint,
				  rain,
				  decin, epochin,
				  obslat, obslong, obsalt,
				  ppm_mu_ra_in, ppm_mu_dec_in,
				  TIMESYSTEM_UTC)), tint,
			  raout, decout, epochout,
			  ppm_mu_ra_out, ppm_mu_dec_out,
			  TIMESYSTEM_UTC) + tint
	      - c->outputsubtractval;
	  }
	}
      }
      else
	error(ERR_CODEERROR);
      break;
    case TIMETYPE_BJD:

      /* BJD -> BJD, no system change */
      if(c->inputsys == c->outputsys) {
	if(rain != raout || decin != decout || epochin != epochout ||
	   ppm_mu_ra_in != ppm_mu_ra_out || ppm_mu_dec_in != ppm_mu_dec_out) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd((bjd2jd(delt, tint,
				  rain, decin, epochin,
				  obslat, obslong, obsalt,
				  ppm_mu_ra_in,
				  ppm_mu_dec_in, c->inputsys)), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out,
			  ppm_mu_dec_out,
			  c->inputsys) + tint - c->outputsubtractval;
	  }
	}
	else {
	  for(i=0; i < N; i++)
	    t[i] = t[i] + c->inputsubtractval - c->outputsubtractval;
	}
      }


      /* BJD -> BJD, UTC -> TDB */
      else if(c->inputsys == TIMESYSTEM_UTC) {
	if(rain != raout || decin != decout || epochin != epochout ||
	   ppm_mu_ra_in != ppm_mu_ra_out || ppm_mu_dec_in != ppm_mu_dec_out) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd((bjd2jd(jdutc2jdtdb(delt,tint), tint,
				  rain, decin, epochin,
				  obslat, obslong, obsalt,
				  ppm_mu_ra_in,
				  ppm_mu_dec_in,
				  TIMESYSTEM_TDB)), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out,
			  ppm_mu_dec_out,
			  TIMESYSTEM_TDB) + tint - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jdutc2jdtdb(delt,tint) + tint
	      - c->outputsubtractval;
	  }
	}
      }

      /* BJD -> BJD, TDB -> UTC */
      else if(c->inputsys == TIMESYSTEM_TDB) {
	if(rain != raout || decin != decout || epochin != epochout ||
	   ppm_mu_ra_in != ppm_mu_ra_out || ppm_mu_dec_in != ppm_mu_dec_out) {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jd2bjd((bjd2jd(jdtdb2jdutc(delt,tint), tint,
				  rain, decin, epochin,
				  obslat, obslong, obsalt,
				  ppm_mu_ra_in,
				  ppm_mu_dec_in,
				  TIMESYSTEM_UTC)), tint,
			  raout, decout, epochout,
			  obslat, obslong, obsalt,
			  ppm_mu_ra_out,
			  ppm_mu_dec_out,
			  TIMESYSTEM_UTC) + tint - c->outputsubtractval;
	  }
	} else {
	  for(i=0; i < N; i++) {
	    tint = floor(t[i] + c->inputsubtractval);
	    delt = (t[i] + c->inputsubtractval) - floor(t[i] + c->inputsubtractval);
	    t[i] = jdtdb2jdutc(delt,tint) + tint
	      - c->outputsubtractval;
	  }
	}
      }
      else
	error(ERR_CODEERROR);
      break;
    default:
      error(ERR_CODEERROR);
    }
    break;
#endif
  default:
    error(ERR_CODEERROR);
  }
}
