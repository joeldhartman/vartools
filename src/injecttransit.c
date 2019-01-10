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

#define RSUNTOAU 0.004649126
#define YEARSTODAYS 365.242199
#define RJUPTORSUN 0.1027922
#define MJUPTOMSUN 0.000954638698
#define TWOPI 6.28318530717958647692528676656

double get_rand_h_given_k(double k)
{
  /* This function generates a random h or k given the other, assuming random e and omega */
  double maxh, u, N, C, lk, h, testd;
  lk = fabs(k);
  lk = log(k);
  maxh = sqrt(1. - k*k);
  u = rand() / (RAND_MAX + 0.0);
  N = log(maxh + sqrt(maxh*maxh + k*k)) - lk;
  C = exp((u*N) + lk);
  h = (C*C - k*k)/2./C;
  testd = rand() / (RAND_MAX + 0.0);
  if(testd > 0.5)
    h = -h;
  return h;
}

double get_maxcosi(double e, double omega, double r1, double r2)
{
  /* This function returns the maximum value of cos(i) that yields transits r1 and r2 are in units of the semi-major axis */
  double v1, v2, v3;
  v1 = 1. + e*sin(omega);
  v2 = 1. - e*e;
  v3 = r1 + r2;
  return v3*v1/v2;
}

void doinjecttransit(int N, double *t, double *mag, double *sig, int lc, int lcreal, _Injecttransit *c, char *modeloutname)
{
  int i, j;
  FILE *outfile;
  double ldcoeffs[4], *delmag, *phase, a, rprstar, dval, dval2, e, omega;

  /* Set the parameters */
  for(i=0;i<c->Nparam;i++)
    {
      switch(c->paramtype[i])
	{
	case PERTYPE_SPECIFIED:
	  c->paraminject[i][lc] = c->paramspec[i][lcreal][0];
	  break;
	case PERTYPE_FIX:
	  c->paraminject[i][lc] = c->paramfix[i];
	  break;
	case PERTYPE_EXPR:
	  c->paraminject[i][lc] = EvaluateExpression(lcreal, lc, 0, c->paramexpr[i]);
	  break;
	case PERTYPE_UNIFORMRAND:
	  switch(i)
	    {
	    case INJECTTR_IDX_PERIOD:
	      c->paraminject[i][lc] = c->minp + (c->maxp - c->minp)*(rand() / (RAND_MAX + 0.0));
	      break;
	    case INJECTTR_IDX_RP:
	      c->paraminject[i][lc] = c->minRp + (c->maxRp - c->minRp)*(rand() / (RAND_MAX + 0.0));
	      break;
	    case INJECTTR_IDX_MP:
	      c->paraminject[i][lc] = c->minMp + (c->maxMp - c->minMp)*(rand() / (RAND_MAX + 0.0));
	      break;
	    default:
	      c->paraminject[i][lc] = (rand() / (RAND_MAX + 0.0));
	      break;
	    }
	  break;
	case PERTYPE_LOGRAND:
	  switch(i)
	    {
	    case INJECTTR_IDX_PERIOD:
	      dval = log(c->minp) + (log(c->maxp) - log(c->minp))*(rand() / (RAND_MAX + 0.0));
	      break;
	    case INJECTTR_IDX_RP:
	      dval = log(c->minRp) + (log(c->maxRp) - log(c->minRp))*(rand() / (RAND_MAX + 0.0));
	      break;
	    case INJECTTR_IDX_MP:
	      dval = log(c->minMp) + (log(c->maxMp) - log(c->minMp))*(rand() / (RAND_MAX + 0.0));
	      break;
	    }
	  c->paraminject[i][lc] = exp(dval);
	  break;
	case PERTYPE_UNIFORMRANDFREQ:
	  dval = c->minf + (c->maxf - c->minf)*(rand() / (RAND_MAX + 0.0));
	  c->paraminject[i][lc] = 1./dval;
	  break;
	case PERTYPE_LOGRANDFREQ:
	  dval = log(c->minf) + (log(c->maxf) - log(c->minf))*(rand() / (RAND_MAX + 0.0));
	  c->paraminject[i][lc] = 1./exp(dval);
	  break;
	default:
	  error(-1);
	}
    }

  /* get the ld-coeffs into the format needed by the mandelagoltransitmodel function */
  for(i=0;i<c->Nld;i++)
    ldcoeffs[i] = c->paraminject[INJECTTR_IDX_LD + i][lc];

  /* Get a and rp/rstar */
  a = (c->paraminject[INJECTTR_IDX_PERIOD][lc] / YEARSTODAYS);
  a = a*a*(c->paraminject[INJECTTR_IDX_MSTAR][lc] + (c->paraminject[INJECTTR_IDX_MP][lc]*MJUPTOMSUN));
  a = pow(a,0.33333333333);
  a = a / (c->paraminject[INJECTTR_IDX_RSTAR][lc]*RSUNTOAU);
  rprstar = (c->paraminject[INJECTTR_IDX_RP][lc] * RJUPTORSUN)/c->paraminject[INJECTTR_IDX_RSTAR][lc];

  /* Get e/omega */
  if(!c->eomegatype)
    {
      if(c->paramtype[INJECTTR_IDX_OMEGA] == PERTYPE_UNIFORMRAND)
	{
	  c->paraminject[INJECTTR_IDX_OMEGA][lc] *= TWOPI;
	}
      else
	{
	  c->paraminject[INJECTTR_IDX_OMEGA][lc] *= (TWOPI / 360.);
	}
      e = c->paraminject[INJECTTR_IDX_E][lc];
      omega = c->paraminject[INJECTTR_IDX_OMEGA][lc];
    }
  else
    {
      if(c->paramtype[INJECTTR_IDX_H] == PERTYPE_UNIFORMRAND &&
	 c->paramtype[INJECTTR_IDX_K] == PERTYPE_UNIFORMRAND)
	{
	  dval = c->paraminject[INJECTTR_IDX_H][lc];
	  dval2 = c->paraminject[INJECTTR_IDX_K][lc] * TWOPI;
	  c->paraminject[INJECTTR_IDX_H][lc] = dval*sin(dval2);
	  c->paraminject[INJECTTR_IDX_K][lc] = dval*cos(dval2);
	}
      else if(c->paramtype[INJECTTR_IDX_H] == PERTYPE_UNIFORMRAND)
	{
	  c->paraminject[INJECTTR_IDX_H][lc] = get_rand_h_given_k(c->paraminject[INJECTTR_IDX_K][lc]);
	}
      else if(c->paramtype[INJECTTR_IDX_K] == PERTYPE_UNIFORMRAND)
	{
	  c->paraminject[INJECTTR_IDX_K][lc] = get_rand_h_given_k(c->paraminject[INJECTTR_IDX_H][lc]);
	}
      e = sqrt(c->paraminject[INJECTTR_IDX_H][lc]*c->paraminject[INJECTTR_IDX_H][lc] + c->paraminject[INJECTTR_IDX_K][lc]*c->paraminject[INJECTTR_IDX_K][lc]);
      omega = atan2(c->paraminject[INJECTTR_IDX_H][lc],c->paraminject[INJECTTR_IDX_K][lc]);
    }

  /* Get sin_i */
  if(c->paramtype[INJECTTR_IDX_SINI] == PERTYPE_UNIFORMRAND)
    {
      dval2 = get_maxcosi(e, omega, (1./a), (rprstar / a));
      dval = dval2 * rand() / (RAND_MAX + 0.0);
      c->paraminject[INJECTTR_IDX_SINI][lc] = sqrt(1. - (dval*dval));
    }

  if((delmag = (double *) malloc(N * sizeof(double))) == NULL ||
     (phase = (double *) malloc(N * sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  /* Calculate the phases */
  for(i=0;i<N;i++)
    {
      phase[i] = t[i]/c->paraminject[INJECTTR_IDX_PERIOD][lc] + c->paraminject[INJECTTR_IDX_PHASE][lc];
      phase[i] -= floor(phase[i]);
    }

  /* Generate the model transit */
  mandelagoltransitmodel(N, phase, delmag, c->ldtype, ldcoeffs, c->paraminject[INJECTTR_IDX_SINI][lc], a, e, rprstar, omega);

  /* Open the file to output the model to if we're doing that */
  if(c->omodel)
    {
      if((outfile = fopen(modeloutname,"w")) == NULL)
	error2(ERR_CANNOTWRITE, modeloutname);
    }

  for(j=0;j<N;j++)
    {
      delmag[j] = -2.5*log(1.0 - (1.0 - delmag[j])*c->paraminject[INJECTTR_IDX_DILUTE][lc])/M_LN10;
      if(c->omodel)
	{
	  fprintf(outfile,"%f %f %f %f\n",t[j],mag[j],delmag[j],sig[j]);
	}
      mag[j] += delmag[j];
    }
  if(!c->eomegatype)
    c->paraminject[INJECTTR_IDX_OMEGA][lc] *= 360./TWOPI;

  if(c->omodel)
    fclose(outfile);
}
