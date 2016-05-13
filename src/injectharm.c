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

#define TWOPI 6.28318530717958647692528676656

void doinjectharm(int N, double *t, double *mag, double *sig, int lc, int lcreal, _Injectharm *c, char *modeloutname)
{
  int i, j;
  double dval, m, ph, per;
  FILE *outfile;

  /* Set the period */
  switch(c->pertype)
    {
    case PERTYPE_SPECIFIED:
      c->periodinject[lc] = c->periods[lcreal][0];
      break;
    case PERTYPE_FIX:
      c->periodinject[lc] = c->fixperiod;
      break;
    case PERTYPE_UNIFORMRAND:
      c->periodinject[lc] = c->minp + (c->maxp - c->minp)*(rand() / (RAND_MAX + 0.0));
      break;
    case PERTYPE_LOGRAND:
      dval = log(c->minp) + (log(c->maxp) - log(c->minp))*(rand() / (RAND_MAX + 0.0));
      c->periodinject[lc] = exp(dval);
      break;
    case PERTYPE_UNIFORMRANDFREQ:
      c->periodinject[lc] = c->minf + (c->maxf - c->minf)*(rand() / (RAND_MAX + 0.0));
      c->periodinject[lc] = 1./c->periodinject[lc];
      break;
    case PERTYPE_LOGRANDFREQ:
      dval = log(c->minp) + (log(c->maxp) - log(c->minp))*(rand() / (RAND_MAX + 0.0));
      c->periodinject[lc] = 1./exp(dval);
      break;
    default:
      error(-1);
    }

  /* Now set each of the harmonics and sub-harmonics */
  for(i=0;i<=c->Nharm;i++)
    {
      switch(c->harm_amptype[i])
	{
	case PERTYPE_SPECIFIED:
	  c->harm_amp[lc][i] = c->harm_ampspec[i][lcreal][0];
	  break;
	case PERTYPE_FIX:
	  c->harm_amp[lc][i] = c->harm_ampfix[i];
	  break;
	case PERTYPE_UNIFORMRAND:
	  c->harm_amp[lc][i] = c->harm_minamp[i] + (c->harm_maxamp[i] - c->harm_minamp[i])*(rand() / (RAND_MAX + 0.0));
	  break;
	case PERTYPE_LOGRAND:
	  c->harm_amp[lc][i] = log(c->harm_minamp[i]) + (log(c->harm_maxamp[i]) - log(c->harm_minamp[i]))*(rand() / (RAND_MAX + 0.0));
	  c->harm_amp[lc][i] = exp(c->harm_amp[lc][i]);
	  break;
	default:
	  error(-1);
	}
      if(i && c->harm_amprel[i])
	c->harm_amp[lc][i] *= c->harm_amp[lc][0];
      switch(c->harm_phasetype[i])
	{
	case PERTYPE_SPECIFIED:
	  c->harm_phase[lc][i] = c->harm_phasespec[i][lcreal][0];
	  break;
	case PERTYPE_FIX:
	  c->harm_phase[lc][i] = c->harm_phasefix[i];
	  break;
	case PERTYPE_UNIFORMRAND:
	  c->harm_phase[lc][i] = (rand() / (RAND_MAX + 0.0));
	  break;
	default:
	  error(-1);
	}
      if(i && c->harm_phaserel[i])
	c->harm_phase[lc][i] += ((double) (i + 1)) * c->harm_phase[lc][0];
    }
  for(i=0;i<c->Nsubharm;i++)
    {
      switch(c->subharm_amptype[i])
	{
	case PERTYPE_SPECIFIED:
	  c->subharm_amp[lc][i] = c->subharm_ampspec[i][lcreal][0];
	  break;
	case PERTYPE_FIX:
	  c->subharm_amp[lc][i] = c->subharm_ampfix[i];
	  break;
	case PERTYPE_UNIFORMRAND:
	  c->subharm_amp[lc][i] = c->subharm_minamp[i] + (c->subharm_maxamp[i] - c->subharm_minamp[i])*(rand() / (RAND_MAX + 0.0));
	  break;
	case PERTYPE_LOGRAND:
	  c->subharm_amp[lc][i] = log(c->subharm_minamp[i]) + (log(c->subharm_maxamp[i]) - log(c->subharm_minamp[i]))*(rand() / (RAND_MAX + 0.0));
	  c->subharm_amp[lc][i] = exp(c->subharm_amp[lc][i]);
	  break;
	default:
	  error(-1);
	}
      if(c->subharm_amprel[i])
	c->subharm_amp[lc][i] *= c->harm_amp[lc][0];
      switch(c->subharm_phasetype[i])
	{
	case PERTYPE_SPECIFIED:
	  c->subharm_phase[lc][i] = c->subharm_phasespec[i][lcreal][0];
	  break;
	case PERTYPE_FIX:
	  c->subharm_phase[lc][i] = c->subharm_phasefix[i];
	  break;
	case PERTYPE_UNIFORMRAND:
	  c->subharm_phase[lc][i] = (rand() / (RAND_MAX + 0.0));
	  break;
	default:
	  error(-1);
	}
      if(c->subharm_phaserel[i])
	c->subharm_phase[lc][i] += c->harm_phase[lc][0] / ((double) (i + 2));
    }

  /* Open the file to output the model to if we're doing that */
  if(c->omodel)
    {
      if((outfile = fopen(modeloutname,"w")) == NULL)
	error2(ERR_CANNOTWRITE, modeloutname);
    }

  for(j=0;j<N;j++)
    {
      m = 0.;
      for(i=0;i<=c->Nharm;i++)
	{
	  per = c->periodinject[lc] / (i + 1.0);
	  ph = TWOPI * (c->harm_phase[lc][i] + (t[j]/per));

	  m += c->harm_amp[lc][i]*cos(ph);
	}
      for(i=0;i<c->Nsubharm;i++)
	{
	  per = c->periodinject[lc] * (i + 1.0);
	  ph = TWOPI * (c->subharm_phase[lc][i] + (t[j]/per));
	  m += c->subharm_amp[lc][i]*cos(ph);
	}
      if(c->omodel)
	{
	  fprintf(outfile,"%f %f %f %f\n",t[j],mag[j], m, sig[j]);
	}
      mag[j] += m;
    }
  if(c->omodel)
    fclose(outfile);
}


