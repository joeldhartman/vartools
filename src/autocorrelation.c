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

void autocorrelation(double *t, double *mag, double *sig, int N, double tmin, double tmax, double tstep, char *outname)
{
  int i, j, u, l;
  double t_;
  long double var1;
  double var2, mag_ave, var3;
  FILE *outfile;

  double *udcf, *eudcf;
  int *Nudcf;

  int Nbins;

  if((outfile = fopen(outname,"w")) == NULL)
    {
      fprintf(stderr,"Cannot write autocorrelation to %s\n",outname);
      exit(3);
    }
  fprintf(outfile,"#Time_separation Autocorrelation Autocorrelation_error Npairs\n");

  /* allocate space for the vectors */
  Nbins = ceil((tmax - tmin)/tstep) + 1;
  if((udcf = (double *) malloc(Nbins * sizeof(double))) == NULL ||
     (Nudcf = (int *) malloc(Nbins * sizeof(int))) == NULL ||
     (eudcf = (double *) malloc(Nbins * sizeof(double))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(3);
    }

  /* Zero the vectors */
  for(i=0;i<Nbins;i++)
    {
      udcf[i] = 0.;
      Nudcf[i] = 0;
    }

  /* Get the error weighted average */
  var2 = 0.;
  var3 = 0.;
  for(i=0;i<N;i++)
    {
      var2 += mag[i]/sig[i]/sig[i];
      var3 += 1./sig[i]/sig[i];
    }
  mag_ave = var2/var3;

  /* calculate all point-wise correlations, and add them to the appropriate ACF bin */
  for(i=0;i<N;i++)
    {
      for(j=i;j<N;j++)
	{
	  t_ = fabs(t[i] - t[j]);
	  /* Find the bin that this point contributes to */
	  l = floor((t_ - tmin + tstep/2.)/tstep);
	  if(l > -1 && l < Nbins)
	    {
	      /* increment the udcf */
	      udcf[l] += (mag[i] - mag_ave)*(mag[j] - mag_ave)/sig[i]/sig[j];
	      eudcf[l] += ((mag[i] - mag_ave)*(mag[i] - mag_ave)/sig[i]/sig[i]) + ((mag[j] - mag_ave)*(mag[j] - mag_ave)/sig[j]/sig[j]);
	      Nudcf[l]++;
	    }
	}
    }

  /* Output the ACF */

  t_ = tmin;
  for(i=0,t_=tmin;i<Nbins;i++,t_ += tstep)
    {
      if(Nudcf[i] > 1)
	{
	  udcf[i] /= (double) Nudcf[i];
	  eudcf[i] = sqrt(eudcf[i]) / (double) Nudcf[i];
	  fprintf(outfile,"%f %f %f %d\n",t_, udcf[i], eudcf[i], Nudcf[i]);
	}
    }
  free(udcf);
  free(eudcf);
  free(Nudcf);
  fclose(outfile);
}

