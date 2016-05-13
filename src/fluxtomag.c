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

void fluxtomag(double *t, double *mag, double *sig, int N, double mag_constant1, double offset)
{
  int i;
  double mag2;
  mag2 = mag_constant1 + offset;
  for(i=0; i < N; i++)
    {
      if(mag[i] > 0)
	{
	  sig[i] = 1.0857*sig[i]/mag[i];
	  mag[i] = mag2 - 2.5*log(mag[i])/M_LN10;
	}
      else
	{
	  sig[i] = sqrt(-1.);
	  mag[i] = sqrt(-1.);
	}
    }
}

void difffluxtomag(double *t, double *mag, double *sig, int N, double mag_star, double mag_constant1, double offset)
{
  int i;
  double fluxstar, mag2;

  /* Note the input flux is stored initially in the "mag" vector, the rest of the vartools program makes no distinction between whether this is a magnitude or a flux, it is up to the user to know what it is */

  mag2 = mag_constant1 + offset;
  fluxstar = pow((double) 10.0,(mag_constant1 - mag_star)/2.5);
  for(i=0;i<N;i++)
    {
      if(mag[i] >= fluxstar)
	{
	  mag[i] = sqrt(-1);
	  sig[i] = sqrt(-1);
	}
      else
	{
	  mag[i] = mag2 - 2.5*log(fluxstar - mag[i])/M_LN10;
	  sig[i] = 1.0857*sig[i]/(fluxstar - mag[i]);
	}
    }
}

