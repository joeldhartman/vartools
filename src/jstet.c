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

/* This file contains functiosn to compute Stetson's J variability statistics for the program vartools by J. Hartman */

#define abs_(val) ((val)>0?(val):-(val))
void w_ave(int ngood, double *data, double *sig, double *ws, double *ave)
{
  int j;
  double weight = 0.0, a=2.0, b=2.0, delta, dumval;
  *ave=0.0;
  for(j=0;j<ngood;j++)
    {
      if(!isnan(data[j]) && !isnan(ws[j]) && !isnan(sig[j]))
	{
	  *ave += ws[j]*data[j]/(sig[j]*sig[j]);
	  weight += ws[j]/(sig[j]*sig[j]);
	}
    }
  *ave /= weight;
  dumval = sqrt((double) ngood/(double) (ngood-1));
  for(j=0;j<ngood;j++)
    {
      if(!isnan(data[j]) && !isnan(sig[j]))
	{
	  delta=dumval*(data[j]-*ave)/sig[j];
	  ws[j]=1.0/(1.0+pow(abs_(delta)/a,b));
	}
    }
}

void getJstet(int ngood_in, double tmin, double wkmax, double *time_in, double *mag_in, double *sig_in, double *wtave, double *jst, double *kur, double *lst, int lcnum, int lclistnum, int usemask, _Variable *maskvar)
{
  int i, j, flag = 0, npair = 0;
  double dt, wave = 0.0, sigma = 0.0, *ws, wtave1, delta, *delt, *pk, *wk, dumval, s=0.0, weight;
  int ngood;
  double *time, *mag, *sig;
  double *time_mask = NULL, *mag_mask = NULL, *sig_mask = NULL;

  if(ngood_in < 1) {
    *wtave = 0.;
    *jst = 0.;
    *kur = 0.;
    *lst = 0.;
    return;
  }

  if(!usemask) {
    ngood = ngood_in;
    time = time_in;
    mag = mag_in;
    sig = sig_in;
  } else {
    if((time_mask = (double *) malloc(ngood_in * sizeof(double))) == NULL ||
       (mag_mask = (double *) malloc(ngood_in * sizeof(double))) == NULL ||
       (sig_mask = (double *) malloc(ngood_in * sizeof(double))) == NULL)
      error(ERR_MEMALLOC);
    ngood = 0;
    for(i = 0; i < ngood_in; i++) {
      if(EvaluateVariable_Double(lclistnum, lcnum, i, maskvar) > VARTOOLS_MASK_TINY) {
	time_mask[ngood] = time_in[i];
	mag_mask[ngood] = mag_in[i];
	sig_mask[ngood] = sig_in[i];
	ngood++;
      }
    }
    if(ngood < 1 ) {
      *wtave = 0.;
      *jst = 0.;
      *kur = 0.;
      *lst = 0.;
      free(time_mask);
      free(mag_mask);
      free(sig_mask);
      return;
    }
    time = time_mask;
    mag = mag_mask;
    sig = sig_mask;
  }


  /* Compute weighted average magnitude */
  for(i=0; i<ngood; i++)
    {
      if(!isnan(mag[i]) && !isnan(sig[i]))
	{
	  sigma += 1.0/(sig[i]*sig[i]);
	  wave += mag[i]/(sig[i]*sig[i]);
	}
    }
  wave /= sigma;

  /* Compute avestet */
  ws = (double *) malloc(2*ngood*sizeof(double));
  for(i=0; i<ngood; i++) ws[i] = 1.0;
  w_ave(ngood,mag,sig,ws,wtave);
  do
    {
      w_ave(ngood,mag,sig,ws,&wtave1);
      delta = abs_(*wtave-wtave1);
      *wtave=wtave1;
    } while(delta > 0.00001);
  free(ws);

  /* Compute j_stet */
  dumval = sqrt((double) ngood/(double) (ngood-1));
  delt = (double *) malloc(2*ngood*sizeof(double));
  pk = (double *) malloc(2*ngood*sizeof(double));
  wk = (double *) malloc(2*ngood*sizeof(double));
  for(j=0;j<ngood;j++)
    if(!isnan(mag[j]) && !isnan(sig[j]))
      delt[j]=dumval*(mag[j]- *wtave)/sig[j];
  *kur = 0.0;
  for(j=0;j<ngood;j++){
    if(!isnan(mag[j]) && !isnan(sig[j]))
      {
	*kur += abs_(delt[j]);
	s += (delt[j]*delt[j]);
      }
  }
  *kur /= (double) ngood;
  s = 0.798*sqrt(s/(double) ngood);
  *kur /= s;

  flag = 0;
  for(j=0;j<ngood-1;j++){
    if(!isnan(mag[j]) && !isnan(mag[j+1]) && !isnan(sig[j]) && !isnan(sig[j+1]))
      {
	dt=abs_(time[j+1]-time[j]);
	if(dt > tmin && !flag){
	  pk[npair]=(delt[j]*delt[j])-1.0;
	  wk[npair]=0.1;
	  npair++;
	  flag = 0;
	}
	if(dt > tmin && flag) flag=0;
	if(dt <= tmin){
	  pk[npair]=delt[j+1]*delt[j];
	  wk[npair]=1.0;
	  npair++;
	  flag = 1;
	}
      }
  }
  if(!flag){
    pk[npair]=(delt[ngood-1]*delt[ngood-1])-1.0;
    wk[npair]=0.1;
    npair++;
  }
  weight = 0.0;
  *jst=0.0;
  for(j=0;j<npair;j++)
    {
      s=pk[j]/abs_(pk[j])*sqrt(abs_(pk[j]));
      *jst += wk[j]*s;
      weight += wk[j];
    }
  *jst /= weight;
  /* Compute lst */
  *lst = (*jst) * (*kur) * weight / wkmax;
  *jst = (*jst) * weight / wkmax;
  free(delt);
  free(pk);
  free(wk);
  if(time_mask != NULL) free(time_mask);
  if(mag_mask != NULL) free(mag_mask);
  if(sig_mask != NULL) free(sig_mask);
}

#undef abs_
