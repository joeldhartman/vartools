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

#ifdef _HAVE_GSL
#include <gsl/gsl_errno.h>
#include <gsl/gsl_fft_complex.h>

#define REAL(z,i) ((z)[2*(i)])
#define IMAG(z,i) ((z)[2*(i)+1])

#endif

void RunFFTCommand(ProgramData *p, int lcindex, int threadindex, _FFT *f)
{
#ifdef _HAVE_GSL
  int i, j, k;
  double *tmpdata = NULL;

  gsl_fft_complex_wavetable * wavetable = NULL;
  gsl_fft_complex_workspace * workspace = NULL;

  if(p->NJD[threadindex] <= 0) return;

  if((tmpdata = (double *) malloc(2*p->NJD[threadindex]*sizeof(double))) == NULL)
    error(ERR_MEMALLOC);

  for(i=0; i < p->NJD[threadindex]; i++) {
    if(f->inputvar_real != NULL) {
      REAL(tmpdata,i) = EvaluateVariable_Double(lcindex, threadindex, i, f->inputvar_real);
    } else {
      REAL(tmpdata,i) = 0.0;
    }
    if(f->inputvar_imag != NULL) {
      IMAG(tmpdata,i) = EvaluateVariable_Double(lcindex, threadindex, i, f->inputvar_imag);
    } else {
      IMAG(tmpdata,i) = 0.0;
    }
  }
  
  wavetable = gsl_fft_complex_wavetable_alloc(p->NJD[threadindex]);
  workspace = gsl_fft_complex_workspace_alloc(p->NJD[threadindex]);

  if(f->isforward) 
    gsl_fft_complex_forward(tmpdata, 1, p->NJD[threadindex], wavetable, workspace);
  else
    gsl_fft_complex_inverse(tmpdata, 1, p->NJD[threadindex], wavetable, workspace);


  for(i=0; i < p->NJD[threadindex]; i++) {
    if(f->outputvar_real != NULL) {
      SetVariable_Value_Double(lcindex, threadindex, i, f->outputvar_real, REAL(tmpdata,i));
    }
    if(f->outputvar_imag != NULL) {
      SetVariable_Value_Double(lcindex, threadindex, i, f->outputvar_imag, IMAG(tmpdata,i));
    }
  }

  gsl_fft_complex_wavetable_free(wavetable);
  gsl_fft_complex_workspace_free(workspace);

  if(tmpdata != NULL)
    free(tmpdata);
#endif
  return;
}

int ParseFFTCommand(int *iret, int argc, char **argv, ProgramData *p, _FFT *f)
{
#ifdef _HAVE_GSL
  int i;
  f->inputvar_real = NULL; f->inputvar_imag = NULL;
  f->outputvar_real = NULL; f->outputvar_imag = NULL;
  i = *iret;
  if(i >= argc)
    return 1;
  if(!strcmp(argv[i],"NULL")) {
    f->inputvarname_real = malloc(1*sizeof(char));
    f->inputvarname_real[0] = '\0';
  } else {
    if((f->inputvarname_real = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    sprintf(f->inputvarname_real,"%s",argv[i]);
  }

  i++;
  if(i >= argc)
    return 1;
  if(!strcmp(argv[i],"NULL")) {
    f->inputvarname_imag = malloc(1*sizeof(char));
    f->inputvarname_imag[0] = '\0';
  } else {
    if((f->inputvarname_imag = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    sprintf(f->inputvarname_imag,"%s",argv[i]);
  }

  i++;
  if(i >= argc)
    return 1;
  if(!strcmp(argv[i],"NULL")) {
    f->outputvarname_real = malloc(1*sizeof(char));
    f->outputvarname_real[0] = '\0';
  } else {
    if((f->outputvarname_real = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    sprintf(f->outputvarname_real,"%s",argv[i]);
  }

  i++;
  if(i >= argc)
    return 1;
  if(!strcmp(argv[i],"NULL")) {
    f->outputvarname_imag = malloc(1*sizeof(char));
    f->outputvarname_imag[0] = '\0';
  } else {
    if((f->outputvarname_imag = (char *) malloc((strlen(argv[i])+1)*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
    sprintf(f->outputvarname_imag,"%s",argv[i]);
  }

  *iret = i;
  return 0;
    
#endif
  return 1;
}
