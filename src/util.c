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

void checkUTCFormat(char *UTCformat, int *UTCindex)
{
  int i = 0, j;
  int nscanned = 0;
  for(j=0; j < 6; j++)
    UTCindex[j] = -1;
  if(UTCformat == NULL)
    error(ERR_INVALID_UTC_FORMAT);
  while(UTCformat[i] != '\0') {
    if(UTCformat[i] == '%') {
      i++;
      if(UTCformat[i] == '\0') {
	error(ERR_INVALID_UTC_FORMAT);
      }
      j = -1;
      switch(UTCformat[i]) {
      case '%':
	break;
      case 'Y':
	j = 0;
	break;
      case 'M':
	j = 1;
	break;
      case 'D':
	j = 2;
	break;
      case 'h':
	j = 3;
	break;
      case 'm':
	j = 4;
	break;
      case 's':
	j = 5;
	break;
      default:
	error(ERR_INVALID_UTC_FORMAT);
      };
      if(j != -1) {
	if(UTCindex[j] != -1)
	  error(ERR_INVALID_UTC_FORMAT);
	if(nscanned == 6)
	  error(ERR_INVALID_UTC_FORMAT);
	UTCformat[i] = 'f';
	UTCindex[j] = nscanned;
	nscanned++;
      }
    }
    i++;
  }
  if(nscanned <= 0)
    error(ERR_INVALID_UTC_FORMAT);
  return;
}

void convertUTCtoJD(char *inputUTC, char *UTCformat, int *UTCindex, double *outJD)
{
  int yr, mo, day, hr, min, A, B, i;
  double jd;
  float sec, flt;

  float arg[6];

  for(i=0; i < 6; i++)
    arg[i] = 0.;

  sscanf(inputUTC,UTCformat,&arg[0],&arg[1],&arg[2],&arg[3],
	 &arg[4],&arg[5],&arg[6]);

  if(UTCindex[0] != -1) {
    flt = arg[UTCindex[0]];
    if(flt - floor(flt) > 0.5) {
      yr = (int) ceil(flt);
    } else {
      yr = (int) floor(flt);
    }
  }
  if(UTCindex[1] != -1) {
    flt = arg[UTCindex[1]];
    if(flt - floor(flt) > 0.5) {
      mo = (int) ceil(flt);
    } else {
      mo = (int) floor(flt);
    }
  }
  if(UTCindex[2] != -1) {
    flt = arg[UTCindex[2]];
    if(flt - floor(flt) > 0.5) {
      day = (int) ceil(flt);
    } else {
      day = (int) floor(flt);
    }
  }
  if(UTCindex[3] != -1) {
    flt = arg[UTCindex[3]];
    if(flt - floor(flt) > 0.5) {
      hr = (int) ceil(flt);
    } else {
      hr = (int) floor(flt);
    }
  }
  if(UTCindex[4] != -1) {
    flt = arg[UTCindex[4]];
    if(flt - floor(flt) > 0.5) {
      min = (int) ceil(flt);
    } else {
      min = (int) floor(flt);
    }
  }
  if(UTCindex[5] != -1) {
    sec = arg[UTCindex[5]];
  }


  if(mo <= 2) {
    yr -= 1;
    mo += 12;
  }
  A = floor(yr / 100.);
  B = 2 - A + floor((A / 4.));
  jd = floor(365.25 * yr) + floor(30.6001 * (mo + 1)) + day + 1720994.5;
  jd = jd + (hr/24.) + (min/1440.) + (sec/86400.);
  if(yr > 1583) jd += B;
  *outJD = jd;
}

void printtostring_indentwrap(OutText *text, const char *stoadd, int Ntab_indent)
{
  int l, lold, j, k, k1, k2, space_indent, m;
  while(Ntab_indent*TAB_SPACE_SIZE >= LINEWRAP_LENGTH-1) {
    Ntab_indent--;
  }
  space_indent = Ntab_indent*TAB_SPACE_SIZE;
  l = strlen(stoadd);
  lold = text->len_s;
  if(!text->space)
    {
      text->Nchar_cur_line = 0;
      text->space = MAXLEN;
      if((text->s = malloc(MAXLEN)) == NULL)
	error(ERR_MEMALLOC);
      text->s[0] = '\0';
    }

  k = 0;
  j = lold;
  do {
    k1 = k;
    while(stoadd[k] != '\0' && stoadd[k] != '\t'
	  && stoadd[k] != '\n' && stoadd[k] != ' ')
      k++;
    while(text->len_s + (k-k1) + 2 >= text->space) {
      text->space = text->space * 2;
      if((text->s = realloc(text->s, text->space)) == NULL)
	error(ERR_MEMALLOC);
    }
    if(text->Nchar_cur_line + (k - k1) > LINEWRAP_LENGTH) {
      text->s[j] = '\n';
      j++;
      while(text->len_s + (k-k1) + 2 + space_indent >= text->space) {
	text->space = text->space * 2;
	if((text->s = realloc(text->s, text->space)) == NULL)
	  error(ERR_MEMALLOC);
      }
      for(m=0; m < space_indent; m++) {
	text->s[j] = ' ';
	j++;
	text->len_s += 1;
      }
      text->Nchar_cur_line = space_indent;
      text->s[j] = '\0';
      text->len_s += 1;
    }
    while(k1 < k) {
      if(stoadd[k1] == '~')
	text->s[j] = ' ';
      else
	text->s[j] = stoadd[k1];
      text->Nchar_cur_line += 1;
      j++;
      k1++;
    }
    text->s[j] = '\0';
    text->len_s = j;
    while(stoadd[k] == '\t' || stoadd[k] == ' ' || stoadd[k] == '\n') {
      if(stoadd[k] == '\t') {
	while(text->len_s + TAB_SPACE_SIZE + 1 >= text->space) {
	  text->space = text->space * 2;
	  if((text->s = realloc(text->s, text->space)) == NULL)
	    error(ERR_MEMALLOC);
	}
	if(text->Nchar_cur_line + TAB_SPACE_SIZE >= LINEWRAP_LENGTH) {
	  text->s[j] = '\n';
	  j++;
	  while(text->len_s + space_indent + 1 >= text->space) {
	    text->space = text->space * 2;
	    if((text->s = realloc(text->s, text->space)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(m=0; m < space_indent; m++) {
	    text->s[j] = ' ';
	    j++;
	    text->len_s += 1;
	  }
	  text->s[j] = '\0';
	  text->len_s += 1;
	  text->Nchar_cur_line = space_indent;
	}
	else {
	  for(k2 = 0; k2 < TAB_SPACE_SIZE; k2++) {
	    text->s[j] = ' ';
	    j++;
	    text->s[j] = '\0';
	    text->len_s += 1;
	    text->Nchar_cur_line += 1;
	  }
	}
	k++;
      }
      else if(stoadd[k] == ' ') {
	while(text->len_s + 2 >= text->space) {
	  text->space = text->space * 2;
	  if((text->s = realloc(text->s, text->space)) == NULL)
	    error(ERR_MEMALLOC);
	}
	if(text->Nchar_cur_line + 1 >= LINEWRAP_LENGTH) {
	  text->s[j] = '\n';
	  j++;
	  while(text->len_s + space_indent + 1 >= text->space) {
	    text->space = text->space * 2;
	    if((text->s = realloc(text->s, text->space)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(m=0; m < space_indent; m++) {
	    text->s[j] = ' ';
	    j++;
	    text->len_s += 1;
	  }
	  text->s[j] = '\0';
	  text->len_s += 1;
	  text->Nchar_cur_line = space_indent;
	}
	if(text->Nchar_cur_line != space_indent) {
	  text->s[j] = ' ';
	  j++;
	  text->s[j] = '\0';
	  text->len_s += 1;
	  text->Nchar_cur_line += 1;
	}
	k++;
      }
      else if(stoadd[k] == '\n') {
	while(text->len_s + 1 >= text->space) {
	  text->space = text->space * 2;
	  if((text->s = realloc(text->s, text->space)) == NULL)
	    error(ERR_MEMALLOC);
	}
	text->s[j] = '\n';
	j++;
	while(text->len_s + space_indent + 1 >= text->space) {
	  text->space = text->space * 2;
	  if((text->s = realloc(text->s, text->space)) == NULL)
	    error(ERR_MEMALLOC);
	}
	for(m=0; m < space_indent; m++) {
	  text->s[j] = ' ';
	  j++;
	  text->len_s += 1;
	}
	text->s[j] = '\0';
	text->len_s += 1;
	text->Nchar_cur_line = space_indent;
	k++;
      }
    }
  } while(stoadd[k] != '\0');
}


void printtostring_nowrap(OutText *text, const char *stoadd)
{
  int l, lold, j, k, k1, k2;
  l = strlen(stoadd);
  lold = text->len_s;
  if(!text->space)
    {
      text->Nchar_cur_line = 0;
      text->space = MAXLEN;
      if((text->s = malloc(MAXLEN)) == NULL)
	error(ERR_MEMALLOC);
      text->s[0] = '\0';
    }

  k = 0;
  j = lold;
  do {
    k1 = k;
    while(stoadd[k] != '\0' && stoadd[k] != '\t'
	  && stoadd[k] != '\n' && stoadd[k] != ' ')
      k++;
    while(text->len_s + (k-k1) + 2 >= text->space) {
      text->space = text->space * 2;
      if((text->s = realloc(text->s, text->space)) == NULL)
	error(ERR_MEMALLOC);
    }
    while(k1 < k) {
      if(stoadd[k1] == '~')
	text->s[j] = ' ';
      else
	text->s[j] = stoadd[k1];
      text->Nchar_cur_line += 1;
      j++;
      k1++;
    }
    text->s[j] = '\0';
    text->len_s = j;
    while(stoadd[k] == '\t' || stoadd[k] == ' ' || stoadd[k] == '\n') {
      if(stoadd[k] == '\t') {
	while(text->len_s + TAB_SPACE_SIZE + 1 >= text->space) {
	  text->space = text->space * 2;
	  if((text->s = realloc(text->s, text->space)) == NULL)
	    error(ERR_MEMALLOC);
	}
	for(k2 = 0; k2 < TAB_SPACE_SIZE; k2++) {
	  text->s[j] = ' ';
	  j++;
	  text->s[j] = '\0';
	  text->len_s += 1;
	  text->Nchar_cur_line += 1;
	}
	k++;
      }
      else if(stoadd[k] == ' ') {
	while(text->len_s + 2 >= text->space) {
	  text->space = text->space * 2;
	  if((text->s = realloc(text->s, text->space)) == NULL)
	    error(ERR_MEMALLOC);
	}
	if(text->Nchar_cur_line != 0) {
	  text->s[j] = ' ';
	  j++;
	  text->s[j] = '\0';
	  text->len_s += 1;
	  text->Nchar_cur_line += 1;
	}
	k++;
      }
      else if(stoadd[k] == '\n') {
	while(text->len_s + 1 >= text->space) {
	  text->space = text->space * 2;
	  if((text->s = realloc(text->s, text->space)) == NULL)
	    error(ERR_MEMALLOC);
	}
	text->s[j] = '\n';
	j++;
	text->s[j] = '\0';
	text->len_s += 1;
	text->Nchar_cur_line = 0;
	k++;
      }
    }
  } while(stoadd[k] != '\0');
}

void printtostring_nowrap_tabindent(OutText *text, const char *stoadd)
{
  int l, lold, j, k, k1, k2;
  l = strlen(stoadd);
  lold = text->len_s;
  if(!text->space)
    {
      text->Nchar_cur_line = 0;
      text->space = MAXLEN;
      if((text->s = malloc(MAXLEN)) == NULL)
	error(ERR_MEMALLOC);
      text->s[0] = '\0';
    }

  k = 0;
  j = lold;
  k1 = k;
  while(stoadd[k] != '\0') {
    if(stoadd[k] != '\n')
      k++;
    else
      k += 2;
  }
  while(text->len_s + (k-k1) + 2 >= text->space) {
    text->space = text->space * 2;
    if((text->s = realloc(text->s, text->space)) == NULL)
      error(ERR_MEMALLOC);
  }
  while(k1 < k) {
    text->s[j] = stoadd[k1];
    text->Nchar_cur_line += 1;
    if(stoadd[k1] == '\n') {
      j++;
      text->s[j] = '\t';
      text->Nchar_cur_line = 1;
    }
    j++;
    k1++;
  }
  text->s[j] = '\0';
  text->len_s = j;
}


void InitOutTextStruct(OutText *text)
{
  text->s = NULL;
  text->space = 0;
  text->len_s = 0;
  text->Nchar_cur_line = 0;
}

void printtostring(OutText *text, const char *stoadd)
{
  int l, lold, j, k, k1, k2;
  l = strlen(stoadd);
  lold = text->len_s;
  if(!text->space)
    {
      text->Nchar_cur_line = 0;
      text->space = MAXLEN;
      if((text->s = malloc(MAXLEN)) == NULL)
	error(ERR_MEMALLOC);
      text->s[0] = '\0';
    }

  k = 0;
  j = lold;
  do {
    k1 = k;
    while(stoadd[k] != '\0' && stoadd[k] != '\t'
	  && stoadd[k] != '\n' && stoadd[k] != ' ')
      k++;
    while(text->len_s + (k-k1) + 2 >= text->space) {
      text->space = text->space * 2;
      if((text->s = realloc(text->s, text->space)) == NULL)
	error(ERR_MEMALLOC);
    }
    if(text->Nchar_cur_line + (k - k1) > LINEWRAP_LENGTH) {
      text->s[j] = '\n';
      j++;
      text->s[j] = '\0';
      text->len_s += 1;
      text->Nchar_cur_line = 0;
    }
    while(k1 < k) {
      if(stoadd[k1] == '~')
	text->s[j] = ' ';
      else
	text->s[j] = stoadd[k1];
      text->Nchar_cur_line += 1;
      j++;
      k1++;
    }
    text->s[j] = '\0';
    text->len_s = j;
    while(stoadd[k] == '\t' || stoadd[k] == ' ' || stoadd[k] == '\n') {
      if(stoadd[k] == '\t') {
	while(text->len_s + TAB_SPACE_SIZE + 1 >= text->space) {
	  text->space = text->space * 2;
	  if((text->s = realloc(text->s, text->space)) == NULL)
	    error(ERR_MEMALLOC);
	}
	if(text->Nchar_cur_line + TAB_SPACE_SIZE >= LINEWRAP_LENGTH) {
	  text->s[j] = '\n';
	  j++;
	  text->s[j] = '\0';
	  text->len_s += 1;
	  text->Nchar_cur_line = 0;
	}
	else {
	  for(k2 = 0; k2 < TAB_SPACE_SIZE; k2++) {
	    text->s[j] = ' ';
	    j++;
	    text->s[j] = '\0';
	    text->len_s += 1;
	    text->Nchar_cur_line += 1;
	  }
	}
	k++;
      }
      else if(stoadd[k] == ' ') {
	while(text->len_s + 2 >= text->space) {
	  text->space = text->space * 2;
	  if((text->s = realloc(text->s, text->space)) == NULL)
	    error(ERR_MEMALLOC);
	}
	if(text->Nchar_cur_line + 1 >= LINEWRAP_LENGTH) {
	  text->s[j] = '\n';
	  j++;
	  text->s[j] = '\0';
	  text->len_s += 1;
	  text->Nchar_cur_line = 0;
	}
	if(text->Nchar_cur_line != 0) {
	  text->s[j] = ' ';
	  j++;
	  text->s[j] = '\0';
	  text->len_s += 1;
	  text->Nchar_cur_line += 1;
	}
	k++;
      }
      else if(stoadd[k] == '\n') {
	while(text->len_s + 1 >= text->space) {
	  text->space = text->space * 2;
	  if((text->s = realloc(text->s, text->space)) == NULL)
	    error(ERR_MEMALLOC);
	}
	text->s[j] = '\n';
	j++;
	text->s[j] = '\0';
	text->len_s += 1;
	text->Nchar_cur_line = 0;
	k++;
      }
    }
  } while(stoadd[k] != '\0');
}



void phaselc(int N, double *t, double *mag, double *sig, double period, int is_T0_given, double T0, char *phasevarname, _Variable *phasevar, int threadid, double startphase)
{
  int i;
  double ph, t0;
  double stopphase = startphase + 1.;
  if(period > 0.)
    {
      if(is_T0_given)
	t0 = T0;
      else
	t0 = t[0];
      for(i=0;i<N;i++)
	{
	  ph = (t[i]-t0)/period;
	  if(ph >= 0) {
	    ph -= (int) ph;
	  }
	  else {
	    ph += (double) ( -((int) ph) + 1);
	  }
	  if(ph > stopphase) {
	    ph = ph - ceil(ph - stopphase);
	  } else if(ph < startphase) {
	    ph = ph + ceil(startphase - ph);
	  }
	  if(phasevarname == NULL)
	    t[i] = ph;
	  else
	    (*((double ***) phasevar->dataptr))[threadid][i] = ph;
	}
      //mysort3(N,t,mag,sig);
    } 
  else if(phasevarname != NULL)
    {
      for(i=0;i<N;i++) {
	(*((double ***) phasevar->dataptr))[threadid][i] = t[i];
      }
    }
}

void spline(double *x,double *y,int n,double yp1,double ypn,
            double *y2,double *u){

  int i,k;
  double p,qn,sig,un;

  if (yp1 > 0.99e30){
    y2[0]=u[0]=0.0;
  }else {
    y2[0] = -0.5;
    u[0]=(3.0/(x[1]-x[0]))*((y[1]-y[0])/(x[1]-x[0])-yp1);
  }

  for (i=1;i<=n-2;i++) {
    sig=(x[i]-x[i-1])/(x[i+1]-x[i-1]);
    p=sig*y2[i-1]+2.0;
    y2[i]=(sig-1.0)/p;
    u[i]=(y[i+1]-y[i])/(x[i+1]-x[i]) - (y[i]-y[i-1])/(x[i]-x[i-1]);
    u[i]=(6.0*u[i]/(x[i+1]-x[i-1])-sig*u[i-1])/p;
  }

  if (ypn > 0.99e30){
    qn=un=0.0;
  }else{
    qn=0.5;
    un=(3.0/(x[n-1]-x[n-2]))*(ypn-(y[n-1]-y[n-2])/(x[n-1]-x[n-2]));
  }

  y2[n-1]=(un-qn*u[n-2])/(qn*y2[n-2]+1.0);

  for (k=n-2;k>=0;k--){
    y2[k]=y2[k]*y2[k+1]+u[k];
  }
  return;
}

void splint(double *xa,double *ya,double *y2a,int n,double x,double *y)
{
	int klo,khi,k;
	double h,b,a;

	/* Check to see if the point needs to be extrapolated */
	if(x > xa[n-1])
	  {
	    h = xa[n-1] - xa[n-2];
	    (*y) = ya[n-1] + (x - xa[n-1])*(((ya[n-1]-ya[n-2])/h) + (2.0*y2a[n-1]+y2a[n-1])*(h/6.0));
	    return;
	  }
	else if(x < xa[0])
	  {
	    h = xa[1] - xa[0];
	    (*y) = ya[0] + (x - xa[0])*(((ya[1]-ya[0])/h) + (-2.0*y2a[0] - y2a[1])*(h/6.0));
	    return;
	  }

	klo=0;
	khi=n-1;


	while (khi-klo > 1) {
		k=(khi+klo) >> 1;
		if (xa[k] > x) khi=k;
		else klo=k;}



	h=xa[khi]-xa[klo];
	if (h == 0.0) printf("Bad XA input to routine SPLINT\n");
	a=(xa[khi]-x)/h;
	b=(x-xa[klo])/h;


	*y=a*ya[klo]+b*ya[khi]+((a*a*a-a)*y2a[klo]+(b*b*b-b)*y2a[khi])*(h*h)/6.0;

  return;
}

#ifndef MIN
#define MIN(A,B) ((A) < (B) ? (A) : (B))
#endif

void spline_monotonic(int N, double *x, double *y, double *yprime)
{
  /* This function is an implementation of the C1 monotonic cubic polynomial interpolation proposed by Steffen, 1990 A&A 239, 443. This procedure calculates the slopes to use. */
  int i,k;
  double p,qn,sig,un, si0, si1, h0, h1, a;

  yprime[0] = 0.0;
  yprime[N-1] = 0.0;

  if(N < 3)
    return;

  si1 = (y[1] - y[0])/(x[1] - x[0]);
  for(i=1;i<=N-2;i++) {
    si0 = si1;
    si1 = (y[i+1] - y[i])/(x[i+1] - x[i]);

    if(si0 * si1 <= 0.)
      yprime[i] = 0.;
    else {
      h0 = x[i] - x[i-1];
      h1 = x[i+1] - x[i];
      p = (si0*h1 + si1*h0)/(h0 + h1);
      if (fabs(p) > 2.*fabs(si0) || fabs(p) > 2.*fabs(si1))
	{
	  a = (si0 > 0. ? 1. : -1.);
	  yprime[i] = 2.0*a*MIN(fabs(si0),fabs(si1));
	}
      else
	yprime[i] = p;
    }
  }
}

double splint_monotonic(int N, double *x, double *y, double *yprime, double xt)
{
  int klo,khi,k;
  double h,b,a,c,d,delx,s,retval,del;

  /* Check to see if the point needs to be extrapolated */
  if(xt > x[N-1])
    {
      klo = N-2;
      khi = N-1;
    }
  else if(xt < x[0])
    {
      klo = 0;
      khi = 1;
    }

  else
    {
      klo=0;
      khi=N-1;


      while (khi-klo > 1) {
	k=(khi+klo) >> 1;
	if (x[k] > xt) khi=k;
	else klo=k;}
    }


  h=x[khi]-x[klo];
  if (h == 0.0) printf("Bad XA input to routine SPLINT\n");
  s = (y[khi] - y[klo])/h;

  a = (yprime[klo] + yprime[khi] - 2.0*s)/h/h;
  b = (3.0*s - 2.0*yprime[klo] - yprime[khi])/h;
  c = yprime[klo];
  d = y[klo];
  del = (xt - x[klo]);

  retval = a*del*del*del + b*del*del + c*del + d;

  return retval;
}

void SwapLightCurvePoints(int i, int j, ProgramData *p, int lc)
{
  double ***dblptr;
  double ****dbl2ptr;
  short ***shortptr;
  short ****short2ptr;
  int ***intptr;
  int ****int2ptr;
  char ***charptr;
  char ****char2ptr;
  char ****stringptr;
  char *****string2ptr;
  float ***floatptr;
  float ****float2ptr;
  long ***longptr;
  long ****long2ptr;
  int k, u, Nc;
  _DataFromLightCurve *d;
  char *ts;
  int sizemax = 0;
  double td;
  short tsh;
  int ti;
  char tc;
  long tl;
  float tf;

  for(k=0; k < p->NDataFromLightCurve; k++) {
    d = &(p->DataFromLightCurve[k]);
    Nc = d->Ncolumns;
    if(Nc == 0) {
      switch(d->datatype)
	{
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = (double ***) d->dataptr;
	  td = (*dblptr)[lc][i];
	  (*dblptr)[lc][i] = (*dblptr)[lc][j];
	  (*dblptr)[lc][j] = td;
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  dblptr = (double ***) d->dataptr;
	  td = (*dblptr)[lc][i];
	  (*dblptr)[lc][i] = (*dblptr)[lc][j];
	  (*dblptr)[lc][j] = td;
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = (char ****) d->dataptr;
	  if(d->maxstringlength > sizemax) {
	    if(!sizemax) {
	      sizemax = d->maxstringlength;
	      if((ts = (char *) malloc(sizemax)) == NULL)
		error(ERR_MEMALLOC);
	    }
	    else {
	      sizemax = d->maxstringlength;
	      if((ts = (char *) realloc(ts,sizemax)) == NULL)
		error(ERR_MEMALLOC);
	    }
	  }
	  memcpy(ts, ((*stringptr)[lc][i]), d->maxstringlength);
	  memcpy(((*stringptr)[lc][i]),((*stringptr)[lc][j]),
		 d->maxstringlength);
	  memcpy(((*stringptr)[lc][j]),ts,d->maxstringlength);
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = (int ***) d->dataptr;
	  ti = (*intptr)[lc][i];
	  (*intptr)[lc][i] = (*intptr)[lc][j];
	  (*intptr)[lc][j] = ti;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = (short ***) d->dataptr;
	  tsh = (*shortptr)[lc][i];
	  (*shortptr)[lc][i] = (*shortptr)[lc][j];
	  (*shortptr)[lc][j] = tsh;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = (float ***) d->dataptr;
	  tf = (*floatptr)[lc][i];
	  (*floatptr)[lc][i] = (*floatptr)[lc][j];
	  (*floatptr)[lc][j] = tf;
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = (long ***) d->dataptr;
	  tl = (*longptr)[lc][i];
	  (*longptr)[lc][i] = (*longptr)[lc][j];
	  (*longptr)[lc][j] = tl;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = (char ***) d->dataptr;
	  tc = (*charptr)[lc][i];
	  (*charptr)[lc][i] = (*charptr)[lc][j];
	  (*charptr)[lc][j] = tc;
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
    }
    else {
      for(u=0; u < Nc; u++) {
	switch(d->datatype)
	  {
	  case VARTOOLS_TYPE_DOUBLE:
	    dbl2ptr = (double ****) d->dataptr;
	    td = (*dbl2ptr)[lc][u][i];
	    (*dbl2ptr)[lc][u][i] = (*dbl2ptr)[lc][u][j];
	    (*dbl2ptr)[lc][u][j] = td;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    dbl2ptr = (double ****) d->dataptr;
	    td = (*dbl2ptr)[lc][u][i];
	    (*dbl2ptr)[lc][u][i] = (*dbl2ptr)[lc][u][j];
	    (*dbl2ptr)[lc][u][j] = td;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    string2ptr = (char *****) d->dataptr;
	    if(d->maxstringlength > sizemax) {
	      if(!sizemax) {
		sizemax = d->maxstringlength;
		if((ts = (char *) malloc(sizemax)) == NULL)
		  error(ERR_MEMALLOC);
	      }
	      else {
		sizemax = d->maxstringlength;
		if((ts = (char *) realloc(ts,sizemax)) == NULL)
		  error(ERR_MEMALLOC);
	      }
	    }
	    memcpy(ts, ((*string2ptr)[lc][u][i]), d->maxstringlength);
	    memcpy(((*string2ptr)[lc][u][i]),((*string2ptr)[lc][u][j]),
		   d->maxstringlength);
	    memcpy(((*string2ptr)[lc][u][j]),ts,d->maxstringlength);
	    break;
	  case VARTOOLS_TYPE_INT:
	    int2ptr = (int ****) d->dataptr;
	    ti = (*int2ptr)[lc][u][i];
	    (*int2ptr)[lc][u][i] = (*int2ptr)[lc][u][j];
	    (*int2ptr)[lc][u][j] = ti;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    short2ptr = (short ****) d->dataptr;
	    tsh = (*short2ptr)[lc][u][i];
	    (*short2ptr)[lc][u][i] = (*short2ptr)[lc][u][j];
	    (*short2ptr)[lc][u][j] = tsh;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    float2ptr = (float ****) d->dataptr;
	    tf = (*float2ptr)[lc][u][i];
	    (*float2ptr)[lc][u][i] = (*float2ptr)[lc][u][j];
	    (*float2ptr)[lc][u][j] = tf;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    long2ptr = (long ****) d->dataptr;
	    tl = (*long2ptr)[lc][u][i];
	    (*long2ptr)[lc][u][i] = (*long2ptr)[lc][u][j];
	    (*long2ptr)[lc][u][j] = tl;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    char2ptr = (char ****) d->dataptr;
	    tc = (*char2ptr)[lc][u][i];
	    (*char2ptr)[lc][u][i] = (*char2ptr)[lc][u][j];
	    (*char2ptr)[lc][u][j] = tc;
	    break;
	  default:
	    error(ERR_BADTYPE);
	    break;
	  }
      }
    }
  }
  if(sizemax > 0)
    free(ts);
}

void CopyLightCurvePoints(int i, int j, ProgramData *p, int lc)
/* In this case all of the light curve data at time index i is copied to
   time index j. No swapping is done. */
{
  double ***dblptr;
  double ****dbl2ptr;
  short ***shortptr;
  short ****short2ptr;
  int ***intptr;
  int ****int2ptr;
  char ***charptr;
  char ****char2ptr;
  char ****stringptr;
  char *****string2ptr;
  float ***floatptr;
  float ****float2ptr;
  long ***longptr;
  long ****long2ptr;
  int k, u, Nc;
  _DataFromLightCurve *d;

  for(k=0; k < p->NDataFromLightCurve; k++) {
    d = &(p->DataFromLightCurve[k]);
    Nc = d->Ncolumns;
    if(Nc == 0) {
      switch(d->datatype)
	{
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = (double ***) d->dataptr;
	  (*dblptr)[lc][j] = (*dblptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  dblptr = (double ***) d->dataptr;
	  (*dblptr)[lc][j] = (*dblptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = (char ****) d->dataptr;
	  memcpy(((*stringptr)[lc][j]),((*stringptr)[lc][i]),
		 d->maxstringlength);
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = (int ***) d->dataptr;
	  (*intptr)[lc][j] = (*intptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = (short ***) d->dataptr;
	  (*shortptr)[lc][j] = (*shortptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = (float ***) d->dataptr;
	  (*floatptr)[lc][j] = (*floatptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = (long ***) d->dataptr;
	  (*longptr)[lc][j] = (*longptr)[lc][i];
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = (char ***) d->dataptr;
	  (*charptr)[lc][j] = (*charptr)[lc][i];
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
    }
    else {
      for(u=0; u < Nc; u++) {
	switch(d->datatype)
	  {
	  case VARTOOLS_TYPE_DOUBLE:
	    dbl2ptr = (double ****) d->dataptr;
	    (*dbl2ptr)[lc][u][j] = (*dbl2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    dbl2ptr = (double ****) d->dataptr;
	    (*dbl2ptr)[lc][u][j] = (*dbl2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_STRING:
	    string2ptr = (char *****) d->dataptr;
	    memcpy(((*string2ptr)[lc][u][j]),((*string2ptr)[lc][u][i]),
		   d->maxstringlength);
	    break;
	  case VARTOOLS_TYPE_INT:
	    int2ptr = (int ****) d->dataptr;
	    (*int2ptr)[lc][u][j] = (*int2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    short2ptr = (short ****) d->dataptr;
	    (*short2ptr)[lc][u][j] = (*short2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    float2ptr = (float ****) d->dataptr;
	    (*float2ptr)[lc][u][j] = (*float2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_LONG:
	    long2ptr = (long ****) d->dataptr;
	    (*long2ptr)[lc][u][j] = (*long2ptr)[lc][u][i];
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    char2ptr = (char ****) d->dataptr;
	    (*char2ptr)[lc][u][j] = (*char2ptr)[lc][u][i];
	    break;
	  default:
	    error(ERR_BADTYPE);
	    break;
	  }
      }
    }
  }
}

void AverageLightCurvePoints(int i, int j, int h, ProgramData *p, int lc)
/* Take all of the light curve data between time indices i and j, average them
   and store the result in time index h. */
{
  double ***dblptr;
  double ****dbl2ptr;
  short ***shortptr;
  short ****short2ptr;
  int ***intptr;
  int ****int2ptr;
  char ***charptr;
  char ****char2ptr;
  char ****stringptr;
  char *****string2ptr;
  float ***floatptr;
  float ****float2ptr;
  long ***longptr;
  long ****long2ptr;
  int k, u, Nc, Nbin;
  _DataFromLightCurve *d;
  double td;
  short tsh;
  int ti;
  char tc;
  long tl;
  float tf;

  /* Only do something if j is larger than i */
  if(j <= i) return;

  Nbin = (j-i)+1;

  for(k=0; k < p->NDataFromLightCurve; k++) {
    d = &(p->DataFromLightCurve[k]);
    Nc = d->Ncolumns;
    if(Nc == 0) {
      switch(d->datatype)
	{
	case VARTOOLS_TYPE_DOUBLE:
	  dblptr = (double ***) d->dataptr;
	  td = 0.0;
	  /* Handle the error vector differently */
	  if((*dblptr)[lc] != p->sig[lc]) {
	    for(k=i; k <= j; k++) {
	      td += (*dblptr)[lc][k];
	    }
	    td = td/Nbin;
	  } else {
	    for(k=i; k <= j; k++) {
	      td += ((*dblptr)[lc][k])*((*dblptr)[lc][k]);
	    }
	    td = sqrt(td/Nbin);
	  }
	  (*dblptr)[lc][h] = td;
	  break;
	case VARTOOLS_TYPE_CONVERTJD:
	  dblptr = (double ***) d->dataptr;
	  td = 0.0;
	  /* Handle the error vector differently */
	  if((*dblptr)[lc] != p->sig[lc]) {
	    for(k=i; k <= j; k++) {
	      td += (*dblptr)[lc][k];
	    }
	    td = td/Nbin;
	  } else {
	    for(k=i; k <= j; k++) {
	      td += ((*dblptr)[lc][k])*((*dblptr)[lc][k]);
	    }
	    td = sqrt(td/Nbin);
	  }
	  (*dblptr)[lc][h] = td;
	  break;
	case VARTOOLS_TYPE_STRING:
	  stringptr = (char ****) d->dataptr;
	  /* For strings, just adopt the middle value */
	  k = (i+j)/2;
	  if(h != k)
	    memcpy(((*stringptr)[lc][h]),((*stringptr)[lc][k]),
		   d->maxstringlength);
	  break;
	case VARTOOLS_TYPE_INT:
	  intptr = (int ***) d->dataptr;
	  ti = 0;
	  for(k=i; k <= j; k++) {
	    ti += (*intptr)[lc][k];
	  }
	  ti = ti/Nbin;
	  (*intptr)[lc][h] = ti;
	  break;
	case VARTOOLS_TYPE_SHORT:
	  shortptr = (short ***) d->dataptr;
	  tsh = 0;
	  for(k=i; k <= j; k++) {
	    tsh += (*shortptr)[lc][k];
	  }
	  tsh = tsh/Nbin;
	  (*shortptr)[lc][h] = tsh;
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  floatptr = (float ***) d->dataptr;
	  tf = 0.;
	  for(k=i; k <= j; k++) {
	    tf += (*floatptr)[lc][k];
	  }
	  tf = tf/Nbin;
	  (*floatptr)[lc][h] = tf;
	  break;
	case VARTOOLS_TYPE_LONG:
	  longptr = (long ***) d->dataptr;
	  tl = 0;
	  for(k=i; k <= j; k++) {
	    tl += (*longptr)[lc][k];
	  }
	  (*longptr)[lc][h] = tl;
	  break;
	case VARTOOLS_TYPE_CHAR:
	  charptr = (char ***) d->dataptr;
	  /* For chars just adopt the middle value */
	  k = (i+j)/2;
	  if(h != k)
	    (*charptr)[lc][h] = (*charptr)[lc][k];
	  break;
	default:
	  error(ERR_BADTYPE);
	  break;
	}
    }
    else {
      for(u=0; u < Nc; u++) {
	switch(d->datatype)
	  {
	  case VARTOOLS_TYPE_DOUBLE:
	    dbl2ptr = (double ****) d->dataptr;
	    /* Handle the error vector differently */
	    if((*dbl2ptr)[lc][u] != p->sig[lc]) {
	      for(k=i; k <= j; k++) {
		td += (*dbl2ptr)[lc][u][k];
	      }
	      td = td/Nbin;
	    } else {
	      for(k=i; k <= j; k++) {
		td += ((*dbl2ptr)[lc][u][k])*((*dbl2ptr)[lc][u][k]);
	      }
	      td = sqrt(td/Nbin);
	    }
	    (*dbl2ptr)[lc][u][h] = td;
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    dbl2ptr = (double ****) d->dataptr;
	    td = 0.0;
	    /* Handle the error vector differently */
	    if((*dbl2ptr)[lc][u] != p->sig[lc]) {
	      for(k=i; k <= j; k++) {
		td += (*dbl2ptr)[lc][u][k];
	      }
	      td = td/Nbin;
	    } else {
	      for(k=i; k <= j; k++) {
		td += ((*dbl2ptr)[lc][u][k])*((*dbl2ptr)[lc][u][k]);
	      }
	      td = sqrt(td/Nbin);
	    }
	    (*dbl2ptr)[lc][u][h] = td;
	    break;
	  case VARTOOLS_TYPE_STRING:
	    string2ptr = (char *****) d->dataptr;
	    /* For strings, just adopt the middle value */
	    k = (i+j)/2;
	    if(h != k)
	      memcpy(((*string2ptr)[lc][u][h]),((*string2ptr)[lc][u][k]),
		     d->maxstringlength);
	    break;
	  case VARTOOLS_TYPE_INT:
	    int2ptr = (int ****) d->dataptr;
	    ti = 0.0;
	    for(k=i; k <= j; k++) {
	      ti += (*int2ptr)[lc][u][k];
	    }
	    ti = ti/Nbin;
	    (*int2ptr)[lc][u][h] = ti;
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    short2ptr = (short ****) d->dataptr;
	    tsh = 0.0;
	    for(k=i; k <= j; k++) {
	      tsh += (*short2ptr)[lc][u][k];
	    }
	    tsh = tsh/Nbin;
	    (*short2ptr)[lc][u][h] = tsh;
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    float2ptr = (float ****) d->dataptr;
	    tf = 0.0;
	    for(k=i; k <= j; k++) {
	      tf += (*float2ptr)[lc][u][k];
	    }
	    tf = tf/Nbin;
	    (*float2ptr)[lc][u][h] = tf;
	    break;
	  case VARTOOLS_TYPE_LONG:
	    long2ptr = (long ****) d->dataptr;
	    tl = 0;
	    for(k=i; k <= j; k++) {
	      tl += (*long2ptr)[lc][u][k];
	    }
	    tl = tl/Nbin;
	    (*long2ptr)[lc][u][h] = tl;
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    char2ptr = (char ****) d->dataptr;
	    /* For chars, just adopt the middle value */
	    k = (i+j)/2;
	    if(i != k)
	      (*char2ptr)[lc][u][h] = (*char2ptr)[lc][u][k];
	    break;
	  default:
	    error(ERR_BADTYPE);
	    break;
	  }
      }
    }
  }
}

void mergeequallctimes(ProgramData *p, int lc)
/* This function goes through the light curve in thread lc finds
   groups of points with equal times and merges them by averaging all
   associated light curve data. This function assumes the light curves
   are sorted in time. */
{
  int i, j, k;
  if(p->NJD[lc] < 2) return;

  /* First check if there are any equal times */
  for(i=1; i < p->NJD[lc]; i++) {
    if(p->t[lc][i-1] == p->t[lc][i])
      break;
  }

  /* All points are distinct, nothing to do */
  if(i >= p->NJD[lc]) return;
    
  for(i=0,j=0; i < p->NJD[lc]; i++) {
    k=i+1;
    while(k < p->NJD[lc] ? p->t[lc][k] == p->t[lc][i] : 0) k++;
    if(k - i > 1) {
      AverageLightCurvePoints(i, (k-1), j, p, lc);
      i = k-1;
    } else {
      if(i != j) {
	CopyLightCurvePoints(i, j, p, lc);
      }
    }
    j++;
  }
  p->NJD[lc] = j;
}

int sortlcbytime_internal(int i0, int iN, double *t, int lc, ProgramData *p)
{
  double v;
  int i, j, k;

  if(iN - i0 <= 1) return 0;

  v = t[i0];
  i = i0;
  j = iN;
  for(;;)
    {
      while(++i < iN ? t[i] < v : 0) { }
      while(t[--j] > v) {}
      if(i >= j) break;
      SwapLightCurvePoints(i, j, p, lc);
    }
  SwapLightCurvePoints(i0,i-1,p,lc);

  sortlcbytime_internal(i0, i-1, t, lc, p);
  sortlcbytime_internal(i, iN, t, lc, p);

}

int sortlcbytime(int size, double *t, int lc, ProgramData *p)
{
  int i, j, k;
  double v;

  if(size<=1) return 0;

  /* Check if the lc is already sorted, return if it is */
  for(i=1; i < size; i++) {
    if(t[i] < t[i-1])
      break;
  }
  if(i == size)
    return 0;

  v = t[0];
  i = 0;
  j = size;
  for(;;)
    {
      while(++i < size ? t[i] < v : 0) { }
      while(t[--j] > v) {}
      if(i >= j) break;
      SwapLightCurvePoints(i, j, p, lc);
    }
  SwapLightCurvePoints(0, i-1, p, lc);
  sortlcbytime_internal(0, i-1, t, lc, p);
  sortlcbytime_internal(i, size, t, lc, p);

  /* Check for any = time points */
  for(i=1; i < size; i++) {
    if(t[i] == t[i-1])
      return 1;
  }

  return 0;
}

int check_isspecialchar(char c){
  if(c == ' ' ||
     c == '\n' ||
     c == '*' ||
     c == '\t' ||
     c == '\\' ||
     c == '?' ||
     c == '(' ||
     c == ')' ||
     c == '[' ||
     c == ']' ||
     c == '#' ||
     c == '}' ||
     c == '{' ||
     c == '\'' ||
     c == '\"' ||
     c == ';' ||
     c == ':' ||
     c == '&' ||
     c == '`' ||
     c == '!' ||
     c == '$' ||
     c == '>' ||
     c == '<' ||
     c == '|' ||
     c == '%' ||
     c == '~' ||
     c == '^' ||
     c == '@')
    return 1;
  else
    return 0;
}
