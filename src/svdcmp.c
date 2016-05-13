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
/* This file contains the pythag, svdcmp and svdfit routines from Numerical Recipes in C, Press et al., 1992. They have been modified to begin array indexing from 0 rather than from 1 and to use doubles instead of floats. */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define ABS_(val) ((val) > 0 ? (val) : -(val))

#define SIGN_(A,B) ((B) >= 0.0 ? ABS_(A) : (-(ABS_(A))))
#define MAX_(A,B) ((A) > (B) ? (A) : (B))
#define MIN_(A,B) ((A) < (B) ? (A) : (B))
#define TINY 1.0e-20
#define TOL 1.0e-9

double pythag(double a, double b)
{
  double absa,absb;
  absa = ABS_(a);
  absb = ABS_(b);
  if(absa > absb) return absa*sqrt(1.0+((absb*absb)/(absa*absa)));
  else return (absb == 0.0 ? 0.0 : absb*sqrt(1.0 + ((absa*absa)/(absb*absb))));
}

void svdcmp(double **a, int m, int n, double *w, double **v)
{
  int flag, i, its, j, jj, k, l, nm;
  double anorm, c, f, g, h, s, scale, x, y, z, *rv1;

  rv1 = (double *) malloc(n*sizeof(double));
  g=scale=anorm=0.0;
  for(i=0;i<n;i++) {
    l=i+1;
    rv1[i]=scale*g;
    g=s=scale=0.0;
    if (i < m) {
      for(k=i;k<m;k++) scale += ABS_(a[k][i]);
      if(scale) {
	for(k=i;k<m;k++) {
	  a[k][i] /= scale;
	  s += a[k][i]*a[k][i];
	}
	f=a[i][i];
	g = -SIGN_(sqrt(s),f);
	h=f*g-s;
	a[i][i]=f-g;
	for (j=l;j<n;j++) {
	  for(s=0.0,k=i;k<m;k++) s += a[k][i]*a[k][j];
	  f=s/h;
	  for (k=i;k<m;k++) a[k][j] += f*a[k][i];
	}
	for (k=i;k<m;k++) a[k][i] *= scale;
      }
    }
    w[i] = scale * g;
    g = s = scale = 0.0;
    if (i < m && i != n - 1) {
      for(k=l;k<n;k++) scale += ABS_(a[i][k]);
      if(scale) {
	for(k=l;k<n;k++) {
	  a[i][k] /= scale;
	  s += a[i][k]*a[i][k];
	}
	f = a[i][l];
	g = -SIGN_(sqrt(s),f);
	h = f*g-s;
	a[i][l] = f-g;
	for (k=l;k<n;k++) rv1[k] = a[i][k]/h;
	for (j=l;j<m;j++) {
	  for (s=0.0,k=l;k<n;k++) s += a[j][k]*a[i][k];
	  for (k=l;k<n;k++) a[j][k] += s*rv1[k];
	}
	for(k=l;k<n;k++) a[i][k] *= scale;
      }
    }
    anorm=MAX_(anorm,(ABS_(w[i])+ABS_(rv1[i])));
  }
  for(i=n-1;i>=0;i--) {
    if (i < n-1) {
      if(g) {
	for (j=l;j<n;j++)
	  v[j][i]=(a[i][j]/a[i][l])/g;
	for(j=l;j<n;j++) {
	  for(s=0.0,k=l;k<n;k++) s += a[i][k]*v[k][j];
	  for(k=l;k<n;k++) v[k][j] += s*v[k][i];
	}
      }
      for(j=l;j<n;j++) v[i][j] = v[j][i] = 0.0;
    }
    v[i][i] = 1.0;
    g = rv1[i];
    l = i;
  }
  for(i=MIN_(m,n)-1;i>=0;i--) {
    l = i+1;
    g = w[i];
    for (j=l;j < n; j++) a[i][j] = 0.0;
    if (g) {
      g = 1.0/g;
      for (j=l;j<n;j++) {
	for (s=0.0,k=l;k<m;k++) s += a[k][i] * a[k][j];
	f = (s/a[i][i]) * g;
	for (k=i;k<m;k++) a[k][j] += f*a[k][i];
      }
      for (j=i;j<m;j++) a[j][i] *= g;
    } else for (j=i;j<m;j++) a[j][i] = 0.0;
    ++a[i][i];
  }
  for (k=n-1;k>=0;k--) {
    for (its=1;its<=30;its++) {
      flag = 1;
      for (l=k;l>=0;l--) {
	nm=l-1;
	if ((double)(ABS_(rv1[l])+anorm) == anorm) {
	  flag = 0;
	  break;
	}
	if ((double)(ABS_(w[nm])+anorm) == anorm) break;
      }
      if(flag) {
	c=0.0;
	s=1.0;
	for (i=l;i<=k;i++) {
	  f = s*rv1[i];
	  rv1[i] = c*rv1[i];
	  if ((double) (ABS_(f)+anorm) == anorm) break;
	  g = w[i];
	  h = pythag(f,g);
	  w[i] = h;
	  h = 1.0/h;
	  c = g*h;
	  s = -f*h;
	  for (j=0;j<m;j++) {
	    y = a[j][nm];
	    z = a[j][i];
	    a[j][nm] = y*c+z*s;
	    a[j][i] = z*c-y*s;
	  }
	}
      }
      z = w[k];
      if (l == k) {
	if (z < 0.0) {
	  w[k] = -z;
	  for (j=0;j<n;j++) v[j][k] = -v[j][k];
	}
	break;
      }
      if (its == 30)
	{
	  fprintf(stderr,"Error: No convergence in 30 svdcmp iterations\n");
	  exit(5);
	}
      x = w[l];
      nm = k-1;
      y=w[nm];
      g=rv1[nm];
      h=rv1[k];
      f = ((y-z)*(y+z)+(g-h)*(g+h))/(2.0*h*y);
      g = pythag(f,1.0);
      f = ((x-z)*(x+z)+h*((y/(f+SIGN_(g,f)))-h))/x;
      c = s = 1.0;
      for (j=l;j<=nm;j++) {
	i = j+1;
	g = rv1[i];
	y=w[i];
	h = s*g;
	g = c*g;
	z = pythag(f,h);
	rv1[j] = z;
	c = f/z;
	s = h/z;
	f = x*c+g*s;
	g = g*c-x*s;
	h=y*s;
	y *= c;
	for (jj=0;jj<n;jj++) {
	  x = v[jj][j];
	  z = v[jj][i];
	  v[jj][j] = x*c + z*s;
	  v[jj][i] = z*c - x*s;
	}
	z = pythag(f,h);
	w[j] = z;
	if (z) {
	  z = 1.0/z;
	  c = f*z;
	  s = h*z;
	}
	f = c*g + s*y;
	x = c*y - s*g;
	for (jj=0;jj<m;jj++) {
	  y = a[jj][j];
	  z = a[jj][i];
	  a[jj][j] = y*c + z*s;
	  a[jj][i] = z*c - y*s;
	}
      }
      rv1[l] = 0.0;
      rv1[k] = f;
      w[k] = x;
    }
  }
  free(rv1);
}


void svbksb(double **u, double *w, double **v, int m, int n, double *b, double *x)
{
  int jj,j,i;
  double s,*tmp;
  tmp = (double *) malloc(n*sizeof(double));
  for(j=0;j<n;j++) {
    s = 0.0;
    if(w[j]) {
      for(i=0;i<m;i++) {if(!isnan(b[i])) s += u[i][j]*b[i];}
      s /= w[j];
    }
    tmp[j] = s;
  }
  for(j=0;j<n;j++) {
    s = 0.0;
    for(jj=0;jj<n;jj++) s += v[j][jj]*tmp[jj];
    x[j]=s;
  }
  free(tmp);
}
