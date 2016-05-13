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
/* This file contains functions that are used by the -findblends VARTOOLS
   command, by Joel Hartman. 18 May 2009.
*/

#include "commands.h"
#include "programdata.h"
#include "functions.h"

#ifndef DEG2RAD
#define DEG2RAD(A) ((A)*0.01745329252)
#endif
#ifndef RAD2DEG
#define RAD2DEG(A) ((A)*57.29577951)
#endif
#ifndef ABS
#define ABS(A) ((A) >= 0 ? (A) : (-(A)))
#endif
#ifndef SIGNONEARG
#define SIGNONEARG(A) ((A) >= 0 ? (1) : (-1))
#endif
#ifndef dmax
#define dmax(A,B) ((A)>=(B)?(A):(B))
#endif
#ifndef dmin
#define dmin(A,B) ((A)<=(B)?(A):(B))
#endif
#ifndef FLT_EPSILON
#define FLT_EPSILON 0.000001
#endif
#ifndef PI
#define PI 3.1415926535897932384626433832795
#endif
#ifndef PI_2
#define PI_2 1.57079632679490
#endif
#ifndef PI2
#define PI2 6.28318530717958
#endif
#ifndef PI3_2
#define PI3_2 4.71238898038469
#endif

/* Bisection with indices used to sort x given in indx */
int findXindx(double *x, int *indx, double xval, int i1, int N)
{
  int u;
  if(x[indx[i1]] >= xval)
    return(i1);
  if(x[indx[N-1]] < xval)
    return(N);
  while(1)
    {
      u = (N+i1)/2;
      if(u == N || u == i1)
	return(u);
      if(x[indx[u]] >= xval)
	{
	  if(x[indx[u-1]] < xval)
	    return(u);
	  N = u;
	}
      if(x[indx[u]] <= xval)
	{
	  if(u == N-1)
	    return(N);
	  if(x[indx[u+1]] >= xval)
	    return(u+1);
	  i1 = u;
	}
    }
}

double distxy(double x1, double y1, double x2, double y2)
{
  return sqrt((x1 - x2)*(x1 - x2) + (y1 - y2)*(y1 - y2));
}

double distradec(double ra1, double dec1, double ra2, double dec2)
{
  return (cos(DEG2RAD(ra1-ra2))*cos(DEG2RAD(dec1))*cos(DEG2RAD(dec2)) + sin(DEG2RAD(dec1))*sin(DEG2RAD(dec2)));
}

/* inverse arc projection (input output in radians) */
void  astr_iarc(double xi, double eta, double rac, double decc, double *ra, double *dec)
{
  double  dist  = sqrt(xi*xi + eta*eta);

  if(dist < FLT_EPSILON)
  {
    *ra  = rac;
    *dec = decc;
  }
  else
  {
    double  uu    = cos(dist);
    double  vv    = xi *sin(dist)/dist;
    double  ww    = eta*sin(dist)/dist;

    double  cdecc = cos(decc);
    double  sdecc = sin(decc);
    double  crac  = cos(rac);
    double  srac  = sin(rac);

    double  uun   = uu*cdecc*crac - vv*srac - ww*sdecc*crac;
    double  vvn   = uu*cdecc*srac + vv*crac - ww*sdecc*srac;
    double  wwn   = uu*sdecc + ww*cdecc;

    *dec = asin(wwn);

    if(ABS(uun) < FLT_EPSILON)
      *ra = SIGNONEARG(xi)*PI_2;
    else
    {
      *ra = atan(vvn/uun);

      if(ABS(decc + eta) > PI_2)
        *ra += 2.0*(PI - ABS(*ra - rac));
      if(ABS(cos(*dec)) > FLT_EPSILON &&
	 rac + xi/cos(*dec) > PI_2 && rac + xi/cos(*dec) < PI3_2)
        *ra += PI;

      if(*ra < 0.0)
        *ra += PI2;
      else if(*ra > PI2)
	*ra -= PI2;
    }
  }

  return;
}  /* endof astr_iarc() */


/* inverse arc projection (input output in degrees) */
void  astr_irarc(double xi, double eta, double rac, double decc, double *ra, double *dec)
{
  double  rxi   = DEG2RAD(xi);
  double  reta  = DEG2RAD(eta);
  double  rrac  = DEG2RAD(rac);
  double  rdecc = DEG2RAD(decc);

  astr_iarc(rxi, reta, rrac, rdecc, ra, dec);

  *ra  = RAD2DEG(*ra);
  *dec = RAD2DEG(*dec);

  return;
}  /* endof astr_irarc() */


/* This function finds stars that are likely to be blends with other real variables. */
void findblends(int Nvars, int *N, double **t, double **mag, double **sig, _FindBlends *c)
{
  /* Procedure is
     1. for each star in the star list find all matches to stars in the variables list.
     2. Calculate the amplitude of each star list light curve at each of the periods from the matches.
     3. For each star in the variables list choose the matching star in the star list with the largest flux amplitude as the probable real variable.
  */

  int i, j, k, jnext, Nstarlist, *sortvarid, *sortstarid, *matchidstarlist, **matchidvarlist, *Nmatchvars, Nmatchstars, foundgood, Npointsinlc, sizeinlc = 0, bestid;
  FILE *starlist, *inlc, *outmatches;
  char *line, **starnames;
  double *starx, *stary, xmin, xmax, mindx, maxdx, mindxnext, maxdxnext, **matchamps, *inJD, *inmag, *influx, *inerr, per, *harmA, *harmB, fundA, fundB, meanval, bestamp, rdeg, cosrdeg, ra_ll, dec_ll, ra_ul, dec_ul, ra_lr, dec_lr, ra_ur, dec_ur;
  size_t line_size = MAXLEN;

  line = malloc(line_size);

  /* Allocate memory for vectors used by the harmonic fitting routine */
  if(c->Nharm > 0)
    {
      if((harmA = (double *) malloc(c->Nharm * sizeof(double))) == NULL ||
	 (harmB = (double *) malloc(c->Nharm * sizeof(double))) == NULL)
	error(ERR_MEMALLOC);
    }

  if(c->sepstarlist)
    {
      /* In this case the star list is distinct from the variables list, we read it in now */
      Nstarlist = 0;
      if((starlist = fopen(c->starlistname,"r")) == NULL)
	{
	  error2(ERR_FILENOTFOUND,c->starlistname);
	}
      while(gnu_getline(&line,&line_size,starlist) >= 0)
	{
	  if(line[0] != '#')
	    Nstarlist++;
	}
      if((starx = (double *) malloc(Nstarlist * sizeof(double))) == NULL ||
	 (stary = (double *) malloc(Nstarlist * sizeof(double))) == NULL ||
	 (starnames = (char **) malloc(Nstarlist * sizeof(char *))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0;i<Nstarlist;i++)
	{
	  if((starnames[i] = (char *) malloc(MAXLEN * sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	}
      rewind(starlist);
      i=0;
      while(gnu_getline(&line,&line_size,starlist) >= 0)
	{
	  if(line[0] != '#')
	    {
	      sscanf(line,"%s %lf %lf",starnames[i],&starx[i],&stary[i]);
	      i++;
	    }
	}
      fclose(starlist);
    }
  else
    {
      /* Otherwise we initialize it to the variables list */
      Nstarlist = Nvars;
      if((starx = (double *) malloc(Nstarlist * sizeof(double))) == NULL ||
	 (stary = (double *) malloc(Nstarlist * sizeof(double))) == NULL ||
	 (starnames = (char **) malloc(Nstarlist * sizeof(char *))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0;i<Nstarlist;i++)
	{
	  if((starnames[i] = (char *) malloc(MAXLEN * sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	  strcpy(starnames[i],c->varnames[i]);
	  starx[i] = c->varx[i];
	  stary[i] = c->vary[i];
	}
    }

  /* Get the x-coordinate sorted ids for the varlist and starlist */
  if((sortvarid = (int *) malloc(Nvars * sizeof(int))) == NULL ||
     (sortstarid = (int *) malloc(Nstarlist * sizeof(int))) == NULL ||
     (matchidstarlist = (int *) malloc(Nvars * sizeof(int))) == NULL ||
     (matchidvarlist = (int **) malloc(Nvars * sizeof(int *))) == NULL ||
     (matchamps = (double **) malloc(Nvars * sizeof(double *))) == NULL ||
     (Nmatchvars = (int *) calloc(Nvars, sizeof(int))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0;i<Nvars;i++)
    sortvarid[i] = i;
  for(i=0;i<Nstarlist;i++)
    sortstarid[i] = i;

  mysort2dblint_id(Nvars, c->varx, sortvarid);
  mysort2dblint_id(Nstarlist, starx, sortstarid);

  /* do the x-y point matching */
  xmin = c->varx[sortvarid[0]];
  xmax = c->varx[sortvarid[Nvars-1]];
  jnext = 0;
  foundgood = 0;
  if(c->radec)
    {
      rdeg = c->matchrad / 3600.0;
      cosrdeg = cos(DEG2RAD(rdeg));
    }

  for(i=0;i<Nstarlist;i++)
    {
      Nmatchstars = 0;
      if(i == 0)
	{
	  if(c->radec)
	    {
	      /* If we're doing ra/dec matching, use inverse arc projection to
		 find ra/dec coordinates of 2*radius X 2*radius box around
		 the star */
	      astr_irarc(-rdeg,-rdeg,starx[sortstarid[i]],stary[sortstarid[i]],&ra_ll,&dec_ll);
	      astr_irarc(-rdeg,rdeg,starx[sortstarid[i]],stary[sortstarid[i]],&ra_ul,&dec_ul);
	      astr_irarc(rdeg,-rdeg,starx[sortstarid[i]],stary[sortstarid[i]],&ra_lr,&dec_lr);
	      astr_irarc(rdeg,rdeg,starx[sortstarid[i]],stary[sortstarid[i]],&ra_ur,&dec_ur);
	      mindx = ra_ul;
	      if(mindx > ra_ll) mindx = ra_ll;
	      maxdx = ra_ur;
	      if(maxdx < ra_lr) maxdx = ra_lr;

	      /* Deal with the case where the pole is included in the search region */
	      if(stary[sortstarid[i]] - rdeg <= -90.0 || stary[sortstarid[i]] + rdeg >= 90.0) {
		mindx = 0.0;
		maxdx = 360.0;
	      }
	    }
	  else
	    {
	      /* Get min/max x */
	      mindx = starx[sortstarid[i]] - c->matchrad;
	      maxdx = starx[sortstarid[i]] + c->matchrad;
	    }
	}
      else
	{
	  mindx = mindxnext;
	  maxdx = maxdxnext;
	}

      /* Now find the min and max X search range for the next star in the list */
      if(i < Nstarlist - 1)
	{
	  if(c->radec)
	    {
	      astr_irarc(-rdeg,-rdeg,starx[sortstarid[i+1]],stary[sortstarid[i+1]],&ra_ll,&dec_ll);
	      astr_irarc(-rdeg,rdeg,starx[sortstarid[i+1]],stary[sortstarid[i+1]],&ra_ul,&dec_ul);
	      astr_irarc(rdeg,-rdeg,starx[sortstarid[i+1]],stary[sortstarid[i+1]],&ra_lr,&dec_lr);
	      astr_irarc(rdeg,rdeg,starx[sortstarid[i+1]],stary[sortstarid[i+1]],&ra_ur,&dec_ur);

	      /* Get min/max ra */
	      mindxnext = ra_ul;
	      if(mindxnext > ra_ll) mindxnext = ra_ll;
	      maxdxnext = ra_ur;
	      if(maxdxnext < ra_lr) maxdxnext = ra_lr;

	      /* Deal with the case where the pole is included in the search region */
	      if(stary[sortstarid[i+1]] - rdeg <= -90.0 || stary[sortstarid[i+1]] + rdeg >= 90.0){
		mindxnext = 0.0;
		maxdxnext = 360.0;
	      }
	    }
	  else
	    {
	      mindxnext = starx[sortstarid[i+1]] - c->matchrad;
	      maxdxnext = starx[sortstarid[i+1]] + c->matchrad;
	    }
	}

      if(!c->radec) /* Do the rectangular coordinate matching */
	{
	  /* Determine the id of the first star above the minimum x */
	  if(!foundgood)
	    {
	      if(xmin < mindx)
		{
		  if(jnext < Nvars ? (c->varx[sortvarid[jnext]] < mindx) : 0)
		    j = findXindx(c->varx,sortvarid,mindx,jnext,Nvars);
		  else
		    j = findXindx(c->varx,sortvarid,mindx,0,Nvars);
		}
	      else
		j= 0;
	    }
	  else
	    j = jnext;
	  foundgood = 0;

	  /* Search over the varlist and find matches */
	  while((j < Nvars ? (c->varx[sortvarid[j]] < maxdx) : 0))
	    {
	      if(fabs(c->vary[sortvarid[j]] - stary[sortstarid[i]]) < c->matchrad)
		{
		  if(distxy(starx[sortstarid[i]],stary[sortstarid[i]],c->varx[sortvarid[j]],c->vary[sortvarid[j]]) < c->matchrad)
		    {
		      matchidstarlist[Nmatchstars] = sortvarid[j];
		      Nmatchstars++;
		      if(!Nmatchvars[sortvarid[j]])
			{
			  if((matchidvarlist[sortvarid[j]] = (int *) malloc(sizeof(int))) == NULL ||
			     (matchamps[sortvarid[j]] = (double *) malloc(sizeof(double))) == NULL)
			    error(ERR_MEMALLOC);
			}
		      else
			{
			  if((matchidvarlist[sortvarid[j]] = (int *) realloc(matchidvarlist[sortvarid[j]], (Nmatchvars[sortvarid[j]] + 1) * sizeof(int))) == NULL ||
			     (matchamps[sortvarid[j]] = (double *) realloc(matchamps[sortvarid[j]], (Nmatchvars[sortvarid[j]] + 1) * sizeof(double))) == NULL)
			    error(ERR_MEMALLOC);
			}
		      matchidvarlist[sortvarid[j]][Nmatchvars[sortvarid[j]]] = sortstarid[i];
		      Nmatchvars[sortvarid[j]]++;
		    }
		}

	      /* Determine if this j index is the appropriate index to begin searching at for the next star in the starlist */
	      if((j > 0 ? (c->varx[sortvarid[j-1]] < mindxnext && c->varx[sortvarid[j]] > mindxnext) : (c->varx[sortvarid[j]] > mindxnext)))
		{
		  jnext = j;
		  foundgood = 1;
		}
	      j++;
	    }
	}
      else
	{
	  /* Do the ra/dec matching */
	  /* Deal first with the case where 360 is not included */
	  if(maxdx > mindx)
	    {
	      /* Determine the id of the first star above the minimum ra */
	      if(!foundgood)
		{
		  if(xmin < mindx)
		    {
		      if(jnext < Nvars ? (c->varx[sortvarid[jnext]] < mindx) : 0)
			j = findXindx(c->varx,sortvarid,mindx,jnext,Nvars);
		      else
			j = findXindx(c->varx,sortvarid,mindx,0,Nvars);
		    }
		  else
		    j = 0;
		}
	      else
		j = jnext;
	      foundgood = 0;
	      /* Now search over the varlist and find matches */
	      while((j < Nvars ? (c->varx[sortvarid[j]] < maxdx) : 0))
		{
		  if(fabs(c->vary[sortvarid[j]] - stary[sortstarid[i]]) < rdeg)
		    {
		      if(distradec(starx[sortstarid[i]],stary[sortstarid[i]],c->varx[sortvarid[j]],c->vary[sortvarid[j]]) > cosrdeg)
			{
			  matchidstarlist[Nmatchstars] = sortvarid[j];
			  Nmatchstars++;
			  if(!Nmatchvars[sortvarid[j]])
			    {
			      if((matchidvarlist[sortvarid[j]] = (int *) malloc(sizeof(int))) == NULL ||
				 (matchamps[sortvarid[j]] = (double *) malloc(sizeof(double))) == NULL)
				error(ERR_MEMALLOC);
			    }
			  else
			    {
			      if((matchidvarlist[sortvarid[j]] = (int *) realloc(matchidvarlist[sortvarid[j]], (Nmatchvars[sortvarid[j]] + 1) * sizeof(int))) == NULL ||
				 (matchamps[sortvarid[j]] = (double *) realloc(matchamps[sortvarid[j]], (Nmatchvars[sortvarid[j]] + 1) * sizeof(double))) == NULL)
				error(ERR_MEMALLOC);
			    }
			  matchidvarlist[sortvarid[j]][Nmatchvars[sortvarid[j]]] = sortstarid[i];
			  Nmatchvars[sortvarid[j]]++;
			}
		    }

		  /* Determine if this j index is the appropriate index to begin searching at for the next star in the starlist */
		  if((j > 0 ? (c->varx[sortvarid[j-1]] < mindxnext && c->varx[sortvarid[j]] > mindxnext) : (c->varx[sortvarid[j]] > mindxnext)))
		    {
		      jnext = j;
		      foundgood = 1;
		    }
		  j++;
		}
	      if(!foundgood) jnext = j;
	    }

	  /* Now deal with the case where 360 is included */
	  else
	    {
	      /* Again determine the starting position - this time we'll search first over stars between the minimum and 360 and then go back and search over stars between 0 and the maximum ra */
	      if(!foundgood)
		{
		  if(xmax > mindx)
		    {
		      if((jnext < Nvars ? (c->varx[sortvarid[jnext]] < mindx) : 0))
			j = findXindx(c->varx,sortvarid,mindx,jnext,Nvars);
		      else
			j = findXindx(c->varx,sortvarid,mindx,0,Nvars);
		    }
		  else
		    j = Nvars;
		}
	      else
		j = jnext;
	      foundgood = 0;
	      while(j < Nvars)
		{
		  if(fabs(c->vary[sortvarid[j]] - stary[sortstarid[i]]) < rdeg)
		    {
		      if(distradec(starx[sortstarid[i]],stary[sortstarid[i]],c->varx[sortvarid[j]],c->vary[sortvarid[j]]) > cosrdeg)
			{
			  matchidstarlist[Nmatchstars] = sortvarid[j];
			  Nmatchstars++;
			  if(!Nmatchvars[sortvarid[j]])
			    {
			      if((matchidvarlist[sortvarid[j]] = (int *) malloc(sizeof(int))) == NULL ||
				 (matchamps[sortvarid[j]] = (double *) malloc(sizeof(double))) == NULL)
				error(ERR_MEMALLOC);
			    }
			  else
			    {
			      if((matchidvarlist[sortvarid[j]] = (int *) realloc(matchidvarlist[sortvarid[j]], (Nmatchvars[sortvarid[j]] + 1) * sizeof(int))) == NULL ||
				 (matchamps[sortvarid[j]] = (double *) realloc(matchamps[sortvarid[j]], (Nmatchvars[sortvarid[j]] + 1) * sizeof(double))) == NULL)
				error(ERR_MEMALLOC);
			    }
			  matchidvarlist[sortvarid[j]][Nmatchvars[sortvarid[j]]] = sortstarid[i];
			  Nmatchvars[sortvarid[j]]++;
			}
		    }

		  /* Determine if this j index is the appropriate index to begin searching at for the next star in the starlist */
		  if((j > 0 ? (c->varx[sortvarid[j-1]] < mindxnext && c->varx[sortvarid[j]] > mindxnext) : (c->varx[sortvarid[j]] > mindxnext)))
		    {
		      jnext = j;
		      foundgood = 1;
		    }
		  j++;
		}
	      j = 0;
	      while((j < Nvars ? (c->varx[j] < maxdx) : 0))
		{
		  if(fabs(c->vary[sortvarid[j]] - stary[sortstarid[i]]) < rdeg)
		    {
		      if(distradec(starx[sortstarid[i]],stary[sortstarid[i]],c->varx[sortvarid[j]],c->vary[sortvarid[j]]) > cosrdeg)
			{
			  matchidstarlist[Nmatchstars] = sortvarid[j];
			  Nmatchstars++;
			  if(!Nmatchvars[sortvarid[j]])
			    {
			      if((matchidvarlist[sortvarid[j]] = (int *) malloc(sizeof(int))) == NULL ||
				 (matchamps[sortvarid[j]] = (double *) malloc(sizeof(double))) == NULL)
				error(ERR_MEMALLOC);
			    }
			  else
			    {
			      if((matchidvarlist[sortvarid[j]] = (int *) realloc(matchidvarlist[sortvarid[j]], (Nmatchvars[sortvarid[j]] + 1) * sizeof(int))) == NULL ||
				 (matchamps[sortvarid[j]] = (double *) realloc(matchamps[sortvarid[j]], (Nmatchvars[sortvarid[j]] + 1) * sizeof(double))) == NULL)
				error(ERR_MEMALLOC);
			    }
			  matchidvarlist[sortvarid[j]][Nmatchvars[sortvarid[j]]] = sortstarid[i];
			  Nmatchvars[sortvarid[j]]++;
			}
		    }

		  /* Determine if this j index is the appropriate index to begin searching at for the next star in the starlist */
		  if((j > 0 ? (c->varx[sortvarid[j-1]] < mindxnext && c->varx[sortvarid[j]] > mindxnext) : (c->varx[sortvarid[j]] > mindxnext)))
		    {
		      jnext = j;
		      foundgood = 1;
		    }
		  j++;
		}
	      if(!foundgood)
		jnext = j;
	    }
	}
      /* Now get the amplitudes for this star if it matches to any variable */
      if(Nmatchstars > 0)
	{
	  /* First we have to get the light curve for the star, we read it in from a file if the starlist is distinct from the varlist, otherwise we just copy it over from the input varlist lcs */
	  if(c->sepstarlist)
	    {
	      if((inlc = fopen(starnames[sortstarid[i]],"r")) == NULL)
		error2(ERR_FILENOTFOUND, starnames[sortstarid[i]]);
	      Npointsinlc = 0;
	      while(gnu_getline(&line,&line_size,inlc) >= 0)
		{
		  if(line[0] != '#')
		    Npointsinlc++;
		}
	      rewind(inlc);
	    }
	  else
	    Npointsinlc = N[sortstarid[i]];
	  if(Npointsinlc > sizeinlc)
	    {
	      if(!sizeinlc)
		{
		  if((inJD = (double *) malloc(Npointsinlc * sizeof(double))) == NULL ||
		     (inmag = (double *) malloc(Npointsinlc * sizeof(double))) == NULL ||
		     (influx = (double *) malloc(Npointsinlc * sizeof(double))) == NULL ||
		     (inerr = (double *) malloc(Npointsinlc * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      else
		{
		  if((inJD = (double *) realloc(inJD, Npointsinlc * sizeof(double))) == NULL ||
		     (inmag = (double *) realloc(inmag, Npointsinlc * sizeof(double))) == NULL ||
		     (influx = (double *) realloc(influx, Npointsinlc * sizeof(double))) == NULL ||
		     (inerr = (double *) realloc(inerr, Npointsinlc * sizeof(double))) == NULL)
		    error(ERR_MEMALLOC);
		}
	      sizeinlc = Npointsinlc;
	    }
	  if(c->sepstarlist)
	    {
	      Npointsinlc = 0;
	      if(!c->converttoflux)
		{
		  while(gnu_getline(&line,&line_size,inlc) >= 0)
		    {
		      if(line[0] != '#')
			{
			  sscanf(line,"%lf %lf %lf",&inJD[Npointsinlc],&influx[Npointsinlc],&inerr[Npointsinlc]);
			  Npointsinlc++;
			}
		    }
		}
	      else
		{
		  while(gnu_getline(&line,&line_size,inlc) >= 0)
		    {
		      if(line[0] != '#')
			{
			  sscanf(line,"%lf %lf %lf",&inJD[Npointsinlc],&inmag[Npointsinlc],&inerr[Npointsinlc]);
			  /* Convert mag to flux */
			  influx[Npointsinlc] = exp(0.4*log(10.)*(c->zeromag - inmag[Npointsinlc]));
			  inerr[Npointsinlc] = 0.4*log(10.)*inerr[Npointsinlc]*influx[Npointsinlc];
			  Npointsinlc++;
			}
		    }
		}
	      fclose(inlc);
	    }
	  else
	    {
	      Npointsinlc=0;
	      for(Npointsinlc=0;Npointsinlc < N[sortstarid[i]];Npointsinlc++)
		{
		  inJD[Npointsinlc] = t[sortstarid[i]][Npointsinlc];
		  inmag[Npointsinlc] = mag[sortstarid[i]][Npointsinlc];
		  inerr[Npointsinlc] = sig[sortstarid[i]][Npointsinlc];
		  if(!c->converttoflux)
		    influx[Npointsinlc] = inmag[Npointsinlc];
		  else
		    {
		      influx[Npointsinlc] = exp(0.4*log(10.)*(c->zeromag - inmag[Npointsinlc]));
		      inerr[Npointsinlc] = 0.4*log(10.)*inerr[Npointsinlc]*influx[Npointsinlc];
		    }
		}
	    }

	  /* TODO: Here is where we should add the option to run arbitrary
	     VARTOOLS processing commands on the starlist lcs */

	  /* Now get the amplitudes for the light curve at each of the
             test periods. */
	  for(k = 0; k < Nmatchstars; k++)
	    {
	      per = c->periods[matchidstarlist[k]][0];
	      if(per <= 0.)
		{
		  per = inJD[Npointsinlc - 1] - inJD[0];
		}
	      /* TODO: use a more sophisticated scheme for getting the amplitudes
		 e.g. automatically determine Nharm */
	      dokillharms(Npointsinlc, inJD, influx, inerr, 1, &per, 0, c->Nharm, NULL, NULL, &harmA, &harmB, &fundA, &fundB, &meanval, 0, NULL, &(matchamps[matchidstarlist[k]][Nmatchvars[matchidstarlist[k]]-1]), 1, KILLHARM_OUTTYPE_DEFAULT, -1.);
	    }
	}
      if(!foundgood) jnext = j;
    }

  /* Output the full list of matches if we're doing that */
  if(c->outputmatches)
    {
      if((outmatches = fopen(c->outmatchesfilename,"w")) == NULL)
	error2(ERR_CANNOTWRITE,c->outmatchesfilename);
      for(i=0;i<Nvars;i++)
	{
	  fprintf(outmatches,"%s",c->varnames[i]);
	  for(j=0;j<Nmatchvars[i];j++)
	    {
	      k = matchidvarlist[i][j];
	      fprintf(outmatches," %s %g",starnames[k],matchamps[i][j]);
	    }
	  fprintf(outmatches,"\n");
	}
      fclose(outmatches);
    }

  /* Now for each lc in the varlist choose the matching lc in the starlist with the highest flux amplitude */
  for(i=0;i<Nvars;i++)
    {
      if(!Nmatchvars[i])
	{
	  sprintf(c->varblendnames[i],"No_match");
	  c->blendamps[i] = -1.;
	}
      else
	{
	  bestid = 0;
	  bestamp = matchamps[i][0];
	  for(k=1;k<Nmatchvars[i];k++)
	    {
	      if(matchamps[i][k] > bestamp)
		{
		  bestid = k;
		  bestamp = matchamps[i][k];
		}
	    }
	  strcpy(c->varblendnames[i],starnames[matchidvarlist[i][bestid]]);
	  c->blendamps[i] = bestamp;
	}
    }

  /* Do clean-up */

  free(line);
  if(c->Nharm > 0)
    {
      free(harmA);
      free(harmB);
    }
  free(starx);
  free(stary);
  for(i=0;i<Nstarlist;i++)
    free(starnames[i]);
  free(starnames);
  free(sortvarid);
  free(sortstarid);
  free(matchidstarlist);
  for(i=0;i<Nvars;i++)
    {
      if(Nmatchvars[i])
	{
	  free(matchidvarlist[i]);
	  free(matchamps[i]);
	}
    }
  free(matchidvarlist);
  free(matchamps);
  free(Nmatchvars);
  free(inJD);
  free(inmag);
  free(influx);
  free(inerr);
}
