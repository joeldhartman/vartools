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
/* Routines to print the tables, print the results, read in and out files, etc. Part of the vartools program by J. Hartman. */

#include "commands.h"
#include "programdata.h"
#include "functions.h"
#include <sys/stat.h>
#include <sys/types.h>
#include <assert.h>
#include <errno.h>


#ifdef USECFITSIO
#include "fitsio.h"
#endif

#define MIN_CHUNK 64

/* This is a copy of the GNU getstr and getline commands, included
   here for portability in case the user is not using a gnu
   compiler. Modifications to the function names to avoid conflicting
   with the existing GNU functions.  The original was written by
   Jan Brittenson.
*/
int
gnu_getstr (lineptr, n, stream, terminator, offset)
     char **lineptr;
     size_t *n;
     FILE *stream;
     char terminator;
     int offset;
{
  int nchars_avail;		/* Allocated but unused chars in *LINEPTR.  */
  char *read_pos;		/* Where we're reading into *LINEPTR. */
  int ret;

  if (!lineptr || !n || !stream)
    {
      errno = EINVAL;
      return -1;
    }

  if (!*lineptr)
    {
      *n = MIN_CHUNK;
      *lineptr = malloc (*n);
      if (!*lineptr)
	{
	  errno = ENOMEM;
	  return -1;
	}
    }

  nchars_avail = *n - offset;
  read_pos = *lineptr + offset;

  for (;;)
    {
      int save_errno;
      register int c = getc (stream);

      save_errno = errno;

      /* We always want at least one char left in the buffer, since we
	 always (unless we get an error while reading the first char)
	 NUL-terminate the line buffer.  */

      assert((*lineptr + *n) == (read_pos + nchars_avail));
      if (nchars_avail < 2)
	{
	  if (*n > MIN_CHUNK)
	    *n *= 2;
	  else
	    *n += MIN_CHUNK;

	  nchars_avail = *n + *lineptr - read_pos;
	  *lineptr = realloc (*lineptr, *n);
	  if (!*lineptr)
	    {
	      errno = ENOMEM;
	      return -1;
	    }
	  read_pos = *n - nchars_avail + *lineptr;
	  assert((*lineptr + *n) == (read_pos + nchars_avail));
	}

      if (ferror (stream))
	{
	  /* Might like to return partial line, but there is no
	     place for us to store errno.  And we don't want to just
	     lose errno.  */
	  errno = save_errno;
	  return -1;
	}

      if (c == EOF)
	{
	  /* Return partial line, if any.  */
	  if (read_pos == *lineptr)
	    return -1;
	  else
	    break;
	}

      *read_pos++ = c;
      nchars_avail--;

      if (c == terminator)
	/* Return the line.  */
	break;
    }

  /* Done - NUL terminate and return the number of chars read.  */
  *read_pos = '\0';

  ret = read_pos - (*lineptr + offset);
  return ret;
}

int
gnu_getline (lineptr, n, stream)
     char **lineptr;
     size_t *n;
     FILE *stream;
{
  return gnu_getstr (lineptr, n, stream, '\n', 0);
}


/* This function returns the index of a string after skipping one column*/
int skipone(char *line)
{
  int i = 0;
  while(line[i] == ' ' || line[i] == '\t')
    i++;
  while(line[i] != ' ' && line[i] != '\t' && line[i] != '\n' && line[i] != '\0')
    i++;
  while(line[i] == ' ' || line[i] == '\t')
    i++;
  return i;
}

/* This function parses one column of a string and returns the index for the start of the next column */
int parseone(char *line, void *val, int vartype)
{
  int i = 0, j = 0;
  char *line2 = NULL;
  if((line2 = (char *) malloc((strlen(line)+1)*sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  while(line[i] == ' ' || line[i] == '\t')
    i++;
  while(line[i] != ' ' && line[i] != '\t' && line[i] != '\n' && line[i] != '\0')
    {
      line2[j] = line[i];
      i++;
      j++;
    }
  line2[j] = '\0';
  switch(vartype)
    {
    case VARTOOLS_TYPE_DOUBLE:
      *((double *) val) = atof(line2);
      break;
    case VARTOOLS_TYPE_STRING:
      sprintf((char *) val, "%s", line2);
      break;
    case VARTOOLS_TYPE_INT:
      *((int *) val) = atoi(line2);
      break;
    case VARTOOLS_TYPE_SHORT:
      *((short *) val) = (short) atoi(line2);
      break;
    case VARTOOLS_TYPE_FLOAT:
      *((float *) val) = (float) atof(line2);
      break;
    case VARTOOLS_TYPE_LONG:
      *((long *) val) = atol(line2);
      break;
    case VARTOOLS_TYPE_CHAR:
      *((char *) val) = line2[0];
      break;
    default:
      error(ERR_BADTYPE);
      break;
    }
  while(line[i] == ' ' || line[i] == '\t')
    i++;
  if(line2 != NULL) free(line2);
  return(i);
}

/* This function parses one column of a string and returns the index for the start of the next column */
int parseone_growstring(char *line, void *val, int vartype, int *colstrlen)
{
  int i = 0, j = 0;
  char *line2 = NULL;
  if((line2 = (char *) malloc((strlen(line)+1)*sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  while(line[i] == ' ' || line[i] == '\t')
    i++;
  while(line[i] != ' ' && line[i] != '\t' && line[i] != '\n' && line[i] != '\0')
    {
      line2[j] = line[i];
      i++;
      j++;
    }
  line2[j] = '\0';
  switch(vartype)
    {
    case VARTOOLS_TYPE_DOUBLE:
      *((double *) val) = atof(line2);
      break;
    case VARTOOLS_TYPE_STRING:
      if(strlen(line2) >= *(colstrlen)) {
	if(!(*colstrlen)) {
	  *colstrlen = strlen(line2) + 1;
	  if((*((char **) val) = (char *) malloc((*colstrlen)*sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	} else {
	  *colstrlen = strlen(line2) + 1;
	  if((*((char **) val) = (char *) realloc(*((char **) val), (*colstrlen)*sizeof(char))) == NULL)
	    error(ERR_MEMALLOC);
	}
      }
      sprintf(*((char **) val), "%s", line2);
      break;
    case VARTOOLS_TYPE_INT:
      *((int *) val) = atoi(line2);
      break;
    case VARTOOLS_TYPE_SHORT:
      *((short *) val) = (short) atoi(line2);
      break;
    case VARTOOLS_TYPE_FLOAT:
      *((float *) val) = (float) atof(line2);
      break;
    case VARTOOLS_TYPE_LONG:
      *((long *) val) = atol(line2);
      break;
    case VARTOOLS_TYPE_CHAR:
      *((char *) val) = line2[0];
      break;
    default:
      error(ERR_BADTYPE);
      break;
    }
  while(line[i] == ' ' || line[i] == '\t')
    i++;
  if(line2 != NULL) free(line2);
  return(i);
}


/* This function skips one column of a string and returns the index for the start of the next column, splitting on the specified delimiter string */
int skiponedelimstring(char *line, char *delim)
{
  int i = 0, k, test;
  while(line[i] != '\n' && line[i] != '\0')
    {
      if(line[i] == delim[0]) {
	test = 1;
	k = 1;
	while(delim[k] != '\0' ? (line[i+k] != '\n' && line[i+k] != '\0' ? line[i+k] == delim[k] : 0 ) : 0) k++;
	if(delim[k] != '\0') { test = 0; }
      } else {
	test = 0;
      }
      if(test) {
	i = i + k;
	break;
      }
      else {
	i++;
      }
    }
  return(i);
}


/* This function parses one column of a string and returns the index for the start of the next column */
int parseonedelimstring(char *line, void *val, int vartype, char *delim)
{
  int i = 0, j = 0, k, test;
  char *line2 = NULL;
  if((line2 = (char *) malloc((strlen(line)+1)*sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  while(line[i] != '\n' && line[i] != '\0')
    {
      if(line[i] == delim[0]) {
	test = 1;
	k = 1;
	while(delim[k] != '\0' ? (line[i+k] != '\n' && line[i+k] != '\0' ? line[i+k] == delim[k] : 0 ) : 0) k++;
	if(delim[k] != '\0') { test = 0; }
      } else {
	test = 0;
      }
      if(test) {
	i = i + k;
	break;
      }
      else {
	line2[j] = line[i];
	i++;
	j++;
      }
    }
  line2[j] = '\0';
  switch(vartype)
    {
    case VARTOOLS_TYPE_DOUBLE:
      *((double *) val) = atof(line2);
      break;
    case VARTOOLS_TYPE_STRING:
      sprintf(*((char **) val), "%s", line2);
      break;
    case VARTOOLS_TYPE_INT:
      *((int *) val) = atoi(line2);
      break;
    case VARTOOLS_TYPE_SHORT:
      *((short *) val) = (short) atoi(line2);
      break;
    case VARTOOLS_TYPE_FLOAT:
      *((float *) val) = (float) atof(line2);
      break;
    case VARTOOLS_TYPE_LONG:
      *((long *) val) = atol(line2);
      break;
    case VARTOOLS_TYPE_CHAR:
      *((char *) val) = line2[0];
      break;
    default:
      error(ERR_BADTYPE);
      break;
    }
  if(line2 != NULL) free(line2);
  return(i);
}

/* This function skips one column of a string and returns the index for the start of the next column, splitting on the specified delimiter character */
int skiponedelimchar(char *line, char delim)
{
  int i = 0, k, test;
  while(line[i] != '\n' && line[i] != '\0')
    {
      if(line[i] == delim) {
	i++;
	break;
      }
      else {
	i++;
      }
    }
  return(i);
}


/* This function parses one column of a string and returns the index for the start of the next column */
int parseonedelimchar(char *line, void *val, int vartype, char delim)
{
  int i = 0, j = 0, k, test;
  char *line2 = NULL;
  if((line2 = (char *) malloc((strlen(line)+1)*sizeof(char))) == NULL)
    error(ERR_MEMALLOC);
  while(line[i] != '\n' && line[i] != '\0')
    {
      if(line[i] == delim) {
	i++;
	break;
      }
      else {
	line2[j] = line[i];
	i++;
	j++;
      }
    }
  line2[j] = '\0';
  switch(vartype)
    {
    case VARTOOLS_TYPE_DOUBLE:
      *((double *) val) = atof(line2);
      break;
    case VARTOOLS_TYPE_STRING:
      sprintf(*((char **) val), "%s", line2);
      break;
    case VARTOOLS_TYPE_INT:
      *((int *) val) = atoi(line2);
      break;
    case VARTOOLS_TYPE_SHORT:
      *((short *) val) = (short) atoi(line2);
      break;
    case VARTOOLS_TYPE_FLOAT:
      *((float *) val) = (float) atof(line2);
      break;
    case VARTOOLS_TYPE_LONG:
      *((long *) val) = atol(line2);
      break;
    case VARTOOLS_TYPE_CHAR:
      *((char *) val) = line2[0];
      break;
    default:
      error(ERR_BADTYPE);
      break;
    }
  if(line2 != NULL) free(line2);
  return(i);
}

void GetOutputFilenameFromCommand(char *lcoutname, char *lcname, char *outdir,
				  int lc_name_num, char *lcnamecommand)
{
  FILE *lc_name_pipe;
  char *execcommand;
  int size_execcommand;
  char *line = NULL;
  size_t line_size = 2048;
  int i1 = 0;

  line = (char *) malloc(line_size*sizeof(char));

  size_execcommand = 10 + strlen(lcname) + strlen(outdir) + ceil(log((double) lc_name_num)) + strlen(lcnamecommand) + 1;

  if((execcommand = (char *) malloc((size_execcommand+1))) == NULL)
    error(ERR_MEMALLOC);

  sprintf(execcommand,"echo %s %s %d | %s", lcname, outdir, lc_name_num, lcnamecommand);

  if((lc_name_pipe = popen(execcommand,"r")) == NULL) {
    error2(ERR_OUTPUTFILENAMECOMMAND,execcommand);
  }

  if(gnu_getline(&line,&line_size,lc_name_pipe) <= 0) {
    error2(ERR_OUTPUTFILENAMECOMMAND,execcommand);
  }
  i1 = 0;
  while(line[i1] != '\n' && line[i1] != '\0') {
    lcoutname[i1] = line[i1];
    i1++;
  }
  if(!i1)
    error2(ERR_OUTPUTFILENAMECOMMAND,execcommand);
  
  pclose(lc_name_pipe);
  free(line);
  free(execcommand);
}


void GetOutputFilename(char *lcoutname, char *lcname, char *outdir,
		       char *suffix, char *format, int lc_name_num)
/* Determines the name for an output file following the options
   associated with, for example, the -o vartools command.

     lcoutname - where the output file name to use is returned.

     lcname - the input file name.

     outdir - the base directory to output the file to.

     suffix - a string which is appended to the end of the file name.

     format - an optional format string to use in determining the
              filename. This is parsed following the rules of the
              "nameformat" option to the vartools -o command. This is
              ignored if it is NULL.

     lc_name_num - the light curve number index (needed in case the
              format contains at %d flag).
*/
{
  int i1, i2, i3, i4, i5, i6;

  char tmpstring[MAXLEN];

  i1 = 0; i2 = 0; i5 = -1;
  while(lcname[i1] != '\0')
    {
      if(lcname[i1] == '/')
	i2 = i1 + 1;
      if(lcname[i1] == '.')
	i5 = i1 - 1;
      i1++;
    }

  if(format == NULL || format[0] == '\0') {
    /* We're using the default name */
    sprintf(lcoutname,"%s/%s.%s",outdir,&(lcname[i2]),suffix);
  }
  else {
    /* Use the given format to determine the output name */
    sprintf(lcoutname,"%s/",outdir);
    i1 = strlen(lcoutname);
    i3 = 0;
    while(format[i3] != '\0')
      {
	if(format[i3] != '%')
	  {
	    lcoutname[i1] = format[i3];
	    i1++;
	    lcoutname[i1] = '\0';
	    i3++;
	  }
	else
	  {
	    i3++;
	    if(format[i3] == 's')
	      {
		i3++;
		sprintf(&lcoutname[i1],"%s",&(lcname[i2]));
		i1 = strlen(lcoutname);
	      }
	    else if(format[i3] == 'b')
	      {
		i3++;
		i6 = i2;
		if(i5 > i2) {
		  for(; i6 <= i5; i6++) {
		    lcoutname[i1] = lcname[i6];
		    i1++;
		  }
		  lcoutname[i1] = '\0';
		}
		else {
		  sprintf(&lcoutname[i1],"%s",&(lcname[i2]));
		  i1 = strlen(lcoutname);
		}
	      }
	    else if(format[i3] == 'd')
	      {
		i3++;
		sprintf(&lcoutname[i1],"%d",lc_name_num+1);
		i1 = strlen(lcoutname);
	      }
	    else if(format[i3] == '0')
	      {
		i3++;
		tmpstring[0] = '%';
		tmpstring[1] = '0';
		i4 = 2;
		while(format[i3] >= '1' && format[i3] <= '9')
		  {
		    tmpstring[i4] = format[i3];
		    i4++;
		    i3++;
		  }
		if(format[i3] != 'd')
		  error(ERR_INVALIDOUTPUTFORMAT);
		i3++;
		tmpstring[i4] = 'd';
		i4++;
		tmpstring[i4] = '\0';
		sprintf(&lcoutname[i1],tmpstring,lc_name_num+1);
		i1 = strlen(lcoutname);
	      }
	    else if(format[i3] == '%')
	      {
		i3++;
		lcoutname[i1] = '%';
		i1++;
		lcoutname[i1] = '\0';
	      }
	    else
	      error(ERR_INVALIDOUTPUTFORMAT);
	  }
      }
  }
}

void ReadGlobalDecorr(ProgramData *p, Command *c)
{
  int i, j, k, Ncommands;
  FILE *global_file;
  char *line;
  size_t line_size = MAXLEN;
  line = malloc(line_size);
  Ncommands = p->Ncommands;
  for(i=0;i<Ncommands;i++)
    {
      if(c[i].cnum == CNUM_DECORR)
	{
	  c[i].Decorr->size_globaldecorrvector = 0;
	  for(j=0;j<c[i].Decorr->N_globalterms;j++)
	    {
	      if((global_file = fopen(c[i].Decorr->global_file_names[j],"r")) == NULL)
		error2(ERR_FILENOTFOUND,c[i].Decorr->global_file_names[j]);
	      c[i].Decorr->N_globaldecorr_JD = 0;
	      while(gnu_getline(&line,&line_size,global_file) >= 0)
		if(line[0] != '#')
		  c[i].Decorr->N_globaldecorr_JD++;
	      if(c[i].Decorr->size_globaldecorrvector == 0)
		{
		  c[i].Decorr->size_globaldecorrvector = c[i].Decorr->N_globaldecorr_JD;
		  if(p->matchstringid)
		    {
		      if((c[i].Decorr->globaldecorr_stringid = (char **) malloc(c[i].Decorr->size_globaldecorrvector * sizeof(char *))) == NULL ||
			 (c[i].Decorr->globaldecorr_stringid_idx = (int *) malloc(c[i].Decorr->size_globaldecorrvector * sizeof(int))) == NULL)
			error(ERR_MEMALLOC);
		      for(k=0;k<c[i].Decorr->size_globaldecorrvector;k++)
			{
			  c[i].Decorr->globaldecorr_stringid_idx[k] = k;
			  if((c[i].Decorr->globaldecorr_stringid[k] = (char *) malloc(MAXIDSTRINGLENGTH * sizeof(char))) == NULL)
			    error(ERR_MEMALLOC);
			}
		    }
		  else
		    {
		      if((c[i].Decorr->globaldecorr_JD = (double *) malloc(c[i].Decorr->size_globaldecorrvector * sizeof(double))) == NULL)
			error(ERR_MEMALLOC);
		    }
		  if((c[i].Decorr->globaldecorr_terms = (double **) malloc(c[i].Decorr->size_globaldecorrvector * sizeof(double *))) == NULL)
		    error(ERR_MEMALLOC);
		  for(k=0;k<c[i].Decorr->size_globaldecorrvector;k++)
		    {
		      if((c[i].Decorr->globaldecorr_terms[k] = (double *) malloc(c[i].Decorr->N_globalterms * sizeof(double))) == NULL)
			error(ERR_MEMALLOC);
		    }
		}
	      else if(c[i].Decorr->N_globaldecorr_JD != c[i].Decorr->size_globaldecorrvector)
		error(ERR_INVALIDGLOBALDECORR);
	      rewind(global_file);
	      k = 0;
	      if(p->matchstringid)
		{
		  while(gnu_getline(&line,&line_size,global_file) >= 0)
		    {
		      if(line[0] != '#')
			{
			  sscanf(line,"%s %lf",c[i].Decorr->globaldecorr_stringid[k],&c[i].Decorr->globaldecorr_terms[k][j]);
			  k++;
			}
		    }
		  mysortstringint(k, MAXIDSTRINGLENGTH, c[i].Decorr->globaldecorr_stringid, c[i].Decorr->globaldecorr_stringid_idx);
		}
	      else
		{
		  while(gnu_getline(&line,&line_size,global_file) >= 0)
		    {
		      if(line[0] != '#')
			{
			  sscanf(line,"%lf %lf",&c[i].Decorr->globaldecorr_JD[k],&c[i].Decorr->globaldecorr_terms[k][j]);
			  k++;
			}
		    }
		}
	      fclose(global_file);
	    }
	}
    }
  free(line);
}


void Filldecorr_matrix(ProgramData *p, Command *c, int lc)
{
  /* This routine copies the data from the lcdecorr_terms_in tensor and the globaldecorr_terms matrix into the decorr_terms tensor. */

  int i, j, k, m, l;
  double tmpterm;
  for(i=0; i < p->Ncommands ; i++)
    {
      if(c[i].cnum == CNUM_DECORR)
	{
	  /* First add the lcdecorr_terms_in tensor */
	  m = c[i].Decorr->N_globalterms;
	  for(j=0;j < p->NJD[lc]; j++)
	    for(k = 0; k < c[i].Decorr->N_lcterms; k++)
	      c[i].Decorr->decorr_terms[lc][j][k + m] = c[i].Decorr->lcdecorr_terms_in[k][lc][j];
	  if(c[i].Decorr->subtractfirstterm)
	    {
	      for(k = 0; k < c[i].Decorr->N_lcterms; k++)
		{
		  tmpterm = c[i].Decorr->decorr_terms[lc][0][k + m];
		  for(j=0; j < p->NJD[lc]; j++)
		    {
		      c[i].Decorr->decorr_terms[lc][j][k + m] -= tmpterm;
		    }
		}
	    }

	  /* Next add the globaldecorr_terms matrix */
	  if(m > 0)
	    {
	      if(p->matchstringid)
		{
		  /* Use string-ids to match the global decorr file to the light curve */
		  k = 0;
		  j = 0;
		  while (k < p->NJD[lc] && j < c[i].Decorr->N_globaldecorr_JD)
		    {
		      while(k < p->NJD[lc] ? strncmp(c[i].Decorr->globaldecorr_stringid[c[i].Decorr->globaldecorr_stringid_idx[j]],p->stringid[lc][p->stringid_idx[lc][k]],MAXIDSTRINGLENGTH) > 0 : 0)
			{
			  for(l=0;l<m;l++)
			    c[i].Decorr->decorr_terms[lc][p->stringid_idx[lc][k]][l] = sqrt(-1);
			  k++;
			}
		      if(k < p->NJD[lc] ? !strncmp(c[i].Decorr->globaldecorr_stringid[c[i].Decorr->globaldecorr_stringid_idx[j]], p->stringid[lc][p->stringid_idx[lc][k]],MAXIDSTRINGLENGTH) : 0)
			{
			  for(l=0;l<m;l++)
			    c[i].Decorr->decorr_terms[lc][p->stringid_idx[lc][k]][l] = c[i].Decorr->globaldecorr_terms[c[i].Decorr->globaldecorr_stringid_idx[j]][l];
			  k++;

			}
		      j++;
		    }
		  for(;k < p->NJD[lc]; k++)
		    for(l=0;l<m;l++)
		      c[i].Decorr->decorr_terms[lc][p->stringid_idx[l][k]][l] = sqrt(-1);

		  if(c[i].Decorr->subtractfirstterm)
		    {
		      for(k = 0; k < c[i].Decorr->N_globalterms; k++)
			{
			  tmpterm = c[i].Decorr->decorr_terms[lc][0][k];
			  if(!isnan(tmpterm))
			    {
			      for(j=0; j < p->NJD[lc]; j++)
				{
				  c[i].Decorr->decorr_terms[lc][j][k] -= tmpterm;
				}
			    }
			}
		    }
		}
	      else
		{
		  /* Use JDs to match the global decorr file to the light curve */
		  k = 0;
		  j = 0;
		  while (k < p->NJD[lc] && j < c[i].Decorr->N_globaldecorr_JD)
		    {
		      while(k < p->NJD[lc] ? (c[i].Decorr->globaldecorr_JD[j] > p->t[lc][k] + p->JDTOL) : 0)
			{
			  for(l=0;l<m;l++)
			    c[i].Decorr->decorr_terms[lc][k][l] = sqrt(-1);
			  k++;
			}
		      if(k < p->NJD[lc] ? (c[i].Decorr->globaldecorr_JD[j] > p->t[lc][k] - p->JDTOL) : 0)
			{
			  for(l=0;l<m;l++)
			    c[i].Decorr->decorr_terms[lc][k][l] = c[i].Decorr->globaldecorr_terms[j][l];
			  k++;

			}
		      j++;
		    }
		  for(;k < p->NJD[lc]; k++)
		    for(l=0;l<m;l++)
		      c[i].Decorr->decorr_terms[lc][k][l] = sqrt(-1);

		  if(c[i].Decorr->subtractfirstterm)
		    {
		      for(k = 0; k < c[i].Decorr->N_globalterms; k++)
			{
			  tmpterm = c[i].Decorr->decorr_terms[lc][0][k];
			  if(!isnan(tmpterm))
			    {
			      for(j=0; j < p->NJD[lc]; j++)
				{
				  c[i].Decorr->decorr_terms[lc][j][k] -= tmpterm;
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

void dotab(FILE *outfile,int tabflag)
{
  if(tabflag)
    fprintf(outfile,"\t");
  else
    fprintf(outfile," ");
}

void dotab_buffer(char *c,int tabflag)
{
  if(tabflag)
    sprintf(c,"\t");
  else
    sprintf(c," ");
}


void docolnumber(int numbercolumns)
{
  static int colnum = 1;
  if(numbercolumns)
    printf("%d_",colnum);
  colnum++;
}

#ifdef _USEBINARY_LC
/* Prepare a new light curve header */
/* stopped working here, write the routines to allow output of binary light curves 
void make_new_binarylc_header(ProgramData *p, char ***lc_header, 
			      int *memsize_lc_header)
{
}

void write_binary_lightcurve_header()
{
}

void write_binary_lightcurve(ProgramData *p, int threadid, int lcid, 
			     char *outname,
			     int usecolumnformat, int Nvars, 
			     _Variable **variables,
			     char **formats, int noclobber)
{
}
*/
#endif

#ifdef USECFITSIO
void write_fits_lightcurve(ProgramData *p, int threadid, int lcid, 
			   char *outname,
			   int usecolumnformat, int Nvars, 
			   _Variable **variables,
			   char **formats, int noclobber, int copyheaderfrominput, int logcommandline, char **description, char **units)
{
  int status = 0, status2 = 0, tfields, idx, strlenval;
  long nrows;
  char **ttype = NULL, **tform = NULL, **tunit = NULL;
  fitsfile *outfile, *input_header_file;
  int i, closefile=1, N, j;
  double *t, *mag, *sig;
  double *outdbl_vec = NULL;
  float *outfloat_vec = NULL;
  int *outint_vec = NULL;
  short *outshort_vec = NULL;
  long *outlong_vec = NULL;
  char *outchar_vec = NULL;
  char **outstring_vec = NULL;
  char tryout[MAXLEN];
  char commentpar[MAXLEN];
  char *fmt;
  double nulldbl;
  float nullflt;

  int *maxlenstringvec = NULL;

  nulldbl = (double) sqrt(-1.0);
  nullflt = (float) nulldbl;

#ifdef PARALLEL
  if(p->Nproc_allow > 1) {
    while(pthread_mutex_trylock(&(p->cfitsio_mutex)));
  }
#endif

  fits_create_file(&outfile,outname,&status);
  if(status) {
    if(!noclobber) {
      status = 0;
      sprintf(tryout,"!%s",outname);
      fits_create_file(&outfile,tryout,&status);
      if(status) {
	fits_report_error(stderr, status);
	error(ERR_FITSERROR);
      }
    }
    else {
      fits_report_error(stderr, status);
      error(ERR_FITSERROR);
    }
  }

  if(p->is_inputlc_fits[lcid] && copyheaderfrominput) {
    /* Copy the header from the input light curve */
    if((fits_open_file(&input_header_file,p->lcnames[lcid],READONLY,&status2)))
      error2(ERR_CANNOTOPEN,p->lcnames[lcid]);
    fits_movabs_hdu(input_header_file, 1, NULL, &status2);
    fits_copy_header(input_header_file,outfile,&status);
    fits_close_file(input_header_file,&status2);
  }
  


  if(logcommandline) {
    fits_write_history(outfile, p->cmdline, &status);
  }

  /* prepare the table for output */
  if(!usecolumnformat) {
    tfields = 3;
  } else {
    tfields = Nvars;
  }
  if((ttype = (char **) malloc(tfields*sizeof(char *))) == NULL ||
     (tform = (char **) malloc(tfields*sizeof(char *))) == NULL ||
     (tunit = (char **) malloc(tfields*sizeof(char *))) == NULL ||
     (maxlenstringvec = (int *) malloc(tfields*sizeof(int))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0; i < tfields; i++) {
    if((ttype[i] = (char *) malloc(MAXLEN*sizeof(char))) == NULL ||
       (tform[i] = (char *) malloc(MAXLEN*sizeof(char))) == NULL ||
       (tunit[i] = (char *) malloc(MAXLEN*sizeof(char))) == NULL)
      error(ERR_MEMALLOC);
  }
  if(!usecolumnformat) {
    sprintf(ttype[0],"time");
    sprintf(ttype[1],"mag");
    sprintf(ttype[2],"err");
    sprintf(tform[0],"D");
    sprintf(tform[1],"D");
    sprintf(tform[2],"D");
    sprintf(tunit[0],"d");
    sprintf(tunit[1],"mag");
    sprintf(tunit[2],"mag");
  }
  else {
    for(j=0; j < Nvars; j++) {
      sprintf(ttype[j],"%s",variables[j]->varname);
      switch(variables[j]->datatype) {
      case VARTOOLS_TYPE_DOUBLE:
      case VARTOOLS_TYPE_CONVERTJD:
	sprintf(tform[j],"D");
	break;
      case VARTOOLS_TYPE_FLOAT:
	sprintf(tform[j],"E");
	break;
      case VARTOOLS_TYPE_INT:
	if(sizeof(int) == 4) {
	  sprintf(tform[j],"J");
	  break;
	} else if(sizeof(int) == 8) {
	  sprintf(tform[j],"K");
	  break;
	} else if(sizeof(int) == 2) {
	  sprintf(tform[j],"I");
	  break;
	} else {
	  error(ERR_BADTYPE);
	}
      case VARTOOLS_TYPE_LONG:
	if(sizeof(long) == 8) {
	  sprintf(tform[j],"K");
	  break;
	} else if(sizeof(long) == 4) {
	  sprintf(tform[j],"J");
	  break;
	} else if(sizeof(long) == 2) {
	  sprintf(tform[j],"I");
	  break;
	} else {
	  error(ERR_BADTYPE);
	}
      case VARTOOLS_TYPE_SHORT:
	if(sizeof(short) == 1) {
	  sprintf(tform[j],"B");
	  break;
	} else if(sizeof(short) == 2) {
	  sprintf(tform[j],"I");
	  break;
	} else if(sizeof(short) == 4) {
	  sprintf(tform[j],"J");
	  break;
	} else if(sizeof(short) == 8) {
	  sprintf(tform[j],"K");
	  break;
	} else {
	  error(ERR_BADTYPE);
	}
      case VARTOOLS_TYPE_CHAR:
	sprintf(tform[j],"A");
	break;
      case VARTOOLS_TYPE_STRING:
	switch(variables[j]->vectortype) {
	case VARTOOLS_VECTORTYPE_CONSTANT:
	  maxlenstringvec[j] = strlen(((char **) variables[j]->dataptr)[0]);
	  break;
	case VARTOOLS_VECTORTYPE_SCALAR:
	case VARTOOLS_VECTORTYPE_INLIST:
	  if(variables[j]->vectortype == VARTOOLS_VECTORTYPE_SCALAR) {
	    idx = threadid;
	  } else {
	    idx = lcid;
	  }
	  maxlenstringvec[j] = strlen((((char ***) variables[j]->dataptr)[0])[idx]);
	  break;
	case VARTOOLS_VECTORTYPE_LC:
	  maxlenstringvec[j] = 0;
	  for(i=0; i < p->NJD[threadid]; i++) {
	    strlenval = strlen((*((char ****) variables[j]->dataptr))[threadid][i]);
	    if(strlenval > maxlenstringvec[j]) maxlenstringvec[j] = strlenval;
	  }
	  break;
	default:
	  maxlenstringvec[j] = 1;
	  break;
	}
	if(maxlenstringvec[j] <= 0) maxlenstringvec[j] = 1;
	sprintf(tform[j],"%dA",maxlenstringvec[j]);
	break;
      default:
	error(ERR_BADTYPE);
      }
      if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	sprintf(tunit[j],"%s",formats[j]);
      }
      else {
	if(variables[j] == p->tvar) {
	  sprintf(tunit[j],"d");
	}
	else if(variables[j] == p->magvar || variables[j] == p->sigvar) {
	  sprintf(tunit[j],"mag");
	}
	else {
	  sprintf(tunit[j],"unknown");
	}
      }
    }
  }
  nrows = (long) p->NJD[threadid];
  fits_create_tbl(outfile, BINARY_TBL, nrows, tfields, ttype, tform,
		  tunit, NULL, &status);

  if(status){
    fits_report_error(stderr, status);
    error(ERR_FITSERROR);
  }

  if(description != NULL) {
    for(j=0; j < Nvars; j++) {
      if(description[j][0] != '\0') break;
    }
    if(j < Nvars) {
      sprintf(commentpar, "--BEGIN-ASTROPY-SERIALIZED-COLUMNS--");
      fits_write_comment(outfile, commentpar, &status);
      if(status){
	fits_report_error(stderr, status);
	error(ERR_FITSERROR);
      }
      sprintf(commentpar, "datatype:");
      fits_write_comment(outfile, commentpar, &status);
      if(status){
	fits_report_error(stderr, status);
	error(ERR_FITSERROR);
      }
      for(j=0; j < Nvars; j++) {
	sprintf(commentpar, "- name: %.*s", strlen(ttype[j]), ttype[j]);
	fits_write_comment(outfile, commentpar, &status);
	if(status){
	  fits_report_error(stderr, status);
	  error(ERR_FITSERROR);
	}
	if(units != NULL && units[j][0] != '\0') {
	  sprintf(commentpar, "  unit: %.*s", strlen(units[j]), units[j]);
	  fits_write_comment(outfile, commentpar, &status);
	  if(status){
	    fits_report_error(stderr, status);
	    error(ERR_FITSERROR);
	  }
	}
/*	else if(formats != NULL && formats[j][0] != '\0') {
	  sprintf(commentpar, "  unit: %.*s", strlen(formats[j]), formats[j]);
	  fits_write_comment(outfile, commentpar, &status);
	  if(status){
	    fits_report_error(stderr, status);
	    error(ERR_FITSERROR);
	  }
	}*/
	if(description != NULL && description[j][0] != '\0') {
	  sprintf(commentpar, "  description: %.*s", 55, description[j]);
	  fits_write_comment(outfile, commentpar, &status);
	  if(status){
	    fits_report_error(stderr, status);
	    error(ERR_FITSERROR);
	  }
	}
      }
      sprintf(commentpar, "meta:");
      fits_write_comment(outfile, commentpar, &status);
      if(status){
	fits_report_error(stderr, status);
	error(ERR_FITSERROR);
      }
      sprintf(commentpar, "  __serialized_columns__: {}");
      fits_write_comment(outfile, commentpar, &status);
      if(status){
	fits_report_error(stderr, status);
	error(ERR_FITSERROR);
      }
      sprintf(commentpar, "--END-ASTROPY-SERIALIZED-COLUMNS--");
      fits_write_comment(outfile, commentpar, &status);
      if(status){
	fits_report_error(stderr, status);
	error(ERR_FITSERROR);
      }
    }
  }

  if(!usecolumnformat) {
    t = p->t[threadid];
    mag = p->mag[threadid];
    sig = p->sig[threadid];
    fits_write_colnull(outfile, TDOUBLE, 1, 1, 1, nrows,
			  (void *) t, (void *) (&nulldbl),
		       &status);
    if(status){
      fits_report_error(stderr, status);
      error(ERR_FITSERROR);
    }
    fits_write_colnull(outfile, TDOUBLE, 2, 1, 1, nrows,
		       (void *) mag, (void *) (&nulldbl),
		       &status);
    if(status){
      fits_report_error(stderr, status);
      error(ERR_FITSERROR);
    }
    fits_write_colnull(outfile, TDOUBLE, 3, 1, 1, nrows,
			  (void *) sig, (void *) (&nulldbl),
		       &status);
    if(status){
      fits_report_error(stderr, status);
      error(ERR_FITSERROR);
    }
  } else {
    N = p->NJD[threadid];
    for(j=0; j < Nvars; j++) {
      switch(variables[j]->vectortype) {
      case VARTOOLS_VECTORTYPE_CONSTANT:
	switch(variables[j]->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  if(outdbl_vec == NULL) {
	    if((outdbl_vec = (double *) malloc(N*sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outdbl_vec[i] = ((double *) variables[j]->dataptr)[0];
	  }
	  if(fits_write_colnull(outfile, TDOUBLE, j+1, 1, 1, nrows,
			     (void *) outdbl_vec, (void *) (&nulldbl),
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  if(outfloat_vec == NULL) {
	    if((outfloat_vec = (float *) malloc(N*sizeof(float))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outfloat_vec[i] = ((float *) variables[j]->dataptr)[0];
	  }
	  if(fits_write_colnull(outfile, TFLOAT, j+1, 1, 1, nrows,
			     (void *) outfloat_vec, (void *) (&nullflt),
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  if(outint_vec == NULL) {
	    if((outint_vec = (int *) malloc(N*sizeof(int))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outint_vec[i] = ((int *) variables[j]->dataptr)[0];
	  }
	  if(fits_write_col(outfile, TINT, j+1, 1, 1, nrows,
			     (void *) outint_vec,
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_LONG:
	  if(outlong_vec == NULL) {
	    if((outlong_vec = (long *) malloc(N*sizeof(long))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outlong_vec[i] = ((long *) variables[j]->dataptr)[0];
	  }
	  if(fits_write_col(outfile, TLONG, j+1, 1, 1, nrows,
			     (void *) outlong_vec,
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_SHORT:
	  if(outshort_vec == NULL) {
	    if((outshort_vec = (short *) malloc(N*sizeof(short))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outshort_vec[i] = ((short *) variables[j]->dataptr)[0];
	  }
	  if(fits_write_col(outfile, TSHORT, j+1, 1, 1, nrows,
			     (void *) outshort_vec,
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_CHAR:
	  if(outchar_vec == NULL) {
	    if((outchar_vec = (char *) malloc(N)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outchar_vec[i] = ((char *) variables[j]->dataptr)[0];
	  }
	  if(fits_write_col(outfile, TBYTE, j+1, 1, 1, nrows,
			     (void *) outchar_vec,
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_STRING:
	  if(outstring_vec != NULL) {
	    for(i=0; i < N; i++)
	      free(outstring_vec[i]);
	    free(outstring_vec);
	  }
	  if((outstring_vec = (char **) malloc(N*sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(i=0; i < N; i++) {
	    if((outstring_vec[i] = (char *) malloc((maxlenstringvec[j]+1))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    sprintf(outstring_vec[i],"%s",((char **) variables[j]->dataptr)[0]);
	  }
	  if(fits_write_col(outfile, TBYTE, j+1, 1, 1, nrows,
			     (void *) outstring_vec,
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	default:
	  error(ERR_BADTYPE);
	}
	break;
      case VARTOOLS_VECTORTYPE_SCALAR:
      case VARTOOLS_VECTORTYPE_INLIST:
	if(variables[j]->vectortype == VARTOOLS_VECTORTYPE_SCALAR) {
	  idx = threadid;
	} else {
	  idx = lcid;
	}
	switch(variables[j]->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  if(outdbl_vec == NULL) {
	    if((outdbl_vec = (double *) malloc(N*sizeof(double))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outdbl_vec[i] = (((double **) variables[j]->dataptr)[0])[idx];
	  }
	  if(fits_write_colnull(outfile, TDOUBLE, j+1, 1, 1, nrows,
			     (void *) outdbl_vec, (void *) (&nulldbl),
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  if(outfloat_vec == NULL) {
	    if((outfloat_vec = (float *) malloc(N*sizeof(float))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outfloat_vec[i] = (((float **) variables[j]->dataptr)[0])[idx];
	  }
	  if(fits_write_colnull(outfile, TFLOAT, j+1, 1, 1, nrows,
			     (void *) outfloat_vec, (void *) (&nullflt),
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  if(outint_vec == NULL) {
	    if((outint_vec = (int *) malloc(N*sizeof(int))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outint_vec[i] = (((int **) variables[j]->dataptr)[0])[lcid];
	  }
	  if(fits_write_col(outfile, TINT, j+1, 1, 1, nrows,
			     (void *) outint_vec,
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_LONG:
	  if(outlong_vec == NULL) {
	    if((outlong_vec = (long *) malloc(N*sizeof(long))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outlong_vec[i] = (((long **) variables[j]->dataptr)[0])[idx];
	  }
	  if(fits_write_col(outfile, TLONG, j+1, 1, 1, nrows,
			     (void *) outlong_vec,
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_SHORT:
	  if(outshort_vec == NULL) {
	    if((outshort_vec = (short *) malloc(N*sizeof(short))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outshort_vec[i] = (((short **) variables[j]->dataptr)[0])[idx];
	  }
	  if(fits_write_col(outfile, TSHORT, j+1, 1, 1, nrows,
			     (void *) outshort_vec,
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_CHAR:
	  if(outchar_vec == NULL) {
	    if((outchar_vec = (char *) malloc(N)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    outchar_vec[i] = (((char **) variables[j]->dataptr)[0])[idx];
	  }
	  if(fits_write_col(outfile, TBYTE, j+1, 1, 1, nrows,
			     (void *) outchar_vec,
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_STRING:
	  if(outstring_vec != NULL) {
	    for(i=0; i < N; i++)
	      free(outstring_vec[i]);
	    free(outstring_vec);
	  }
	  if((outstring_vec = (char **) malloc(N*sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(i=0; i < N; i++) {
	    if((outstring_vec[i] = (char *) malloc((maxlenstringvec[j]+1))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    sprintf(outstring_vec[i],"%s",(((char ***) variables[j]->dataptr)[0])[idx]);
	  }
	  if(fits_write_col(outfile, TSTRING, j+1, 1, 1, nrows,
			     (void *) outstring_vec,
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	default:
	  error(ERR_BADTYPE);
	}
	break;
      case VARTOOLS_VECTORTYPE_LC:
	switch(variables[j]->datatype) {
	case VARTOOLS_TYPE_DOUBLE:
	case VARTOOLS_TYPE_CONVERTJD:
	  if(fits_write_colnull(outfile, TDOUBLE, j+1, 1, 1, nrows,
				(void *) (*((double ***) variables[j]->dataptr))[threadid], 
				(void *) (&nulldbl),
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_FLOAT:
	  if(fits_write_colnull(outfile, TFLOAT, j+1, 1, 1, nrows,
				(void *) (*((float ***) variables[j]->dataptr))[threadid], 
				(void *) (&nullflt),
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_INT:
	  if(fits_write_col(outfile, TINT, j+1, 1, 1, nrows,
				(void *) (*((int ***) variables[j]->dataptr))[threadid], 
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_LONG:
	  if(fits_write_col(outfile, TLONG, j+1, 1, 1, nrows,
				(void *) (*((long ***) variables[j]->dataptr))[threadid], 
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_SHORT:
	  if(fits_write_col(outfile, TSHORT, j+1, 1, 1, nrows,
				(void *) (*((short ***) variables[j]->dataptr))[threadid], 
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_CHAR:
	  if(fits_write_col(outfile, TBYTE, j+1, 1, 1, nrows,
				(void *) (*((char ***) variables[j]->dataptr))[threadid], 
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	case VARTOOLS_TYPE_STRING:
	  if(outstring_vec != NULL) {
	    for(i=0; i < N; i++)
	      free(outstring_vec[i]);
	    free(outstring_vec);
	  }
	  if((outstring_vec = (char **) malloc(N*sizeof(char *))) == NULL)
	    error(ERR_MEMALLOC);
	  for(i=0; i < N; i++) {
	    if((outstring_vec[i] = (char *) malloc((maxlenstringvec[j]+1))) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  for(i=0; i < N; i++) {
	    sprintf(outstring_vec[i],"%s",(*((char ****) variables[j]->dataptr))[threadid][i]);
	  }
	  if(fits_write_col(outfile, TSTRING, j+1, 1, 1, nrows,
			     (void *) outstring_vec,
				&status)) {
	    error(ERR_FITS_WRITECOLUMN);
	  }
	  break;
	default:
	  error(ERR_BADTYPE);
	}
	break;
      default:
	error(ERR_CODEERROR);
      }
    }
  }

  if(p->fits_header_adds != NULL) {
    int (*fitskeyfunctocall)(fitsfile *, int, char *, ...);
    void *fitskeyvalptr;
    for(i=0; i < p->fits_header_adds[threadid].N_added_keywords; i++) {
      if(p->fits_header_adds[threadid].hdrterms[i].hdutouse == 0)
	continue;
      if(p->fits_header_adds[threadid].hdrterms[i].updateexisting) {
	fitskeyfunctocall = &(fits_update_key);
      } else {
	fitskeyfunctocall = &(fits_write_key);
      }
      switch(p->fits_header_adds[threadid].hdrterms[i].datatype) {
      case TDOUBLE:
	fitskeyvalptr = &(p->fits_header_adds[threadid].hdrterms[i].dbl_val);
	break;
      case TINT:
	fitskeyvalptr = &(p->fits_header_adds[threadid].hdrterms[i].int_val);
	break;
      case TLONG:
	fitskeyvalptr = &(p->fits_header_adds[threadid].hdrterms[i].long_val);
	break;
      case TSTRING:
	fitskeyvalptr = p->fits_header_adds[threadid].hdrterms[i].string_val;
	break;
      default:
	error(ERR_CODEERROR);
	break;
      }
      fitskeyfunctocall(outfile,
			p->fits_header_adds[threadid].hdrterms[i].datatype, 
			p->fits_header_adds[threadid].hdrterms[i].keyname, 
			fitskeyvalptr,
			p->fits_header_adds[threadid].hdrterms[i].comment,
			&status);
    }

    fits_movabs_hdu(outfile, 1, NULL, &status2);
    for(i=0; i < p->fits_header_adds[threadid].N_added_keywords; i++) {
      if(p->fits_header_adds[threadid].hdrterms[i].hdutouse == 1)
	continue;
      if(p->fits_header_adds[threadid].hdrterms[i].updateexisting) {
	fitskeyfunctocall = &(fits_update_key);
      } else {
	fitskeyfunctocall = &(fits_write_key);
      }
      switch(p->fits_header_adds[threadid].hdrterms[i].datatype) {
      case TDOUBLE:
	fitskeyvalptr = &(p->fits_header_adds[threadid].hdrterms[i].dbl_val);
	break;
      case TINT:
	fitskeyvalptr = &(p->fits_header_adds[threadid].hdrterms[i].int_val);
	break;
      case TLONG:
	fitskeyvalptr = &(p->fits_header_adds[threadid].hdrterms[i].long_val);
	break;
      case TSTRING:
	fitskeyvalptr = p->fits_header_adds[threadid].hdrterms[i].string_val;
	break;
      default:
	error(ERR_CODEERROR);
	break;
      }
      fitskeyfunctocall(outfile,
			p->fits_header_adds[threadid].hdrterms[i].datatype, 
			p->fits_header_adds[threadid].hdrterms[i].keyname, 
			fitskeyvalptr,
			p->fits_header_adds[threadid].hdrterms[i].comment,
			&status);
    }
  }


  fits_close_file(outfile,&status);
  if(status) {
    fits_report_error(stderr, status);
    error(ERR_FITSERROR);
  }

#ifdef PARALLEL
  if(p->Nproc_allow > 1) {
    pthread_mutex_unlock(&(p->cfitsio_mutex));
  }
#endif

  for(i=0; i < tfields; i++) {
    free(ttype[i]); free(tform[i]); free(tunit[i]);
  }
  if(ttype != NULL) free(ttype);
  if(tform != NULL) free(tform);
  if(tunit != NULL) free(tunit);
  if(maxlenstringvec != NULL) free(maxlenstringvec);
  if(outdbl_vec != NULL) free(outdbl_vec);
  if(outfloat_vec != NULL) free(outfloat_vec);
  if(outint_vec != NULL) free(outint_vec);
  if(outshort_vec != NULL) free(outshort_vec);
  if(outlong_vec != NULL) free(outlong_vec);
  if(outchar_vec != NULL) free(outchar_vec);
  if(outstring_vec != NULL) {
    for(i=0; i < N; i++) {
      free(outstring_vec[i]);
    }
    free(outstring_vec);
  }
}
#endif

void writelightcurves(ProgramData *p, int threadid, int lcid, char *outname,
		      int usecolumnformat, int Nvars, _Variable **variables,
		      char **formats, int noclobber, char sepchar, int logcommandline)
{
  FILE *out;
  int i, closefile=1, N, j, idx;
  double *t, *mag, *sig;
  double outdbl;
  float outfloat;
  int outint;
  short outshort;
  long outlong;
  char outchar;
  char *outstring;
  char *fmt;
  char fmtint[] = "%d";
  char fmtdouble[] = "%.17g";
  char fmtfloat[] = "%f";
  char fmtstring[] = "%s";
  char fmtchar[] = "%c";
  char fmtlong[] = "%d";
  char fmtshort[] = "%d";
  struct stat st;
  if(!noclobber) {
    if(!strncmp(outname,"-",1) && strlen(outname) == 1)
      {
	out = stdout;
	closefile = 0;
      }
    else if((out = fopen(outname,"w")) == NULL)
      error2(ERR_CANNOTWRITE, outname);
  }
  else {
    if(!strncmp(outname,"-",1) && strlen(outname) == 1)
      {
	out = stdout;
	closefile = 0;
      }
    else {
      if(stat(outname,&st)) {
	if((out = fopen(outname,"w")) == NULL)
	  error2(ERR_CANNOTWRITE, outname);
      } else {
	error2(ERR_FILEEXISTS_NOCLOBBER, outname);
      }
    }
  }

  if(logcommandline) {
    fprintf(out,"#%s\n",p->cmdline);
  }

  if(!usecolumnformat) {
    t = p->t[threadid];
    mag = p->mag[threadid];
    sig = p->sig[threadid];
    N = p->NJD[threadid];
    for(i=0;i<N;i++)
      if(!isnan(mag[i]))
	fprintf(out,"%17.9f%c%9.5f%c%9.5f\n",t[i],sepchar,mag[i],sepchar,sig[i]);
  } else {
    for(i=0;i<p->NJD[threadid];i++) {
      for(j=0; j < Nvars; j++) {
	if(j) fprintf(out,"%c",sepchar);
	switch(variables[j]->vectortype) {
	case VARTOOLS_VECTORTYPE_CONSTANT:
	  switch(variables[j]->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    outdbl = *((double *) variables[j]->dataptr);
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outdbl);
	    } else {
	      fprintf(out,fmtdouble,outdbl);
	    }
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    outdbl = *((double *) variables[j]->dataptr);
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outdbl);
	    } else {
	      fprintf(out,fmtdouble,outdbl);
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    outfloat = *((float *) variables[j]->dataptr);
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outfloat);
	    } else {
	      fprintf(out,fmtfloat,outfloat);
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    outint = *((int *) variables[j]->dataptr);
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outint);
	    } else {
	      fprintf(out,fmtint,outint);
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    outlong = *((long *) variables[j]->dataptr);
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outlong);
	    } else {
	      fprintf(out,fmtlong,outlong);
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    outshort = *((short *) variables[j]->dataptr);
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outshort);
	    } else {
	      fprintf(out,fmtshort,outshort);
	    }
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    outchar = *((char *) variables[j]->dataptr);
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outchar);
	    } else {
	      fprintf(out,fmtchar,outchar);
	    }
	    break;
	  case VARTOOLS_TYPE_STRING:
	    outstring = *((char **) variables[j]->dataptr);
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outstring);
	    } else {
	      fprintf(out,fmtstring,outstring);
	    }
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	  break;
	case VARTOOLS_VECTORTYPE_SCALAR:
	case VARTOOLS_VECTORTYPE_INLIST:
	  if(variables[j]->vectortype == VARTOOLS_VECTORTYPE_SCALAR)
	    idx = threadid;
	  else
	    idx = lcid;
	  switch(variables[j]->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    outdbl = (*((double **) variables[j]->dataptr))[idx];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outdbl);
	    } else {
	      fprintf(out,fmtdouble,outdbl);
	    }
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    outdbl = (*((double **) variables[j]->dataptr))[idx];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outdbl);
	    } else {
	      fprintf(out,fmtdouble,outdbl);
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    outfloat = (*((float **) variables[j]->dataptr))[idx];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outfloat);
	    } else {
	      fprintf(out,fmtfloat,outfloat);
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    outint = (*((int **) variables[j]->dataptr))[idx];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outint);
	    } else {
	      fprintf(out,fmtint,outint);
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    outlong = (*((long **) variables[j]->dataptr))[idx];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outlong);
	    } else {
	      fprintf(out,fmtlong,outlong);
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    outshort = (*((short **) variables[j]->dataptr))[idx];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outshort);
	    } else {
	      fprintf(out,fmtshort,outshort);
	    }
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    outchar = (*((char **) variables[j]->dataptr))[idx];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outchar);
	    } else {
	      fprintf(out,fmtchar,outchar);
	    }
	    break;
	  case VARTOOLS_TYPE_STRING:
	    outstring = (*((char ***) variables[j]->dataptr))[idx];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outstring);
	    } else {
	      fprintf(out,fmtstring,outstring);
	    }
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	  break;
	case VARTOOLS_VECTORTYPE_LC:
	  switch(variables[j]->datatype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    outdbl = (*((double ***) variables[j]->dataptr))[threadid][i];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outdbl);
	    } else {
	      fprintf(out,fmtdouble,outdbl);
	    }
	    break;
	  case VARTOOLS_TYPE_CONVERTJD:
	    outdbl = (*((double ***) variables[j]->dataptr))[threadid][i];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outdbl);
	    } else {
	      fprintf(out,fmtdouble,outdbl);
	    }
	    break;
	  case VARTOOLS_TYPE_FLOAT:
	    outfloat = (*((float ***) variables[j]->dataptr))[threadid][i];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outfloat);
	    } else {
	      fprintf(out,fmtfloat,outfloat);
	    }
	    break;
	  case VARTOOLS_TYPE_INT:
	    outint = (*((int ***) variables[j]->dataptr))[threadid][i];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outint);
	    } else {
	      fprintf(out,fmtint,outint);
	    }
	    break;
	  case VARTOOLS_TYPE_LONG:
	    outlong = (*((long ***) variables[j]->dataptr))[threadid][i];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outlong);
	    } else {
	      fprintf(out,fmtlong,outlong);
	    }
	    break;
	  case VARTOOLS_TYPE_SHORT:
	    outshort = (*((short ***) variables[j]->dataptr))[threadid][i];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outshort);
	    } else {
	      fprintf(out,fmtshort,outshort);
	    }
	    break;
	  case VARTOOLS_TYPE_CHAR:
	    outchar = (*((char ***) variables[j]->dataptr))[threadid][i];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outchar);
	    } else {
	      fprintf(out,fmtchar,outchar);
	    }
	    break;
	  case VARTOOLS_TYPE_STRING:
	    outstring = (*((char ****) variables[j]->dataptr))[threadid][i];
	    if(formats[j] != NULL ? formats[j][0] != '\0' : 0) {
	      fprintf(out,formats[j],outstring);
	    } else {
	      fprintf(out,fmtstring,outstring);
	    }
	    break;
	  default:
	    error(ERR_BADTYPE);
	  }
	  break;
	default:
	  error(ERR_CODEERROR);
	}
      }
      fprintf(out,"\n");
    }
  }


  if(closefile)
    fclose(out);
}

void DoOutputLightCurve(ProgramData *p, _Outputlcs *c, int lcid, int threadid)
{
  int i1, i2, i3, i4, i5, i6;
  char outname[MAXLEN], tmpstring[MAXLEN];
  /* First determine the output name of the light curve */
  if(p->listflag || p->Ncopycommands > 0)
    {
      if(c->useoutnamecommand) {
	for(i1=0; i1 < MAXLEN; i1++) outname[i1] = '\0';
	GetOutputFilenameFromCommand(outname, p->lcnames[lcid], c->outdir,
				     lcid+1, c->outnamecommand);
      }
      else if(c->namefrominlist) {
	sprintf(outname,"%s/%s",c->outdir,c->inputlistoutnames[lcid]);
      }
      else {
	i1 = 0;
	i2 = 0;
	i5 = -1;
	while(p->lcnames[lcid][i1] != '\0')
	  {
	    if(p->lcnames[lcid][i1] == '/')
	      i2 = i1 + 1;
	    if(p->lcnames[lcid][i1] == '.')
	      i5 = i1 - 1;
	    i1++;
	  }
	if(!c->useformat)
	  sprintf(outname,"%s/%s",c->outdir,&p->lcnames[lcid][i2]);
	else
	  {
	    sprintf(outname,"%s/",c->outdir);
	    i1=strlen(outname);
	    i3=0;
	    while(c->format[i3] != '\0')
	      {
		if(c->format[i3] != '%')
		  {
		    outname[i1] = c->format[i3];
		    i1++;
		    outname[i1] = '\0';
		    i3++;
		  }
		else
		  {
		    i3++;
		    if(c->format[i3] == 's')
		      {
			i3++;
			sprintf(&outname[i1],"%s",&p->lcnames[lcid][i2]);
			i1 = strlen(outname);
		      }
		    else if(c->format[i3] == 'b')
		      {
			i3++;
			i6 = i2;
			if(i5 >= i2) {
			  for(; i6 <= i5; i6++) {
			    outname[i1] = p->lcnames[lcid][i6];
			    i1++;
			  }
			  outname[i1] = '\0';
			}
			else {
			  sprintf(&outname[i1],"%s",&p->lcnames[lcid][i2]);
			  i1 = strlen(outname);
			}
		      }
		    else if(c->format[i3] == 'd')
		      {
			i3++;
			sprintf(&outname[i1],"%d",lcid+1);
			i1 = strlen(outname);
		      }
		    else if(c->format[i3] == '0')
		      {
			i3++;
			tmpstring[0] = '%';
			tmpstring[1] = '0';
			i4 = 2;
			while(c->format[i3] >= '1' && c->format[i3] <= '9')
			{
			  tmpstring[i4] = c->format[i3];
			  i4++;
			  i3++;
			}
			if(c->format[i3] != 'd')
			  error(ERR_INVALIDOUTPUTFORMAT);
			i3++;
			tmpstring[i4] = 'd';
			i4++;
			tmpstring[i4] = '\0';
			sprintf(&outname[i1],tmpstring,lcid+1);
			i1 = strlen(outname);
		      }
		    else if(c->format[i3] == '%')
		      {
			i3++;
			outname[i1] = '%';
			i1++;
			outname[i1] = '\0';
		      }
		    else
		      error(ERR_INVALIDOUTPUTFORMAT);
		  }
	      }
	  }
	}
#ifdef USECFITSIO
      if(c->outfits) {
	/* Check if the name has .fits at the end, if not, append it */
	i4 = strlen(outname);
	if(i4 > 5 ? !strcmp(&(outname[i4-5]),".fits") : 0) {
	  write_fits_lightcurve(p, threadid, lcid, outname, c->usecolumnformat,
				c->Nvar, c->variables, c->printfformats,
				c->noclobber, c->copyheaderfrominput,
				c->logcommandline, c->descriptions, c->units);
	} else {
	  sprintf(outname,"%s.fits",outname);
	  write_fits_lightcurve(p, threadid, lcid, outname, c->usecolumnformat,
				c->Nvar, c->variables, c->printfformats,
				c->noclobber, c->copyheaderfrominput,
				c->logcommandline, c->descriptions, c->units);
	}
      }
      else
#endif
	writelightcurves(p, threadid, lcid, outname, c->usecolumnformat, c->Nvar, c->variables, c->printfformats, c->noclobber, c->sepchar, c->logcommandline);
    }
  else if(p->fileflag && !p->Ncopycommands)
    {
#ifdef USECFITSIO
      if(c->outfits) {
	/* Check if the name has .fits at the end, if not, append it */
	i4 = strlen(c->outdir);
	if(i4 > 5 ? !strcmp(&(c->outdir[i4-5]),".fits") : 0) {
	  write_fits_lightcurve(p, threadid, lcid, c->outdir, c->usecolumnformat, c->Nvar, c->variables, c->printfformats, c->noclobber,c->copyheaderfrominput,c->logcommandline,c->descriptions,c->units);
	} else {
	  sprintf(outname,"%s.fits",c->outdir);
	  write_fits_lightcurve(p, threadid, lcid, outname, c->usecolumnformat, c->Nvar, c->variables, c->printfformats, c->noclobber,c->copyheaderfrominput,c->logcommandline,c->descriptions,c->units);
	}
      }
      else
#endif
	writelightcurves(p, threadid, lcid, c->outdir, c->usecolumnformat, c->Nvar, c->variables, c->printfformats, c->noclobber, c->sepchar, c->logcommandline);
    }
}

void ReadDatesFiles(ProgramData *p, Command *c)
{
  int i;
  for(i=0;i<p->Ncommands;i++)
    if(c[i].cnum == CNUM_JSTET)
      c[i].Jstet->wkmax = readdates(c[i].Jstet->datesname,c[i].Jstet->Jstet_time);
}

double readdates(char *datesname,double Jstet_time)
{
  FILE *dates;
  char *line;
  size_t line_size = MAXLEN;
  double *JD, dt, wkmax = 0.0;
  int Njd, i, flag = 0;

  line = malloc(line_size);

  if((dates = fopen(datesname,"r")) == NULL)
    {
      fprintf(stderr,"Cannot Open %s\n",datesname);
      exit(4);
    }
  Njd = 0;
  while(gnu_getline(&line,&line_size,dates) >= 0)
    if(line[0] != '#')
      Njd++;
  rewind(dates);

  if((JD = (double *) malloc(Njd * sizeof(double))) == NULL)
    {
      fprintf(stderr,"Memory Allocation Error\n");
      exit(3);
    }
  Njd = 0;
  while(gnu_getline(&line,&line_size,dates) >= 0)
    if(line[0] != '#')
      {
	sscanf(line,"%lf ",&(JD[Njd]));
	Njd++;
      }
  fclose(dates);
  for(i=0;i<Njd - 1; i++)
    {
      dt = fabs(JD[i+1] - JD[i]);
      if(dt > Jstet_time && !flag)
	{
	  flag = 0;
	  wkmax += 0.1;
	}
      if(dt > Jstet_time && flag) flag = 0;
      if(dt <= Jstet_time)
	{
	  flag = 1;
	  wkmax += 1.0;
	}
    }
  if(!flag) wkmax += 0.1;
  free(line);
  free(JD);
  return(wkmax);
}

void Switchtobasename(ProgramData *p, int lc)
{
  int i1, i2, i;
  char tempstring[MAXLEN];
  if(lc < 0)
    {
      for(i=0;i<p->Nlcs;i++)
	{
	  i1 = 0;
	  i2 = 0;
	  while(p->lcnames[i][i1] != '\n' && p->lcnames[i][i1] != '\0')
	    {
	      if(p->lcnames[i][i1] == '/')
		i2 = i1+1;
	      i1++;
	    }
	  sprintf(tempstring,"%s",&p->lcnames[i][i2]);
	  sprintf(p->lcnames[i],"%s",tempstring);
	}
    }
  else
    {
      i1 = 0;
      i2 = 0;
      while(p->lcnames[lc][i1] != '\n' && p->lcnames[lc][i1] != '\0')
	{
	  if(p->lcnames[lc][i1] == '/')
	    i2 = i1+1;
	  i1++;
	}
      sprintf(tempstring,"%s",&p->lcnames[lc][i2]);
      sprintf(p->lcnames[lc],"%s",tempstring);
    }
}

void vAdd_Keyword_To_OutputLC_FitsHeader(ProgramData *p, int lcnum, char *keyname,
					char *comment, int hdutouse,
					int updateexisting,
					int dtype, va_list varlist)
{
  int int_val;
  long long_val;
  double dbl_val;
  char *string_val;
  int i;

#ifdef USECFITSIO
  switch(dtype) {
  case VARTOOLS_TYPE_DOUBLE:
    dbl_val = va_arg(varlist,double);
    break;
  case VARTOOLS_TYPE_INT:
    int_val = va_arg(varlist,int);
    break;
  case VARTOOLS_TYPE_LONG:
    long_val = va_arg(varlist,long);
    break;
  case VARTOOLS_TYPE_STRING:
    string_val = va_arg(varlist,char *);
    break;
  default:
    error(ERR_CODEERROR);
    break;
  }
  if(p->fits_header_adds == NULL) {
    error(ERR_CODEERROR);
  }
  if((p->fits_header_adds[lcnum].N_added_keywords + 1) > p->fits_header_adds[lcnum].size_added_keywords_vec) {
    if(!p->fits_header_adds[lcnum].size_added_keywords_vec) {
      p->fits_header_adds[lcnum].size_added_keywords_vec = p->fits_header_adds[lcnum].N_added_keywords + 1;
      if((p->fits_header_adds[lcnum].hdrterms = (_vartools_header_entry *) malloc(p->fits_header_adds[lcnum].size_added_keywords_vec*sizeof(_vartools_header_entry))) == NULL)
	error(ERR_MEMALLOC);
      for(i=0; i < p->fits_header_adds[lcnum].size_added_keywords_vec; i++) {
	p->fits_header_adds[lcnum].hdrterms[i].datatype = -1;
	p->fits_header_adds[lcnum].hdrterms[i].keyname = NULL;
	p->fits_header_adds[lcnum].hdrterms[i].keyname_veclen = 0;
	p->fits_header_adds[lcnum].hdrterms[i].dbl_val = 0.0;
	p->fits_header_adds[lcnum].hdrterms[i].int_val = 0;
	p->fits_header_adds[lcnum].hdrterms[i].long_val = 0;
	p->fits_header_adds[lcnum].hdrterms[i].string_val = NULL;
	p->fits_header_adds[lcnum].hdrterms[i].string_val_veclen = 0;
	p->fits_header_adds[lcnum].hdrterms[i].comment = NULL;
	p->fits_header_adds[lcnum].hdrterms[i].comment_veclen = 0;
	p->fits_header_adds[lcnum].hdrterms[i].hdutouse = 0;
	p->fits_header_adds[lcnum].hdrterms[i].updateexisting = 1;
      }
    } else {
      if((p->fits_header_adds[lcnum].hdrterms = (_vartools_header_entry *) realloc(p->fits_header_adds[lcnum].hdrterms, (p->fits_header_adds[lcnum].N_added_keywords + 1)*sizeof(_vartools_header_entry))) == NULL)
	error(ERR_MEMALLOC);
      for(i=p->fits_header_adds[lcnum].size_added_keywords_vec; i < p->fits_header_adds[lcnum].N_added_keywords + 1; i++) {
	p->fits_header_adds[lcnum].hdrterms[i].datatype = -1;
	p->fits_header_adds[lcnum].hdrterms[i].keyname = NULL;
	p->fits_header_adds[lcnum].hdrterms[i].keyname_veclen = 0;
	p->fits_header_adds[lcnum].hdrterms[i].dbl_val = 0.0;
	p->fits_header_adds[lcnum].hdrterms[i].int_val = 0;
	p->fits_header_adds[lcnum].hdrterms[i].long_val = 0;
	p->fits_header_adds[lcnum].hdrterms[i].string_val = NULL;
	p->fits_header_adds[lcnum].hdrterms[i].string_val_veclen = 0;
	p->fits_header_adds[lcnum].hdrterms[i].comment = NULL;
	p->fits_header_adds[lcnum].hdrterms[i].comment_veclen = 0;
	p->fits_header_adds[lcnum].hdrterms[i].hdutouse = 0;
	p->fits_header_adds[lcnum].hdrterms[i].updateexisting = 1;
      }
      p->fits_header_adds[lcnum].size_added_keywords_vec = p->fits_header_adds[lcnum].N_added_keywords + 1;
    }
  }
  i = p->fits_header_adds[lcnum].N_added_keywords;
  p->fits_header_adds[lcnum].hdrterms[i].hdutouse = hdutouse;
  p->fits_header_adds[lcnum].hdrterms[i].updateexisting = updateexisting;
  switch(dtype) {
  case VARTOOLS_TYPE_DOUBLE:
    p->fits_header_adds[lcnum].hdrterms[i].datatype = TDOUBLE;
    p->fits_header_adds[lcnum].hdrterms[i].dbl_val = dbl_val;
    break;
  case VARTOOLS_TYPE_INT:
    p->fits_header_adds[lcnum].hdrterms[i].datatype = TINT;
    p->fits_header_adds[lcnum].hdrterms[i].int_val = int_val;
    break;
  case VARTOOLS_TYPE_LONG:
    p->fits_header_adds[lcnum].hdrterms[i].datatype = TLONG;
    p->fits_header_adds[lcnum].hdrterms[i].int_val = long_val;
    break;
  case VARTOOLS_TYPE_STRING:
    p->fits_header_adds[lcnum].hdrterms[i].datatype = TSTRING;
    if(strlen(string_val)+1 > p->fits_header_adds[lcnum].hdrterms[i].string_val_veclen) {
      if(!p->fits_header_adds[lcnum].hdrterms[i].string_val_veclen) {
	p->fits_header_adds[lcnum].hdrterms[i].string_val_veclen = strlen(string_val)+1;
	if((p->fits_header_adds[lcnum].hdrterms[i].string_val = (char *) malloc(p->fits_header_adds[lcnum].hdrterms[i].string_val_veclen*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	p->fits_header_adds[lcnum].hdrterms[i].string_val_veclen = strlen(string_val)+1;
	if((p->fits_header_adds[lcnum].hdrterms[i].string_val = (char *) realloc(p->fits_header_adds[lcnum].hdrterms[i].string_val, p->fits_header_adds[lcnum].hdrterms[i].string_val_veclen*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
      }
    }
    sprintf(p->fits_header_adds[lcnum].hdrterms[i].string_val,"%s",string_val);
    break;
  default:
    error(ERR_CODEERROR);
    break;
  }
  
  if(keyname != NULL) {
    if(strlen(keyname)+1 > p->fits_header_adds[lcnum].hdrterms[i].keyname_veclen) {
      if(!p->fits_header_adds[lcnum].hdrterms[i].keyname_veclen) {
	p->fits_header_adds[lcnum].hdrterms[i].keyname_veclen = strlen(keyname)+1;
	if((p->fits_header_adds[lcnum].hdrterms[i].keyname = (char *) malloc(p->fits_header_adds[lcnum].hdrterms[i].keyname_veclen*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	p->fits_header_adds[lcnum].hdrterms[i].keyname_veclen = strlen(keyname)+1;
	if((p->fits_header_adds[lcnum].hdrterms[i].keyname = (char *) realloc(p->fits_header_adds[lcnum].hdrterms[i].keyname, p->fits_header_adds[lcnum].hdrterms[i].keyname_veclen*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
      }
    }
    sprintf(p->fits_header_adds[lcnum].hdrterms[i].keyname,"%s",keyname);
  }

  if(comment != NULL) {
    if(strlen(comment)+1 > p->fits_header_adds[lcnum].hdrterms[i].comment_veclen) {
      if(!p->fits_header_adds[lcnum].hdrterms[i].comment_veclen) {
	p->fits_header_adds[lcnum].hdrterms[i].comment_veclen = strlen(comment)+1;
	if((p->fits_header_adds[lcnum].hdrterms[i].comment = (char *) malloc(p->fits_header_adds[lcnum].hdrterms[i].comment_veclen*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
      } else {
	p->fits_header_adds[lcnum].hdrterms[i].comment_veclen = strlen(comment)+1;
	if((p->fits_header_adds[lcnum].hdrterms[i].comment = (char *) realloc(p->fits_header_adds[lcnum].hdrterms[i].comment, p->fits_header_adds[lcnum].hdrterms[i].comment_veclen*sizeof(char))) == NULL)
	  error(ERR_MEMALLOC);
      }
    }
    sprintf(p->fits_header_adds[lcnum].hdrterms[i].comment,"%s",comment);
  }

  p->fits_header_adds[lcnum].N_added_keywords++;

#else
  return;
#endif

}

void Add_Keyword_To_OutputLC_FitsHeader(ProgramData *p, int lcnum, char *keyname,
					char *comment, int hdutouse,
					int updateexisting,
					int dtype, ...)
{
  va_list varlist;
  va_start(varlist, dtype);
  vAdd_Keyword_To_OutputLC_FitsHeader(p, lcnum, keyname, comment, hdutouse,
				      updateexisting, dtype, varlist);
  va_end(varlist);
  return;
}

void Run_AddFitsKeyword_Command(ProgramData *p, _AddFitsKeyword *addfitskeyword,
				int lcnum, int lc_name_num) 
{
  char tmpstring[MAXLEN];
  if(addfitskeyword->combinelckeyword) {
    char tmpstring_keyname[9];
    int Nlcgroup, i, i1, i3;
    int *tmpindx = NULL;
    int *tmpindx2 = NULL;
    if(p->NJD[lcnum] <= 0) return;
    if((tmpindx = (int *) malloc(p->NJD[lcnum]*sizeof(int))) == NULL ||
       (tmpindx2 = (int *) malloc(p->NJD[lcnum]*sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < p->NJD[lcnum]; i++) {
      tmpindx[i] = EvaluateVariable_Int(lc_name_num, 
					lcnum, i,
					addfitskeyword->lcnumvar);
      tmpindx2[i] = i;
    }
    sort_generic(p->NJD[lcnum], 0, tmpindx2, 1, VARTOOLS_TYPE_INT, 
		 (void *) tmpindx, 1);
    for(i=0; i < p->NJD[lcnum]; i++) {
      if(i == 0 || (i > 0 && (tmpindx[tmpindx2[i]] != tmpindx[tmpindx2[i-1]]))) {
	i1=0;
	i3=0;
	while(i1 < 8 && addfitskeyword->keyname[i3] != '\0') {
	  if(addfitskeyword->keyname[i3] != '%') {
	    tmpstring_keyname[i1] = addfitskeyword->keyname[i3];
	    i1++;
	    tmpstring_keyname[i1] = '\0';
	    i3++;
	  } else {
	    i3++;
	    if(addfitskeyword->keyname[i3] == 'd') {
	      i3++;
	      snprintf(&(tmpstring_keyname[i1]), (9-i1), "%d", tmpindx[tmpindx2[i]]);
	      i1 = strlen(tmpstring_keyname);
	      tmpstring_keyname[i1] = '\0';
	    } else if(addfitskeyword->keyname[i3] == '%') {
	      tmpstring_keyname[i1] = '%';
	      i1++;
	      tmpstring_keyname[i1] = '\0';
	      i3++;
	    }
	    else
	      error(ERR_INVALIDKEYNAMEFORMAT);
	  }
	}
	if(addfitskeyword->keyval_source == VARTOOLS_SOURCE_FIXED) {
	  switch(addfitskeyword->dtype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, tmpstring_keyname,
					       addfitskeyword->comment_string,
					       addfitskeyword->hdutouse,
					       addfitskeyword->updateexisting,
					       addfitskeyword->dtype,
					       addfitskeyword->dbl_fixval);
	    break;
	  case VARTOOLS_TYPE_INT:
	    Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, tmpstring_keyname,
					       addfitskeyword->comment_string,
					       addfitskeyword->hdutouse,
					       addfitskeyword->updateexisting,
					       addfitskeyword->dtype,
					       addfitskeyword->int_fixval);
	    break;
	  case VARTOOLS_TYPE_LONG:
	    Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, tmpstring_keyname,
					       addfitskeyword->comment_string,
					       addfitskeyword->hdutouse,
					       addfitskeyword->updateexisting,
					       addfitskeyword->dtype,
					       addfitskeyword->long_fixval);
	    break;
	  case VARTOOLS_TYPE_STRING:
	    Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, tmpstring_keyname,
					       addfitskeyword->comment_string,
					       addfitskeyword->hdutouse,
					       addfitskeyword->updateexisting,
					       addfitskeyword->dtype,
					       addfitskeyword->string_fixval);
	    break;
	  default:
	    error(ERR_CODEERROR);
	    break;
	  }
	}
	else if(addfitskeyword->keyval_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
	  switch(addfitskeyword->dtype) {
	  case VARTOOLS_TYPE_DOUBLE:
	    Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, tmpstring_keyname,
					       addfitskeyword->comment_string,
					       addfitskeyword->hdutouse,
					       addfitskeyword->updateexisting,
					       addfitskeyword->dtype,
					       EvaluateVariable_Double(lc_name_num,
								       lcnum, tmpindx2[i],
								       addfitskeyword->keyval_var));
	    break;
	  case VARTOOLS_TYPE_INT:
	    Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, tmpstring_keyname,
					       addfitskeyword->comment_string,
					       addfitskeyword->hdutouse,
					       addfitskeyword->updateexisting,
					       addfitskeyword->dtype,
					       EvaluateVariable_Int(lc_name_num, 
								    lcnum, tmpindx2[i],
								    addfitskeyword->keyval_var));
	    break;
	  case VARTOOLS_TYPE_LONG:
	    Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, tmpstring_keyname,
					       addfitskeyword->comment_string,
					       addfitskeyword->hdutouse,
					       addfitskeyword->updateexisting,
					       addfitskeyword->dtype,
					       EvaluateVariable_Long(lc_name_num, 
								     lcnum, tmpindx2[i],
								     addfitskeyword->keyval_var));
	    break;
	  case VARTOOLS_TYPE_STRING:
	    EvaluateVariable_String(lc_name_num, lcnum, tmpindx2[i],
				    addfitskeyword->keyval_var, tmpstring);
	    Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, tmpstring_keyname,
					       addfitskeyword->comment_string,
					       addfitskeyword->hdutouse,
					       addfitskeyword->updateexisting,
					       addfitskeyword->dtype,
					       tmpstring);
	    break;
	  default:
	    error(ERR_CODEERROR);
	    break;
	  }
	}
      }
    }
    if(tmpindx != NULL) free(tmpindx);
    if(tmpindx2 != NULL) free(tmpindx2);
  } else {
    if(addfitskeyword->keyval_source == VARTOOLS_SOURCE_FIXED) {
      switch(addfitskeyword->dtype) {
      case VARTOOLS_TYPE_DOUBLE:
	Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, addfitskeyword->keyname,
					   addfitskeyword->comment_string,
					   addfitskeyword->hdutouse,
					   addfitskeyword->updateexisting,
					   addfitskeyword->dtype,
					   addfitskeyword->dbl_fixval);
	break;
      case VARTOOLS_TYPE_INT:
	Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, addfitskeyword->keyname,
					   addfitskeyword->comment_string,
					   addfitskeyword->hdutouse,
					   addfitskeyword->updateexisting,
					   addfitskeyword->dtype,
					   addfitskeyword->int_fixval);
	break;
      case VARTOOLS_TYPE_LONG:
	Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, addfitskeyword->keyname,
					   addfitskeyword->comment_string,
					   addfitskeyword->hdutouse,
					   addfitskeyword->updateexisting,
					   addfitskeyword->dtype,
					   addfitskeyword->long_fixval);
	break;
      case VARTOOLS_TYPE_STRING:
	Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, addfitskeyword->keyname,
					   addfitskeyword->comment_string,
					   addfitskeyword->hdutouse,
					   addfitskeyword->updateexisting,
					   addfitskeyword->dtype,
					   addfitskeyword->string_fixval);
	break;
      default:
	error(ERR_CODEERROR);
	break;
      }
    }
    else if(addfitskeyword->keyval_source == VARTOOLS_SOURCE_EXISTINGVARIABLE) {
      switch(addfitskeyword->dtype) {
      case VARTOOLS_TYPE_DOUBLE:
	Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, addfitskeyword->keyname,
					   addfitskeyword->comment_string,
					   addfitskeyword->hdutouse,
					   addfitskeyword->updateexisting,
					   addfitskeyword->dtype,
					   EvaluateVariable_Double(lc_name_num,
								   lcnum, 0,
								   addfitskeyword->keyval_var));
	break;
      case VARTOOLS_TYPE_INT:
	Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, addfitskeyword->keyname,
					   addfitskeyword->comment_string,
					   addfitskeyword->hdutouse,
					   addfitskeyword->updateexisting,
					   addfitskeyword->dtype,
					   EvaluateVariable_Int(lc_name_num, 
								lcnum, 0,
								addfitskeyword->keyval_var));
	break;
      case VARTOOLS_TYPE_LONG:
	Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, addfitskeyword->keyname,
					   addfitskeyword->comment_string,
					   addfitskeyword->hdutouse,
					   addfitskeyword->updateexisting,
					   addfitskeyword->dtype,
					   EvaluateVariable_Long(lc_name_num, 
								 lcnum, 0,
								 addfitskeyword->keyval_var));
	  break;
      case VARTOOLS_TYPE_STRING:
	EvaluateVariable_String(lc_name_num, lcnum, 0,
				addfitskeyword->keyval_var, tmpstring);
	Add_Keyword_To_OutputLC_FitsHeader(p, lcnum, addfitskeyword->keyname,
					   addfitskeyword->comment_string,
					   addfitskeyword->hdutouse,
					   addfitskeyword->updateexisting,
					   addfitskeyword->dtype,
					   tmpstring);
	break;
      default:
	error(ERR_CODEERROR);
	break;
      }
    }
  }
}


void Reset_outlc_fitsheader_additions(ProgramData *p, int lcnum)
{
  if(p->fits_header_adds != NULL) {
    if(p->fits_header_adds[lcnum].N_added_keywords > 0) {
      int i;
      for(i = 0; i < p->fits_header_adds[lcnum].N_added_keywords; i++) {
	p->fits_header_adds[lcnum].hdrterms[i].datatype = -1;
	if(p->fits_header_adds[lcnum].hdrterms[i].keyname_veclen > 0)
	  p->fits_header_adds[lcnum].hdrterms[i].keyname[0] = '\0';
	p->fits_header_adds[lcnum].hdrterms[i].dbl_val = 0.0;
	p->fits_header_adds[lcnum].hdrterms[i].int_val = 0;
	p->fits_header_adds[lcnum].hdrterms[i].long_val = 0;
	p->fits_header_adds[lcnum].hdrterms[i].hdutouse = 0;
	p->fits_header_adds[lcnum].hdrterms[i].updateexisting = 1;
	if(p->fits_header_adds[lcnum].hdrterms[i].string_val_veclen > 0)
	  p->fits_header_adds[lcnum].hdrterms[i].string_val[0] = '\0';
	if(p->fits_header_adds[lcnum].hdrterms[i].comment_veclen > 0)
	  p->fits_header_adds[lcnum].hdrterms[i].comment[0] = '\0';
      }
      p->fits_header_adds[lcnum].N_added_keywords = 0;
    }
  }
}

	
	

