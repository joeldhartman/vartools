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

/* This file contains sorting functions for the program vartools by J. Hartman */

void sort_generic(int N, int isreverse, int *index, int Nms, ...)
/* This function can be used to sort one or more data vectors of various types.
   The syntax is as follows:
      N - number of datapoints in each of the vectors.
      isreverse - 0 = sort from low to high.
                  1 = sort from high to low.
      index - = A array of integer indices to be provided if we don't want to
                transpose the data for at least one of the vectors. Set this
                to NULL if all of the vectors will be transposed. If this is
                not null, it should initially store the values 0 to N-1 in
                order. The values of index will then be transposed such that
                dataN[index[i=0 ... N-1]] will be sorted based on the key.
      Nms - The number of data vectors. This must be >= 1 or else VARTOOLS
            will abort with an error.
      ...
          - For each data vector provide the additional input arguments
            given below. Note that the first vector described will be the
            key used for the sorting. Additional vectors will be re-ordered
            with the key.

            int datatype = the datatype for this vector. Allowed values are
                           VARTOOLS_TYPE_DOUBLE, VARTOOLS_TYPE_STRING,
                           VARTOOLS_TYPE_INT, VARTOOLS_TYPE_FLOAT,
                           VARTOOLS_TYPE_LONG, VARTOOLS_TYPE_CHAR,
                           VARTOOLS_TYPE_SHORT, VARTOOLS_TYPE_USERDEF (for
                               an array of user-defined structures).
            void *dataptr = the data vector (e.g. if the vector is
                           s of type (double *), one would pass (void *)s
                           for this argument).
            int useindex = 0 - the data in the vector will be transposed.
                           1 - the data will not be transposed and instead
                               will be indexed. Note that if this is one,
                               then index above cannot be NULL.

          - If datatype is VARTOOLS_TYPE_STRING the following additional
            argument is needed.

            size_t sizeobj = the size of single string element.

          - If datatype is VARTOOLS_TYPE_USERDEF the following additional
            two arguments are needed.

            size_t sizeobj = the size of a single element.

            int (*comp_func)((void *),(void *)) = a pointer to a function
                used to compare two elements. This function is of the form
                used by the qsort function (e.g. strcmp for strings).
                If this vector is not the key then this can be NULL.
*/
{
  mysort_generic_struct *s;
  int i;
  va_list varlist;
  int datatype;
  void *dataptr;
  comp_func_type comp_func;
  int useindex;
  size_t sizeobj;

  va_start(varlist, Nms);
  if(Nms <= 0) {
    error(ERR_MYSORT_GENERIC_BADCALL);
  }
  if((s = (mysort_generic_struct *) malloc(Nms * sizeof(mysort_generic_struct))) == NULL)
    error(ERR_MEMALLOC);
  for(i=0; i < Nms; i++) {
    datatype = va_arg(varlist,int);
    dataptr = va_arg(varlist,void *);
    useindex = va_arg(varlist,int);
    if(useindex && index == NULL)
      error(ERR_MYSORT_GENERIC_BADCALL);
    if(datatype == VARTOOLS_TYPE_STRING) {
      sizeobj = va_arg(varlist,size_t);
      comp_func = NULL;
    }
    else if(datatype == VARTOOLS_TYPE_USERDEF) {
      sizeobj = va_arg(varlist,size_t);
      comp_func = va_arg(varlist,comp_func_type);
    }
    else {
      sizeobj = 0;
      comp_func = NULL;
    }
    s[i].datatype = datatype;
    s[i].dataptr = dataptr;
    s[i].comp_func = comp_func;
    s[i].useindex = useindex;
    s[i].sizeobj = sizeobj;
  }
  va_end(varlist);

  vsort_generic(N, isreverse, index, Nms, s);

  free(s);
}

void create_mysort_swapstruct(mysort_swapstruct **s)
{
  (*s) = (mysort_swapstruct *) malloc(sizeof(mysort_swapstruct));
  (*s)->sizeswaps = 1024;
  (*s)->i1 = (int *) malloc((*s)->sizeswaps * sizeof(int));
  (*s)->i2 = (int *) malloc((*s)->sizeswaps * sizeof(int));
  (*s)->Nswaps = 0;
}

void destroy_mysort_swapstruct(mysort_swapstruct *s)
{
  free(s->i1);
  free(s->i2);
  free(s);
}

void add_swap(mysort_swapstruct *s, int i1, int i2)
{
  if(s->Nswaps == s->sizeswaps) {
    s->sizeswaps *= 2;
    if((s->i1 = (int *) realloc(s->i1, s->sizeswaps * sizeof(int))) == NULL ||
       (s->i2 = (int *) realloc(s->i2, s->sizeswaps * sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
  }
  s->i1[s->Nswaps] = i1;
  s->i2[s->Nswaps] = i2;
  s->Nswaps++;
}

void do_swaps(mysort_swapstruct *sw, int Nms, mysort_generic_struct *s) {
  int i, j, *i1, *i2;
  double *dblptr;
  double dblval;
  char **stringptr;
  char *stringval = NULL;
  int sizestring = 0;
  float *fltptr;
  float fltval;
  int *intptr;
  int intval;
  long *longptr;
  long longval;
  short *shortptr;
  short shortval;
  char *charptr;
  char charval;
  if(!sw->Nswaps)
    return;
  i1 = sw->i1;
  i2 = sw->i2;
  for(i=0; i < Nms; i++) {
    if(!s[i].useindex) {
      switch(s[i].datatype) {
      case VARTOOLS_TYPE_DOUBLE:
	dblptr = (double *) s[i].dataptr;
	for(j=0; j < sw->Nswaps; j++) {
	  dblval = dblptr[i1[j]]; dblptr[i1[j]] = dblptr[i2[j]];
	  dblptr[i2[j]] = dblval;
	}
	break;
      case VARTOOLS_TYPE_FLOAT:
	fltptr = (float *) s[i].dataptr;
	for(j=0; j < sw->Nswaps; j++) {
	  fltval = fltptr[i1[j]]; fltptr[i1[j]] = fltptr[i2[j]];
	  fltptr[i2[j]] = fltval;
	}
	break;
      case VARTOOLS_TYPE_INT:
	intptr = (int *) s[i].dataptr;
	for(j=0; j < sw->Nswaps; j++) {
	  intval = intptr[i1[j]]; intptr[i1[j]] = intptr[i2[j]];
	  intptr[i2[j]] = intval;
	}
	break;
      case VARTOOLS_TYPE_LONG:
	longptr = (long *) s[i].dataptr;
	for(j=0; j < sw->Nswaps; j++) {
	  longval = longptr[i1[j]]; longptr[i1[j]] = longptr[i2[j]];
	  longptr[i2[j]] = longval;
	}
	break;
      case VARTOOLS_TYPE_SHORT:
	shortptr = (short *) s[i].dataptr;
	for(j=0; j < sw->Nswaps; j++) {
	  shortval = shortptr[i1[j]]; shortptr[i1[j]] = shortptr[i2[j]];
	  shortptr[i2[j]] = shortval;
	}
	break;
      case VARTOOLS_TYPE_CHAR:
	charptr = (char *) s[i].dataptr;
	for(j=0; j < sw->Nswaps; j++) {
	  charval = charptr[i1[j]]; charptr[i1[j]] = charptr[i2[j]];
	  charptr[i2[j]] = charval;
	}
	break;
      case VARTOOLS_TYPE_STRING:
      case VARTOOLS_TYPE_USERDEF:
	stringptr = (char **) s[i].dataptr;
	if(s[i].sizeobj > sizestring) {
	  if(!sizestring) {
	    if((stringval = (char *) malloc(s[i].sizeobj)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  else {
	    if((stringval = (char *) realloc(stringval, s[i].sizeobj)) == NULL)
	      error(ERR_MEMALLOC);
	  }
	  sizestring = s[i].sizeobj;
	}
	for(j=0; j < sw->Nswaps; j++) {
	  memcpy(stringval, stringptr[i1[j]], s[i].sizeobj);
	  memcpy(stringptr[i1[j]], stringptr[i2[j]], s[i].sizeobj);
	  memcpy(stringptr[i2[j]], stringval, s[i].sizeobj);
	}
	break;
      default:
	error(ERR_BADTYPE);
      }
    }
  }
  if(stringval != NULL)
    free(stringval);
}

void vsort_generic(int N, int isreverse, int *index, int Nms, mysort_generic_struct *s)
{
  int i, isstring;
  int *tmpindex;
  mysort_swapstruct *swaps = NULL;

  void sort_index_double(int N, int i1, int *index, double *data, mysort_swapstruct *s);
  void reverse_sort_index_double(int N, int i1, int *index, double *data, mysort_swapstruct *s);
  void sort_index_float(int N, int i1, int *index, float *data, mysort_swapstruct *s);
  void reverse_sort_index_float(int N, int i1, int *index, float *data, mysort_swapstruct *s);
  void sort_index_int(int N, int i1, int *index, int *data, mysort_swapstruct *s);
  void reverse_sort_index_int(int N, int i1, int *index, int *data, mysort_swapstruct *s);
  void sort_index_short(int N, int i1, int *index, short *data, mysort_swapstruct *s);
  void reverse_sort_index_short(int N, int i1, int *index, short *data, mysort_swapstruct *s);
  void sort_index_long(int N, int i1, int *index, long *data, mysort_swapstruct *s);
  void reverse_sort_index_long(int N, int i1, int *index, long *data, mysort_swapstruct *s);
  void sort_index_char(int N, int i1, int *index, char *data, mysort_swapstruct *s);
  void reverse_sort_index_char(int N, int i1, int *index, char *data, mysort_swapstruct *s);
  void sort_index_string(int N, int i1, int *index, char **data, mysort_swapstruct *s);
  void reverse_sort_index_string(int N, int i1, int *index, char **data, mysort_swapstruct *s);
  void sort_index_userdef(int N, int i1, int *index, char **data, mysort_swapstruct *s, comp_func_type comp_func);
  void reverse_sort_index_userdef(int N, int i1, int *index, char **data, mysort_swapstruct *s, comp_func_type comp_func);

  /* Check for bad calls to this function */
  if(Nms < 1) {
    error(ERR_MYSORT_GENERIC_BADCALL);
  }

  if(!N)
    return;

  /* Check to see if the call is a special case, if so
     we can use one of the pre-defined sorting routines */
  if(s[0].datatype == VARTOOLS_TYPE_DOUBLE){
    if(Nms == 1) {
      if(index == NULL) {
	if(!isreverse) {
	  mysort1(N, (double *) s[0].dataptr);
	  return;
	} else {
	  mysort1_rev(N, (double *) s[0].dataptr);
	  return;
	}
      }
      else if(index != NULL) {
	if(s[0].useindex && !isreverse) {
	  mysort2dblint_id(N, (double *) s[0].dataptr, index);
	  return;
	}
      }
    }
    else if(Nms == 2) {
      if(index == NULL) {
	if(s[1].datatype == VARTOOLS_TYPE_DOUBLE) {
	  if(!isreverse) {
	    mysort2(N, (double *) s[0].dataptr, (double *) s[1].dataptr);
	    return;
	  }
	  else {
	    mysort2_rev(N, (double *) s[0].dataptr, (double *) s[1].dataptr);
	    return;
	  }
	}
	else if(s[1].datatype == VARTOOLS_TYPE_INT) {
	  if(!isreverse) {
	    mysort2dblint(N, (double *) s[0].dataptr, (int *) s[1].dataptr);
	    return;
	  }
	}
      }
    }
    else if(Nms == 3) {
      if(index == NULL) {
	if(s[1].datatype == VARTOOLS_TYPE_DOUBLE &&
	   s[2].datatype == VARTOOLS_TYPE_DOUBLE) {
	  if(!isreverse) {
	    mysort3(N, (double *) s[0].dataptr, (double *) s[1].dataptr,
		    (double *) s[2].dataptr);
	    return;
	  } else {
	    mysort3_rev(N, (double *) s[0].dataptr, (double *) s[1].dataptr,
			(double *) s[2].dataptr);
	    return;
	  }
	}
	else if(s[1].datatype == VARTOOLS_TYPE_DOUBLE &&
		s[2].datatype == VARTOOLS_TYPE_INT) {
	  if(!isreverse) {
	    mysort3_int(N, (double *) s[0].dataptr, (double *) s[1].dataptr,
			(int *) s[2].dataptr);
	    return;
	  }
	}
      }
    }
    else if(Nms == 4){
      if(index == NULL) {
	if(s[1].datatype == VARTOOLS_TYPE_DOUBLE &&
	   s[2].datatype == VARTOOLS_TYPE_DOUBLE &&
	   s[3].datatype == VARTOOLS_TYPE_DOUBLE) {
	  if(!isreverse) {
	    mysort4(N, (double *) s[0].dataptr, (double *) s[1].dataptr,
		    (double *) s[2].dataptr, (double *) s[3].dataptr);
	    return;
	  } else {
	    mysort4_rev(N, (double *) s[0].dataptr, (double *) s[1].dataptr,
			(double *) s[2].dataptr, (double *) s[3].dataptr);
	    return;
	  }
	}
      }
    }
  }
  else if(s[0].datatype == VARTOOLS_TYPE_STRING) {
    if(Nms == 1 && index != NULL) {
      if(!isreverse && s[0].useindex) {
	mysortstringint(N, (int)(s[0].sizeobj), (char **) s[0].dataptr, index);
	return;
      }
    }
  }

  /* If we got here then one of the predefined sort functions does not
     fit the template, so use one of the generic sorting functions.
     The idea here is to first sort an array of indices and keep track of
     all swaps that are needed. Then go back and swap the data as needed.
  */
  if(index == NULL) {
    if((tmpindex = (int *) malloc(N * sizeof(int))) == NULL)
      error(ERR_MEMALLOC);
    for(i=0; i < N; i++)
      tmpindex[i] = i;
  }
  else
    tmpindex = index;


  create_mysort_swapstruct(&swaps);

  /* Call the appropriate sorting function based on the data type of the
     key and whether or not the sort is to be reversed */
  switch(s[0].datatype) {
  case VARTOOLS_TYPE_DOUBLE:
    if(!isreverse) {
      sort_index_double(N, 0, tmpindex, (double *) s[0].dataptr, swaps);
    }
    else {
      reverse_sort_index_double(N, 0, tmpindex, (double *) s[0].dataptr, swaps);
    }
    break;
  case VARTOOLS_TYPE_FLOAT:
    if(!isreverse) {
      sort_index_float(N, 0, tmpindex, (float *) s[0].dataptr, swaps);
    }
    else {
      reverse_sort_index_float(N, 0, tmpindex, (float *) s[0].dataptr, swaps);
    }
    break;
  case VARTOOLS_TYPE_INT:
    if(!isreverse) {
      sort_index_int(N, 0, tmpindex, (int *) s[0].dataptr, swaps);
    }
    else {
      reverse_sort_index_int(N, 0, tmpindex, (int *) s[0].dataptr, swaps);
    }
    break;
  case VARTOOLS_TYPE_SHORT:
    if(!isreverse) {
      sort_index_short(N, 0, tmpindex, (short *) s[0].dataptr, swaps);
    }
    else {
      reverse_sort_index_short(N, 0, tmpindex, (short *) s[0].dataptr, swaps);
    }
    break;
  case VARTOOLS_TYPE_LONG:
    if(!isreverse) {
      sort_index_long(N, 0, tmpindex, (long *) s[0].dataptr, swaps);
    }
    else {
      reverse_sort_index_long(N, 0, tmpindex, (long *) s[0].dataptr, swaps);
    }
    break;
  case VARTOOLS_TYPE_CHAR:
    if(!isreverse) {
      sort_index_char(N, 0, tmpindex, (char *) s[0].dataptr, swaps);
    }
    else {
      reverse_sort_index_char(N, 0, tmpindex, (char *) s[0].dataptr, swaps);
    }
    break;
  case VARTOOLS_TYPE_STRING:
    if(!isreverse) {
      sort_index_string(N, 0, tmpindex, (char **) s[0].dataptr, swaps);
    }
    else {
      reverse_sort_index_string(N, 0, tmpindex, (char **) s[0].dataptr, swaps);
    }
    break;
  case VARTOOLS_TYPE_USERDEF:
    if(!isreverse) {
      sort_index_userdef(N, 0, tmpindex, (char **) s[0].dataptr, swaps, s[0].comp_func);
    }
    else {
      reverse_sort_index_userdef(N, 0, tmpindex, (char **) s[0].dataptr, swaps, s[0].comp_func);
    }
    break;
  default:
    error(ERR_BADTYPE);
  }

  /* Carry out the swaps as needed */
  do_swaps(swaps, Nms, s);

  destroy_mysort_swapstruct(swaps);

  if(index == NULL && tmpindex != NULL)
    free(tmpindex);

  return;
}

void sort_index_double(int N, int i1, int *index, double *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  double v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] < v : 0) { }
      while(data[index[--j]] > v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  sort_index_double(i-i1-1, i1, index, data, s);
  sort_index_double(Nmax-i, i, index, data, s);
}

void reverse_sort_index_double(int N, int i1, int *index, double *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  double v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] > v : 0) { }
      while(data[index[--j]] < v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  reverse_sort_index_double(i-i1-1, i1, index, data, s);
  reverse_sort_index_double(Nmax-i, i, index, data, s);
}

void sort_index_float(int N, int i1, int *index, float *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  float v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] < v : 0) { }
      while(data[index[--j]] > v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  sort_index_float(i-i1-1, i1, index, data, s);
  sort_index_float(Nmax-i, i, index, data, s);
}

void reverse_sort_index_float(int N, int i1, int *index, float *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  float v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] > v : 0) { }
      while(data[index[--j]] < v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  reverse_sort_index_float(i-i1-1, i1, index, data, s);
  reverse_sort_index_float(Nmax-i, i, index, data, s);
}

void sort_index_int(int N, int i1, int *index, int *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  int v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] < v : 0) { }
      while(data[index[--j]] > v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  sort_index_int(i-i1-1, i1, index, data, s);
  sort_index_int(Nmax-i, i, index, data, s);
}

void reverse_sort_index_int(int N, int i1, int *index, int *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  int v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] > v : 0) { }
      while(data[index[--j]] < v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  reverse_sort_index_int(i-i1-1, i1, index, data, s);
  reverse_sort_index_int(Nmax-i, i, index, data, s);
}

void sort_index_long(int N, int i1, int *index, long *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  long v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] < v : 0) { }
      while(data[index[--j]] > v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  sort_index_long(i-i1-1, i1, index, data, s);
  sort_index_long(Nmax-i, i, index, data, s);
}

void reverse_sort_index_long(int N, int i1, int *index, long *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  long v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] > v : 0) { }
      while(data[index[--j]] < v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  reverse_sort_index_long(i-i1-1, i1, index, data, s);
  reverse_sort_index_long(Nmax-i, i, index, data, s);
}

void sort_index_short(int N, int i1, int *index, short *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  short v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] < v : 0) { }
      while(data[index[--j]] > v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  sort_index_short(i-i1-1, i1, index, data, s);
  sort_index_short(Nmax-i, i, index, data, s);
}

void reverse_sort_index_short(int N, int i1, int *index, short *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  float v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] > v : 0) { }
      while(data[index[--j]] < v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  reverse_sort_index_short(i-i1-1, i1, index, data, s);
  reverse_sort_index_short(Nmax-i, i, index, data, s);
}

void sort_index_char(int N, int i1, int *index, char *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  char v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] < v : 0) { }
      while(data[index[--j]] > v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  sort_index_char(i-i1-1, i1, index, data, s);
  sort_index_char(Nmax-i, i, index, data, s);
}

void reverse_sort_index_char(int N, int i1, int *index, char *data, mysort_swapstruct *s) {
  int i, j, Nmax, k;
  char v;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  v = data[index[i]];
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? data[index[i]] > v : 0) { }
      while(data[index[--j]] < v) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  reverse_sort_index_char(i-i1-1, i1, index, data, s);
  reverse_sort_index_char(Nmax-i, i, index, data, s);
}

void sort_index_string(int N, int i1, int *index, char **data, mysort_swapstruct *s) {
  int i, j, Nmax, k;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? strcmp(data[index[i]], data[index[i1]]) < 0 : 0) { }
      while(strcmp(data[index[--j]], data[index[i1]]) > 0) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  sort_index_string(i-i1-1, i1, index, data, s);
  sort_index_string(Nmax-i, i, index, data, s);
}

void reverse_sort_index_string(int N, int i1, int *index, char **data, mysort_swapstruct *s) {
  int i, j, Nmax, k;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? strcmp(data[index[i]], data[index[i1]]) > 0 : 0) { }
      while(strcmp(data[index[--j]], data[index[i1]]) < 0) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  reverse_sort_index_string(i-i1-1, i1, index, data, s);
  reverse_sort_index_string(Nmax-i, i, index, data, s);
}

 void sort_index_userdef(int N, int i1, int *index, char **data, mysort_swapstruct *s, comp_func_type comp_func) {
  int i, j, Nmax, k;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? comp_func((void *) data[index[i]], (void *) data[index[i1]]) < 0 : 0) { }
      while(comp_func((void *) data[index[--j]], (void *) data[index[i1]]) > 0) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  sort_index_userdef(i-i1-1, i1, index, data, s, comp_func);
  sort_index_userdef(Nmax-i, i, index, data, s, comp_func);
}

 void reverse_sort_index_userdef(int N, int i1, int *index, char **data, mysort_swapstruct *s, comp_func_type comp_func){
  int i, j, Nmax, k;

  if(N<=1) return;

  i = i1;
  j = i1 + N;
  Nmax = j;
  for(;;)
    {
      while(++i < Nmax ? comp_func((void *)data[index[i]], (void *)data[index[i1]]) > 0 : 0) { }
      while(comp_func((void *)data[index[--j]], (void *)data[index[i1]]) < 0) { }
      if(i >= j) break;
      add_swap(s, i, j);
      k = index[i]; index[i] = index[j]; index[j] = k;
    }
  add_swap(s, i-1, i1);
  k = index[i-1]; index[i-1] = index[i1]; index[i1] = k;
  reverse_sort_index_userdef(i-i1-1, i1, index, data, s, comp_func);
  reverse_sort_index_userdef(Nmax-i, i, index, data, s, comp_func);
}


void mysort2dblint_id(int N, double* data1, int* data2)
{
  /* In this sorting routine only the indexes held in data2 are changed, the
     result will be that data1[data2[i]] is sorted where i goes from 0 to N */
  int i, j;
  double v, t1;
  int t2;

  if(N<=1) return;

  v = data1[data2[0]];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[data2[i]] < v : 0) { }
    while(data1[data2[--j]] > v) { }
    if(i >= j) break;
    /*t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;*/
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
  }
  /*t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;*/
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  mysort2dblint_id(i-1,data1,data2);
  mysort2dblint_id(N-i,data1,data2+i);
}


#define SSORT_SWAP(a,b) {char *t=x[a]; \
			 x[a]=x[b]; x[b]=t; }
#define SSORT_I2C(i) x[i][depth]

void ssort_vecswap(int i, int j, int n, char **x)
{
  while(n-- > 0) {
    SSORT_SWAP(i, j);
    i++;
    j++;
  }
}

void ssort1(char **x, int n, int depth) {
  int a, b, c, d, r, v;
  if(n<= 1)
    return;
  a = rand() % n;
  SSORT_SWAP(0, a);
  v = SSORT_I2C(0);
  a = b = 1;
  c = d = n-1;
  for(;;) {
    while(b <= c && (r = SSORT_I2C(b)-v) <= 0) {
      if(r == 0) {SSORT_SWAP(a, b); a++;}
      b++;
    }
    while(b <= c && (r = SSORT_I2C(c)-v) >= 0) {
      if (r == 0) { SSORT_SWAP(c, d); d--; }
      c--;
    }
    if(b > c) break;
    SSORT_SWAP(b, c);
    b++;
    c++;
  }
  r = (a < (b-a) ? a : (b-a)); ssort_vecswap(0, b-r, r, x);
  r = (d-c < (n-d-1) ? (d-c) : (n-d-1)); ssort_vecswap(b, n-r, r, x);
  r = b-a; ssort1(x, r, depth);
  if (SSORT_I2C(r) != 0)
    ssort1(x+r, a+n-d-1, depth+1);
  r = d-c; ssort1(x + n-r, r, depth);
}

void tmpvecswap(int i, int j, int n, int *x) {
  int t;
  while(n-- > 0) {
    t = x[i]; x[i] = x[j]; x[j] = t;
    i++;
    j++;
  }
}

void tmpmysortstringint(int N, int depth, char **data1, int *data2)
{
  int a, b, c, d, r, v, t;
  if (N <= 1)
    return;
  a = rand() % N;
  t = data2[0];
  data2[0] = data2[a];
  data2[a] = t;
  v = data1[data2[0]][depth];
  a = b = 1;
  c = d = N - 1;
  for(;;) {
    while(b <= c && (r = data1[data2[b]][depth] - v) <= 0) {
      if (r == 0) {t = data2[a]; data2[a] = data2[b]; data2[b] = t; a++;}
      b++;
    }
    while(b <= c && (r = data1[data2[c]][depth] - v) >= 0) {
      if (r == 0) {t = data2[c]; data2[c] = data2[d]; data2[d] = t; d--;}
      c--;
    }
    if(b > c) break;
    t = data2[b]; data2[b] = data2[c]; data2[c] = t;
    b++;
    c--;
  }
  r = (a < (b-a) ? a : (b-a)); tmpvecswap(0, b-r, r, data2);
  r = ((d-c) < (N-d-1) ? (d-c) : (N-d-1)); tmpvecswap(b, N-r, r, data2);
  r = b-a; tmpmysortstringint(r, depth, data1, data2);
  if (data1[data2[r]][depth] != 0)
    tmpmysortstringint(a+N-d-1, depth+1, data1, data2+r);
  r = d - c;
  tmpmysortstringint(r, depth, data1, data2 + N - r);
}

void mysortstringint(int N, int sizestr, char **data1, int *data2)
{
  tmpmysortstringint(N, 0, data1, data2);
}

void mysortstringint_old(int N, int sizestr, char **data1, int *data2)
{
  /* For this sort we only change the indexes held in data2, the strings in data1 do not need to be moved, the output will be such that data1[data2[i]] will be sorted when i goes from 0...N-1 */
  int i, j;
  /*char v[MAXIDSTRINGLENGTH], t1[MAXIDSTRINGLENGTH];*/
  int t2;

  if(N<=1) return;

  /*strncpy(v,data1[0],sizestr);*/
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? strcmp(data1[data2[i]], data1[data2[0]]) < 0 : 0) { }
    while(strcmp(data1[data2[--j]], data1[data2[0]]) > 0) { }
    if(i >= j) break;
    /*    strncpy(t1, data1[i], sizestr);
    if(i != j)
      strncpy(data1[i], data1[j], sizestr);
      strncpy(data1[j], t1, sizestr);*/
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
  }
/*  strncpy(t1, data1[i-1], sizestr);
  if(i-1 != 0)
    strncpy(data1[i-1], data1[0], sizestr);
    strncpy(data1[0], t1, sizestr);*/
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  mysortstringint(i-1,sizestr,data1,data2);
  mysortstringint(N-i,sizestr,data1,data2+i);
}


void mysort4(int N, double* data1, double* data2, double* data3, double *data4)
{
  int i, j;
  double v, t1, t2, t3, t4;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] < v : 0) { }
    while(data1[--j] > v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
    t3 = data3[i]; data3[i] = data3[j]; data3[j] = t3;
    t4 = data4[i]; data4[i] = data4[j]; data4[j] = t4;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  t3 = data3[i-1]; data3[i-1] = data3[0]; data3[0] = t3;
  t4 = data4[i-1]; data4[i-1] = data4[0]; data4[0] = t4;
  mysort4(i-1,data1,data2, data3, data4);
  mysort4(N-i,data1+i,data2+i, data3+i, data4+i);
}


void mysort4_rev(int N, double* data1, double* data2, double* data3, double *data4)
{
  int i, j;
  double v, t1, t2, t3, t4;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] > v : 0) { }
    while(data1[--j] < v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
    t3 = data3[i]; data3[i] = data3[j]; data3[j] = t3;
    t4 = data4[i]; data4[i] = data4[j]; data4[j] = t4;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  t3 = data3[i-1]; data3[i-1] = data3[0]; data3[0] = t3;
  t4 = data4[i-1]; data4[i-1] = data4[0]; data4[0] = t4;
  mysort4_rev(i-1,data1,data2, data3, data4);
  mysort4_rev(N-i,data1+i,data2+i, data3+i, data4+i);
}

void mysort3_int(int N, double* data1, double* data2, int* data3)
{
  int i, j;
  double v, t1, t2;
  int t3;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] < v : 0) { }
    while(data1[--j] > v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
    t3 = data3[i]; data3[i] = data3[j]; data3[j] = t3;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  t3 = data3[i-1]; data3[i-1] = data3[0]; data3[0] = t3;
  mysort3_int(i-1,data1,data2, data3);
  mysort3_int(N-i,data1+i,data2+i, data3+i);
}


void mysort3(int N, double* data1, double* data2, double* data3)
{
  int i, j;
  double v, t1, t2, t3;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] < v : 0) { }
    while(data1[--j] > v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
    t3 = data3[i]; data3[i] = data3[j]; data3[j] = t3;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  t3 = data3[i-1]; data3[i-1] = data3[0]; data3[0] = t3;
  mysort3(i-1,data1,data2, data3);
  mysort3(N-i,data1+i,data2+i, data3+i);
}

void mysort3_rev(int N, double* data1, double* data2, double* data3)
{
  int i, j;
  double v, t1, t2, t3;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] > v : 0) { }
    while(data1[--j] < v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
    t3 = data3[i]; data3[i] = data3[j]; data3[j] = t3;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  t3 = data3[i-1]; data3[i-1] = data3[0]; data3[0] = t3;
  mysort3_rev(i-1,data1,data2, data3);
  mysort3_rev(N-i,data1+i,data2+i, data3+i);
}

void mysort2(int N, double* data1, double* data2)
{
  int i, j;
  double v, t1, t2;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] < v : 0) { }
    while(data1[--j] > v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  mysort2(i-1,data1,data2);
  mysort2(N-i,data1+i,data2+i);
}

void mysort2_rev(int N, double* data1, double* data2)
{
  int i, j;
  double v, t1, t2;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] > v : 0) { }
    while(data1[--j] < v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  mysort2_rev(i-1,data1,data2);
  mysort2_rev(N-i,data1+i,data2+i);
}

void mysort2ptr(int N, int* data1, double*** data2)
{
  int i, j;
  int v, t1;
  double **t2;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] < v : 0) { }
    while(data1[--j] > v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  mysort2ptr(i-1,data1,data2);
  mysort2ptr(N-i,data1+i,data2+i);
}

void mysort3ptrint(int N, int* data1, void*** data2, int *data3)
{
  int i, j;
  int v, t1;
  void **t2;
  int t3;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] < v : 0) { }
    while(data1[--j] > v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
    t3 = data3[i]; data3[i] = data3[j]; data3[j] = t3;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  t3 = data3[i-1]; data3[i-1] = data3[0]; data3[0] = t3;
  mysort3ptrint(i-1,data1,data2,data3);
  mysort3ptrint(N-i,data1+i,data2+i,data3+i);
}

void mysort4ptrint(int N, int* data1, void*** data2, int *data3, int *data4)
{
  int i, j;
  int v, t1;
  void **t2;
  int t3;
  int t4;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] < v : 0) { }
    while(data1[--j] > v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
    t3 = data3[i]; data3[i] = data3[j]; data3[j] = t3;
    t4 = data4[i]; data4[i] = data4[j]; data4[j] = t4;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  t3 = data3[i-1]; data3[i-1] = data3[0]; data3[0] = t3;
  t4 = data4[i-1]; data4[i-1] = data4[0]; data4[0] = t4;
  mysort4ptrint(i-1,data1,data2,data3,data4);
  mysort4ptrint(N-i,data1+i,data2+i,data3+i,data4+i);
}

void mysort1int(int N, int* data1)
{
  int i, j;
  int v, t1;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] < v : 0) { }
    while(data1[--j] > v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  mysort1int(i-1,data1);
  mysort1int(N-i,data1+i);
}


void mysort1(int N, double* data1)
{
  int i, j;
  double v, t1;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] < v : 0) { }
    while(data1[--j] > v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  mysort1(i-1,data1);
  mysort1(N-i,data1+i);
}

void mysort1_rev(int N, double* data1)
{
  int i, j;
  double v, t1;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] > v : 0) { }
    while(data1[--j] < v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  mysort1_rev(i-1,data1);
  mysort1_rev(N-i,data1+i);
}


void mysort2dblint(int N, double* data1, int* data2)
{
  int i, j;
  double v, t1;
  int t2;

  if(N<=1) return;

  v = data1[0];
  i = 0;
  j = N;
  for(;;)
  {
    while(++i < N ? data1[i] < v : 0) { }
    while(data1[--j] > v) { }
    if(i >= j) break;
    t1 = data1[i]; data1[i] = data1[j]; data1[j] = t1;
    t2 = data2[i]; data2[i] = data2[j]; data2[j] = t2;
  }
  t1 = data1[i-1]; data1[i-1] = data1[0]; data1[0] = t1;
  t2 = data2[i-1]; data2[i-1] = data2[0]; data2[0] = t2;
  mysort2dblint(i-1,data1,data2);
  mysort2dblint(N-i,data1+i,data2+i);
}

