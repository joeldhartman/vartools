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
#define Initialize_function_type void (*)(char *, int *, int *, int *, size_t *)
#define ShowSyntax_function_type void (*)(FILE *)
#define ShowHelp_function_type void (*)(FILE *)
#define ShowExample_function_type void (*)(FILE *)

#ifndef _USERLIB_STRUCT_DEFINE
typedef struct {
  void *ParseCL_function_ptr;
  void *RunCommand_function_ptr;

  void (*ShowSyntax_function)(FILE *);
  void (*ShowHelp_function)(FILE *);
  void (*ShowExample_function)(FILE *);
  void (*Initialize_function)(char *, int *, int *, int *, size_t *);

  char commandname[256];
  void *userdata;
  size_t sizeuserdata;
  int RequireReadAll;
  int RequireSortLC;
  int RequireDistinctTimes;
} _UserLib;

typedef struct {
  int datatype;
  void *dataptr;
  int Ncolumns;
  int output;
  char name[MAXLEN];
} _ParseFixSpecFixcolumnStruct;

typedef struct {
  int doinitialize;
  void *fixptr;
} _ParseParameter_InitializeStruct;

#ifdef _HAVE_PYTHON
typedef struct {
  void *ParsePythonCommand_ptr;
  void *RunPythonCommand_ptr;
  void *InitPythonCommand_ptr;
} _VartoolsPythonLibStruct;
#endif

#define _USERLIB_STRUCT_DEFINE
#endif
