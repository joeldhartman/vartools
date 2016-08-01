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
#ifndef _USERFUNC_STRUCT_DEFINE
#define MAXFUNCNAMELENGTH 256
typedef struct {
  char funcname[MAXFUNCNAMELENGTH];
  double (*EvalFunction_ptr)(double *);
  int Nargs;
} _UserFunc;
#define _USERFUNC_STRUCT_DEFINE
#endif

#ifndef _ANALYTIC_USERFUNC_STRUCT_DEFINE
typedef struct {
  char funcname[MAXFUNCNAMELENGTH];
  int Nargs;
  _Variable **input_argvars;
  _Expression *func_expression;
} _AnalyticUserFunc;
#define _ANALYTIC_USERFUNC_STRUCT_DEFINE
#endif
