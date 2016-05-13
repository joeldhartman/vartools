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
#ifndef VARTOOLS_VECTORTYPE_CONSTANT
#include "analytic.h"
#endif

#ifndef _IFSTRUCTAREDEFINED
/* Structure used to describe an if, elif, else, fi structure. */
typedef struct {
  _Expression **expressions;
  int nterms;
  int sizearray;
  char *wasfoundtrue;
} _IfStruct;

typedef struct {
  int sizearray;
  _IfStruct **IfPtrs;
  int curpos;
  char *istrue;
} _IfStack;
#define _IFSTRUCTAREDEFINED
#endif
