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
#include "../config.h"
#ifndef DYNAMICLIB
#define DYNAMICLIB 1
#endif


#ifndef _ANALYTIC_HEADER_INCLUDE

#define VARTOOLS_VECTORTYPE_CONSTANT 0
#define VARTOOLS_VECTORTYPE_SCALAR 1
#define VARTOOLS_VECTORTYPE_INLIST 2
#define VARTOOLS_VECTORTYPE_LC 3
#define VARTOOLS_VECTORTYPE_OUTCOLUMN 4
#define VARTOOLS_VECTORTYPE_PERSTARDATA 5
#define VARTOOLS_VECTORTYPE_ANY 6

#define VARTOOLS_OPERANDTYPE_CONSTANT 0
#define VARTOOLS_OPERANDTYPE_VARIABLE 1
#define VARTOOLS_OPERANDTYPE_EXPRESSION 2
#define VARTOOLS_OPERANDTYPE_FUNCTION 3
#define VARTOOLS_OPERANDTYPE_ITERATORNR 4
#define VARTOOLS_OPERANDTYPE_ITERATORNF 5
#define VARTOOLS_OPERANDTYPE_ARRAYINDEX 6

#define VARTOOLS_OPERATORTYPE_CONSTANT 0
#define VARTOOLS_OPERATORTYPE_ADD 1
#define VARTOOLS_OPERATORTYPE_SUBTRACT 2
#define VARTOOLS_OPERATORTYPE_MULTIPLY 3
#define VARTOOLS_OPERATORTYPE_DIVIDE 4
#define VARTOOLS_OPERATORTYPE_MODULO 5
#define VARTOOLS_OPERATORTYPE_POWER 6
#define VARTOOLS_OPERATORTYPE_GREATERTHAN 7
#define VARTOOLS_OPERATORTYPE_GREATERTHANEQUAL 8
#define VARTOOLS_OPERATORTYPE_LESSTHAN 9
#define VARTOOLS_OPERATORTYPE_LESSTHANEQUAL 10
#define VARTOOLS_OPERATORTYPE_ISEQUAL 11
#define VARTOOLS_OPERATORTYPE_NOTEQUAL 12
#define VARTOOLS_OPERATORTYPE_NOT 13
#define VARTOOLS_OPERATORTYPE_AND 14
#define VARTOOLS_OPERATORTYPE_OR 15
#define VARTOOLS_OPERATORTYPE_BITWISECOMPLEMENT 16
#define VARTOOLS_OPERATORTYPE_BITWISEAND 17
#define VARTOOLS_OPERATORTYPE_BITWISEOR 18

#define VARTOOLS_FUNCTIONCALL_ARRAYINDEXEVAL -1;
#define VARTOOLS_FUNCTIONCALL_EXP 0
#define VARTOOLS_FUNCTIONCALL_LOG 1
#define VARTOOLS_FUNCTIONCALL_LOG10 2
#define VARTOOLS_FUNCTIONCALL_SQRT 3
#define VARTOOLS_FUNCTIONCALL_ABS 4
#define VARTOOLS_FUNCTIONCALL_MAX 5
#define VARTOOLS_FUNCTIONCALL_MIN 6
#define VARTOOLS_FUNCTIONCALL_HYPOT 7
#define VARTOOLS_FUNCTIONCALL_SIN 8
#define VARTOOLS_FUNCTIONCALL_SINDEGR 9
#define VARTOOLS_FUNCTIONCALL_COS 10
#define VARTOOLS_FUNCTIONCALL_COSDEGR 11
#define VARTOOLS_FUNCTIONCALL_TAN 12
#define VARTOOLS_FUNCTIONCALL_TANDEGR 13
#define VARTOOLS_FUNCTIONCALL_ASIN 14
#define VARTOOLS_FUNCTIONCALL_ASINDEGR 15
#define VARTOOLS_FUNCTIONCALL_ACOS 16
#define VARTOOLS_FUNCTIONCALL_ACOSDEGR 17
#define VARTOOLS_FUNCTIONCALL_ATAN2 18
#define VARTOOLS_FUNCTIONCALL_ATAN2DEGR 19
#define VARTOOLS_FUNCTIONCALL_CEIL 20
#define VARTOOLS_FUNCTIONCALL_FLOOR 21
#define VARTOOLS_FUNCTIONCALL_COSH 22
#define VARTOOLS_FUNCTIONCALL_SINH 23
#define VARTOOLS_FUNCTIONCALL_TANH 24
#define VARTOOLS_FUNCTIONCALL_ERF 25
#define VARTOOLS_FUNCTIONCALL_ERFC 26
#define VARTOOLS_FUNCTIONCALL_LGAMMA 27
#define VARTOOLS_FUNCTIONCALL_GAMMA 28
#define VARTOOLS_FUNCTIONCALL_ROUND 29
#define VARTOOLS_FUNCTIONCALL_THETA 30
#define VARTOOLS_FUNCTIONCALL_ACOSH 31
#define VARTOOLS_FUNCTIONCALL_ASINH 32
#define VARTOOLS_FUNCTIONCALL_ATANH 33
#define VARTOOLS_FUNCTIONCALL_RAND 34
#define VARTOOLS_FUNCTIONCALL_GAUSS 35
#define VARTOOLS_FUNCTIONCALL_LEN 36
#define VARTOOLS_FUNCTIONCALL_USERFUNC 37
#define VARTOOLS_FUNCTIONCALL_ISNAN 38

#define VARTOOLS_EXPRESSIONCOMMAND_INDEXTYPE_NOINDEX 0
#define VARTOOLS_EXPRESSIONCOMMAND_INDEXTYPE_SINGLEINDEX 1
#define VARTOOLS_EXPRESSIONCOMMAND_INDEXTYPE_INDEXRANGE 2
#define VARTOOLS_EXPRESSIONCOMMAND_INDEXTYPE_VECTOREXPRESSION 3

typedef struct {
  char *varname;
  char datatype;
  void *dataptr;
  char vectortype;
  OutColumn *outc;
  int vectorlength;
  void *vectorlengthptr;
} _Variable;

typedef struct {
  char op1type;
  char op2type;
  char operatortype;
  void *op1_expression;
  void *op2_expression;
  _Variable *op1_variable;
  _Variable *op2_variable;
  void *op1_functioncall;
  void *op2_functioncall;
  double op1_constant;
  double op2_constant;
} _Expression;

#ifdef DYNAMICLIB
#ifndef _USERLIB_STRUCT_DEFINE
#include "userfunc.h"
#endif
#endif

typedef struct {
  int functionid;
  int Nexpr;
#ifdef DYNAMICLIB
  _UserFunc *UserFunc;
  _AnalyticUserFunc *AnalyticUserFunc;
#endif
  _Expression **arguments;
} _FunctionCall;

#define _ANALYTIC_HEADER_INCLUDE
#endif
