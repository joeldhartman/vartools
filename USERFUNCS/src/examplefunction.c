#include "../../src/vartools.h"

/* This is the source code for a user-defined function to be used with
   vartools.

   This library defines a number of functions that can be included in
   vartools analytic expressions.
*/

void examplefunction_Initialize(ProgramData *p) 
/* This function is used to initialize the library. Every library loadable
   with the name $LIBNAME.so that is loaded with the VARTOOLS -F option must 
   contain a function with the name $LIBNAME_Initialize which is used to
   register the new analytic functions with VARTOOLS.
*/
{
  double addvals(double *);

  /* Use the VARTOOLS_RegisterUserFunction procedure to register each
     function that this library
     provides. VARTOOLS_RegisterUserFunction takes as input a pointer
     to the ProgramData structure p, the name of the function to use
     in the analytic expression evaluator, the number of arguments
     required by the function, and a pointer to the function */
  VARTOOLS_RegisterUserFunction(p, "addvals", 2, &addvals);

}

double addvals(double *param) {
  return param[0]+param[1];
}
