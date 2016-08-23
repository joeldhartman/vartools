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
     required by the function, a pointer to the function, and then
     either 0 if no additional help will be provided for this
     function, or 1 if you will provide some text to describe the
     purpose of the function and its arguments. If you give one, then
     you will need to provide 1 + 2*Narg additional arguments. The
     first is a string giving a brief statement of the purpose of the
     function. Following this you should provide strings giving a
     variable name for each argument, and a brief description of it.*/
  VARTOOLS_RegisterUserFunction(p, "addvals", 2, &addvals,1,"add two numbers","a", "first number to add","b","second number to add");

}

double addvals(double *param) {
  return param[0]+param[1];
}
