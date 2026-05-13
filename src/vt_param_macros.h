/*
 * vt_param_macros.h — macros for per-LC variable/expression parameter support
 *
 * Three macro families are provided:
 *
 *   VT_PARAM_COMPANIONS(name)
 *     Declares the three companion fields that must be added ALONGSIDE the
 *     existing fixed-value field in a command struct (commands.h).  The
 *     existing field (e.g. "double period;") is left untouched.
 *
 *   VT_PARSE_DOUBLE(cmdptr, field, argv, i)   [and _INT, _FLOAT, _LONG]
 *     Parses one command-line token.  Recognises the "var varname" and
 *     "expr expression" keywords in addition to a plain numeric literal.
 *     Uses c[] and cn from the enclosing scope (as in parsecommandline.c).
 *     Increments i as needed.
 *
 *   VT_EVAL_DOUBLE(cmdptr, field, lcindex, threadindex)   [and _INT, _FLOAT, _LONG]
 *     Evaluates the parameter for light curve lcindex.  Returns the fixed
 *     value when source==VARTOOLS_SOURCE_FIXED, looks up the per-star
 *     variable when source==VARTOOLS_SOURCE_EXISTINGVARIABLE, and evaluates
 *     the expression when source==VARTOOLS_SOURCE_EVALEXPRESSION.
 *
 * All constants (VARTOOLS_SOURCE_*, VARTOOLS_TYPE_*, VARTOOLS_VECTORTYPE_*)
 * are defined in programdata.h and analytic.h, which must be included before
 * this header.
 *
 * Usage example (commands.h):
 *
 *   typedef struct {
 *     double period;              // existing fixed-value field — keep as-is
 *     VT_PARAM_COMPANIONS(period); // adds period_source, period_var, period_expr
 *     int    nharm;
 *     VT_PARAM_COMPANIONS(nharm);
 *   } _Killharm;
 *
 * Usage example (parsecommandline.c):
 *
 *   VT_PARSE_DOUBLE(c[cn].Killharm, period, argv, i);
 *
 * Usage example (killharm.c execution function):
 *
 *   double p = VT_EVAL_DOUBLE(kh, period, lcindex, threadindex);
 */

#ifndef VT_PARAM_MACROS_H
#define VT_PARAM_MACROS_H

/* -----------------------------------------------------------------------
 * Struct companion fields (add after the existing fixed-value field)
 * ----------------------------------------------------------------------- */

/* Three companion fields for a double, float, int, or long parameter. */
#define VT_PARAM_COMPANIONS(name) \
    int         name##_source; \
    _Variable  *name##_var; \
    _Expression *name##_expr

/* Initialise companions to the fixed-value defaults (use after malloc). */
#define VT_INIT_PARAM(ptr, name) \
    do { (ptr)->name##_source = VARTOOLS_SOURCE_FIXED; \
         (ptr)->name##_var    = NULL; \
         (ptr)->name##_expr   = NULL; } while (0)

/* -----------------------------------------------------------------------
 * Parse-time macros (for use in parsecommandline.c where c[] and cn are
 * in scope)
 * ----------------------------------------------------------------------- */

/* _CS variants: same as above but take an explicit Command * (cs) instead
 * of using c[cn] from scope.  Use these inside external parse functions
 * (e.g. ParseWWZCommand, ParseNonlinfitCommand) that receive a Command *. */

#define VT_PARSE_DOUBLE_CS(cs, cmdptr, field, argv, i_) \
    do { \
        if (!strcmp((argv)[(i_)], "var")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EXISTINGVARIABLE; \
            (i_)++; \
            parse_setparam_existingvariable((cs), (argv)[(i_)], \
                &((cmdptr)->field##_var), \
                VARTOOLS_VECTORTYPE_PERSTARDATA, VARTOOLS_TYPE_DOUBLE); \
        } else if (!strcmp((argv)[(i_)], "expr")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EVALEXPRESSION; \
            (i_)++; \
            parse_setparam_expr((cs), (argv)[(i_)], &((cmdptr)->field##_expr)); \
        } else { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_FIXED; \
            (cmdptr)->field = atof((argv)[(i_)]); \
        } \
    } while (0)

#define VT_PARSE_INT_CS(cs, cmdptr, field, argv, i_) \
    do { \
        if (!strcmp((argv)[(i_)], "var")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EXISTINGVARIABLE; \
            (i_)++; \
            parse_setparam_existingvariable((cs), (argv)[(i_)], \
                &((cmdptr)->field##_var), \
                VARTOOLS_VECTORTYPE_PERSTARDATA, VARTOOLS_TYPE_INT); \
        } else if (!strcmp((argv)[(i_)], "expr")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EVALEXPRESSION; \
            (i_)++; \
            parse_setparam_expr((cs), (argv)[(i_)], &((cmdptr)->field##_expr)); \
        } else { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_FIXED; \
            (cmdptr)->field = atoi((argv)[(i_)]); \
        } \
    } while (0)

#define VT_PARSE_LONG_CS(cs, cmdptr, field, argv, i_) \
    do { \
        if (!strcmp((argv)[(i_)], "var")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EXISTINGVARIABLE; \
            (i_)++; \
            parse_setparam_existingvariable((cs), (argv)[(i_)], \
                &((cmdptr)->field##_var), \
                VARTOOLS_VECTORTYPE_PERSTARDATA, VARTOOLS_TYPE_LONG); \
        } else if (!strcmp((argv)[(i_)], "expr")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EVALEXPRESSION; \
            (i_)++; \
            parse_setparam_expr((cs), (argv)[(i_)], &((cmdptr)->field##_expr)); \
        } else { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_FIXED; \
            (cmdptr)->field = atol((argv)[(i_)]); \
        } \
    } while (0)

#define VT_PARSE_DOUBLE(cmdptr, field, argv, i_) \
    do { \
        if (!strcmp((argv)[(i_)], "var")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EXISTINGVARIABLE; \
            (i_)++; \
            parse_setparam_existingvariable(&(c[cn]), (argv)[(i_)], \
                &((cmdptr)->field##_var), \
                VARTOOLS_VECTORTYPE_PERSTARDATA, VARTOOLS_TYPE_DOUBLE); \
        } else if (!strcmp((argv)[(i_)], "expr")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EVALEXPRESSION; \
            (i_)++; \
            parse_setparam_expr(&(c[cn]), (argv)[(i_)], &((cmdptr)->field##_expr)); \
        } else { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_FIXED; \
            (cmdptr)->field = atof((argv)[(i_)]); \
        } \
    } while (0)

#define VT_PARSE_INT(cmdptr, field, argv, i_) \
    do { \
        if (!strcmp((argv)[(i_)], "var")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EXISTINGVARIABLE; \
            (i_)++; \
            parse_setparam_existingvariable(&(c[cn]), (argv)[(i_)], \
                &((cmdptr)->field##_var), \
                VARTOOLS_VECTORTYPE_PERSTARDATA, VARTOOLS_TYPE_INT); \
        } else if (!strcmp((argv)[(i_)], "expr")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EVALEXPRESSION; \
            (i_)++; \
            parse_setparam_expr(&(c[cn]), (argv)[(i_)], &((cmdptr)->field##_expr)); \
        } else { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_FIXED; \
            (cmdptr)->field = atoi((argv)[(i_)]); \
        } \
    } while (0)

#define VT_PARSE_FLOAT(cmdptr, field, argv, i_) \
    do { \
        if (!strcmp((argv)[(i_)], "var")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EXISTINGVARIABLE; \
            (i_)++; \
            parse_setparam_existingvariable(&(c[cn]), (argv)[(i_)], \
                &((cmdptr)->field##_var), \
                VARTOOLS_VECTORTYPE_PERSTARDATA, VARTOOLS_TYPE_FLOAT); \
        } else if (!strcmp((argv)[(i_)], "expr")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EVALEXPRESSION; \
            (i_)++; \
            parse_setparam_expr(&(c[cn]), (argv)[(i_)], &((cmdptr)->field##_expr)); \
        } else { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_FIXED; \
            (cmdptr)->field = (float)atof((argv)[(i_)]); \
        } \
    } while (0)

#define VT_PARSE_LONG(cmdptr, field, argv, i_) \
    do { \
        if (!strcmp((argv)[(i_)], "var")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EXISTINGVARIABLE; \
            (i_)++; \
            parse_setparam_existingvariable(&(c[cn]), (argv)[(i_)], \
                &((cmdptr)->field##_var), \
                VARTOOLS_VECTORTYPE_PERSTARDATA, VARTOOLS_TYPE_LONG); \
        } else if (!strcmp((argv)[(i_)], "expr")) { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_EVALEXPRESSION; \
            (i_)++; \
            parse_setparam_expr(&(c[cn]), (argv)[(i_)], &((cmdptr)->field##_expr)); \
        } else { \
            (cmdptr)->field##_source = VARTOOLS_SOURCE_FIXED; \
            (cmdptr)->field = atol((argv)[(i_)]); \
        } \
    } while (0)

/* -----------------------------------------------------------------------
 * Evaluation macros (for use in command execution functions)
 *
 * lcindex     : 0-based light-curve index in the input list
 * threadindex : thread index (pass lcindex when not using threads)
 * ----------------------------------------------------------------------- */

#define VT_EVAL_DOUBLE(cmdptr, field, lcindex, threadindex) \
    ((cmdptr)->field##_source == VARTOOLS_SOURCE_EVALEXPRESSION \
        ? EvaluateExpression((lcindex), (threadindex), 0, (cmdptr)->field##_expr) \
        : (cmdptr)->field##_source == VARTOOLS_SOURCE_EXISTINGVARIABLE \
        ? EvaluateVariable_Double((lcindex), (threadindex), 0, (cmdptr)->field##_var) \
        : (double)(cmdptr)->field)

#define VT_EVAL_INT(cmdptr, field, lcindex, threadindex) \
    ((cmdptr)->field##_source == VARTOOLS_SOURCE_EVALEXPRESSION \
        ? (int)EvaluateExpression((lcindex), (threadindex), 0, (cmdptr)->field##_expr) \
        : (cmdptr)->field##_source == VARTOOLS_SOURCE_EXISTINGVARIABLE \
        ? EvaluateVariable_Int((lcindex), (threadindex), 0, (cmdptr)->field##_var) \
        : (cmdptr)->field)

#define VT_EVAL_FLOAT(cmdptr, field, lcindex, threadindex) \
    ((cmdptr)->field##_source == VARTOOLS_SOURCE_EVALEXPRESSION \
        ? (float)EvaluateExpression((lcindex), (threadindex), 0, (cmdptr)->field##_expr) \
        : (cmdptr)->field##_source == VARTOOLS_SOURCE_EXISTINGVARIABLE \
        ? EvaluateVariable_Float((lcindex), (threadindex), 0, (cmdptr)->field##_var) \
        : (cmdptr)->field)

#define VT_EVAL_LONG(cmdptr, field, lcindex, threadindex) \
    ((cmdptr)->field##_source == VARTOOLS_SOURCE_EVALEXPRESSION \
        ? (long)EvaluateExpression((lcindex), (threadindex), 0, (cmdptr)->field##_expr) \
        : (cmdptr)->field##_source == VARTOOLS_SOURCE_EXISTINGVARIABLE \
        ? EvaluateVariable_Long((lcindex), (threadindex), 0, (cmdptr)->field##_var) \
        : (cmdptr)->field)

/* -----------------------------------------------------------------------
 * PERTYPE-based var/expr evaluation (for commands using the PERTYPE_*
 * enum system rather than VARTOOLS_SOURCE_* companions).
 *
 * These evaluate a field that has separate _var and _expr pointer fields
 * and a PERTYPE-valued type field.  Use after the existing PERTYPE_FIX
 * case in the runtime switch.
 * ----------------------------------------------------------------------- */

#define VT_PERTYPE_EVAL_DOUBLE(cmdptr, typefield, field, lcindex, threadindex) \
    ((typefield) == PERTYPE_EXPR \
        ? EvaluateExpression((lcindex), (threadindex), 0, (cmdptr)->field##_expr) \
        : (typefield) == PERTYPE_VAR \
        ? EvaluateVariable_Double((lcindex), (threadindex), 0, (cmdptr)->field##_var) \
        : (double)(cmdptr)->field)

#endif /* VT_PARAM_MACROS_H */
