/* SCCS $Id: survS.h,v 5.2 1998-09-01 09:46:28 therneau Exp $
/*
** The next line is needed on Sun: S.h includes time.h, which defines
**   the variable "time", a variable name that I use often.
** But on Linux, you must leave the line out: time.h does not define "time",
**   and "select.h" (also pulled in by S.h) depends on time.h.
** Real solution: offer the user something other than S.h, which defines
**   what he needs without including the zillion .h files a compile of S needs!
*/

/* #if (defined(SUNOS))   ok- this didn't work, how do we do it? */
#define _TIME_H   /* a hack to stop inclusion of time.h */
/*  #endif  */

#include "S.h"
#define ALLOC(a,b) S_alloc(a,b,S_evaluator)
