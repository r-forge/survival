/* SCCS $Id: survS.h,v 5.1 1998-08-30 14:46:35 therneau Exp $
/*
** some macros
*/
#define _TIME_H   /* a hack to stop inclusion of time.h */
#include "S.h"
#define ALLOC(a,b) S_alloc(a,b,S_evaluator)
