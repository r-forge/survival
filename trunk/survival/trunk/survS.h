/* SCCS $Id: survS.h,v 5.4 1998-12-22 09:43:12 therneau Exp $
/*
**   The S.h file defines a few things that I need, and hundreds that I don't.
** In particular, on some architectures, it defines a variable "time"
** which of course conflicts with lots of my C-code, 'time' being a natural
** variable name for survival models.
**   Thanks to Brian Ripley for suggesting a machine independent way of
** fixing this.
**
** The S_alloc function changed it's argument list from version 4 to 5, and
**   the ALLOC macro allows me to have common C code for the two versions,
**   with only this file "survS.h" changed.
*/
#define time timexxx
#include "S.h"
#undef time

#define ALLOC(a,b) S_alloc(a,b,S_evaluator)
