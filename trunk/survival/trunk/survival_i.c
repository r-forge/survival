/* SCCS $Id: survival_i.c,v 4.4 1993-01-13 00:26:55 therneau Exp $  */
#include "S.h"
#include "cdefs.h"
#if defined(IRIS4D) || defined(DECSTATION) /* ( */
#define void long
/* Because Mips C compiler doesn't like void to long coercion. */
#endif /* ) */
extern vector *model_frame();
extern void agexact();
extern void agfit2();
extern void agfit_null();
extern void aghaz2();
extern void agres12();
extern void coxdetail();
extern void agsurv1();
extern void agsurv2();
extern void coxfit2();
extern void coxfit_null();
extern void coxhaz2();
extern void coxres12();
extern void coxscho();
extern void survdiff2();
extern void survexp2();
extern void survfit2();
extern void survindex2();
extern void survreg();
extern void survreg_g();
x_h survival_init[]  = {
{SYMBOL(agexact),(long)agexact,NULL},
{SYMBOL(agfit2),(long)agfit2,NULL},
{SYMBOL(agfit_null),(long)agfit_null,NULL},
{SYMBOL(aghaz2),(long)aghaz2,NULL},
{SYMBOL(agres12),(long)agres12,NULL},
{SYMBOL(coxdetail),(long)coxdetail,NULL},
{SYMBOL(agsurv1),(long)agsurv1,NULL},
{SYMBOL(agsurv2),(long)agsurv2,NULL},
{SYMBOL(coxfit2),(long)coxfit2,NULL},
{SYMBOL(coxfit_null),(long)coxfit_null,NULL},
{SYMBOL(coxhaz2),(long)coxhaz2,NULL},
{SYMBOL(coxres12),(long)coxres12,NULL},
{SYMBOL(coxscho),(long)coxscho,NULL},
{SYMBOL(survdiff2),(long)survdiff2,NULL},
{SYMBOL(survexp2),(long)survexp2,NULL},
{SYMBOL(survfit2),(long)survfit2,NULL},
{SYMBOL(survindex2),(long)survindex2,NULL},
{SYMBOL(survreg),(long)survreg,NULL},
{SYMBOL(survreg_g),(long)survreg_g,NULL},
{NULL,0L,NULL} };
