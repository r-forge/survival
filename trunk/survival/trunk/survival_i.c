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
extern void agres1();
extern void agres2();
extern void agsurv1();
extern void agsurv2();
extern void coxfit2();
extern void coxfit_null();
extern void coxhaz2();
extern void coxres1();
extern void coxres2();
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
{SYMBOL(agres1),(long)aghaz2,NULL},
{SYMBOL(agres2),(long)aghaz2,NULL},
{SYMBOL(agsurv1),(long)agsurv1,NULL},
{SYMBOL(agsurv2),(long)agsurv2,NULL},
{SYMBOL(coxfit2),(long)coxfit2,NULL},
{SYMBOL(coxfit_null),(long)coxfit_null,NULL},
{SYMBOL(coxhaz2),(long)coxhaz2,NULL},
{SYMBOL(coxres1),(long)coxhaz2,NULL},
{SYMBOL(coxres2),(long)coxhaz2,NULL},
{SYMBOL(survdiff2),(long)survdiff2,NULL},
{SYMBOL(survexp2),(long)survexp2,NULL},
{SYMBOL(survfit2),(long)survfit2,NULL},
{SYMBOL(survindex2),(long)survindex2,NULL},
{SYMBOL(survreg),(long)survreg,NULL},
{SYMBOL(survreg_g),(long)survreg_g,NULL},
{NULL,0L,NULL} };
