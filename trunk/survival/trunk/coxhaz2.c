/*  SCCS $Id: coxhaz2.c,v 5.1 1998-08-30 14:52:41 therneau Exp $
/*
** Compute the hazard and cumulative hazard functions.
**
** Input
**      n       number of subjects
**      score   the vector of subject scores, i.e., exp(beta*z)
**      strata  is =1 for the last obs of a strata
**      mark    carried forward from the coxfit routine
**
** Output
**      hazard  for each subject, the increment in the cumulative hazard
**                 computed at that subject's observation time.  If two
**                 subjects in the same strata have a tied time, then the
**                 hazard is set to 0 for all but the first of the ties.
**      cumhaz  The cumulative hazard within each strata.  If there are no
**                 strata then cumhaz = cumsum(hazard).
**
** The martingale residual will be status[i] - score[i]*cumhaz[i]
*/
#include <stdio.h>
#include "survproto.h"

void coxhaz2(long   *n,      double *score,   long   *mark, 
	     long   *strata, double *hazard,  double *cumhaz)
    {
    register int i;
    register double temp;

    temp=0;
    for (i= *n-1; i>=0; i--) {
	if (strata[i]==1) temp=0;
	temp += score[i];
	score[i] = temp;
	}

    temp=0;
    for (i=0; i<*n; i++) {
	hazard[i] = mark[i]/score[i];
	temp += hazard[i];
	cumhaz[i] = temp;
	if (strata[i]==1) temp=0;
	}
    }
