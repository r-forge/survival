/* SCCS $Id: agfit_null.c,v 4.2 1992-08-10 13:15:01 grill Exp $  */
/*
** Fit a "null" model.  We just need the loglik and cumulative hazard
**
** Input
**      n       number of subjects
**      start   start times
**      stop    stop times
**      event   =1 if there was an event at time 'stop'
**      offset  the vector of linear predictors
**      strata  is =1 for the last obs of a strata
**
** Output
**      loglik  (Breslow approx)
**      cumhaz  The cumulative hazard experienced by each subject.
**
** Scratch
**      hazard(n)
**
** The martingale residual will be event[i] - score[i]*cumhaz[i]
*/
#include <math.h>

void agfit_null2(n, start, stop, event, offset, strata, loglik, hazard, cumhaz)
double  offset[],
	start[],
	stop[],
	loglik[],
	hazard[],
	cumhaz[];
long    n[1],
	strata[],
	event[];

    {
    register int i,j,k;
    register double denom;
    double temp;
    double time;
    int deaths;

    loglik[0]=0;
    for (i=0; i<*n; ) {
	if (event[i] ==1) {
	    /*
	    ** compute the sum of weights over the risk set
	    **   and count the deaths
	    */
	    denom =0;
	    deaths =0;
	    time = stop[i];
	    for (k=i; k<*n; k++) {
		if (start[k] < time) denom += exp(offset[k]);
		if (stop[k]==time) {
		    deaths += event[k];
		    loglik[0] += offset[k];
		    }
		if (strata[k]==1) break;
		}

	    loglik[0] -= deaths *log(denom);
	    hazard[i] = deaths/denom;
	    for (k=i; k<*n && stop[k]==time; k++) {
		i++;
		if (strata[k]==1) break;
		}
	    }
	else i++;
	}

    j=0;
    for (i=0; i<*n; i++) {
	temp =0;
	for (k=j; k<=i; k++)
	    if (start[i] < stop[k]) temp += hazard[k];
	cumhaz[i] = temp;
	if (strata[i]==1) j=i+1;
	}
    }
