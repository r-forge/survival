/* SCCS $Id: agfit_null.c,v 4.4 1993-01-30 19:34:57 therneau Exp $  */
/*
** Fit a "null" model.  We just need the loglik
**
** Input
**      n       number of subjects
**      method  ==1 for efron
**      start   start times
**      stop    stop times
**      event   =1 if there was an event at time 'stop'
**      offset  the vector of linear predictors
**      strata  is =1 for the last obs of a strata
**
** Output
**      loglik  (Breslow approx)
**
*/
#include <math.h>

void agfit_null(n, method, start, stop, event, offset, strata, loglik)
double  offset[],
	start[],
	stop[],
	loglik[];
long    n[1],
	method[1],
	strata[],
	event[];

    {
    register int i,j,k;
    register double denom;
    double e_denom;
    double temp;
    double time;
    int deaths;
    double itemp;

    loglik[0]=0;
    for (i=0; i<*n; ) {
	if (event[i] ==1) {
	    /*
	    ** compute the sum of weights over the risk set
	    **   and count the deaths
	    */
	    denom =0;
	    e_denom =0;
	    deaths =0;
	    time = stop[i];
	    for (k=i; k<*n; k++) {
		if (start[k] < time) denom += exp(offset[k]);
		if (stop[k]==time && event[k]==1) {
		    deaths ++;
		    e_denom += exp(offset[k]);
		    loglik[0] += offset[k];
		    }
		if (strata[k]==1) break;
		}

	    itemp =0;
	    for (k=i; k<*n && stop[k]==time; k++) {
		if (event[k]==1) {
		    temp = *method * itemp / deaths;
		    loglik[0] -= log(denom - temp*e_denom);
		    itemp++;
		    }
		i++;
		if (strata[k]==1) break;
		}
	    }
	else i++;
	}
    }
