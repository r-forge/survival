/* SCCS  $Id: agmart.c,v 4.2 1993-01-30 19:38:50 therneau Exp $
/*
** Compute the martingale residual for a counting process model
**
** Input
**      n       number of subjects
**      method  will be ==1 for the Efron method
**      start
**      stop    vector of (start, stop] times for the subjects
**      event   vector of status values
**      score   the vector of subject scores, i.e., exp(beta*z)
**      strata  is =1 for the last obs of a strata
**
** Output
**      resid   martingale residual
**
** The martingale residual is more of a nuisance for the Efron method
**
**  Data must be sorted by time within strata, with events first
*/
#include <stdio.h>

void agmart(n, method, start, stop, event, score, strata, resid)
double  score[],
	resid[],
	start[],
	stop[];
long    n[1],
	method[1],
	event[],
	strata[];

    {
    register int i,k;
    double deaths, denom, e_denom;
    double hazard, e_hazard;
    double temp, time;
    int nused;
    int person;

    nused = *n;
    strata[nused-1] =1;  /* Failsafe */

    for (i=0; i<nused; i++)  resid[i]=event[i];
    for (person=0; person<nused;) {
	if (event[person]==0) person++;
	else {
	    denom =0;
	    e_denom =0;
	    time = stop[person];
	    deaths=0;
	    for (k=person; k<nused; k++) {
		if (start[k] < time) {
		    denom += score[k];
		    if (stop[k]==time && event[k]==1) {
			deaths++;
			e_denom += score[k];
			}
		     }
		if (strata[k]==1) break;
		}

	    /*
	    ** Do "expected" for the risk set
	    */
	    hazard =0;
	    e_hazard=0;
	    for (k=0; k<deaths; k++) {
		temp = *method *(k/deaths);
		hazard += 1/(denom - temp*e_denom);
		e_hazard += (1-temp)/(denom - temp*e_denom);
		}
	    for (k=person; k<nused; k++) {
		if (start[k] < time) {
		    if (stop[k]==time && event[k]==1)
			    resid[k] -= score[k]*e_hazard;
		    else    resid[k] -= score[k]*hazard;
		    }
		if (stop[k]==time) person++;
		if (strata[k]==1) break;
		}
	    }
	}
    }
