/* SCCS $Id: coxfit_null.c,v 4.1 1992-03-04 16:51:49 therneau Exp $  */
/*
** Special case: fit the "Null" model.  All that is needed are the loglik
**     and the cumulative hazard
**  the input parameters are
**
**       nused        :number of people
**       time(n)      :time of event or censoring for person i
**       status(n)    :status for the ith person    1=dead , 0=censored
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       offset(n)    :the linear predictor for each subject
**
**  returned parameters
**       loglik       : the log-likelihood for the data
**       cumhaz       : the cumulative hazard experienced by each subject
**
**  the data must be sorted by ascending time within strata
*/
#include <math.h>

void coxfit_null( nusedx, time, status, offset, strata,
		  loglik, cumhaz)

long    *nusedx,
	strata[],
	status[];
double  offset[],
	time[];
double  loglik[],  /* returned values */
	cumhaz[];
{
    register int i,j, person;
    int deaths;
    int n;
    double  denom;
    double  hazard, ttime;

    n = *nusedx;

    /*
    ** pass 1- cumhaz will contain the increment in hazard at that time
    */
    for (i=0; i<n; i++) cumhaz[i]=0;
    strata[n-1] =1;  /* just in case */
    for (person=n-1; person>=0; ) {
	if (strata[person] == 1) {
	    denom = 0;
	    }

	denom += exp(offset[person]);
	ttime = time[person];
	deaths= status[person];
	for (i=person-1; i >=0 && time[i]==ttime && strata[i]==0; i--) {
	    denom += exp(offset[i]);
	    deaths += status[i];
	    }

	for (j=person; j>i; j--) {
	    if (status[j]==1) loglik[0] += offset[j] - log(denom);
	    }

	cumhaz[j+1] = deaths/denom;
	person=j;
	}

    /*
    ** pass2, do a cumulative sum
    */
    hazard=0;
    for (i=0; i<n; i++) {
	hazard += cumhaz[i];
	cumhaz[i] = exp(offset[i])* hazard;
	if (strata[i]==1) hazard=0;
	}
    }
