/* SCCS  $Id: coxmart.c,v 4.1 1993-01-21 09:40:16 therneau Exp $      */
/*
** Compute the martingale residual for a Cox model
**
** Input
**      n       number of subjects
**      method  will be ==1 for the Efron method
**      time    vector of times
**      status  vector of status values
**      score   the vector of subject scores, i.e., exp(beta*z)
**      strata  is =1 for the last obs of a strata
**      mark    carried forward from the coxfit routine
**
** Output
**      expected the expected number of events for the subject
**
** The martingale residual is more of a nuisance for the Efron method
**
*/
#include <stdio.h>

void coxmart(n, method, time, status, strata, score, expect)
double  score[],
	expect[],
	time[];
long    n[1],
	method[1],
	status[],
	strata[];

    {
    register int i,j;
    int lastone;
    double deaths, denom, e_denom;
    double hazard;
    double temp;
    double downwt;

    strata[*n-1] =1;  /* Failsafe */

    /* Pass 1-- store the risk denominator in 'expect' */
    for (i= *n -1; i>=0; i--) {
	if (strata[i]==1) denom =0;
	denom += score[i];
	if (deaths>0 && (i==0 || strata[i-1]==1 ||  time[i-1]!=time[i]))
		expect[i] = denom;
	}

    /* Pass 2-- now do the work */
    deaths=0;
    e_denom=0;
    hazard =0;
    lastone = 0;
    for (i= 0; i<*n; i++) {
	if (expect[i]!=0) denom = expect[i];
	expect[i] = status[i];
	deaths += status[i];
	e_denom += score[i]*status[i];
	if (strata[i]==1 ||  time[i+1]!=time[i]) {
	    /*last subject of a set of tied times */
	    if (deaths<2 || *method==0) {
		hazard += deaths/denom;
		for (j=lastone; j<=i; j++) {
		    expect[j] -= score[j]*hazard;
		    }
		}
	    else {
		temp = hazard;
		for (j=0; j<deaths; j++) {
		    downwt = j /deaths;
		    hazard +=  1/(denom - e_denom* downwt);
		    temp   +=  (1-downwt)/(denom - e_denom* downwt);
		    }
		for (j=lastone; j<=i; j++) {
		    if (status[j]==0) expect[j] = -score[j]*hazard;
		    else  expect[j] -=  score[j]* temp;
		    }
		}
	    lastone =i +1;
	    deaths =0;
	    e_denom =0;
	    }
	}
    }
