/* SCCS $Id: survfit2.c,v 4.7 1993-05-03 14:39:12 therneau Exp $
/*
** Fit the survival curve
**  Input
**    n=# of subjects
**    y[ny,n]    - matrix of time and status values
**    ny   - number of columns of y
**    wt[n] - vector of case weights
**    strata[n] - ==1 at the last obs of each strata
**    method- 1= km  2= fleming-harrington
**    error  -1= Greenwood, 2=Tsiatis
**    mark[n], risksum[n], wtsum[n] -- work arrays
** Output
**    surv  - the survival
**    varh  - the variance of the hazard function
**    nsurv - returned, number of survival time points
**    y[,1] - contains the survival times
**    risksum-the weighted N at that time
**    method - # of strata
**    strata[0: (n-1)]= last obs strata 1,2, etc
*/
#include <math.h>

void survfit2(sn, y, ny, wt, strata, method, error,mark,surv,
		  varh, risksum, snsurv)
long *sn;
long *snsurv;
long *method, *error;
long *ny;
long strata[];
double mark[];
double wt[], y[];
double varh[];
double surv[];
double risksum[];
{
    register int i,j;
    double hazard, varhaz;
    double sum, km;
    double *time, *status;
    double temp;
    int n;
    int nsurv, nstrat;

    n = *sn;
    time =y;
    status = y+n;
    strata[n-1] =1;   /*just in case the parent routine forgot */
    /*
    **  initialize a couple of arrays
    **    mark(i) contains the number of deaths at this particular time point
    **    risksum contains the running # at risk
    */
    for (i=0; i<(n-1); i++)  /*first pass, just mark the ties */
	if (time[i]==time[i+1] && strata[i]==0)  mark[i+1]= -1;
	else                                     mark[i+1]=  1;
    mark[0]=1;

    temp =0;
    for (i=n-1; i>=0; i--) {
	if (strata[i]==1) sum =0;
	sum  += wt[i];
	temp += status[i] * wt[i];

	if (mark[i] == 1) {
	    mark[i] = temp;
	    risksum[i] = sum;
	    temp =0;
	    }
	}

    /*
    ** the hazard starts out at zero;
    */
    nsurv=0;
    nstrat=0;
    km =1;
    hazard  =0;
    varhaz  =0;
    surv[0]=1.0;
    varh[0] = 0.0;
    for(i=0; i<n; i++) {
	if (mark[i] <0) continue;    /*extra data line */
	if (mark[i] >0) {
	    if (*method==1) {
	       km *= (risksum[i] - mark[i]) / risksum[i];
	       surv[nsurv] = km;;

		if (*error==1 )
		     varhaz += mark[i]/(risksum[i]*(risksum[i]-mark[i]));
		else varhaz += mark[i]/(risksum[i]*risksum[i]);
		varh[nsurv] = varhaz;
		}
	    else  if (*method==2) {
		hazard += mark[i]/risksum[i];
		surv[nsurv] = exp(-hazard);
		if (*error==1 )
		     varhaz += mark[i]/(risksum[i]*(risksum[i]-mark[i]));
		else varhaz += mark[i]/(risksum[i]*risksum[i]);
		varh[nsurv] = varhaz;
		}

	    else  if (*method==3) {
		for (j=0; j<mark[i]; j++) {
		    hazard += 1/(risksum[i] -j);
		    if (*error==1 )
			 varhaz += 1/(risksum[i]*(risksum[i]-1));
		    else varhaz += 1/(risksum[i]*risksum[i]);
		    }
		surv[nsurv] = exp(-hazard);
		varh[nsurv] = varhaz;
		}
	    }

	time[nsurv] = time[i];
	mark[nsurv] = mark[i];
	risksum[nsurv] = risksum[i];
	nsurv++;

	if (strata[i]==1) {
	    strata[nstrat]= nsurv;
	    nstrat++;
	    if (nsurv<n) {
		surv[nsurv] =1;
		varh[nsurv] = 0;
		}
	    km=1;
	    hazard  =0;
	    varhaz  =0;
	    }
	else if (nsurv < n) {
	    surv[nsurv] = surv[nsurv-1];
	    varh[nsurv] = varh[nsurv-1];
	    }
	}

    *method = nstrat;
    *snsurv = nsurv;
    }

