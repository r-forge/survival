/* SCCS $Id: survfit2.c,v 4.2 1992-03-24 09:23:30 therneau Exp $  */
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
**    strata[0]= # of strata, strata[1:n]= last obs strata 1,2, etc
*/
#include <math.h>

void survfit(sn, y, ny, wt, strata, method, error,mark,surv,
		  varh, risksum, snsurv)
long *sn;
long mark[], *snsurv;
long *method, *error;
long *ny;
long strata[];
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
    /*
    **  initialize a couple of arrays
    **    mark(i) contains the number of deaths at this particular time point
    **    risksum contains the running # at risk
    */
    strata[n-1] =1;   /*just in case the parent routine forgot */
    for (i=n-1; i>0; i--) {
	if (strata[i]==1) {
	    sum=0;
	    j=0;
	    temp=0;
	    }
	sum += wt[i];
	risksum[i] = sum;

	if ( time[i]==time[i-1] && strata[i-1] ==0) {
	    j +=status[i] * wt[i];
	    mark[i] =0;
	    }
	else {
	    mark[i] = j +status[i]* wt[i];
	    temp=0;
	    j = 0;
	    }
	}

    mark[0] = j+ status[0]*wt[0];
    risksum[0] = sum + wt[0];

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
    for(i=0; i<n; ) {
	if (mark[i] >0) {
	    if (*error==1 )
		 varhaz += mark[i]/(risksum[i]*(risksum[i]-mark[i]));
	    else varhaz += mark[i]/(risksum[i]*risksum[i]);
	    varh[nsurv] = varhaz;

	    if (*method==2) {
		/* fh  estimator is easy */
		hazard += mark[i]/risksum[i];
		surv[nsurv] = exp(-hazard);
		}

	    else  {
	       km *= (risksum[i] - mark[i]) / risksum[i];
	       surv[nsurv] = km;;
	       }
	    }
	time[nsurv] = time[i];
	mark[nsurv] = mark[i];
	risksum[nsurv] = risksum[i];
	nsurv++;

	if (strata[i]==1) {
	    nstrat++;
	    strata[nstrat]= nsurv;
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

	/* walk past any tied survival times */
	temp = time[i];
	for (; i<n && time[i]==temp; i++);
	}

    strata[0] = nstrat;
    *snsurv = nsurv;
    }

