/* SCCS  $Id: coxscore.c,v 4.1 1993-01-30 19:51:50 therneau Exp $
/*
** Compute the score residuals for a Cox model
**
** Input
**      nx      number of subjects
**      nvarx   number of variables in the covariance matrix
**      y       matrix of time and status values
**      strata  =1 for the last obs of each strata
**      covar2  the matrix of covariates, rows=variables, columns=subjects
**                (the S executive stores matrices in the Fortran ordering)
**      score   the vector of subject scores, i.e., exp(beta*z)
**      method  ==1 for efron method
**
** Output
**      resid   a matrix of the same shape as x
**
** Scratch
**      scratch,  from which a and a2 are carved
*/
#include <stdio.h>
extern double **dmatrix();

void coxscore(nx, nvarx, y, covar2, strata, score, method, resid2, scratch)
long    nx[1],
	nvarx[1],
	*method,
	strata[];
double  y[],
	*covar2,
	*scratch,
	*resid2,
	score[];

    {
    register int i,j, k;
    register double temp;
    int n, nvar;
    double deaths;
    int dd;
    double *time, *status;
    double *a, *a2;
    double denom, e_denom;
    double **covar;
    double **resid;
    double hazard;
    double downwt, temp2;

    n = *nx;
    nvar  = *nvarx;
    time = y;
    status = y+n;
    a = scratch;
    a2 = a+nvar;
    /*
    **  Set up the ragged array
    */
    covar=  dmatrix(covar2, n, nvar);
    resid=  dmatrix(resid2, n, nvar);

    e_denom=0;
    deaths=0;
    for (i=0; i<nvar; i++) a2[i] =0;
    strata[n-1] =1;  /*failsafe */
    for (i=n-1; i >=0; i--) {
	if (strata[i]==1) {
	    denom =0;
	    for (j=0; j<nvar; j++) a[j] =0;
	    }

	denom += score[i];
	if (status[i]==1) {
	    deaths++;
	    e_denom += score[i];
	    for (j=0; j<nvar; j++) a2[j] += score[i]*covar[j][i];
	    }
	for (j=0; j<nvar; j++) {
	    a[j] += score[i] * covar[j][i];
	    resid[j][i] =0;
	    }

	if (deaths>0 && (i==0 || strata[i-1]==1 || time[i]!=time[i-1])){
	    /* last obs of a set of tied death times */
	    if (deaths <2 || *method==0) {
		hazard = deaths/denom;
		for (j=0; j<nvar; j++)  {
		    temp = (a[j]/denom);     /* xbar */
		    for (k=i; k<n; k++) {
			temp2 = covar[j][k] - temp;
			if (time[k]==time[i] && status[k]==1)
				resid[j][k] += temp2;
			resid[j][k] -= temp2* score[k] * hazard;
			if (strata[k]==1) break;
			}
		    }
		}
	    else {  /* the harder case */
		for (dd=0; dd<deaths; dd++) {
		    downwt = dd/deaths;
		    temp = denom - downwt* e_denom;
		    hazard = 1/temp;
		    for (j=0; j<nvar; j++) {
			temp = (a[j] - downwt*a2[j])/ temp;
			for (k=i; k<n; k++) {
			    temp2 = covar[j][k] - temp;
			    if (time[k]==time[i] && status[k]==1) {
				resid[j][k] += temp2/deaths;
				resid[j][k] -= temp2 * score[k] * hazard *
						    (1 - downwt);
				}
			    else resid[j][k]-= temp2*score[k] * hazard;
			    if (strata[k]==1) break;
			    }
			}
		    }
		}
	    e_denom =0;
	    deaths =0;
	    for (j=0; j<nvar; j++)  a2[j] =0;
	    }
	}
    }
