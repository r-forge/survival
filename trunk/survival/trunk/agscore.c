/* SCCS: $Id: agscore.c,v 1.1 1993-01-30 19:51:22 therneau Exp $
/*
** Do the score residuals
**
** Input
**      nx      number of subjects
**      nvarx   number of variables in the covariance matrix
**      y       matrix of start, stop, and event
**      strata  =1 for the last obs of each strata
**      covar2  the matrix of covariates, rows=variables, columns=subjects
**                (the S executive stores matrices in the Fortran ordering)
**      score   the vector of subject scores, i.e., exp(beta*z)
**      method  ==1 for efron approx
**
** Output
**      resid2  matrix of score residuals, same "shape" of matrix as covar2
**
** Scratch
**      a       vector of length 3*nvar
*/
#include <stdio.h>
extern double **dmatrix();

void agscore(nx, nvarx, y, covar2, strata, score,
		method, resid2, a)
long    nx[1],
	nvarx[1],
	method[1],
	strata[];
double  y[],
	*covar2,
	score[],
	*a,
	*resid2;

    {
    register int i,k;
    register double temp;
    int n, nvar;
    int person;
    double denom, time;
    double *a2, *mean;
    double e_denom;
    double hazard;
    double  deaths, downwt;
    int dd;
    double *start, *stop, *event;
    double **covar,
	   **resid;
    double temp1, temp2;
    double *mh1, *mh2, *mh3;

    n = *nx;
    nvar  = *nvarx;
    start =y;
    stop  = y+n;
    event = y+(n+n);
    /*
    **  Set up the ragged arrays
    */
    covar=  dmatrix(covar2, n, nvar);
    resid = dmatrix(resid2, n, nvar);
    a2  = a+nvar;
    mean= a2 + nvar;
    mh1 = mean + nvar;
    mh2 = mh1 + nvar;
    mh3 = mh2 + nvar;

    for (person=0; person<n; ) {
	if (event[person]==0) person++;
	else {
	    /*
	    ** compute the mean over the risk set, also hazard at this time
	    */
	    denom =0;
	    e_denom =0;
	    deaths =0;
	    for (i=0; i<nvar; i++) {
		a[i] =0;
		a2[i]=0;
		}
	    time = stop[person];
	    for (k=person; k<n; k++) {
		if (start[k] < time) {
		    denom += score[k];
		    for (i=0; i<nvar; i++) {
			a[i] = a[i] + score[k]*covar[i][k];
			}
		     if (stop[k]==time && event[k]==1) {
			deaths++;
			e_denom += score[k];
			for (i=0; i<nvar; i++)
			    a2[i] = a2[i] + score[k]*covar[i][k];
			}
		     }
		if (strata[k]==1) break;
		}

	    /* add things in for everyone in the risk set*/
	    if (deaths <2 || *method==0) {
		/* easier case */
		hazard = deaths/denom;
		for (i=0; i<nvar; i++) mean[i] = a[i]/denom;
		for (k=person; k<n; k++) {
		    if (start[k] < time) {
			for (i=0; i<nvar; i++)
			    resid[i][k] -= (covar[i][k] -mean[i])*score[k]*hazard;
			if (stop[k]==time) {
			    person++;
			    if (event[k]==1)
				for (i=0; i<nvar; i++)
				    resid[i][k] += (covar[i][k] -mean[i]);
			    }
			}
		    if (strata[k]==1) break;
		    }
		}

	    else {
		temp1=0;
		temp2=0;
		for (i=0; i<nvar; i++) {
		    mh1[i] =0;
		    mh2[i] =0;
		    mh3[i] =0;
		    }
		for (dd=0; dd<deaths; dd++){
		    downwt = dd/deaths;
		    hazard = 1/(denom - downwt*e_denom);
		    temp1 += hazard;
		    temp2 += (1-downwt) * hazard;
		    for (i=0; i<nvar; i++) {
			mean[i] = (a[i] - downwt*a2[i])*hazard;
			mh1[i]  += mean[i] * hazard;
			mh2[i]  += mean[i] * (1-downwt) * hazard;
			mh3[i]  += mean[i]/deaths;
			}
		    }
		for (k=person; k<n; k++) {
		    if (start[k] < time) {
			if (stop[k]==time && event[k]==1) {
			    for (i=0; i<nvar; i++) {
				resid[i][k] += covar[i][k] - mh3[i];
				resid[i][k] -= score[k]*covar[i][k]*temp2;
				resid[i][k] += score[k]* mh2[i];
				}
			    }
			else {
			    for (i=0; i<nvar; i++)
				resid[i][k] -= score[k]*(covar[i][k]*temp1 - mh1[i]);
			    }
			}
		    if (strata[k]==1) break;
		    }
		for (person; stop[person]==time; person++)
		    if (strata[person]==1) break;
		}
	    }
	}
    }
