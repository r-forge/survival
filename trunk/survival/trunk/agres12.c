/* SCCS: $Id: agres12.c,v 1.3 1992-05-01 16:28:44 splus Exp $
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
**      hazard  carried forward from the cumhaz routine
**      cumhaz  carried forward from the cumhaz routine
**
** Output
**      resid2  matrix of score residuals, same "shape" of matrix as covar2
**
** Scratch
**      a       vector of length nvar
*/
#include <stdio.h>
extern double **dmatrix();

void agres1(nx, nvarx, y, covar2, strata, score,
		hazard, cumhaz, resid2)
long    nx[1],
	nvarx[1],
	strata[];
double  y[],
	*covar2,
	score[],
	hazard[],
	cumhaz[],
	*resid2;

    {
    register int i,j, k;
    register double temp;
    int n, nvar;
    int person;
    double denom, time;
    double *a;
    double *start, *stop, *event;
    double **covar,
	   **resid,
	   **wmean;

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
    a = (double *)S_alloc(nvar, sizeof(double));

    for (person=0; person<n; person++) {
	/* first piece for each person is covar * martingale resid */
	for (j=0; j<nvar; j++) resid[j][person] = covar[j][person] *
				   (event[person] - score[person]*cumhaz[person]);
	}
    for (person=0; person<n; ) {
	if (hazard[person]==0) person++;
	else {
	    /*
	    ** compute the mean over the risk set
	    */
	    denom =0;
	    for (i=0; i<nvar; i++) {
		a[i] =0;
		}
	    time = stop[person];
	    for (k=person; k<n; k++) {
		if (start[k] < time) {
		    denom += score[k];
		    for (i=0; i<nvar; i++) {
			a[i] = a[i] + score[k]*covar[i][k];
			}
		     }
		if (strata[k]==1) break;
		}

	    /* add mean times increment in hazard */
	    for (k=person; k<n; k++) {
		if (start[k] < time)
		    for (i=0; i<nvar; i++) {
			resid[i][k] += (a[i]/denom)*score[k]*hazard[person];
			}
		if (strata[k]==1) break;
		}

	    /* for events at exactly this time, subtract the mean */
	    for (k=person; k<n & stop[k]==time; k++) {
		/* walk through any tied death times */
		if (event[k]==1) {
		    for (i=0; i<nvar; i++) resid[i][k] -= a[i]/denom;
		    }
		person++;
		if (strata[k]==1) break;
		}
	    }
	}
    }
