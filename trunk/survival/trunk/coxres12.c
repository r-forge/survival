/* SCCS: @(#)coxres1.c	2.1  6/12/91 */
/*
** Do the score residuals
**
** Input
**      nx      number of subjects
**      nvarx   number of variables in the covariance matrix
**      y       matrix of time and status values
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
**      wmean2  matrix of the same size as resid2
*/
#include <stdio.h>
extern double **matrix();

void coxres1(nx, nvarx, y, covar2, strata, score,
		hazard, cumhaz, resid2, wmean2)
long    nx[1],
	nvarx[1],
	strata[];
double  y[],
	*covar2,
	score[],
	hazard[],
	cumhaz[],
	*wmean2,
	*resid2;

    {
    register int i,j, k;
    register double temp;
    int n, nvar;
    double *time, *status;

    double **covar,
	   **resid,
	   **wmean;

    n = *nx;
    nvar  = *nvarx;
    time = y;
    status = y+n;
    /*
    **  Set up the ragged arrays
    */
    covar=  matrix(covar2, n, nvar);
    wmean = matrix(wmean2, n, nvar);
    resid = matrix(resid2, n, nvar);

    /*
    ** set up the matrix of weighted means
    **  use the first row of resid as a temp variable
    */
    for (i=n-1; i>=0; i--) {
	if (strata[i]==1) {
	    temp =0;
	    for (j=0; j<nvar; j++)  resid2[j] =0;
	    }
	temp += score[i];
	for (j=0; j<nvar; j++) {
	    resid2[j] += covar[j][i]* score[i];
	    wmean[j][i] = resid2[j]/temp;
	    }
	}

    /* correct for tied times */
    for (i=0; i<n; ) {
	for (k=i+1; k<n && time[k]==time[i]; k++) {
	    for (j=0; j<nvar; j++)
		wmean[j][k] = wmean[j][i];
	    if (strata[k]==1) break;
	    }
	i=k;
	}

    /*
    ** Now do the actual residuals
    */
    for (i=0; i<n; i++)
	for (j=0; j<nvar; j++)
	    resid[j][i]= status[i]*(covar[j][i]-wmean[j][i]) -
			      covar[j][i] * cumhaz[i] * score[i];

    for (i=0; i<n; i++) {
	temp = hazard[i];     /*size of the hazard function's jump */
	if (temp>0) {
	    for (k=i; k<n; k++) {
		for (j=0; j<nvar; j++)
		    resid[j][k] += wmean[j][i] * temp * score[k];
		if (strata[k]==1) break;
		}
	    }
	}
   }
