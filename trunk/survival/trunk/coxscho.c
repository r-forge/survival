/* SCCS $Id: coxscho.c,v 4.1 1993-01-12 23:38:19 therneau Exp $
/*
** Return the Schoenfeld residuals.
**
**  the input parameters are
**
**       nused        :number of people
**       nvar         :number of covariates
**       y(3,n)       :start, stop, and event for each subject
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       score(n)     :the risk score for the subject
**       method       : =1 if the Efron method was used
**
**  returned parameters
**       covar        :for each death, the row is replaced with the residuals
**
**  work arrays
**       a(nvar)
**       a2(nvar)
**
**  the 2 arrays a and a2 are passed as a single
**    vector of storage, and then broken out.
**
**  the data must be sorted by ascending time within strata, deaths before
**          living within tied times.
*/
#include <math.h>
#include <stdio.h>

double **dmatrix();

void coxscho(nusedx, nvarx, y, covar2, score, strata,  method2, work)

long    *nusedx,
	*nvarx,
	*method2,
	strata[];
double  *covar2,
	*work,
	*score,
	*y;
{
    register int i,j,k,person;
    int     nused, nvar;
    double **covar;
    double *a;
    double *a2;
    double  denom, weight;
    double  time;
    double  temp, temp2;
    double     method;
    int     deaths;
    double efron_wt;
    double  *start,
	    *stop,
	    *event;

    nused = *nusedx;
    nvar  = *nvarx;
    method= *method2;
    /*
    **  Set up the ragged arrays
    */
    covar= dmatrix(covar2, nused, nvar);
    a = work + 2*nvar*nvar;
    a2= a+nvar;
    start =y;
    stop  =y + nused;
    event =y + nused +nused;

    /*
    ** Now walk through the data
    */
    for (person=0; person<nused;) {
	if (event[person]==0) person++;
	else {
	    /*
	    ** compute the mean over the risk set and over the deaths (a & a2)
	    */
	    denom =0;
	    efron_wt =0;
	    for (i=0; i<nvar; i++) {
		a[i] =0;
		a2[i]=0;
		}
	    time = stop[person];
	    deaths=0;
	    for (k=person; k<nused; k++) {
		if (start[k] < time) {
		    weight = score[k];
		    denom += weight;
		    for (i=0; i<nvar; i++) {
			a[i] += weight*covar[i][k];
			}
		    if (stop[k]==time && event[k]==1) {
			deaths += 1;
			efron_wt += weight*event[k];
			for (i=0; i<nvar; i++) a2[i]+= weight*covar[i][k];
			}
		     }
		if (strata[k]==1) break;
		}

	    /*
	    ** Compute the residual for this time point
	    */
	    temp = method * (deaths -1) / 2.0;
	    for (k=person; k<nused && stop[k]==time; k++) {
		if (event[k]==1) {
		    for (i=0; i<nvar; i++) {
			temp2 = (a[i] - temp*a2[i])/(denom - temp*efron_wt);
			covar[i][k] -= temp2;
			}
		    }
		person++;
		if (strata[k]==1) break;
		}
	    }
	}
    }
