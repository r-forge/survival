/* SCCS $Id: agfit2.c,v 4.2 1992-03-30 02:21:55 therneau Exp $  */
/*
** Anderson-Gill formulation of the cox Model
**
**  the input parameters are
**
**       maxiter      :number of iterations
**       nused        :number of people
**       nvar         :number of covariates
**       start(n)     :each row covers the time interval (start,stop]
**       stop(n)      :
**       event(n)     :was there an event at 'stop':1=dead , 0=censored
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       offset(n)    :linear offset
**       eps          :tolerance for convergence.  Iteration continues until
**                       the percent change in loglikelihood is <= eps.
**       maxbe        :the 'infinite beta' check values
**
**  returned parameters
**       means(nv)    :column means of the X matrix
**       beta(nv)     :the vector of answers (at start contains initial est)
**       imat(nv,nv)  :the variance matrix at beta=final, also a ragged array
**                      if flag<0, imat is undefined upon return
**       loglik(2)    :loglik at beta=initial values, at beta=final
**       sctest       :the score test at beta=initial
**       flag         :success flag  0 - ok
**                                  -1 to -nvar: variable __ caused a singular
**                                             matrix failure
**                                   1 to nvar: var __ is going to infinity
**                                 1000 did not converge
**       maxiter      :actual number of iterations used
**
**  work arrays
**       score(n)              the score exp(beta*z)
**       a(nvar)
**       a2(nvar)
**       cmat(nvar,nvar)       ragged array
**       cmat2(nvar,nvar)
**       newbeta(nvar)         always contains the "next iteration"
**       u(nvar)
**
**  the 7 arrays score, a, cmat, newbeta, and u are passed as a single
**    vector of storage, and then broken out.
**
**  calls functions:  cholesky, chsolve, chinv
**
**  the data must be sorted by ascending time within strata, deaths before
**          living within tied times.
*/
#include <math.h>
#include <stdio.h>

double **dmatrix();

void agfit(maxiter, nusedx, nvarx, start, stop, event, covar2, offset, strata,
		 means, beta, maxbe, imat2, loglik, flag, work,
		 eps, sctest)

long    *nusedx,
	*nvarx,
	*maxiter,
	*flag,
	strata[],
	event[];
double  *covar2,
	*imat2;
double  maxbe[],
	means[],
	*work,
	beta[],
	offset[],
	start[],
	stop[];
double  loglik[2],  /* returned values */
	*sctest;
double  *eps;
{
    register int i,j,k,person;
    int     ierr,iter;
    int     nused, nvar;

    double **covar, **cmat, **imat;  /*ragged array versions*/
    double *a, *newbeta, *u;
    double *a2, **cmat2;
    double *score;
    double  denom, zbeta, weight;
    double  time;
    double  temp, temp2;
    double  newlk;
    int     halving;    /*are we doing step halving at the moment? */
    double     method;
    int itemp, deaths;
    double efron_wt, d2;

    nused = *nusedx;
    nvar  = *nvarx;
    method= *sctest;
    /*
    **  Set up the ragged arrays
    */
    covar= dmatrix(covar2, nused, nvar);
    imat = dmatrix(imat2,  nvar, nvar);
    cmat = dmatrix(work,   nvar, nvar);
    cmat2= dmatrix(work + nvar*nvar, nvar, nvar);
    a = work + 2*nvar*nvar;
    a2= a+nvar;
    u = a2+ nvar;
    newbeta = u + nvar;
    score   = newbeta + nvar;

    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable
    ** At the same time, determine maxbeta (input contains the constant).
    */
    for (i=0; i<nvar; i++) {
	temp=0;
	for (person=0; person<nused; person++) temp += covar[i][person];
	temp /= nused;
	means[i] = temp;
	denom=0;
	for (person=0; person<nused; person++) {
	    covar[i][person] -=temp;
	    denom += fabs(covar[i][person]);
	    }
	maxbe[i] /= (denom/nused);
	}

    /*
    ** do the initial iteration step
    */
    loglik[1] =0;
    for (i=0; i<nvar; i++) {
	u[i] =0;
	for (j=0; j<nvar; j++)
	    imat[i][j] =0 ;
	}

    for (person=0; person<nused; person++) {
	zbeta = 0;      /* form the term beta*z   (vector mult) */
	for (i=0; i<nvar; i++)
	    zbeta += beta[i]*covar[i][person];
	score[person] = exp(zbeta + offset[person]);
        }

    for (person=0; person<nused;) {
	if (event[person]==0) person++;
	else {
	    /*
	    ** compute the mean and covariance over the risk set (a and c)
	    */
	    denom =0;
	    efron_wt =0;
	    for (i=0; i<nvar; i++) {
		a[i] =0;
		a2[i]=0;
		for (j=0; j<nvar; j++) {
		    cmat[i][j]=0;
		    cmat2[i][j]=0;
		    }
		}
	    time = stop[person];
	    deaths=0;
	    for (k=person; k<nused; k++) {
		if (start[k] < time) {
		    weight = score[k];
		    denom += weight;
		    for (i=0; i<nvar; i++) {
			a[i] += weight*covar[i][k];
			for (j=0; j<=i; j++)
			    cmat[i][j] += weight*covar[i][k]*covar[j][k];
			}
		    if (stop[k]==time && event[k]==1) {
			deaths += event[k];
			efron_wt += weight*event[k];
			for (i=0; i<nvar; i++) {
			    a2[i]+= weight*covar[i][k];
			    for (j=0; j<=i; j++)
				cmat2[i][j] += weight*covar[i][k]*covar[j][k];
			    }
			}
		     }
		if (strata[k]==1) break;
		}

	    /*
	    ** Add results into u and imat for all events at this time point
	    */
	    itemp = -1;
	    for (k=person; k<nused && stop[k]==time; k++) {
		if (event[k]==1) {
		    itemp++;
		    temp = itemp*method/deaths;
		    d2 = denom - temp*efron_wt;
		    loglik[1] +=  log(score[k]/d2);
		    for (i=0; i<nvar; i++) {
			temp2 = (a[i] - temp*a2[i])/d2;
			u[i] += covar[i][k] - temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] += (cmat[i][j] - temp*cmat2[i][j])/d2-
					   temp2*(a[j]-temp*a2[j])/d2;
			}
		    }
		person++;
		if (strata[k]==1) break;
		}
	    }
	}   /* end  of accumulation loop */

    loglik[0] = loglik[1];   /* save the loglik for iteration zero  */

    /* am I done?
    **   update the betas and test for convergence
    */

    for (i=0; i<nvar; i++) /*use 'a' as a temp to save u0, for the score test*/
	a[i] = u[i];

    ierr= cholesky(imat, nvar);
    if (ierr != 0) {
	*flag= -ierr;
	*maxiter=0;
	return;
	}

    chsolve(imat,nvar,a);        /* a replaced by  a *inverse(i) */

    *sctest=0;
    for (i=0; i<nvar; i++)
	*sctest +=  u[i]*a[i];

    /*
    **  Never, never complain about convergence on the first step.  That way,
    **  if someone HAS to they can force one iter at a time.
    */
    for (i=0; i<nvar; i++) {
	newbeta[i] = beta[i] + a[i];
	}
    if (*maxiter==0) {
	chinv(imat,nvar);
	for (i=1; i<nvar; i++)
	    for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	*flag=0;
	return;   /* and we leave the old beta in peace */
	}

    /*
    ** here is the main loop
    */
    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (iter=1; iter<=*maxiter; iter++) {
	newlk =0;
	for (i=0; i<nvar; i++) {
	    u[i] =0;
	    for (j=0; j<nvar; j++)
		imat[i][j] =0;
	    }

	for (person=0; person<nused; person++) {
	    zbeta = 0;      /* form the term beta*z   (vector mult) */
	    for (i=0; i<nvar; i++)
		zbeta += newbeta[i]*covar[i][person];
	    score[person] = exp(zbeta + offset[person]);
	    }

	for (person=0; person<nused; ) {
	    if (event[person]==0) person++;
	    else {
		/*
		** compute the mean and covariance over the risk set (a and c)
		*/
		efron_wt =0;
		denom =0;
		for (i=0; i<nvar; i++) {
		    a[i] =0;
		    a2[i]=0;
		    for (j=0; j<nvar; j++) {
			cmat[i][j]=0;
			cmat2[i][j]=0;
			}
		    }
		time = stop[person];
		deaths=0;
		for (k=person; k<nused; k++) {
		    if (start[k] < time) {
			weight = score[k];
			denom += weight;
			for (i=0; i<nvar; i++) {
			    a[i] += weight*covar[i][k];
			    for (j=0; j<=i; j++)
				cmat[i][j] += weight*covar[i][k]*covar[j][k];
			    }
			if (stop[k]==time && event[k]==1) {
			    deaths += event[k];
			    efron_wt += weight*event[k];
			    for (i=0; i<nvar; i++) {
				a2[i]+= weight*covar[i][k];
				for (j=0; j<=i; j++)
				    cmat2[i][j] += weight*covar[i][k]*covar[j][k];
				}
			    }
			}
		    if (strata[k]==1) break;
		    }

		itemp = -1;
		for (k=person; k<nused && stop[k]==time; k++) {
		    if (event[k]==1) {
			itemp++;
			temp = itemp*method/deaths;
			d2 = denom - temp*efron_wt;
			newlk +=  log(score[k]/d2);
			for (i=0; i<nvar; i++) {
			    temp2 = (a[i] - temp*a2[i])/d2;
			    u[i] += covar[i][k] - temp2;
			    for (j=0; j<=i; j++)
				imat[j][i] += (cmat[i][j] - temp*cmat2[i][j])/d2-
					       temp2*(a[j]-temp*a2[j])/d2;
			    }
			}

		    person++;
		    if (strata[k]==1) break;
		    }
		}
	    }   /* end  of accumulation loop */

	/* am I done?
	**   update the betas and test for convergence
	*/
	ierr = cholesky(imat, nvar);
	if (ierr != 0) {
	    *flag= -ierr;
	    *maxiter=iter;
	    return;
	    }

	if (fabs(1-(loglik[1]/newlk))<=*eps ) { /* all done */
	    loglik[1] = newlk;
	    chinv(imat, nvar);     /* invert the information matrix */
	    for (i=1; i<nvar; i++)
		for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
	    for (i=0; i<nvar; i++)
		beta[i] = newbeta[i];
	    if (halving==1) *flag= 1000; /*didn't converge after all */
	    else            *flag=0;
	    *maxiter = iter;
	    return;
	    }

	if (iter==*maxiter) break;  /*skip the step halving and etc */

	if (newlk < loglik[1])   {    /*it is not converging ! */
		halving =1;
		for (i=0; i<nvar; i++)
		    newbeta[i] = (newbeta[i] + beta[i]) /2; /*half of old increment */
		}
	    else {
		halving=0;
		loglik[1] = newlk;
		chsolve(imat,nvar,u);

		j=0;
		for (i=0; i<nvar; i++) {
		    beta[i] = newbeta[i];
		    newbeta[i] = newbeta[i] +  u[i];

		    if (fabs(newbeta[i]) > maxbe[i]) {   /*failing ! */
			*flag = 1+i;
			*maxiter = iter;
			return;
			}
		    }
		}
	}   /* return for another iteration */

    loglik[1] = newlk;
    chinv(imat, nvar);
    for (i=1; i<nvar; i++)
	for (j=0; j<i; j++)  imat[i][j] = imat[j][i];
    for (i=0; i<nvar; i++)
	beta[i] = newbeta[i];
    *flag= 1000;
    }
