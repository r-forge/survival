/* SCCS $Id: coxfit2.c,v 4.3 1992-08-10 13:36:48 grill Exp $  */
/*
** here is a cox regression program, written in c
**     uses Efron's approximation for ties
**  the input parameters are
**
**       maxiter      :number of iterations
**       nused        :number of people
**       nvar         :number of covariates
**       time(n)      :time of event or censoring for person i
**       status(n)    :status for the ith person    1=dead , 0=censored
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       offset(n)    :offset for the linear predictor
**       eps          :tolerance for convergence.  Iteration continues until
**                       the percent change in loglikelihood is <= eps.
**       sctest       : on input contains the method 0=Breslow, 1=Efron
**
**  returned parameters
**       means(nv)    : vector of column means of X
**       beta(nv)     :the vector of answers (at start contains initial est)
**       u(nv)        :score vector
**       imat(nv,nv)  :the variance matrix at beta=final, also a ragged array
**                      if flag<0, imat is undefined upon return
**       loglik(2)    :loglik at beta=initial values, at beta=final
**       sctest       :the score test at beta=initial
**       flag         :success flag  0 - ok
**                                  1 to nvar: variable __ caused a singular
**                                             matrix failure
**                                 1000 did not converge
**       maxiter      :actual number of iterations used
**
**  work arrays
**       mark(n)
**       a(nvar)
**       cmat(nvar,nvar)       ragged array
**       newbeta(nvar)         always contains the "next iteration"
**
**  the 5 arrays a, a2, cmat, cmat2, and newbeta are passed as a single
**    vector of storage, and then broken out.
**
**  calls functions:  cholesky, chsolve, chinv
**
**  the data must be sorted by ascending time within strata
*/
#include <math.h>
#include <stdio.h>

double **dmatrix();

void coxfit2(maxiter, nusedx, nvarx, time, status, covar2, offset, strata,
		 means, beta, u, imat2, loglik, flag, mark, work,
		 eps, sctest)

long    *nusedx,
	*nvarx,
	*maxiter,
	*flag,
	strata[],
	status[],
	mark[];
double  *covar2,
	*imat2;
double  u[],
	means[],
	*work,
	beta[],
	offset[],
	time[];
double  loglik[2],  /* returned values */
	*sctest;
double  *eps;
{
    register int i,j,k, person;
    int     ierr,iter;
    int     nused, nvar;

    double **covar, **cmat, **imat;  /*ragged array versions*/
    double *a, *newbeta;
    double *a2, **cmat2;
    double  denom, zbeta, weight;
    double  temp, temp2;
    double  newlk;
    double  d2, efron_wt;
    int     halving;    /*are we doing step halving at the moment? */
    double     method;

    nused = *nusedx;
    nvar  = *nvarx;
    method= *sctest;
    /*
    **  Set up the ragged arrays
    */
    covar= dmatrix(covar2, nused, nvar);
    imat = dmatrix(imat2,  nvar, nvar);
    cmat = dmatrix(work,   nvar, nvar);
    cmat2= dmatrix(work+nvar*nvar, nvar, nvar);
    a = work + 2*nvar*nvar;
    newbeta = a + nvar;
    a2 = newbeta + nvar;

    /*
    **   Mark(i) contains the number of tied deaths at this point,
    **    for the first person of several tied times. It is zero for
    **    the second and etc of a group of tied times.
    */
    j=0;
    for (i=nused-1; i>0; i--) {
	if ((time[i]==time[i-1]) & (strata[i-1] != 1)) {
	    j=j+status[i];
	    mark[i]=0;
	    }
	else  {
	    mark[i]=j+status[i];
	    j=0;
	    }
	}
    mark[0] = j+status[0];

    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable
    */
    for (i=0; i<nvar; i++) {
	temp=0;
	for (person=0; person<nused; person++) temp += covar[i][person];
	temp /= nused;
	means[i] = temp;
	}

    /*
    ** do the initial iteration step
    */
    strata[nused-1] =1;
    loglik[1] =0;
    for (i=0; i<nvar; i++) {
	u[i] =0;
	for (j=0; j<nvar; j++)
	    imat[i][j] =0 ;
	}

    for (person=nused-1; person>=0; person--) {
	if (strata[person] == 1) {
	    efron_wt =0;
	    denom = 0;
	    for (i=0; i<nvar; i++) {
		a[i] = 0;
		a2[i]=0 ;
		for (j=0; j<nvar; j++) {
		    cmat[i][j] = 0;
		    cmat2[i][j]= 0;
		    }
		}
	    }

	zbeta = 0;      /* form the term beta*z   (vector mult) */
	for (i=0; i<nvar; i++)
	    zbeta += beta[i]*covar[i][person];
	weight = exp(zbeta +offset[person]);

	denom += weight;
	efron_wt += status[person] * weight;  /* sum of wts for tied deaths*/

	for (i=0; i<nvar; i++) {
	    a[i] += weight*covar[i][person];
	    for (j=0; j<=i; j++)
		cmat[i][j] += weight*covar[i][person]*covar[j][person];
	    }

	if (status[person]==1) {
	    loglik[1] = loglik[1] + zbeta;
	    for (i=0; i<nvar; i++) {
		    u[i] += covar[i][person];
		a2[i] +=  weight*covar[i][person];
		for (j=0; j<=i; j++)
		    cmat2[i][j] += weight*covar[i][person]*covar[j][person];
		}
	    }
	if (mark[person] >0) {  /* once per unique death time */
	    for (k=0; k<mark[person]; k++) {
		/*
		** trick of the week: if method=0, then temp always =0,
		**  giving the Breslow method.  This was simpler than lots
		**  of if-then-else in the code.
		*/
		temp = (double)k *method / mark[person];
		d2= denom - temp*efron_wt;
		loglik[1] -= log(d2);
		for (i=0; i<nvar; i++) {
		    temp2 = (a[i] - temp*a2[i])/ d2;
		    u[i] -= temp2;
		    for (j=0; j<=i; j++)
			imat[j][i] +=  (cmat[i][j] - temp*cmat2[i][j]) /d2 -
					  temp2*(a[j]-temp*a2[j])/d2;
		    }
		}
	    efron_wt =0;
	    for (i=0; i<nvar; i++) {
		a2[i]=0;
		for (j=0; j<nvar; j++)  cmat2[i][j]=0;
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
	*flag= ierr;
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

	for (person=nused-1; person>=0; person--) {
	    if (strata[person] == 1) {
		efron_wt =0;
		denom = 0;
		for (i=0; i<nvar; i++) {
		    a[i] = 0;
		    a2[i]=0 ;
		    for (j=0; j<nvar; j++) {
			cmat[i][j] = 0;
			cmat2[i][j]= 0;
			}
		    }
		}

	    zbeta = 0;
	    for (i=0; i<nvar; i++)
		zbeta += newbeta[i]*covar[i][person];
	    weight = exp(zbeta +offset[person]);

	    denom += weight;
	    efron_wt += status[person] * weight;  /* sum of wts for tied deaths*/

	    for (i=0; i<nvar; i++) {
		a[i] += weight*covar[i][person];
		for (j=0; j<=i; j++)
		    cmat[i][j] += weight*covar[i][person]*covar[j][person];
		}

	    if (status[person]==1) {
		newlk += zbeta;
		for (i=0; i<nvar; i++) {
		    u[i] += covar[i][person];
		    a2[i] +=  weight*covar[i][person];
		    for (j=0; j<=i; j++)
			cmat2[i][j] += weight*covar[i][person]*covar[j][person];
		    }
		}

	    if (mark[person] >0) {  /* once per unique death time */
		for (k=0; k<mark[person]; k++) {
		    temp = (double)k* method /mark[person];
		    d2= denom - temp*efron_wt;
		    newlk -= log(d2);
		    for (i=0; i<nvar; i++) {
			temp2 = (a[i] - temp*a2[i])/ d2;
			u[i] -= temp2;
			for (j=0; j<=i; j++)
			    imat[j][i] +=  (cmat[i][j] - temp*cmat2[i][j]) /d2 -
					      temp2*(a[j]-temp*a2[j])/d2;
			}
		    }
		efron_wt =0;
		for (i=0; i<nvar; i++) {
		    a2[i]=0;
		    for (j=0; j<nvar; j++)  cmat2[i][j]=0;
		    }
		}
	    }   /* end  of accumulation loop  */

	/* am I done?
	**   update the betas and test for convergence
	*/
	ierr = cholesky(imat, nvar);
	if (ierr != 0) {
	    *flag= ierr;
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
    return;
    }
