/* SCCS  $Id: agfit4.c,v 1.1 1998-10-28 08:42:18 therneau Exp $  */
/* A reentrant version of the agfit program, for random effects modeling
**   with reasonable efficiency (I hope).  The important arrays are saved
**   from call to call so as to speed up the process.  The x-matrix itself
**   is the most important of these.
**
** agfit4_a: Entry and intial iteration step for beta=initial, theta=0
**              (no frailty)
**            Most of the same arguments as agfit2.
**            Allocate and save arrays in static locations.
** agfit4_b: Iterate to convergence given an initial value.
** agfit4_c: Compute residuals and release the saved memory.
**
**     McGilchrist's method for frailty with a fixed theta, but for
**     space savings I assume that many elements of imat are zero
**
**  the input parameters are
**
**       maxiter      :number of iterations
**       nused        :number of people
**       nvar         :number of covariates
**       yy[3,n]      :row 1: start time of event or censoring for person i
**		      :row 2: stop time
**                    :row 3: status for the ith person    1=dead , 0=censored
**       covar(nv,n)  :covariates for person i.
**                        Note that S sends this in column major order.
**       strata(n)    :marks the strata.  Will be 1 if this person is the
**                       last one in a strata.  If there are no strata, the
**                       vector can be identically zero, since the nth person's
**                       value is always assumed to be = to 1.
**       offset(n)    :offset for the linear predictor
**       weights(n)   :case weights
**       eps          :tolerance for convergence.  Iteration continues until
**                       the percent change in loglikelihood is <= eps.
**       tolerch      :tolerance for the Cholesky routines
**       method       : Method 0=Breslow, 1=Efron
**       ptype        : 1 or 3 -- there is a sparse term
**                    : 2 or 3 -- there is a non-sparse term in the model
**       nfrail       : number of frailty groups (sparse terms), 0 if there are
**                        none
**       frail        : a vector containing the frailty groups
**       fbeta        : initial frailty estimates
**       pdiag        : if 0, then for the non-sparse terms only the diagonal
**                        of the variance matrix is penalized, otherwise the
**                        full matrix is used.
**
**  returned parameters
**       means(nv)    : vector of column means of X
**       beta(nv)     : the vector of answers (at start contains initial est)
**       u(nv)        : score vector
**       imat(nv,nv)  : the variance matrix at beta=final
**                      if flag<0, imat is undefined upon return
**       loglik       :loglik at beta=final
**       flag         :success flag  1000  did not converge
**                                   1 to nvar: rank of the solution
**       maxiter      :actual number of iterations used
**       fbeta(nfrail): fitted frailty values
**       fdiag(nfrail + nvar): diagonal of cholesky of the full inverse
**       jmat         : inverse of the cholesky
**       imat         : cholesky of the information matrix
**       expect       : contains the "expected" for each subject
** 
**  work arrays
**       end(n)      :how far to look
**       score(n)     
**       a(nvar+ nfrail), a2(nvar+nfrail)
**       cmat(nvar,nvar+nfrail)       ragged array
**       cmat2(nvar,nvar+nfrail)
**       fdiag                         the diagonal of the sparse information
**       oldbeta(nvar + nfrail)         always contains the "last iteration"
**
**  the work arrays are passed as a single
**    vector of storage, and then broken out.
**
**  calls functions:  cholesky3, chsolve3, chinv2
**
**  the data must be sorted by ascending time within strata
*/
#include <math.h>
#include <stdio.h>
#include "survS.h"
#include "survproto.h"

static double **covar, **cmat, **cmat2;
static int *end;
static double *a, *oldbeta, *a2;
static double *offset, *weights;
static int    *event, *frail, *strata;
static double *score, *start, *stop;
static double *tmean;
static int    ptype, pdiag;
static double *ipen, *upen, logpen;
static long   *zflag;

static double **cmatrix(double *, int, int);

void agfit4_a(long *nusedx, long *nvarx, double *yy, 
	       double *covar2, double *offset2,
	       double *weights2, long *strata2,
	       double *means, double *beta, double *u, 
	       double *loglik, 
	       long *methodx, long *ptype2, long *pdiag2,
	       long *nfrail,  long *frail2)
{
S_EVALUATOR
    int i,j,k, person;
    int     nused, nvar;
    int    nf, nvar2;
    int  deaths, itemp, endp;

    double  denom, zbeta, risk;
    double  temp;
    double  d2, efron_wt;
    double  method;
    double  meanwt, time;

    nused = *nusedx;
    nvar  = *nvarx;
    nf= *nfrail;
    method= *methodx;
    nvar2 = nvar + nf;
    ptype = *ptype2;
    pdiag = *pdiag2;

    /*
    **  Allocate storage for the arrays and vectors
    **  Since they will be used later, sizes are based on what will be
    **    needed with the frailty terms.
    */
    if (nvar >0) {
	covar= cmatrix(covar2, nused, nvar);
	cmat = cmatrix(0, nvar2, nvar+1);
	cmat2= cmatrix(0, nvar2, nvar+1);
        }

    a = Calloc(4*nvar2 + 5*nused, double);
    oldbeta = a + nvar2;
    a2 =  oldbeta + nvar2;
    weights = a2+ nvar2;
    offset  = weights + nused;
    score   = offset + nused;
    tmean   = score + nused;
    start   = tmean + nvar2;
    stop    = start + nused;
    event  = Calloc(3*nused, int);
    strata  = event + nused;
    end     = strata+ nused;
    for (i=0; i<nused; i++) {
	weights[i] = weights2[i];
	offset[i]  = offset2[i];
	event[i]  =  yy[nused + nused +i];
	strata[i] = strata2[i];
	start[i]  = yy[i];
	stop[i]   = yy[nused+i];
        }
    strata[nused-1] =1;  /* failsafe */   

    /* scratch space for penalty 
    **    upen needs to be max(nvar, nfrail), 
    **    ipen max(nfrail, nvar(if pdiag=0) or nvar^2 )
    */
    if (nf > nvar) i=nf; else i=nvar;
    if (nf > nvar*nvar) j=nf; else j=nvar*nvar;
    if (pdiag==0)  upen = Calloc(2*i, double);
    else           upen = Calloc(i+j, double);
    ipen = upen + i;
    if (ptype>1)  zflag = Calloc(nvar, long);
    else          zflag = Calloc(2, long);

    if (nf>0) {
	frail = Calloc(nused, int);
	for (i=0; i<nused; i++) frail[i] = frail2[i];
        }

    /*
    ** Subtract the mean from each covar, as this makes the regression
    **  much more stable
    */
    for (i=0; i<nvar; i++) {
	temp=0;
	for (person=0; person<nused; person++) temp += covar[i][person];
	temp /= nused;
	means[i] = temp;
	for (person=0; person<nused; person++) covar[i][person] -=temp;
	}
	
    /*
    ** Find the loglik of the initial model
    **   (actually, just a no-sparse-terms model) -- loglik only
    */
    strata[nused-1] =1;
    *loglik = 0;

    for (person=0; person<nused; person++) {
        zbeta = 0;      /* form the term beta*z   (vector mult) */
        for (i=0; i<nvar; i++)
            zbeta += beta[i]*covar[i][person];
        score[person] = zbeta + offset[person];  /* save this away */
        }

    for (person=0; person<nused;) {
	if (event[person] ==0) person++;
	else {
            /*
            ** compute the mean over the risk set (ac)
            */
            denom =0;
            efron_wt =0;
            meanwt =0;

            time = stop[person];
            deaths=0;
            for (k=person; k<nused; k++) {
                if (start[k] < time) {
                    end[person] = k;     /*speed up -- last obs in risk set */
                    risk = exp(score[k]) * weights[k];
                    denom += risk;
                    if (stop[k]==time && event[k]==1) {
                        deaths += event[k];
                        efron_wt += risk*event[k];
                        meanwt += weights[k];
                        }
                     }
                if (strata[k]==1) break;
                }
 
            /*
            ** Add up the loglik
            */
            meanwt /= deaths;
            itemp = -1;
            endp = end[person];
            for (k=person; k<= endp && stop[k]==time; k++) {
                if (event[k]==1) {
                    itemp++;
                    temp = itemp*method/deaths;
                    d2 = denom - temp*efron_wt;
                    *loglik +=  weights[k]*score[k] -meanwt *log(d2);
                    }
                person++;
                }
	    }
	}   /* end  of accumulation loop */

    /*
    ** add in the penalty terms
    */
    if (ptype==2 || ptype==3) {
	/* there are non-sparse terms */
	cox_callback(2, beta, upen, ipen, &logpen, zflag);
	*loglik += logpen;
        }
    }

/********************************************************************/

/*
** This call is used for iteration
*/

void agfit4_b(long *maxiter, long *nusedx, long *nvarx, 
	       double *beta, double *u,
	       double *imat2,  double *jmat2, double *loglik, 
	       long *flag,  double *eps, double *tolerch, long *methodx, 
	       long *nfrail, double *fbeta, double *fdiag)
{
S_EVALUATOR
    int i,j,k, person;
    int ii;
    int     iter;
    int     nused, nvar;
    int    nf, nvar2;
    int    fgrp;
    int    halving;
    int    itemp, endp, deaths;

    double  denom, zbeta, risk;
    double  temp, temp2;
    double  newlk;
    double  d2, efron_wt;
    double  meanwt, time;
    double  method;
    double  **imat, **jmat;

    nused = *nusedx;
    nvar  = *nvarx;
    nf= *nfrail;
    method= *methodx;
    nvar2 = nvar + nf;
    if (nvar >0) {
        imat  = dmatrix(imat2, nvar2, nvar);
	jmat  = dmatrix(jmat2, nvar2, nvar);
        }
    else {
	imat = 0;   /*never used, but passed as dummy to chol */
	jmat = 0;
        }

    for (i=0; i<nf; i++) oldbeta[i] = fbeta[i];
    for (i=0; i<nvar; i++) oldbeta[i+nf] = beta[i];

    halving =0 ;             /* =1 when in the midst of "step halving" */
    for (iter=0; iter<=*maxiter; iter++) {
	newlk = 0;
	for (i=0; i<nf; i++) fdiag[i] =0;
	for (i=0; i<nvar2; i++) {
	    u[i] =0;
	    for (j=0; j<nvar; j++)
		    jmat[j][i] =0 ;
            }

        for (person=0; person<nused; person++) {
	    if (nf>0) {
		fgrp = frail[person] -1;
		zbeta = offset[person] + fbeta[fgrp];
	        }
	    else zbeta = offset[person];
	    for (i=0; i<nvar; i++)
		zbeta += beta[i]*covar[i][person];
	    score[person] = zbeta;
	    }
  
	for (person=0; person<nused; ) {
	    if (event[person]==0) person++;
	    else {
		endp = end[person];  /* shorter loops than 1 to n */
                /*
                ** compute the mean and covariance over the risk set (a and c)
                */
                efron_wt =0;
                denom =0;
                meanwt =0;
                for (i=0; i<nvar2; i++) {
                    a[i] =0;
                    a2[i]=0;
		    for (j=0; j<nvar; j++) {
			cmat[j][i] = 0;
			cmat2[j][i]= 0;
                        }
		    }
                time = stop[person];
                deaths=0;
                for (k=person; k<=endp; k++) {
                    if (start[k] < time) {
                        risk = exp(score[k]) * weights[k];
                        denom += risk;

			if (nf>0) {
			    fgrp = frail[k] -1;
			    a[fgrp] += risk;
			    }
			for (i=0; i<nvar; i++) {
			    a[i+nf] += risk*covar[i][k];
			    if (nf>0) cmat[i][fgrp] += risk*covar[i][k];
			    for (j=0; j<=i; j++)
				    cmat[i][j+nf] += covar[i][k]*covar[j][k] *
					    risk;
			    }
                        if (stop[k]==time && event[k]==1) {
                            deaths += event[k];
                            efron_wt += risk;
                            meanwt += weights[k];
			    if (nf>0) {
				u[fgrp] += weights[k];
				a2[fgrp] += risk;
			        }
			    for (i=0; i<nvar; i++) {
				u[i+nf] += weights[k] *covar[i][k];
				a2[i+nf] +=  risk*covar[i][k];
				if (nf>0) cmat2[i][fgrp] += risk*covar[i][k];
				for (j=0; j<=i; j++)
				    cmat2[i][j+nf] +=covar[i][k]*covar[j][k] *
					                  risk;
			        }
			    }
		        }
		    }	
 
                itemp = -1;
                meanwt /= deaths;
                for (k=person; k<=endp && stop[k]==time; k++) {
		    person++;
                    if (event[k]==1) {
                        itemp++;
                        temp = itemp*method/deaths;
                        d2 = denom - temp*efron_wt;
                        newlk +=  weights[k]*score[k] -meanwt *log(d2);
			for (i=0; i<nvar2; i++) {  /* by row of full matrix */
			    temp2 = (a[i] - temp*a2[i])/d2;
			    tmean[i] = temp2;
			    u[i] -= meanwt*temp2;
			    if (i<nf) fdiag[i] += temp2 * (1-temp2);
			    else {
				ii = i-nf;  /*actual row in c/j storage space*/
				for (j=0; j<=i; j++) 
					jmat[ii][j] +=  meanwt*(
					(cmat[ii][j] - temp*cmat2[ii][j]) /d2 -
                                          temp2*tmean[j]);
			        }
			    }
                        }
                    }
                }
            }   /* end  of accumulation loop */

	/*
	** Add in the penalty terms
        */
	if (ptype==1 || ptype==3) {
	    /* there are sparse terms */
    	    cox_callback(1, fbeta, upen, ipen, &logpen, zflag); 
	    if (zflag[0] ==1) {  /* force terms to zero */
		for (i=0; i<nf; i++) {
		    u[i]=0;
		    fdiag[i] =1;
		    for (j=0; j<nvar; j++) jmat[j][i]=0;
		    }
	        }
	    else {
		for (i=0; i<nf; i++) {
		    u[i] += upen[i];
		    fdiag[i] += ipen[i];
		    }
		newlk += logpen;
	        }
	    }

	if (ptype==2 || ptype==3) {
	    /* there are non-sparse terms */
	    cox_callback(2, beta, upen, ipen, &logpen, zflag);
	    newlk += logpen;
	    if (pdiag==0) {
		for (i=0; i<nvar; i++) {
		    u[i+nf] += upen[i];
		    jmat[i][i+nf] += ipen[i];
		    }
	        }
	    else {
		k =0;
		for (i=0; i<nvar; i++) {
		    u[i+nf] += upen[i];
		    for (j=nf; j<nvar2; j++) jmat[i][j] += ipen[k++];
		    }
	        }
	    for (i=0; i<nvar; i++) {
		if (zflag[i] ==1) {
		    u[i+nf]=0;
		    for (j=0; j<i; j++) jmat[i][j+nf]=0;
		    jmat[i+nf][i] =1;
		    }
	        }
	    }

	/* am I done?
	**   update the betas and test for convergence
	*/
	*flag = cholesky3(jmat, nvar2, nf, fdiag, *tolerch);
	if (fabs(1-(*loglik/newlk))<=*eps ) { /* all done */
	    *loglik = newlk;
	    for (i=0; i<nvar; i++) {
	        for (j=0; j<nvar2; j++)  imat[i][j] = jmat[i][j];
	        }
	    chinv3(jmat, nvar2, nf, fdiag);
	    for (i=nf; i<nvar2; i++) {       /*nicer output for S user */
	        fdiag[i] = jmat[i-nf][i];
	        jmat[i-nf][i] =1;
		imat[i-nf][i] =1;
		for (j=i+1; j<nvar2; j++) {
		    jmat[i-nf][j] = 0;
		    imat[i-nf][j] = 0;
		    }
	        }

	    if (halving==1) *flag= 1000; /*didn't converge after all */
	    *maxiter = iter;
	    return;
	    }

	if (iter==*maxiter) break;  /*skip the step halving and etc */

	if (iter>0 && newlk < *loglik)   {    /*it is not converging ! */
		halving =1;
		for (i=0; i<nvar; i++)
		    beta[i] = (oldbeta[i+nf] + beta[i]) /2; 
		for (i=0; i<nf; i++)
		    fbeta[i] = (oldbeta[i] + fbeta[i])/2;
		}
	    else {
		halving=0;
		*loglik = newlk;
		chsolve3(jmat,nvar2, nf, fdiag, u);

		j=0;
		for (i=0; i<nvar; i++) {
		    oldbeta[i] = beta[i];
		    beta[i] += u[i+nf];
		    }
		for (i=0; i<nf; i++) {
		    oldbeta[i] = fbeta[i];
		    fbeta[i] += u[i];
		    }
		}
	}   /* return for another iteration */

    *loglik = newlk;
    for (i=0; i<nvar; i++) 
	for (j=0; j<nvar2; j++) {
            imat[i][j] = jmat[i][j];
        }
    chinv3(jmat, nvar2, nf, fdiag);
    for (i=nf; i<nvar2; i++) {       /*nicer output for S user */
	fdiag[i] = jmat[i-nf][i];
	jmat[i-nf][i] =1;
  	imat[i-nf][i] =1;
  	for (j=i+1; j<nvar2; j++) {
	    jmat[i-nf][j] = 0;
  	    imat[i-nf][j] = 0;
  	    }
        }
    *flag= 1000;
    return;
    }


static double **cmatrix(double *data, int ncol, int nrow)
    {
S_EVALUATOR
    int i,j;
    double **pointer;
    double *temp;
 
    pointer = Calloc(nrow, double *);
    temp =    Calloc(nrow*ncol, double);
    if (data==0){
	for (i=0; i<nrow; i++) {
	    pointer[i] = temp;
	    temp += ncol;
            }
        }
    else {
	for (i=0; i<nrow; i++) {
	    pointer[i] = temp;
	    for (j=0; j<ncol; j++) *temp++ = *data++;
	    }
        }
    return(pointer);
    }

static void cmatrix_free(double **data) 
{
    Free(*data);
    Free(data);
    }


void agfit4_c (long *nusedx, long *nvar, long *methodx, double *expect) {
S_EVALUATOR
    double hazard, 
           denom,
           temp;
    double time, e_denom, wtsum;
    double e_hazard;
    double deaths;

    int    person,
           nused,
           method,
           k;
    int    endp;

    nused = *nusedx;
    method= *methodx;

    for (person=0; person<nused;) {
        if (event[person]==0) person++;
        else {
            denom =0;
            e_denom =0;
            wtsum =0;
            time = stop[person];
            deaths=0;
	    endp = end[person];
            for (k=person; k<=endp; k++) {
                if (start[k] < time) {
                    denom += score[k]*weights[k];
                    if (stop[k]==time && event[k]==1) {
                        deaths++;
                        wtsum += weights[k];
                        e_denom += score[k]*weights[k];
                        }
                     }
                }

            /*
            ** Do "expected" for the risk set
            */
            hazard =0;
            e_hazard=0;
            wtsum /=deaths;
            for (k=0; k<deaths; k++) {
                temp = method *(k/deaths);
                hazard += wtsum/(denom - temp*e_denom);
                e_hazard += wtsum*(1-temp)/(denom - temp*e_denom);
                }
            for (k=person; k<=endp; k++) {
                if (start[k] < time) {
                    if (stop[k]==time && event[k]==1)
                            expect[k] += score[k]*e_hazard;
                    else    expect[k] += score[k]*hazard;
                    }
                if (stop[k]==time) person++;
                }
            }
        }
    
    /*
    ** Free up the extra memory
    */
    Free(zflag);
    Free(upen);
    Free(event);
    Free(a);
    if (*nvar > 0) {
	cmatrix_free(cmat2);
	cmatrix_free(cmat);
	cmatrix_free(covar);
        }
    }
