/*
**  SCCS $Id: survproto.h,v 5.1 1998-08-30 14:46:43 therneau Exp $
** prototypes of all the survival functions
**  along with a few macros
*/


void agexact(long *maxiter,  long *nusedx,   long *nvarx,   double *start, 
	     double *stop,   long *event,    double *covar2,double *offset, 
	     long   *strata, double *means,  double *beta,  double *u, 
	     double *imat2,  double loglik[2], long *flag,  double *work, 
	     long   *work2,  double *eps,     double *sctest);

void agfit2( long   *maxiter,  long   *nusedx,  long   *nvarx, 
	     double *start,    double *stop,    long   *event, 
	     double *covar2,   double *offset,  double *weights,
	     long   *strata,   double *means,   double *beta, 
	     double *u,        double *imat2,   double loglik[2], 
	     long   *flag,     double *work,    long   *end,
	     double *eps,      double *sctest);

void agfit_null(long   *n,      long   *method,   double *start, double *stop, 
		long   *event,  double * offset,  double *weights,
		long   *strata, double loglik[2]);

void aghaz2(long   *n,     double *start,   double *stop,   long   *event, 
	    double *score, long   * strata, double *hazard, double * cumhaz);

void agmart(long   *n,     long   *method,  double *start,   double *stop, 
	    long   *event, double *score,   double *wt,      long   *strata, 
	    double *resid);

void agres12(long   *nx,     long   *nvarx,   double *y,    double *covar2, 
	     long   *strata, double *score,   long *method, double *resid2, 
	     double *a);

void agscore(long   *nx,       long   *nvarx,      double *y,
	     double *covar2,   long   *strata,     double *score,
	     double *weights,  long   *method,     double *resid2, double *a);

void agsurv1(long   *sn,     long   *snvar,  double *y,      double *score, 
	     long   *strata, double *surv,   double *varh,   long   *snsurv,
	     double *xmat,   double *d,      double *varcov, double *yy,
	     long   *shisn,  double *hisy,   double *hisxmat,double *hisrisk, 
	     long   *hisstrat);

void agsurv2(long   *sn,      long   *snvar,    double *y, 
	     double *score,   long   *strata,   double *surv, 
	     double *varh,    double *xmat,     double *varcov, 
	     long   *snsurv,  double *d,        long   *sncurve,
             double *newx,    double *newrisk);

void agsurv3(long   *sn,    long   *snvar,    long   *sncurve, 
	     long   *snpt,  long   *sse,      double *score, 
	     double *sy,    double *r,        double *coef, 
	     double *var,   double *cmean,    long   *scn, 
	     double *cy,    double *cx,       double *ssurv,
	     double *varh,  double *sused,    long   *smethod);

void chinv2  (double **matrix, int n);
int cholesky2(double **matrix, int n);
void chsolve2(double **matrix, int n, double *y);

void coxdetail(long   *nusedx,   long   *nvarx,    long   *ndeadx, 
	       double *y,        double *covar2,   long   *strata,  
	       double *score,    double *weights,  double *means2, 
	       double *u2,       double *var,      double *work);

void coxfit2(long   *maxiter,   long   *nusedx,    long   *nvarx, 
	     double *time,      long   *status,    double *covar2, 
	     double *offset,	double *weights,   long   *strata,
	     double *means,     double *beta,      double *u, 
	     double *imat2,     double loglik[2],  long   *flag, 
	     double *work,	double *eps,       double *sctest);

void coxfit_null(long   *nusedx,    long   *method,   double *time, 
		 long   *status,    double *score,    double *weights, 
		 long   *strata,    double *loglik, double *resid);

void coxhaz2(long   *n,      double *score,   long   *mark, 
	     long   *strata, double *hazard,  double *cumhaz);

void coxmart(long   *sn,     long   *method,    double *time, 
	     long   *status, long   * strata,   double *score, 
	     double *wt,     double *expect);

void coxres12(long   *nx,     long   *nvarx,    double *y, 
	      double *covar2, long   *strata,   double *score, 
	      long   *method, double *scratch);

void coxscho(long   *nusedx,    long   *nvarx,    double *y, 
	     double *covar2,    double *score,    long   *strata,  
	     long   *method2,   double *work);

void coxscore(long   *nx,      long   *nvarx,    double *y, 
	      double *covar2,  long   *strata,   double *score, 
	      double *weights, long   *method,   double *resid2,
	      double *scratch);

double **dmatrix(double *array, int ncol, int nrow);

void init_doloop(int min, int max);
int doloop      (int nloops, int *index);

void pyears1(long   *sn,      long   *sny,      long   *sdoevent, 
	     double *sy,      long   *sedim,    long   *efac, 
	     long   *edims,   double *secut,    double *expect, 
	     double *sedata,  long   *sodim,    long   *ofac, 
	     long   *odims,   double *socut,    long   *smethod, 
	     double *sodata,  double *pyears,   double *pn, 
	     double *pcount,  double *pexpect,  double *offtable);

void pyears2(long   *sn,      long   *sny,   long   *sdoevent, 
	     double *sy,      long   *sodim, long   *ofac, 
	     long   *odims,   double *socut, double *sodata,
	     double *pyears,  double *pn,    double *pcount, 
	     double *offtable);

void pyears3(long   *sdeath,    long   *sn,    long   *sedim, 
	     long   *efac,      long   *edims, double *secut, 
	     double *expect,    double *sx,    double *y, 
	     long   *sntime,    long   *sngrp, double *times,
	     double *esurv,     long   *nsurv);

double pystep(int nc,        int  *index,  int  *index2,   double *wt, 
	      double *data,  long *fac,    long *dims,     double **cuts, 
	      double step,   int  edge);

int rnewton(int    *maxiter,   int  n,        int  nvar,        double *beta, 
	    double *u,         double **imat, double loglik[2], double eps,
	    void (*dolk)(),    void (*doimat)(), 
	    double *newbeta,   double *savediag,  int debug);

void survdiff2(long   *nn,     long   *nngroup,    long   *nstrat, 
	       double *rho,    double *time,       long   *status, 
	       long   *group,  long   *strata,	   double *obs, 
	       double *exp,    double *var,        double *risk, 
	       double *kaplan);

void survfit2(long   *sn,     double *y,        long   *ny, 
	      double *wt,     long   *strata,   long   *method, 
	      long   *error,  double *mark,     double *surv,
	      double *varh,   double *risksum,  long   *snsurv);

void survindex2(long   *n,     double *stime,   long   *strata, 
		long   *ntime, double *time,    long   *nstrat, 
		long   *indx,  long   *indx2);

void survreg(long   *maxiter,    long   *nx,    long   *nvarx, 
	     double *y,          long   *ny,    double *covar2, 
	     double *offset2,    double *beta,  long   *npx, 
	     double *parmsx,     double *u,     double *imatx, 
	     double *loglik,     long   *flag,  double *eps,
	     double *deriv,      long   *dist);
