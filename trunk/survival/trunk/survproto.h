/*
**  SCCS @(#)survproto.h	5.5 02/21/99
** prototypes of all the survival functions
**  along with a few macros
*/

void agexact(long *maxiter,  long *nusedx,   long *nvarx,   double *start, 
	     double *stop,   long *event,    double *covar2,double *offset, 
	     long   *strata, double *means,  double *beta,  double *u, 
	     double *imat2,  double loglik[2], long *flag,  double *work, 
	     long   *work2,  double *eps,    double *tol_chol, double *sctest);

void agfit2( long   *maxiter,  long   *nusedx,  long   *nvarx, 
	     double *start,    double *stop,    long   *event, 
	     double *covar2,   double *offset,  double *weights,
	     long   *strata,   double *means,   double *beta, 
	     double *u,        double *imat2,   double loglik[2], 
	     long   *flag,     double *work,    long   *end,
	     double *eps,      double *tol_chol,double *sctest);

void agfit4_a(long *nusedx, long *nvarx, double *yy, 
	       double *covar2, double *offset2,
	       double *weights2, long *strata2,
	       double *means, double *beta, double *u, 
	       double *loglik, 
	       long *methodx, long *ptype2, long *pdiag2,
	       long *nfrail,  long *frail2);

void agfit4_b(long *maxiter, long *nusedx, long *nvarx, 
	       double *beta, double *u,
	       double *imat2,  double *jmat2, double *loglik, 
	       long *flag,  double *eps, double *tolerch, long *methodx, 
	       long *nfrail, double *fbeta, double *fdiag);

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
int cholesky2(double **matrix, int n, double toler);
void chsolve2(double **matrix, int n, double *y);
void chinv3(double **matrix , int n, int m, double *fdiag);
int cholesky3(double **matrix, int n, int m, double *diag, double toler);
void chsolve3(double **matrix, int n, int m, double *diag, double *y);

void coxdetail(long   *nusedx,   long   *nvarx,    long   *ndeadx, 
	       double *y,        double *covar2,   long   *strata,  
	       double *score,    double *weights,  double *means2, 
	       double *u2,       double *var,      double *work);

void coxfit2(long   *maxiter,   long   *nusedx,    long   *nvarx, 
	     double *time,      long   *status,    double *covar2, 
	     double *offset,	double *weights,   long   *strata,
	     double *means,     double *beta,      double *u, 
	     double *imat2,     double loglik[2],  long   *flag, 
	     double *work,	double *eps,       double *tol_chol,
	     double *sctest);

void coxfit4_a(long *nusedx, long *nvarx, double *yy, 
               double *covar2, double *offset2,
               double *weights2, long *strata2,
               double *means, double *beta, double *u, 
               double *loglik, 
               long *methodx, long *ptype2, long *pdiag2,
               long *nfrail,  long *frail2);

void coxfit4_b(long *maxiter, long *nusedx, long *nvarx, 
               double *beta, double *u,
               double *imat2,  double *jmat2, double *loglik, 
               long *flag,  double *eps, double *tolerch, long *methodx, 
               long *nfrail, double *fbeta, double *fdiag);

void coxfit4_c (long *nusedx, long *nvar, long *methodx, double *expect);

void coxfit_null(long   *nusedx,    long   *method,   double *time, 
		 long   *status,    double *score,    double *weights, 
		 long   *strata,    double *loglik, double *resid);

void coxhaz2(long   *n,      double *score,   long   *mark, 
	     long   *strata, double *hazard,  double *cumhaz);

void coxmart(long   *sn,     long   *method,    double *time, 
	     long   *status, long   * strata,   double *score, 
	     double *wt,     double *expect);

void coxph_wtest(long *nvar2, long *ntest, double *var, double *b,
                 double *scratch, double *tolerch);

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

void init_coxcall1(long *ptr1, s_object **ptr2);
void init_coxcall2(long *ptr1, s_object **ptr2);
void cox_callback (int which, double *coef, double *first, 
                   double *second, double *penalty, long *flag);

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

void survdiff2(long   *nn,     long   *nngroup,    long   *nstrat, 
	       double *rho,    double *time,       long   *status, 
	       long   *group,  long   *strata,	   double *obs, 
	       double *exp,    double *var,        double *risk, 
	       double *kaplan);

void survfit2(long   *sn,      double *y,       double *wt,
	      long   *strata,  long   *method,  long   *error, 
	      double *mark,    double *surv,	double *varh,
	      double *risksum);

void survfit3(long   *sn,        double *y,               double *wt,
	      long   *strata,    long   *method,          long   *error, 
	      long   *nstrat,    double *ntimes_strata,  
	      double *timelist,  double *weighted_event,  double *surv,
	      double *varh,	 double *risksum,         double *enter,
	      double *exit_censored);

void survindex2(long   *n,          double *stime,      long   *strata,
		long   *ntime,      double *time,       long   *nstrat,
		long   *o_n_risk,   long   *o_n_event,  double *o_surv,
		double *o_std_err,  double *o_upper,    double *o_lower, 
		long   *n_risk,     long   *n_event,    double *surv,
		double *std_err,    double *upper,      double *lower,
		double *new_start,  long   *num_extend, long   *times_strata,
		double *temp_times);

void survindex3(long   *n,          double *stime,        long   *strata,
		long   *ntime,      double *time,         long   *nstrat, 
		long   *o_n_risk,   long   *o_n_entered,  long   *o_n_censored,
		long   *o_n_event,  double *o_surv,       double *o_std_err,
		double *o_upper,    double *o_lower,      long   *n_risk, 
		long   *n_entered,  long   *n_censored,   long   *n_event,
		double *surv,       double *std_err,      double *upper,
		double *lower,      double *new_start, 	  long   *num_extend,
		long   *times_strata,                     double *temp_times);

void survreg2(long   *maxiter,   long   *nx,    long   *nvarx, 
	     double *y,          long   *ny,    double *covar2, double *wtx,
	     double *offset2,    double *beta,  long   *nstratx, 
	     long   *stratax,    double *ux,    double *imatx, 
	     double *loglik,     long   *flag,  double *eps,
	     double *tol_chol,   long   *dist,  long   *ddebug);

void survreg4(long   *maxiter,   long   *nx,       long   *nvarx, 
	      double *y,         long   *ny,       double *covar2, 
	      double *wt2,       double *offset2,  double *beta,  
	     long   *nstratx,    long   *stratax,  double *ux,    
	     double *imatx,      double *jmatx,
	     double *loglik,     long   *flag,     double *eps,
	     double *tol_chol,   long   *dist,     long   *ddebug,
             long *ptype2,  	 long   *pdiag2,
	     long *nfrail2,      long   *frail2,   double *fdiag2);

void survreg3(long   *maxiter,   long   *nx,    long   *nvarx, 
	     double *y,          long   *ny,    double *covar2, double *wtx,
	     double *offset2,    double *beta,  long   *nstratx, 
	     long   *stratax,    double *ux,    double *imatx, 
	     double *loglik,     long   *flag,  double *eps,
	     double *tol_chol,   long   *dist,  long   *ddebug);

void survreg5(long   *maxiter,   long   *nx,       long   *nvarx, 
	      double *y,         long   *ny,       double *covar2, 
	      double *wt2,       double *offset2,  double *beta,  
	     long   *nstratx,    long   *stratax,  double *ux,    
	     double *imatx,      double *jmatx,
	     double *loglik,     long   *flag,     double *eps,
	     double *tol_chol,   long   *dist,     long   *ddebug,
             long *ptype2,  	 long   *pdiag2,
	     long *nfrail2,      long   *frail2,   double *fdiag2);
