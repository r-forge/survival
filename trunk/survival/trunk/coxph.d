.BG
.FN coxph
.TL
Proportional Hazards Regression
.DN
Fits either a Cox proportional hazards model or an Anderson-Gill
formulation of a proportional hazards model for multiple events and/or
time dependent covariates.
.CS
coxph(formula=formula(data), data=sys.parent(), subset, 
       na.action, eps=0.0001, inf.ratio=200, init, 
       iter.max=10, method=c("breslow","efron","exact"),
       model=F, x=F, y=T)
.RA
.AG formula
A formula expression as for other survival functions, of the form
Surv( time, status) ~ predictors for the regular Cox proportional
hazards model or Surv( start, stop, status) ~ predictors for the
Anderson-Gill formulation.  For a stratified analysis put a call
to the strata function among the predictors.
.OA
.AG data
a data.frame in which to interpret the variables named in
the formula, or in the subset and the weights argument.
.AG subset
expression saying that only a subset of the rows of the data
should be used in the fit.
.AG na.action
a missing-data filter function, applied to the model.frame, after any
subset argument has been used.  Default is options()$na.action.
.AG eps
convergence criteria.  Iteration will continue until relative change
in log-likelihood is less than eps.  Default is .0001.
.AG inf.ratio
warning criterion.  In certain data cases the actual MLE estimate of a
coefficient is infinity, e.g., a dichotomous variable where one of the
groups has no events.  When this happens the associated coefficient
grows at a steady pace and a race condition will exist in the fitting
routine, either the log likelihood converges, the information matrix
becomes effectively singular, an argument to exp becomes too large for
the computer hardware, the maximum number of interactions is exceeded,
or "average" risk ratio associated with the covariate is deemed to be
greater then inf.ratio. In the latter case, a warning message is
printed.  Default is 200.
.AG init
vector of initial values of the iteration.  Default initial
value is zero for all variables.
.AG iter.max
maximum number of iterations to perform.  Default is 10.
.AG method
method for tie handling.  "breslow" is the default for historical
reasons.  "efron" is a more accurate method of handling ties.  "exact"
refers to the exact partial likelihood (conditional logistic).
.AG model,x,y
flags to control what is returned.  If these are true, then the model
frame, the model matrix, and/or the response is returned as components
of the fitted model, with the same names as the flag arguments. 
.RT 
an object of class "coxph"
.SE
The basic problem is that it is inefficient to return _everything_
one might ever need from a model fit.  The routines for a Kaplan-Meier
survival, differences in survival, and for expected survival return
all necessary information, and have no side effects.  For a Cox
model however ....
.DT
The proportional hazards model is usually expressed in terms of a
single survival time value for each person, with possible censoring.
Anderson and Gill reformulated the same problem as a counting process;
as time marches onward we observe the events for a subject, rather
like watching a Geiger counter.  Extensions to the Cox model are
readily apparent: time dependent covariates, multiple events, and
discontinuous observation periods for a single subject.  A martingale
formulation of the counting process leads to direct proofs of
asymptotic properties, and to a useful set of residuals.
.SH REFERENCES
Terry Therneau, author of local function.

P. Andersen and R. Gill. "Cox's regression model for
counting processes, a large sample study", Annals of Statistics, 
10:1100-1120, 1982.  

T.Therneau, P. Grambsch, and T.Fleming. "Martingale based residuals
for survival models", Biometrika, March 1990.
.SA
survfit, Surv, strata.
.EX
# Create the simplest test data set
#
> test1 <- list(time=  c(4, 3,1,1,2,2,3),
                status=c(1,NA,1,0,1,1,0),
                x=     c(0, 2,1,1,1,0,0),
                sex=   c(0, 0,0,0,1,1,1))

# fit the cox proportional hazards model
# a call to strata in the terms of the formula tells coxph to fit
# a separate hazard function to the different strata.
#
> coxph( Surv(time, status) ~ x + strata(sex), test1)

#
# Create a simple data set for a time-dependent model
#
> test2 <- list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
                event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
                x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) )

# Create a Surv object, using names off of the data.frame test2
# summary will print out a bit more information
#
> summary( coxph( Surv(start, stop, event) ~ x, test2))

.KW ~keyword
.WR
