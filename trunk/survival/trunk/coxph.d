.BG
.FN coxph
.TL
Fit Proportional Hazards Regression Model
.DN
Fits a Cox proportional hazards regression model.
Time dependent variables, time dependent strata, multiple events per subject,
and other extensions are incorporated using the counting process formulation
of Andersen and Gill.
.CS
coxph(formula, data=sys.parent(), subset,
       na.action, weights, eps=0.0001, init,
       iter.max=10, method=c("efron","breslow","exact"),
       singular.ok=T, robust,
       model=F, x=F, y=T)
.RA
.AG formula
a formula object, with the response on the left of a `~' operator, and
the terms on the right.  The response must be a survival object as
returned by the `Surv' function.
.OA
.AG data
a data.frame in which to interpret the variables named in
the `formula', or in the `subset' and the `weights' argument.
.AG subset
expression saying that only a subset of the rows of the data
should be used in the fit.
.AG na.action
a missing-data filter function, applied to the model.frame, after any
subset argument has been used.  Default is `options()$na.action'.
.AG weights
case weights.
.AG eps
convergence threshold.  Iteration will continue until the relative change
in the log-likelihood is less than eps.  Default is .0001.
.AG init
vector of initial values of the iteration.  Default initial
value is zero for all variables.
.AG iter.max
maximum number of iterations to perform.  Default is 10.
.AG method
a character string specifying the method for tie handling.  If there 
are no tied death times all the methods are equivalent.
Nearly all Cox regression programs use the Breslow method by default,
but
not this one.
The Efron approximation is used as the default here, as it is much more
accurate when dealing with tied death times, and is as efficient
computationally.
The exact method computes the exact partial likelihood, which is equivalent
to a conditional logistic model.  If there are a large number of ties the
computational time will be excessive.
.AG singular.ok
logical value indicating how to handle collinearity in the model matrix.
If `TRUE', the program will automatically skip over columns of the X matrix
that are linear combinations of earlier columns.  In this case the
coefficients for such columns will be NA, and the variance matrix will contain
zeros.  For ancillary calculations, such as the linear predictor, the missing
coefficients are treated as zeros.
.AG robust
if `TRUE' a robust variance estimate is returned.  Default is `TRUE' if the
model includes a `cluster()' operative, `FALSE' otherwise.
.AG model, x, y
flags to control what is returned.  If these are true, then the model
frame, the model matrix, and/or the response is returned as components
of the fitted model, with the same names as the flag arguments. 
.RT 
an object of class `"coxph"'. See `coxph.object' for details.
.SE
Depending on the call, the `predict', `residuals', and `survfit' routines may
need to reconstruct the x matrix created by `coxph'.  Differences in the
environment, such as which data frames are attached or the value of
`options()$contrasts', may cause this computation to fail or worse, to be
incorrect.  See the survival overview document for details.
.DT
The proportional hazards model is usually expressed in terms of a
single survival time value for each person, with possible censoring.
Andersen and Gill reformulated the same problem as a counting process;
as time marches onward we observe the events for a subject, rather
like watching a Geiger counter.
The data for a subject is presented as multiple rows or "observations", each
of which applies to an interval of observation (start, stop].
.SH SPECIAL TERMS
There are two special terms that may be used in the model equation.
A 'strata' term identifies a stratified Cox model; separate baseline hazard
functions are fit for each strata.
The `cluster' term is used to compute a robust variance for the model.
The term `+ cluster(id)', where `id == unique(id)', is equivalent to
specifying the `robust=T' argument, and produces an approximate jackknife
estimate of the variance.  If the `id' variable were not unique, but instead
identifies clusters of correlated observations, then the variance estimate
is based on a grouped jackknife.
.SH CONVERGENCE
In certain data cases the actual MLE estimate of a
coefficient is infinity, e.g., a dichotomous variable where one of the
groups has no events.  When this happens the associated coefficient
grows at a steady pace and a race condition will exist in the fitting
routine: either the log likelihood converges, the information matrix
becomes effectively singular, an argument to exp becomes too large for
the computer hardware, or the maximum number of interactions is exceeded.
The routine attempts to detect when this has happened, not always
successfully.
.SH REFERENCES
P. Andersen and R. Gill. "Cox's regression model for
counting processes, a large sample study", 
.ul
Annals of Statistics, 
10:1100-1120, 1982.  
.PP
T. Therneau, P. Grambsch, and T. Fleming. "Martingale based residuals
for survival models", 
.ul
Biometrika, 
March 1990.
.SA
`cluster,' `survfit', `Surv', `strata'.
.EX
# Create the simplest test data set
#
> test1 <- list(time=  c(4, 3,1,1,2,2,3),
                status=c(1,NA,1,0,1,1,0),
                x=     c(0, 2,1,1,1,0,0),
                sex=   c(0, 0,0,0,1,1,1))
> coxph( Surv(time, status) ~ x + strata(sex), test1)  #stratified model

#
# Create a simple data set for a time-dependent model
#
> test2 <- list(start=c(1, 2, 5, 2, 1, 7, 3, 4, 8, 8),
                stop =c(2, 3, 6, 7, 8, 9, 9, 9,14,17),
                event=c(1, 1, 1, 1, 1, 1, 1, 0, 0, 0),
                x    =c(1, 0, 0, 1, 0, 1, 1, 1, 0, 0) )

> summary( coxph( Surv(start, stop, event) ~ x, test2))

.KW survival
.WR
