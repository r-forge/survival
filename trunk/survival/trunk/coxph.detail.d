.BG
.FN coxph.detail
.TL
Details of a Cox Model Fit
.DN
Details of a Cox model fit.
Returns the individual contributions to the first and second derivative
matrix, at each unique event time.
.CS
coxph.detail(fit)
.RA
.AG fit
a Cox model object, i.e., the result of `coxph'.
.RT
a list with components
.RC time
the vector of unique event times
.RC nevent
the number of events at each of these time points.
.RC means
a matrix with one row for each event time and one column for each variable
in the Cox model, containing the weighted mean of the variable at that time,
over all subjects still at risk at that time.  The weights are the risk
weights `exp(x %*% fit$coef)'.
.RC nrisk
number of subjects at risk.
.RC hazard
the hazard increment.
.RC score
the contribution to the score vector (first derivative of the log
partial likelihood) at each time point.
.RC imat
the contribution to the information matrix (second derivative of the
log partial likelihood) at each time point.
.RC varhaz
the variance of the hazard increment.
.RC "x, y, weights"
copies of the input data.
.RC strata
only present for a stratified Cox model, this is
a table giving the number of time points of component `time' that
were contributed by each of the strata.
.DT
This function may be useful for those who wish to investigate new methods or
extensions to the Cox model.  The example below shows one way to calculate
the Schoenfeld residuals.
.SA
`coxph', `residuals.coxph'
.EX
fit   <- coxph(Surv(futime,fustat) ~ age + rx + ecog.ps, fleming, x=T)
fitd  <- coxph.details(fit)
events <- fit$y[,2]==1
etime  <- fit$y[events,1]   #the event times --- may have duplicates
indx   <- match(etime, fitd$time)
sresid <- fit$x[events,] - fitd$means[indx,]
.KW survival
.WR
