.BG
.FN summary.survfit
.TL
Print a Survival Curve
.DN
Returns a list containing the survival curve, confidence limits for the
curve, and other information.
.CS
summary.survfit(fit, times, censored=F, scale=1)
.RA
.AG fit
output from a call to `survfit'.
.OA
.AG times
vector of the times;
the returned matrix will contain 1 row for each time.
This must be in increasing order and missing values are not allowed.
If censored=T, the default time vector contains all the unique times in
`fit.list',
otherwise the default time vector uses only the event (death) times.
.AG censored
logical flag: should the censoring times be included in the output?
This is ignored if the `times' argument is present.
.AG scale
rescale the survival time, e.g., if the input data to survfit were in
days, "scale=365" would scale the printout to years.
.RT
a list with the following components
.AG time
the timepoint on the curve.
.AG surv
the value of the survival curve at time t+0.
.AG n.risk
the number of subjects at risk at time t-0
(but see the comments on weights in the `survfit' help file).
.AG n.event
if the `times' argument is missing, then this column is the number of
events that occurred at time t.
Otherwise, it is the cumulative number of events that have occurred
since the last time listed until time t+0.
.AG std.dev
the standard deviation of the survival value.
.AG lower CI
lower confidence limits for the curve.
.AG upper CI
upper confidence limits for the curve.
.AG call
the statement used to create the `fit' object.
.AG na.action
passed through from `fit', if present.
.SA
`survfit', 'print.summary.survfit'.
.EX
summary( survfit( futime, fustat))
.KW survival
.WR
