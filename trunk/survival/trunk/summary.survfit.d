.BG
.FN summary.survfit
.TL
Summary of a Survival Curve
.DN
Returns a list containing the survival curve, confidence limits for the
curve, and other information.
.CS
summary.survfit(fit, times=<<see below>>, censored=F, scale=1, extend=F)
.RA
.AG fit
the result of a call to the `survfit' function.
.OA
.AG times
vector of times;
the returned matrix will contain 1 row for each time.
This must be in increasing order and missing values are not allowed.
If `censored=T', the default `times' vector contains all the unique times in
`fit',
otherwise the default `times' vector uses only the event (death) times.
.AG censored
logical value: if TRUE, the censoring times are included in the output.
This is ignored if the `times' argument is present.
.AG scale
numeric value to rescale the survival time, e.g., if the input data to
`survfit' were in
days, `scale = 365.25' would scale the output to years.
.AG extend
logical value: if TRUE, prints information for all specified 'times',
even if there are no subjects left at the end of the specified 'times'.
This is only valid if the 'times'
argument is present.
.RT
a list with the following components:
.AG surv
the estimate of survival at time t+0.
.AG time
the timepoints at which the curve has a step.
.AG n.risk
the number of subjects at risk at time t-0
(but see the comments on weights in the `survfit' help file).
.AG n.event
if the `times' argument is missing, this column is the number of
events that occurred at time t.
Otherwise, it is the cumulative number of events that have occurred
since the last time listed until time t+0.
.AG n.entered
if the 'type' argument is 'counting' and the `times' argument is
missing, this column is the number of subjects that entered at time t.
Otherwise, it is the cumulative number of subjects that have entered
since the last time listed until time t+0.  This is only valid if the
'type' argument is 'counting'.
.AG n.exit.censored
if the 'type' argument is 'counting' and the `times' argument is
missing, this column is the number of subjects that left without an
event at time t.
Otherwise, it is the cumulative number of subjects that have left
without an event
since the last time listed until time t+0.  This is only valid if the
'type' argument is 'counting'.
.AG std.err
the standard error of the survival value.
.AG conf.int
level of confidence for the confidence intervals of survival.
.AG lower CI
lower confidence limits for the curve.
.AG upper CI
upper confidence limits for the curve.
.AG strata
indicates stratification of curve estimation.  If `strata' is not `NULL',
there are multiple curves in the result and the `surv', `time', `n.risk', etc. 
vectors will contain multiple curves, pasted end to end. 
The levels of `strata' (a factor) are the labels for the curves.
.AG call
the statement used to create the `fit' object.
.AG na.action
same as for `fit', if present.
.AG table
table of information that is returned from 'print.survfit' function.
.AG type
type of data censoring.  Passed from 'survfit' function.
.SA
`survfit', `print.summary.survfit'.
.EX
summary(survfit(Surv(time, status) ~ group, data = leukemia)) 
.KW survival4
.WR
