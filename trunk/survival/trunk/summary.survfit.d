.BG
.FN summary.survfit
.TL
Print a Survival Curve
.DN
Returns a matrix containing the survival curve, confidence limits for the 
curve, and other information.
.CS
summary.survfit( fit.list, times, censored=F, print.it=T, scale=1,
	digits=3, ...)
.RA
.AG fit.list
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
.AG print.it
logical flag: should output be printed?  Default is true.
.AG scale
rescale the survival time, e.g., if the input data to survfit were in
days, "scale=365" would scale the printout to years.
.AG digits
How many significant digits are desired?  Default is 3.
.RT
a matrix with the following columns
.AG time
the timepoint on the curve.
.AG n.risk
the number of subjects at risk at time t-0
(but see the comments on weights in the `survfit' help file).
.AG n.event
if the `times' argument is missing, then this column is the number of
events that occurred at time t.
Otherwise, it is the cumulative number of events that have occurred
since the last time listed until time t+0.
.AG surv
the value of the survival curve at time t+0.
.AG std.dev
the standard deviation of the survival value.
.AG lower CI
lower confidence limits for the curve.
.AG upper CI
upper confidence limits for the curve.
.SA
`survfit'.
.EX
summary( survfit( futime, fustat))
.KW survival
.WR
