.BG D
.FN survfit.object
.TL
Parametric Survival Model Object
.PP
This class of objects is returned by the `survfit' class of functions
to represent a fitted survival curve.

Objects of this class have methods for the functions `print',
`summary', `plot', `points' and `lines'.
.SH COMPONENTS
The following components must be included in a legitimate `survfit' object.
.AG time
the time points at which the curve has a step.
.AG n.risk
the number of subjects at risk at t.
.AG n.event
the number of events that occur at time t.
.AG surv
the estimate of survival at time t+0.
This may be a vector or a matrix.
.AG strata
if there are multiple curves, this component gives the number of elements of
of the `time' etc. vectors corresponding to the first curve, the second curve,
and so on.  The names of the elements are labels for the curves.
.AG std.err
the standard error of the cumulative hazard or -log(survival).
.AG uppper
upper confidence limit for the survival curve.
.AG lower
lower confidence limit for the survival curve.
.AG conf.type
the approximation used to compute the confidence limits.
.AG conf.int
the level of the confident limits, e.g. 90 or 95%.
.AG na.action
the returned value from the na.action function, if any.  It will be used
in the printout of the curve, e.g., the number of observations deleted due
to missing values.
.AG call
an image of the call that produced the object.
.SA
`survfit', `plot.survfit', `summary.survfit'.
.KW survival
.WR
