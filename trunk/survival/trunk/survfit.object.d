.BG D
.FN survfit.object
.TL
Survival Curve Object
.DN
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
the `time' etc. vectors corresponding to the first curve, the second curve,
and so on.  The names of the elements are labels for the curves.
.AG std.err
the standard error of the cumulative hazard or -log(survival).
.AG upper
upper confidence limit for the survival curve.
.AG lower
lower confidence limit for the survival curve.
.AG conf.type
the approximation used to compute the confidence limits.
.AG conf.int
the level of the confidence limits, e.g. 90 or 95%.
.AG na.action
the returned value from the na.action function, if any.  It will be used
in the printout of the curve, e.g., the number of observations deleted due
to missing values.
.AG call
the call that produced the object.
.SH SUBSCRIPTS
Survfit objects that contain multiple survival curves can be subscripted.
This is most often used to plot a subset of the curves.
Usually a single subscript will be used.  In one particular case \(em
survival curves for multiple covariate values, from a Cox model that includes
a `strata' statement \(em there is a matrix of curves and 2 subscripts may
be used.
(In this case `summary.survfit' will also print the data as a matrix).
.SA
`survfit', `plot.survfit', `summary.survfit'.
.KW survival
.WR
