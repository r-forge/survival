.BG
.FN plot.survfit
.TL
Plot Method for `survfit'
.CS
plot.survfit(survfit, conf.int=<<see below>>, mark.time=T,
 mark=3, col=1, lty=1, lwd=1, cex=1, log=F, xscale=1, yscale=1, 
 firstx=0, firsty=1, xmax, ymin=0, fun,
 xlab="", ylab="", xaxs=<<see below>>, ...)
.RA
.AG survfit
an object of class `survfit', usually returned by the `survfit' function.
.OA
.AG conf.int
determines whether confidence intervals will be plotted.  The default is to
do so if there is only 1 curve, i.e., no strata.
.AG mark.time
controls the labeling of the curves.  If set to `FALSE', no labeling is done.
If `TRUE', then curves are marked at each censoring time which is not also 
a death time.  If `mark.time' is a
numeric vector, then curves are marked at the specified time points.
.AG mark
vector of mark parameters, which will be used to label the curves.
The `lines' help file contains examples of the possible marks.
The vector is reused cyclically if it is shorter than the number of curves.
.AG col
a vector of integers specifying colors for each curve.
The default value is 1.
.AG lty
a vector of integers specifying line types for each curve.
The default value is 1.
.AG lwd
a vector of numeric values for line widths. The default value is 1.
.AG cex
a numeric value specifying the size of the marks.
Not a vector; all marks have the same size.
.AG log
a logical value, if TRUE the y axis wll be on a log scale.
.AG yscale
a numeric value used to multiply the labels on the y axis.
A value of 100, for instance, would be used to give a percent scale.
Only the labels are
changed, not the actual plot coordinates, so that adding a curve with
"`lines(surv.exp(...))'", say, 
will perform as it did without the `yscale' argument.
.AG xscale
a numeric value used like 'yscale' for labels on the x axis. 
A value of 365.25 will give labels in years instead of the original days. 
.AG firstx, firsty
the starting point for the survival curves.  If either of these is set to
`NA' the plot will start at the first time point of the curve.
.PP
If 'new.start' argument is used in 'survfit', 'firstx' is set to the
value of 'new.start'.
.AG xmax
the maximum horizontal plot coordinate.  This can be used to shrink the range
of a plot.  It shortens the curve before plotting it, so that unlike using the
`xlim' graphical parameter, warning messages about out of bounds points are
not generated.
.AG ymin
lower boundary for y values.  Survival curves are most often drawn in the
range of 0-1, even if none of the curves approach zero.
The parameter is ignored if the `fun' argument is present, or if it has been
set to `NA'.
.AG fun
an arbitrary function defining a transformation of the survival curve.
For example `fun=log' is an alternative way to draw a log-survival curve
(but with the axis labeled with log(S) values),
and `fun=sqrt' would generate a curve on square root scale.
Four often used transformations can be specified with a character
argument instead: "log" is the same as using the `log=T' option,
"event" plots cumulative events (f(y) = 1-y), 
"cumhaz" plots the cumulative hazard function (f(y) = -log(y)), and
"cloglog" creates a complimentary log-log survival plot (f(y) =
log(-log(y)) along with log scale for the x-axis). 
.AG xlab
label given to the x-axis.
.AG ylab
label given to the y-axis.
.AG xaxs
the x axis style, as listed in 'par'. Survival curves are traditionally 
drawn with the curve touching the bounding box on the left edge, but 
not touching it on the right edge.  This corresponds to neither of the 
two standard S axis styles of "e" (neither touches) or "i" (both 
touch).  If 'xaxis' is missing or NULL the internal axis style is used 
(xaxs='i') but only after the right endpoint has been extended by 5%. 
.RT
a list with components `x' and `y', containing the coordinates of the last point
on each of the curves (but not the confidence limits).  
This may be useful for labeling.
.SE
A plot of survival curves is produced, one curve for each strata.
.SH DETAILS
The `log=T' option does extra work to avoid log(0), and to try to create a
pleasing result.  If there are zeros, they are plotted by default at
0.8 times the smallest non-zero value on the curve(s).
.SH BUGS
Specifying 'mark.time=T' does not work if 'survfit' was created from  
a 'coxph' model.  In this case, censoring times are not marked on the  
curve(s).  If you want to mark censoring times on the curve(s) resulting  
from a `coxph' fit, provide a vector of times as the `mark.time' argument 
in the call to `plot.survfit' or `lines.survfit'. 
.SA
`par', `survfit', `lines.survfit'.
.EX
leukemia.surv <- survfit(Surv(time, status) ~ group, data = leukemia)
plot(leukemia.surv, lty = 2:3)
legend(100, .9, c("Maintenance", "No Maintenance"), lty = 2:3)
title("Kaplan-Meier Curves\nfor AML Maintenance Study")

lsurv2 <- survfit(Surv(time, status) ~ group, leukemia, type='fleming')
plot(lsurv2, lty=2:3, fun="cumhaz",
	xlab="Months", ylab="Cumulative Hazard")
.KW survival4
.KW hplot
.WR


