.BG
.FN plot.survfit
.TL
Plot method for survfit.
.CS
plot.survfit(survfit, conf.int=<<see below>>, mark.time=T,
 mark=3, col=1, lty=1, lwd=1, cex=1, log=F, yscale=1, xscale=1, xlab="",
 ylab="", xaxs='i', ...)
.RA
.AG survfit
structure returned by survfit.
.OA
.AG conf.int
determines whether confidence intervals will be plotted.  The default is to
do so if there is only 1 curve, i.e., no strata.
.AG mark.time
controls the labeling of the curves.  If set to False, no labeling is done.
If True, then curves are marked at each censoring time.  If mark.time is a
numeric vector, then curves are marked at the specified time points.
.AG mark
vector of mark parameters, which will be used to label the curves.
The `lines' help file contains examples of the possible marks.
The vector is reused cyclically if it is shorter than the number of curves.
.AG col
vector of colors.  The default value is 1.
.AG lty
vector of line types. The default value is 1.
.AG lwd
vector of line widths. The default value is 1.
.AG cex
parameter available to change the size of "mark".
Not a vector; all marks have the same size.
.AG log
logical value: should the y axis be on a log scale?
.AG yscale
will be used to multiply the labels on the y axis.
A value of 100, for instance, would be used to give a percent scale.
Only the labels are
changed, not the actual plot coordinates, so that adding a curve with
"lines(surv.exp(...))", say, will perform as it did without the yscale arg.
.AG yscale
will be used in a similar manner for labels on the x axis.  A value of
365.25 will give labels in years instead of the original days.
.AG xlab
label given to the x-axis.
.AG ylab
label given to the y-axis.
.AG xaxs
the x axis style, as listed in `par'.  The S default option of "r" leads to
curves that are "too far" from the y axis.  This is, of course, just a matter
of esthetic opinion.
.RT
a list with components x and y, containing the coordinates of the last point
on each of the curves.  This may be useful for labeling.
.SE
A plot of survival curves is produced, one curve for each strata.
.EX
plot.survfit(survfit(pgtime,pgstat,g2group),col=c(3,6),log=T)
.KW survival
.WR

