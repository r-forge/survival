.BG
.FN lines.survfit
.TL
Add lines to a survival plot
.DN
Usually used to add the expected survival curve(s) to a Kaplan-Meier plot
generated with plot.survfit.
.CS
lines.survfit(x, type="s", mark=3, col=1, lty=1, lwd=1, mark.time=T,
xscale=1, yscale=1, ...)
.RA
.AG x
a survival object, generated either from the survfit or survexp functions.
.OA
.AG type
the line type, as described in `lines'.  The default is a step function
for survfit objects, and a connected line for survexp objects.
.AG mark, col, lty, lwd
vectors giving the mark symbol, color, line type and line width for the
added curves.
.AG mark.time
controls the labeling of the curves.  If False, no labeling is done.  If
True, then curves are marked at each censoing time.  If mark.time is a numeric
vector, then curves are marked at the specified time points.
.AG xscale
a number used to divide the x values.  If time was originally in days, a
value of 365.24 would give a plotted scale in years.
.AG yscale
a number used to multiply the y values.  A value of 100, for instance, will
give a y axis in percent.
.SE
one or more curves are added to the current plot.
.SH NOTE
Does not yet handle confidence intervals.
.SA
plot.survfit, survfit, survexp
.EX
fit <- survfit(Surv(futime, fustat) ~ surgery, jasa)
ptime <- ifelse(jasa$fustat==0, jasa$futime, mdy.date(4,1,74)-jasa$accept.dt)
age <- jasa$accept.dt - jasa$birth.dt
efit <- survexp(ptime~ratetable(age=age, year=accept.dt, sex=1) + surgery,
		data=jasa, ratetable=survexp.us, conditional=F,times=0:5*300)
plot(fit, col=1:2, lty=1)
lines(efit, col=1:2, lty=2, mark='E', mark.time=1000, cex=1.5)
# Note: the 2 expected curves are overlap
.KW survival
.WR
