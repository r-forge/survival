.BG
.FN lines.survfit
.TL
Add Lines to a Survival Plot
.DN
Often used to add the expected survival curve(s) to a Kaplan-Meier plot
generated with `plot.survfit'.
.CS
lines.survfit(x, type="s", mark=3, col=1, lty=1, lwd=1, mark.time=T,
	xscale=1,  firstx=0, firsty=1, xmax, fun, conf.int=F,  ...)
.RA
.AG x
a survival object, generated from the `survfit' or `survexp' functions.
.OA
.AG type
the line type, as described in `lines'.  The default is a step function
for `survfit' objects, and a connected line for `survexp' objects.
.AG mark, col, lty, lwd
vectors giving the mark symbol, color, line type and line width for the
added curves.
.AG mark.time
controls the labeling of the curves.  
If `FALSE', no labeling is done.  
If `TRUE', then curves are marked at each censoring time.  
If `mark.time' is a numeric vector, then curves are marked at 
the specified time points.
.AG xscale
a number used to divide the x values.  If time was originally in days, a
value of 365.24 would give a plotted scale in years.
.AG firstx, firsty
the starting point for the survival curves.  If either of these is set to
`NA' or <blank> the plot will start at the first time point of the curve.
.AG xmax
the maximum horizontal plot coordinate.    
This shortens the curve before plotting it, so unlike using the
`xlim' graphical parameter, warning messages about out of bounds points are
not generated.
.AG fun
an arbitrary function defining a transformation of the survival curve.
For example `fun=log' is an alternative way to draw a log-survival curve
(but with the axis labeled with log(S) values).
Four often used transformations can be specified with a character
argument instead: `"log"' is the same as using the `log=T' option,
`"event"' plots cumulative events (f(y) =1-y),
`"cumhaz"' plots the cumulative hazard function (f(y) = -log(y))
and `"cloglog"' creates a complementary log-log survival plot 
(f(y) = log(-log(y) along with log scale for the x-axis).
.AG conf.int
if TRUE, confidence bands for the curves are also plotted.
If set to `"only"', then only the CI bands are plotted, and the curve
itself is left off.  
This can be useful for fine control over the colors or line types of a plot.
.RT
a list with components `x' and `y', containing the coordinates of the
last point on each of the curves (but not of the confidence limits).
This may be useful for labeling.
.SE
one or more curves are added to the current plot.
.SA
`lines', `par', `plot.survfit', `survfit', `survexp'.
.EX
fit <- survfit(Surv(futime, status==2) ~ sex, pbc)
plot(fit, mark.time=F, xscale=365.24,
	xlab='Years', ylab='Survival')
lines(fit[1], lwd=2, xscale=365.24)    #darken the first curve and add marks

# add expected survival curves for the two groups,
#   based on the US expected
tdata <- data.frame(age=pbc$age*365.24, sex=pbc$sex +1, 
		    year= rep(mdy.date(1,1,1976), nrow(pbc)))
efit <- survexp(~ sex, data=tdata, ratetable=survexp.us, 
		    times=(0:24)*182)
temp <- lines(efit, lty=2, xscale=365.24, lwd=2:1)
text(temp, c("Male", "Female"), adj= -.1) #labels just past the ends

title(main="Primary Biliary Cirrhosis, Observed and Expected")
.KW survival
.WR
