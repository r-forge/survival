.BG
.FN plot.cox.zph
.TL
Graphical test of proportional hazards
.DN
Displays a graph of the scaled Shoenfeld residuals, along with a smooth.
.CS
plot.cox.zph(x, resid=T, se=T, df=4, nsmo=40)
.RA
.AG x
result of the `cox.zph' function.
.OA
.AG resid
include the residuals on the plot, as well as the smooth fit.
.AG se
if true, confidence bands at 2 standard errors will be added.
.AG df
the degrees of freedom for the natural spline, used as a smooth.
.AG nsmo
number of points used to plot the fitted spline.
.SE
a plot is produced on the current graphics device.
.SA
cox.zph, coxph
.KW survival
.WR
