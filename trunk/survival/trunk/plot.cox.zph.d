.BG
.FN plot.cox.zph
.TL
Graphical test of proportional hazards
.DN
Displays a graph of the scaled Shoenfeld residuals, along with a smooth.
.CS
plot.cox.zph(x, resid=T, se=T, df=4, nsmo=40, var)
.RA
.AG x
result of the `cox.zph' function.
.OA
.AG resid
include the residuals on the plot, as well as the smooth fit.
.AG se
if true, confidence bands at 2 standard errors will be added.
.AG df
the degrees of freedom for the fitted natural spline.  A df value of
2 leads to a linear fit.
.AG nsmo
number of points used to plot the fitted spline.
.AG var
the set of variables for which plots are desired.  By default, plots are
produced in turn for each variable of a model.  Selection of a single variable
allows other features to be added to the plot, e.g., a horizontal line at
zero or a main title.
.SE
a plot is produced on the current graphics device.
.SA
cox.zph, coxph
.EX
vfit <- coxph(Surv(futime,fustat) ~ rx + factor(celltype) + karno +age,
		   data=veteran, x=T)
temp <- cox.zph(vfit)
plot(temp, var=5)      #Look at Karnofsy score
abline(0,0, lty=3)
lines( lm( temp$y[,5] ~ temp$x), lty=4)   #Add the linear fit as well
title(main="VA Lung Study")
.KW survival
.WR
