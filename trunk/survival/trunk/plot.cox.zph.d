.BG
.FN plot.cox.zph
.TL
Graphical Test of Proportional Hazards
.DN
Displays a graph of the scaled Schoenfeld residuals, along with a smooth curve.
.CS
plot.cox.zph(x, resid=T, se=T, df=4, nsmo=40, var, ...)
.RA
.AG x
result of the `cox.zph' function.
.OA
.AG resid
a logical value, if `TRUE' the residuals are included on the plot, as well as the smooth fit.
.AG se
a logical value, if `TRUE', confidence bands at two standard errors
will be added.
.AG df
the degrees of freedom for the fitted natural spline, `df=2' leads
to a linear fit.
.AG nsmo
number of points used to plot the fitted spline.
.AG var
the set of variables for which plots are desired.  By default, plots are
produced in turn for each variable of a model.  Selection of a single variable
allows other features to be added to the plot, e.g., a horizontal line at
zero or a main title.
.PP
This has been superseded by a subscripting method; see the example below.
.AG ...
additional arguments passed to the `plot' function.
.SE
a plot is produced on the current graphics device.
.SA
`cox.zph', `coxph'
.EX
veteran <- data.frame(cancer.vet)
vfit <- coxph(Surv(survival,status) ~ therapy + factor(celltype) +
              Karn..score + age, data=veteran, x=T)
temp <- cox.zph(vfit)
plot(temp, var=5)      # Look at Karnofsy score, old way of doing plot
plot(temp[5])     # New way with subscripting
abline(0, 0, lty=3)
# Add the linear fit as well 
abline(lm(temp$y[,5] ~ temp$x)$coefficients, lty=4, col=3) 
title(main="VA Lung Study")
.KW survival
.WR
