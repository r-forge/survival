.BG
.FN print.survfit
.TL
Short summary of a survival curve
.DN
Print n, number of event, mean and median survival.
.CS
print.survfit(fit, scale=1)
.RA
.AG fit
the result of a call to the survfit function.
.OA
.AG scale
rescale the survival time, e.g., if the input data to survfit were in days,
`scale=365' would scale the printout to years.
.RT
x, with the invisible flag set to prevent printing.
.SE
the number of observations, the number of events, the mean survival and its
standard error, and the median survival with its confidence interval are
printed.  If there are multiple curves, there is one line of output for each.
.DT
The mean and its variance are based on a truncated estimator.  That is, if the
last observation(s) is not a death, then the survival curve estimate does not
go to zero and the mean is undefined.  In such a case, the estimator is based
on an assumption that the true curve goes to zero just beyond the last
observed follow up time; it will systematically underestimate the true mean.
.pp
The median and its confidence interval are defined by drawing a horizontal
line at 0.5 on the plot of the survival curve and it's confidence bands.
The intersection of the line with the lower CI band defines the lower limit
for the median's interval, and similarly for the upper band.  If any of the
intersections is not a point, then we use the smallest point of intersection,
e.g., if the survival curve were exactly equal to 0.5 over an interval.
.SH REFERENCES
Miller, Rupert G Jr. (1981) Survival Analysis, Wiley, New York, p 71.
.SA
summary.survfit
.KW survival
.WR
