.BG
.FN cox.zph
.TL
Test the Proportional Hazards Assumption of a Cox Regression
.DN
Test the proportional hazards assumption for a Cox regression model fit
(`coxph').
.CS
cox.zph(fit, transform="km", global=T)
.RA
.AG fit
the result of fitting a Cox regression model, using the `coxph' function.
.OA
.AG transform
a character string specifying how the survival times should be transformed
before the test is performed.
Possible values are `"km"', `"rank"', `"identity"' or  a
function of one argument.
.AG global
should a global chi-square test be done, in addition to the
per-variable tests.
.AG x
if true, then the result will be a list containing the test table (a matrix),
`x' and `y'.  If false then only the test table is returned.
.RT
an object of class `"cox.zph"', with components:
.AG table
a matrix with one row for each variable, and optionally a last row for
the global test.
Columns of the matrix contain the correlation coefficient between transformed
survival time and the scaled Schoenfeld residuals, a chi-square,
and the two-sided p-value.
For the global test there is no appropriate correlation, so an NA is
entered into the matrix as a placeholder.
.AG x
the transformed time axis.
.AG y
the matrix of scaled Schoenfeld residuals.  There will be one column per
variable and one row per event.  The row labels contain the original event
times (for the identity transform, these will be the same as `x').
.AG call
the calling sequence for the routine.
.PP
The computations require the original `x' matrix of the Cox model fit.
Thus it saves time if the `x=T' option is used in `coxph'.
This function would usually be followed by both a plot and a print of the
result.
The plot gives an estimate of the time-dependent coefficient `beta(t)'.
If the proportional hazards assumption is true, `beta(t)' will be a horizontal
line.  The printout gives a test for `slope=0'.
.SH REFERENCES
P. Grambsch and T. Therneau (1994),
Proportional hazards tests and diagnostics based on weighted residuals.
.ul
Biometrika,
\fB81\fR, 515-26.
.SA
`coxph', `Surv'.
.EX
fit <- coxph( Surv(futime, fustat) ~ age + ecog.ps, ovarian, x=T)
temp<- cox.zph(fit)
print(temp)                  #display the results
plot(temp)                   #plot curves
.KW survival
.WR
