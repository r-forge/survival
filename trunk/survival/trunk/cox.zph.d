.BG
.FN cox.zph
.TL
Test the proportional hazards assumption
.CS
cox.zph(fit, transform='km', global=T, x=F)
.AG fit
the result of fitting a Cox regression model, using the coxph function.
.AG transform
survival times transformed before the test is performed.  Possible values
are the character strings 'km', 'rank', 'identity', or a function of one argument.
.AG global
should a global chisquare test be done, in additoin to the
 per-variable tests.
.AG x
if true, then the result will be a list containing the test table (a matrix),
x and y.  If false then only the test table is returned.
.RT
.AG table
a matrix with one row for each variable, and optionally a last row for
the global test.
Columns of the matrix contain the correlation coefficient between transformed
survival time and the scaled Schoenfeld residuals, a chisquare,
and the two-sided p value.
For the global test there is no appropriate correlation, so an NA is
entered into the matrix as a placeholder.
.AG x
the transformed time axis.
.AG y
the matrix of scaled Schoenfeld residuals.  There will be one column per
variable and one row per event.  The row labels contain the original event
times, which will be the same as `x' for the identity transform.
.pp
The computations require the original `x' matrix of the Cox model fit.
Thus it saves time if the `x=T' option is used in coxph.

This function would usually be executed in conjunction with a plot of the
scaled Schoenfeld residuals versus survival time.  The plot gives an
estimate of the time-dependent coefficient beta(t).
.SH REFERENCES
Terry Therneau, author of local function.
.EX
fit <- coxph( Surv(futime, fustat) ~ age + surgery, jasa, x=T)
temp<- cox.zph(fit, x=T)
print(temp)                  #display the results
plot(temp)                   #plot curves
.KW survival
.WR
