.BG
.FN cox.zph
.TL
Test the proportional hazards assumption
.CS
cox.zph(fit, ranks=T, global=T)
.AG fit
the result of fitting a Cox regression model, using the coxph function.
.AG ranks
if true, survival times are replaced by thier ranks before the test is
performed.
.AG global
do the global test as well as the per-variable tests.
.RT
a matrix with one column for each variable and (optionally) a column for
the global test.
Rows of the matrix contain the correlation coefficent between survival time
or rank survival time and the Shoenfeld residuals, a normal score (based on
Fisher's Z transform), and the two-sided p value.
.pp
This function would usually be executed in conjunction with a plot of the
Shoenfeld residuals versus survival time.
.SH METHOD
The global test is based on the linear predictor lp <- X %*% coef.  It is
equivalent to the univariate test on a new fit with lp as the
covariate.  If the original model has only a single predictor, then the
global test is just a repeat of the univariate test for that predictor.
.EX
fit <- coxph(Surv(futime, fustat) ~ age + surgery, jasa)
cox.zph(fit)
.KW survival
.WR
