.BG
.FN cox.zph
.TL
Test the proportional hazards assumption
.CS
cox.zph(fit, transform='rank', global=F)
.AG fit
the result of fitting a Cox regression model, using the coxph function.
.AG transform
survival times transformed before the test is performed.  Possible values
are the character strings 'rank', 'identity', or a function of one argument.
.AG global
either a global chisquare test, or per-variable tests will be done.
.RT
a matrix with one row for each variable or a vector if the
global test.
Columns of the matrix contain the correlation coefficient between transformed
survival time and the scaled Schoenfeld residuals, a chisquare,
and the two-sided p value.
.pp
This function would usually be executed in conjunction with a plot of the
scaled Schoenfeld residuals versus survival time.
.SH REFERENCES
Terry Therneau, author of local function.
.EX
fit <- coxph( Surv(futime, fustat) ~ age + surgery, jasa)
cox.zph(fit)
.KW survival
.WR
