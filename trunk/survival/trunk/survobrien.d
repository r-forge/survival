.BG
.FN surv.obrien
.TL
Peter O'Brien's test for association of a single variable with survival
.DN
This test is proposed in Biometrics, June 1978.
.CS
surv.obrien(formula, data)
.RA
.AG formula
a valid formula for a cox model, without time dependent covariates.
.OA
.AG data
a data frame.
.RT
a new data frame.  The original time and status variables are removed,
and have been replaced with `start', `stop', and `event'.
If a predictor variable is a factor, it is retained as is.
Other predictor variables have been replaced with time-dependent logit
scores.
.pp
Because of the time dependent variables, the new data frame will have many
more rows that the original data, approximately #rows * #deaths /2.
.SH METHOD
A time-dependent cox model can now be fit to the new data.
The univariate statistic, as originally proposed, is equivalent to
single variable score tests from the time-dependent model.
This equivalence is the rationale for using the time dependent model as a
multivariate extension of the original paper.
.pp
In O'Brien's method, the x variables are re-ranked at each death time.  A
simpler method, proposed by Prentice, ranks the data only once at the
start.
.SH REFERENCES
O'Brien, Peter, "A Nonparametric Test for Association with Censored Data",
Biometrics 34: 243-250, 1978.
.SA
surv.diff
.KW survival
.EX
xx <- surv.obrien(Surv(time, status) ~ age + factor(rx) + ecog.ps,
			       data=fleming)
attach(xx)
coxreg(Surv(start, stop, event) ~ age)
coxreg(Surv(start, stop, event) ~ age + rx + ecog.ps)
.WR
