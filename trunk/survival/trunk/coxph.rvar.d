.BG
.FN coxph.rvar
.TL
Robust variance for a Cox model
.DN
Computes the robust "sandwich" estimator of variance for a proportional
hazards model, and returns a copy of the first argument with the robust
variance added on.
.CS
coxph.rvar(fit, collapse)
.RA
.AG fit
a coxph object, i.e., the result of fitting a Cox model.
.OA
.AG collapse
if the original data contained correlated observations, e.g., multiple
data rows per subject, then this argument contains the id vector that
identifies the subgroups.
.RT
a copy of the input, with two components added
.AG robust.var
the robust variance estimate.
.AG rcall
the call to this function.
.SE
the print and summary methods for coxph recognize and use the robust
variance. The global likelihood ratio and score statistics are
uneffected, but the global Wald test will now be based on the robust
estimator.
.DT
Let r be the matrix of infinitesimal influence functions, i.e.,
r <- residuals(fit, type='dbeta').  Then the robust variance is
v <- t(r) %*% r.  If there are correlated observations, the appropriate rows
or r are first summed, and v is based on the reduced r matrix.  There is
an obvious connection with the ordinary and group jackknife estimates.
.SA
coxph
.EX
fit <- coxph(Surv(futime, fustat) ~ age + rx +ecog.ps, data=fleming)
fit2 <- coxph.rvar(fit)
summary(fit2)
.KW survival
.WR
