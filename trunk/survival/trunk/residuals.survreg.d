.BG
.FN residuals.survreg
.TL
Compute Residuals for `survreg' Objects
.DN
This is a method for the function `residuals' for objects inheriting from
class `survreg'.  
.CS
residuals.survreg(object, type, rsigma=T, collapse=F, weighted=F)
.RA
.AG object
an object inheriting from class `survreg'.
.AG type
type of residuals, with choices of `"response"', `"deviance"',
`"dfbeta"', `"dfbetas"', `"working"', `"ldcase"', `"lsresp"',
`"ldshape"', and `"matrix"'.  See the LaTeX documentation for more
detail.
.OA
.AG rsigma
include the scale parameters in the variance matrix, when doing computations.
(I can think of no good reason not to).
.AG collapse
optional vector of subject groups.  If given, this must be of the same
length as the residuals, and causes the result to be per group residuals.
.AG weighted
give weighted residuals?  Normally residuals are unweighted.
.RT
A vector or matrix of residuals is returned.
Response residuals are on the scale of the original data,
working residuals are on the scale of the linear predictor,
and deviance residuals are on log-likelihood scale.
The dfbeta residuals are a matrix, where the ith row gives the
approximate change in the coefficients due to the addition of subject i.
The dfbetas matrix contains the dfbeta residuals, with each column
scaled by the standard deviation of that coefficient.
.PP
The matrix type produces a matrix based on derivatives of the log-likelihood
function.  Let L be the log-likelihood, p be the linear predictor X %*% coef,
and s be log(sigma).  Then the 6 columns of the matrix are L, dL/dp,
ddL/(dp dp), dL/ds, ddL/(ds ds) and ddL/(dp ds), where d stands for the
derivative and dd the second derivative.  Diagnostics based on these quantities
are discussed in an article by Escobar and Meeker.
The main ones are the likelihood displacement residuals for perturbation
of a case weight (ldcase), the response value (ldresp), and the shape.
.SH REFERENCES
Escobar, L. A. and Meeker, W. Q. (1992).
Assessing influence in regression analysis with censored data.
.ul
Biometrics
\fB48\fR, 507-528.
.EX
# plot figures 2 and 6 from Escobar and Meeker
fit <- survreg(Surv(time,status) ~ age + age^2, data=stanford2,
	dist='lognormal')
rr <- residuals(fit, type='ldcase')
tsplot(rr, xlab="Case Number", ylab='1/2 A')
title(main="Case weight perturbations")

rr <- residuals(fit, type='ldresp')
tsplot(rr, xlab="Case Number", ylab='1/2 A')
title(main="Response pertubations")
.KW survival
.WR
