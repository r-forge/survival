.BG
.FN residuals.glm
.TL
Compute Residuals for `survreg' Objects
.CS
residuals.survreg(object, type)
.PP
This is a method for the function `residuals()' for objects inheriting from
class `glm'.  Several types of residuals are available for `glm' objects,
hence the additional argument:
.AG type
type of residuals, with choices `"deviance"', `"pearson"', `"working"' or
`"matrix"; the first is the default.
.RT
A vector of residuals is returned.
The sum of squared deviance residuals
add up to the deviance.  The pearson residuals are standardized residuals
on the scale of the response.  The working residuals reside on the object,
and are the residuals from the final IRLS fit.
.pp
The matrix type produces a matrix based on derivatives of the log-likelihood
function.  Let L be the log-likelihood, p be the linear preditor X %*% coef,
and s be log(sigma).  Then the 6 columns of the matrix are L, dL/dp,
ddL/(dp dp), dL/ds, ddL/(ds ds) and ddL/(dp ds), where d stands for the
derivative and dd the second derivative.  Diagnstics based on these quantities
are dicussed in an article by Escobar and Meeker.
.SH REFERENCES
Escobar and Meeker (1992). Assessing influence in regression analysis with
censored data. \fIBiometrics,\fP 48, 507-528.
.EX
fit <- survreg(Surv(time,status) ~x, aml)
rr  <- residuals(fit, type='matrix')
.KW survival
.WR
