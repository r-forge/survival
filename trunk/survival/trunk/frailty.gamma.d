.BG
.FN frailty.gamma
.TL
Random effect term for a survival model
.DN
Fit a random effect term for a variable, where the random effect comes
from the gamma distribution.
This function is specific to the coxph and survreg models.
.CS
frailty.gamma(x, sparse=(ngroup>5), theta, df, eps=1e-05, 
method=c("em", "aic", "df", "fixed"), ...)
.RA
.AG x
variable describing the groups.  A separate indicator variable is fit
for each group.
.OA
.AG sparse
use a sparse matrix method of solution.
This parameter saves considerable time and memory when the number of
groups is large, e.g., a problem with one group per family and 200 families.
.AG theta
variance of the random effect.
.AG df
degrees of freedom for the random effect.  
.AG eps
convergence criteria for the outer loop of the algorithm.
.AG method
method of choosing the variance.
The fixed method requirese that theta be given explicitly (and is assumed
when theta is specified), the df method is used when degrees of freedom
is explicitly set.
The em method computes the maximum likelihood estimate, and achieves
the same solution as the EM algorithm found in several references.
This is an iterative method that searches for the maximum of
a profile likelihood.
The aic method chooses the variance (and hence the degrees of freedom) based
on Akaike's information criteria.
.AG ...
optional arguments for the control function, of which trace=T is the most
usual.  
(The set of available arguments depends on the specific control function).
For the AIC method the "caic=T" argument may be used to choose the corrected
AIC criteria.
.RT
an object with class coxph.penal.
.SE
If sparse computation is chosen, then the coefficients for the sparse terms
are not printed by the default print and summary methods.
.SH REFERENCES
Therneau technical report
.SA
coxph, survreg, frailty
.EX
> fit1 <- coxph(Surv(time, status) ~ age + sex + frailty(inst, df=4), 
	           data=lung)
> fit1
Call:
coxph(formula = Surv(time, status) ~ age + sex + frailty(inst, df = 4), data = 
        lung)
 
                         coef se(coef)     se2 Chisq   DF      p 
                  age  0.0176 0.00933  0.00926 3.55  1.00 0.0600
                  sex -0.5133 0.16880  0.16792 9.25  1.00 0.0024
frailty(inst, df = 4)                          3.26  3.99 0.5100
 
Iterations: 3 outer, 8 Newton-Raphson
     Variance of random effect= 0.0383   EM likelihood = -738.7 
Degrees of freedom for terms= 1 1 4 
Likelihood ratio test=19.6  on 5.96 df, p=0.00314
  n=227 (1 observations deleted due to missing)

> round(fit1$frail,2)
 [1]  0.07  0.08 -0.05 -0.06  0.00  0.08 -0.01  0.06 -0.06  0.02 -0.10 -0.05
[13] -0.03  0.16 -0.13 -0.05 -0.02  0.02
.KW survival
.WR
