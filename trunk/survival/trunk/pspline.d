.BG
.FN pspline
.TL
fit a smoothing spline
.DN
This is a modeling term for the survival functions coxph and survreg.
It fits a smoothing spline, using the p-spline basis
.CS
pspline(x, df=4, theta, nterm=2.5 * df, degree=3, eps=0.1, method, ...)
.RA
.AG x
the variable to be fit.  The function does not apply to factor variables.
.OA
.AG df
the desired degrees of freedom.  
One of the arguments `df' or `theta' must be given, but not both.
If `df=0', then the AIC for the term (loglik -df) is used to choose an
"optimal" degrees of freedom.  If AIC is chosen, then an optional
argument `caic=T' can be used to specify the corrected AIC of
Hurvich et. al.
.AG theta
the tuning parameter in the penalized fit.  
It is a monotone function of the degrees of freedom, with theta=1
corresponding to a linear fit and theta=0 to an unconstrained fit
of nterm degrees of freedom.
.AG nterm
the number of basis functions. 
Traditional smoothing splines use one basis per observation, but several
authors have pointed out that the final results of the fit are 
indistinguishable for any number of bases greater than about 
2-3 times the degrees
of freedom. 
.AG degree
the degree of the spline, with cubic splines as the default
.AG eps
tolerance for convergence.
The default states that if 4 degrees of freedom is requested, then a solution
of 3.9 is close enough.\
.AG method
the method for choosing the tuning parameter theta.
If theta is given, then 'fixed' is assumed.
If the degrees of freedom is given, then 'df' is assumed.
If method='aic' then the degrees of freedom is chosen automatically using
Akaike's information criterion.
.AG ...
other arguments to the control function.
.RT
a matrix of basis functions, with the appropriate attributes to be
recognized as a penalized term by the coxph or survreg functions.
.SH REFERENCES
Eilers, Paul H. and Marx, Brian D. (1996).
Flexible smoothing with B-splines and penalties.
    Statistical Science, 11, 89-121.
.pp
Hurvich, C.M. and Simonoff, J.S. and Tsai, Chih-Ling (1998).
Smoothing parameter selection in nonparametric regression using
        an improved Akaike information criterion,
JRSSB, volume 60, 271--293.
 .SA
coxph, survreg, frailty, ridge
.EX
# No evidence for a non-linear age effect
> fit <- coxph(Surv(time, status) ~ pspline(age) + sex + strata(inst), lung)
> fit
                        coef se(coef)     se2 Chisq   DF      p 
pspline(age), linear  0.0189 0.00976  0.00976 3.73  1.00 0.0530
pspline(age), nonlin                          2.14  3.06 0.5500
                 sex -0.4995 0.18304  0.18247 7.45  1.00 0.0064
 
Iterations: 4 outer, 10 Newton-Raphson
     Theta= 0.771 
Degrees of freedom for terms= 4.1 1.0 
Likelihood ratio test=15.6  on 5.05 df, p=0.00852
  n=227 (1 observations deleted due to missing)
.KW survival
.WR
