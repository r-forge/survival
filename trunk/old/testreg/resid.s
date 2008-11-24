fit1 <- survreg(Surv(futime, fustat) ~ age + ecog.ps, ovarian)
fit2 <- censorReg(censor(futime, fustat) ~ age + ecog.ps, ovarian)
fit3 <- survreg(Surv(futime, fustat) ~ age + ecog.ps, ovarian,
		iter=0, init=c(fit2$coef,   log(fit2$scale)))
fit4 <- survreg(Surv(log(futime), fustat) ~age + ecog.ps, ovarian,
		dist='extreme')

print(fit1)
summary(fit4)		


# Hypothesis (and I'm fairly sure): censorReg shares the fault of many
#  iterative codes -- it returns the loglik and variance for iteration k
#  but the coef vector of iteration k+1.  Hence the "all.equal" tests
#  below don't come out perfect.
#

aeq(resid(fit2, type='working')[,1], resid(fit3, type='working'))
aeq(resid(fit2, type='response')[,1], resid(fit3, type='response'))

temp <- sign(resid(fit3, type='working'))
aeq(resid(fit2, type='deviance')[,1], temp*abs(resid(fit3, type='deviance')))
aeq(resid(fit2, type='deviance')[,1], resid(fit3, type='deviance'))

#
# Now check fit1 and fit4, which should follow identical iteration paths
#   These tests should all be true
#
aeq(fit1$coef, fit3$coef)
aeq(fit1$coef, fit4$coef)
 
resid(fit1, type='working')
resid(fit1, type='response')
resid(fit1, type='deviance')
resid(fit1, type='dfbeta')
resid(fit1, type='dfbetas')
resid(fit1, type='ldcase')
resid(fit1, type='ldresp')
resid(fit1, type='ldshape')
resid(fit1, type='matrix')

aeq(resid(fit1, type='working'),resid(fit4, type='working'))
#aeq(resid(fit1, type='response'), resid(fit4, type='response'))#should differ
aeq(resid(fit1, type='deviance'), resid(fit4, type='deviance'))
aeq(resid(fit1, type='dfbeta'),   resid(fit4, type='dfbeta'))
aeq(resid(fit1, type='dfbetas'),  resid(fit4, type='dfbetas'))
aeq(resid(fit1, type='ldcase'),   resid(fit4, type='ldcase'))
aeq(resid(fit1, type='ldresp'),   resid(fit4, type='ldresp'))
aeq(resid(fit1, type='ldshape'),  resid(fit4, type='ldshape'))
aeq(resid(fit1, type='matrix'),   resid(fit4, type='matrix'))
