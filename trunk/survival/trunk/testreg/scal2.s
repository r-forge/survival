#
# debug scale
#
fit1 <- survreg(Surv(time, status) ~ ridge(x, theta=1), test1, debug=1)
fit2 <- survreg(Surv(time, status) ~ ridge(x, theta=1), test1,
		 scale=fit1$scale)

all.equal(fit1$coef, fit2$coef)
all.equal(fit1$log[2], fit2$loglik[2])
