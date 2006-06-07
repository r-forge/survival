#
# Check out the impact of weights on the dfbetas
#    Am I computing them correctly?
#
wtemp <- rep(1,26)
wtemp[c(5,10,15)] <- 2:4
fit <- coxph(Surv(futime, fustat) ~ age + ecog.ps, ovarian, weights=wtemp)
rr <- resid(fit, 'dfbeta')

fit1 <- coxph(Surv(futime, fustat) ~ age + ecog.ps, ovarian, weights=wtemp,
	         subset=(-5))
fit2 <- coxph(Surv(futime, fustat) ~ age + ecog.ps, ovarian, weights=wtemp,
	         subset=(-10))
fit3 <- coxph(Surv(futime, fustat) ~ age + ecog.ps, ovarian, weights=wtemp,
	         subset=(-15))

