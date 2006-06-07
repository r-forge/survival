#
#   Do the test on the simple data set
#
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
fit <-coxph(Surv(time, status) ~x, test1, method='breslow')
fit
fit0 <-coxph(Surv(time, status) ~x, test1, iter=0, method='breslow')
truth0 <- byhand1(0,1)
aeq(truth0$loglik, fit0$loglik[1])
aeq(1/truth0$imat, fit0$var)
aeq(truth0$mart, fit0$resid[c(2:6,1)])
aeq(truth0$scho, resid(fit0, 'schoen'))


fit0 <- coxph(Surv(time, status) ~x, test1, iter=0)
fit0$coef
coxph(Surv(time, status) ~x, test1, iter=1, method='breslow')$coef
coxph(Surv(time, status) ~x, test1, iter=2, method='breslow')$coef
coxph(Surv(time, status) ~x, test1, iter=3, method='breslow')$coef

coxph(Surv(time, status) ~ x, test1, method='efron')
coxph(Surv(time, status) ~ x, test1, method='exact')

resid(fit0)
resid(coxph(Surv(time, status) ~1, test1))
resid(fit0, 'scor')
resid(fit0, 'scho')

resid(fit)
resid(fit, 'scor')
resid(fit, 'scho')

predict(fit, type='lp')
predict(fit, type='risk')
predict(fit, type='expected')
predict(fit, type='terms')
predict(fit, type='lp', se.fit=T)
predict(fit, type='risk', se.fit=T)
predict(fit, type='expected', se.fit=T)
predict(fit, type='terms', se.fit=T)

summary(survfit(fit))
summary(survfit(fit, list(x=2)))
