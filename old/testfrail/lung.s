#
# A test with the lung data
#  This caused problems at one point in the beta release

#
#   First, get rid of some missings
#
lung2 <- na.omit(lung[c('time', 'status', 'wt.loss')])

#
# Test the logliklihoods
#
fit <- coxph(Surv(time, status) ~ pspline(wt.loss,3), lung2, x=T)
fit0<- coxph(Surv(time, status) ~ 1, lung2)
fit1<- coxph(Surv(time, status) ~ fit$x, lung2, iter=0, init=fit$coef)

all.equal(fit$loglik[1], fit0$loglik)
all.equal(fit$loglik[2], fit1$loglik[2])

#
# Check variances
#
imat <- solve(fit1$var)
var2 <- fit$var %*% imat %*% fit$var
all.equal(fit$var2, var2)

rm(lung2, imat, var2, fit, fit0, fit1)
