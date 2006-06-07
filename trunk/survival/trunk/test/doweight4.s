#
# Effect of case weights on expected survival curves post Cox model
#
fit0  <- coxph(Surv(time, status) ~x, testw1, weights=wt, method='breslow',
	       iter=0)
fit0b <- coxph(Surv(time, status) ~x, testw2, method='breslow', iter=0)

surv1 <- survfit(fit0, newdata=list(x=0))
surv2 <- survfit(fit0b, newdata=list(x=0))
aeq(surv1$surv, surv2$surv)
