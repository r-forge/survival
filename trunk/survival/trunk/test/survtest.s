#
# Simple test of (start, stop] Kaplan-Meier curves, using the test2 data
#   set
#
fit1 <- survfit(Surv(start, stop, event) ~1, test2, type='fh2',
                error='tsiatis')
fit2 <- survfit(Surv(start, stop, event) ~x, test2, start.time=3,
		type='fh2')

cfit1<- survfit(coxph(Surv(start, stop, event)~1, test2))
cfit2<- survfit(coxph(Surv(start, stop, event) ~ strata(x), test2, subset=-1))

deaths <- fit1$n.event>0
aeq(fit1$time[deaths], cfit1$time)
aeq(fit1$n.risk[deaths], cfit1$n.risk)
aeq(fit1$n.event[deaths], cfit1$n.event)
aeq(fit1$surv[deaths], cfit1$surv)
aeq(fit1$std.err[deaths], cfit1$std.err)

deaths <- fit2$n.event>0
aeq(fit2$time[deaths], cfit2$time)
aeq(fit2$n.risk[deaths], cfit2$n.risk)
aeq(fit2$n.event[deaths], cfit2$n.event)
aeq(fit2$surv[deaths], cfit2$surv)

fit3 <- survfit(Surv(start, stop, event) ~1, test2) #Kaplan-Meier
aeq(fit3$n, 10)
aeq(fit3$time, c(1:9,14,17))
aeq(fit3$n.risk, c(0,2,3,3,4,5,4,4,5,2,1))
aeq(fit3$n.event,c(0,1,1,0,0,1,1,1,2,0,0))
aeq(fit3$surv[fit3$n.event>0], c(.5, 1/3, 4/15, 1/5, 3/20, 9/100))
