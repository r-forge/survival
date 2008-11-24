#
#  Verify that both surv AND n.risk are right between time points.
#
fit <- survfit(Surv(time, status) ~1, test1)
temp <- summary(fit, time=c(.5,1, 1.5, 2, 2.5, 3, 3.5, 4, 4.5), extend=T)

aeq(temp$n.risk, c(6,6,4,4,2,2,1,1,0))
aeq(temp$surv, c(1, fit$surv[c(1,1,2,2,3,3,4,4)]))
aeq(temp$n.event, c(0,1,0,2,0,0,0,1,0))
aeq(temp$std.err, c(0, (fit$surv*fit$std.err)[c(1,1,2,2,3,3,4,4)]))


fit <- survfit(Surv(start, stop, event) ~1, test2)
temp <- summary(fit, times=c(.5, 1.5, 2.5, 3, 6.5, 14.5, 16.5))
aeq(temp$surv, c(1, fit$surv[c(1,2,3,6, 10,10)]))
aeq(temp$n.risk, c(0, 2, 3, 3, 4, 1,1))
