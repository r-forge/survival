#
# Check that the test statistic computed within aareg and
#  the one recomputed within summary.aareg are the same.
# Of course, they could both be wrong, but at least they'll agree!
# If the maxtime argument is used in summary, it recomputes the test,
#  even if we know that it wouldn't have had to.
#
# Because the 1-variable and >1 variable case have different code, test
#  them both.
#
afit <- aareg(Surv(time, status) ~ age, lung, dfbeta=T)
asum <- summary(afit, maxtime=max(afit$times))
aeq(afit$test.stat, asum$test.stat)
aeq(afit$test.var,  asum$test.var)
aeq(afit$test.var2, asum$test.var2)

print(afit)

afit <- aareg(Surv(time, status) ~ age, lung, dfbeta=T, test='nrisk')
asum <- summary(afit, maxtime=max(afit$times))
aeq(afit$test.stat, asum$test.stat)
aeq(afit$test.var,  asum$test.var)
aeq(afit$test.var2, asum$test.var2)

summary(afit)

#
# Mulitvariate
#
afit <- aareg(Surv(time, status) ~ age + sex + ph.karno + pat.karno, lung,
	      dfbeta=T)
asum <- summary(afit, maxtime=max(afit$times))
aeq(afit$test.stat, asum$test.stat)
aeq(afit$test.var,  asum$test.var)
aeq(afit$test.var2, asum$test.var2)

print(afit)

afit <- aareg(Surv(time, status) ~ age + sex + ph.karno + pat.karno, lung,
	      dfbeta=T, test='nrisk')
asum <- summary(afit, maxtime=max(afit$times))
aeq(afit$test.stat, asum$test.stat)
aeq(afit$test.var,  asum$test.var)
aeq(afit$test.var2, asum$test.var2)

summary(afit)

# Weights play no role in the final computation of the test statistic, given
#  the coefficient matrix, nrisk, and dfbeta as inputs.  (Weights do
#  change the inputs).  So there is no need to reprise the above with
#  case weights.
