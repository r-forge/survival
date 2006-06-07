#  A second test of survfit.ci
# compare results to Bob Gray's functions
#   This test is in the directory, but won't work unless you have
#   installed the cmprsk library!
# An improved file would check for it's existence first
library(cmprsk)
tstat <- ifelse(mgus1$status==0, 0, pgstat+1)
gray1 <- cuminc(mgus1$time, tstat)

plot(gray1[[1]]$time, gray1[[1]]$est, type='l',
     ylim=range(c(gray1[[1]]$est, gray1[[2]]$est)),
     xlab="Time")
lines(gray1[[2]]$time, gray1[[2]]$est, lty=2)

fit2 <- survfit.ci(Surv(time, status) ~ strata(pgstat), mgus1)
matlines(fit2$time, 1-fit2$surv, col=2, lty=1:2, type='s')

# To formally match these is a bit of a nuisance.
#  The cuminc function returns full step function, and survfit.ci only
# the bottoms of the steps.
#  The survfit.ci function returns all time points, cuminc only the jumps.
temp1 <- tapply(gray1[[1]]$est, gray1[[1]]$time, max)
indx1 <- match(names(temp1), fit2$time)
aeq(temp1, 1-fit2$surv[indx1,1])
