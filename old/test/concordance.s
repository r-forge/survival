#
# Simple tests of concordance.  These numbers were derived in multiple
#   codes.
#
tdata <- leukemia[leukemia$group=='Maintained',]
y <- c(1,6,2,7,3,7,3,8,4,4,5)
fit <- survConcordance(Surv(time, status) ~y, tdata)
aeq(fit$stats, c(14,24,2,0,15))

# Lots of ties
tempx <- Surv(c(1,2,2,2,3,4,4,4,5,2), c(1,0,1,0,1,0,1,1,0,1))
tempy <- c(5,5,4,4,3,3,7,6,5,4)
fit2 <- survConcordance(tempx ~ tempy)
aeq(fit2$stats, c(13,13,5,2,12))

# Bigger data
fit3 <- survConcordance(Surv(time, status) ~ age, lung)
aeq(fit3$stats, c(10717, 8706, 596, 28, 5836))

# More ties
fit4 <- survConcordance(Surv(time, status) ~ ph.ecog, lung)
aeq(fit4$stats, c(8392, 4258, 7137, 28, 5836))
