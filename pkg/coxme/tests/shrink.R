library(coxme)
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))
#
# Variable shrinkage
#
ecog0 <- 1*(lung$ph.ecog==0)
ecog1 <- 1*(lung$ph.ecog==1)
ecog2 <- 1*(lung$ph.ecog==2)
ecog3 <- 1*(lung$ph.ecog==3)

fit1 <- coxph(Surv(time, status) ~ age + ridge(ecog0, ecog1, ecog2, ecog3,
                                               scale=FALSE, theta=2), lung)

fit2 <- coxme(Surv(time, status) ~ age + (factor(ph.ecog) |1), lung,
              variance=.5)

aeq(fit1$coef, c(coef(fit2)$fixed, unlist(fit2$frail)))
indx <- c(5,1,2,3,4) #in coxme, shrinkage variables are first
all.equal(fit1$var, as.matrix(fit2$var)[indx, indx])

fit3 <- coxme(Surv(time, status) ~ age + (1|ph.ecog), lung, variance=.5)
all.equal(fit2$var, fit3$var)
all.equal(fit2$loglik, fit3$loglik)

fit4 <- coxme(Surv(time, status) ~ age + (1|ph.ecog), lung)


#fit4 <- coxme(Surv(time, status) ~ age + (factor(ph.ecog) |1), lung,
#              varlist=bdsI)

