# Tests of the weighted Cox model
#
# Similar data set to test1, but add weights,
#                                    a double-death/censor tied time
#                                    a censored last subject
# The latter two are cases covered only feebly elsewhere.
# 
# The data set testw2 has the same data, but done via replication
#
aeq <- function(x,y) all.equal(as.vector(x), as.vector(y))

testw1 <- data.frame(time=  c(1,1,2,2,2,2,3,4,5),
		    status= c(1,0,1,1,1,0,0,1,0),
		    x=      c(2,0,1,1,0,1,0,1,0),
		    wt =    c(1,2,3,4,3,2,1,2,1))
xx <- c(1,2,3,4,3,2,1,2,1)
testw2 <- data.frame(time=   rep(c(1,1,2,2,2,2,3,4,5), xx),
		     status= rep(c(1,0,1,1,1,0,0,1,0), xx),
		     x=      rep(c(2,0,1,1,0,1,0,1,0), xx),
		     id=     rep(1:9, xx))
indx <- match(1:9, testw2$id)
testw2 <- data.frame(time=   rep(c(1,1,2,2,2,2,3,4,5), xx),
		     status= rep(c(1,0,1,1,1,0,0,1,0), xx),
		     x=      rep(c(2,0,1,1,0,1,0,1,0), xx),
		     id=     rep(1:9, xx))
indx <- match(1:9, testw2$id)

fit0 <- coxph(Surv(time, status) ~x, testw1, weights=wt,
		    method='breslow', iter=0)
fit0b <- coxph(Surv(time, status) ~x, testw2, method='breslow', iter=0)
fit  <- coxph(Surv(time, status) ~x, testw1, weights=wt, method='breslow')
fitb <- coxph(Surv(time, status) ~x, testw2, method='breslow')

texp <- function(beta) {  # expected, Breslow estimate
    r <- exp(beta)
    temp <- cumsum(c(1/(r^2 + 11*r +7), 10/(11*r +5), 2/(2*r+1)))
    c(r^2, 1,r,r,1,r,1,r,1)* temp[c(1,1,2,2,2,2,2,3,3)]
    }
aeq(texp(0),  c(1/19, 1/19, rep(103/152, 5), rep(613/456,2))) #verify texp()

xbar <- function(beta) { # xbar, Breslow estimate
    r <- exp(beta)
    temp <- r* rep(c(2*r + 11, 11/10, 1), c(2, 5, 2))
    temp * texp(beta)
    }

fit0
summary(fit)
aeq(resid(fit0), testw1$status - texp(0))
resid(fit0, type='score')
resid(fit0, type='scho')

aeq(resid(fit0, type='mart'), (resid(fit0b, type='mart'))[indx])
aeq(resid(fit0, type='scor'), (resid(fit0b, type='scor'))[indx])
aeq(unique(resid(fit0, type='scho')), unique(resid(fit0b, type='scho')))


aeq(resid(fit, type='mart'), testw1$status - texp(fit$coef))
resid(fit, type='score')
resid(fit, type='scho')
aeq(resid(fit, type='mart'), (resid(fitb, type='mart'))[indx])
aeq(resid(fit, type='scor'), (resid(fitb, type='scor'))[indx])
aeq(unique(resid(fit, type='scho')), unique(resid(fitb, type='scho')))
rr1 <- resid(fit, type='mart')
rr2 <- resid(fit, type='mart', weighted=T)
aeq(rr2/rr1, testw1$wt)

rr1 <- resid(fit, type='score')
rr2 <- resid(fit, type='score', weighted=T)
aeq(rr2/rr1, testw1$wt)

fit  <- coxph(Surv(time, status) ~x, testw1, weights=wt, method='efron')
fit
resid(fit, type='mart')
resid(fit, type='score')
resid(fit, type='scho')

