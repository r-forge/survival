#
# Check out the Efron approx. 
#

fit0 <- coxph(Surv(time, status) ~x,testw1, weights=wt, iter=0)
fit  <- coxph(Surv(time, status) ~x,testw1, weights=wt)
resid(fit0, 'mart')
resid(coxph(Surv(time, status) ~1, testw1, weights=wt))  #Null model

# lfun is the known log-likelihood for this data set, worked out in the
#   appendix of Therneau and Grambsch
# ufun is the score vector and ifun the information matrix
lfun <- function(beta) {
    r <- exp(beta)
    a <- 7*r +3
    b <- 4*r +2
    11*beta - ( log(r^2 + 11*r +7) + 
	(10/3)*(log(a+b) + log(2*a/3 +b) + log(a/3 +b)) + 2*log(2*r +1))
    }
aeq(fit0$log[1], lfun(0))
aeq(fit$log[2], lfun(fit$coef))

ufun <- function(beta, efron=T) {  #score statistic
    r <- exp(beta)
    xbar1 <- (2*r^2+11*r)/(r^2+11*r +7)
    xbar2 <- 11*r/(11*r +5)
    xbar3 <-  2*r/(2*r +1)
    xbar2b<- 26*r/(26*r+12)
    xbar2c<- 19*r/(19*r + 9)
    temp <- 11 - (xbar1 + 2*xbar3)
    if (efron) temp - (10/3)*(xbar2 + xbar2b + xbar2c)
    else       temp - 10*xbar2
    }
print(ufun(fit$coef) < 1e-4)  # Should be true

ifun <- function(beta, efron=T) {  # information matrix
    r <- exp(beta)
    xbar1 <- (2*r^2+11*r)/(r^2+11*r +7)
    xbar2 <- 11*r/(11*r +5)
    xbar3 <-  2*r/(2*r +1)
    xbar2b<- 26*r/(26*r+12)
    xbar2c<- 19*r/(19*r + 9)
    temp <- ((4*r^2 + 11*r)/(r^2+11*r +7) - xbar1^2) +
	    2*(xbar3 - xbar3^2)
    if (efron) temp + (10/3)*((xbar2- xbar2^2) + (xbar2b - xbar2b^2) +
			      (xbar2c -xbar2c^2))
    else       temp + 10 * (xbar2- xbar2^2)
    }

aeq(fit0$var, 1/ifun(0))
aeq(fit$var, 1/ifun(fit$coef))


      
# Make sure that the weights pass through the residuals correctly
rr1 <- resid(fit, type='mart')
rr2 <- resid(fit, type='mart', weighted=T)
aeq(rr2/rr1, testw1$wt)
rr1 <- resid(fit, type='score')
rr2 <- resid(fit, type='score', weighted=T)
aeq(rr2/rr1, testw1$wt)

#
# Look at the individual components
#
dt0 <- coxph.detail(fit0)
dt <- coxph.detail(fit)
aeq(sum(dt$score), ufun(fit$coef))  #score statistic
aeq(sum(dt0$score), ufun(0))
aeq(dt0$hazard, c(1/19, (10/3)*(1/16 + 1/(6+20/3) + 1/(6+10/3)), 2/3))



rm(fit, fit0, rr1, rr2, dt, dt0)
