# First check, use the data from Turnbull, JASA 1974, 169-173.

tdata <- data.frame(time  =c(1,1,1,2,2,2,3,3,3,4,4,4),
		    status=rep(c(1,0,2),4),
		    n     =c(12,3,2,6,2,4,2,0,2,3,3,5))

tfit <- survfit(Surv(time, time, status, type='interval') ~1, tdata, weight=n)
all.equal(round(tfit$surv,3), c(.538, .295, .210, .095))


# Second check, compare to a reversed survival curve
# This is not as simple a test as one might think, because left and right
#  censored observations are not treated symmetrically by the routine:
#  time <= y for left and time> y for right (this is to make the routine
#  correct for the common situation of panel data).
# To get equivalence, make the left censoreds happen just a little bit
#  earlier.  The left-continuous/right-continuous shift is also a bother.
#
fit1 <- survfit(Surv(time, status) ~1, test1)
temp <- ifelse(test1$status==0, 4.99,5) - test1$time
fit2 <- survfit(Surv(temp, status, type='left') ~1, test1)

all.equal(round(fit1$surv[1:2],5), round(1-fit2$surv[3:2],5))

rm(tdata, tfit, fit1, temp, fit2)
