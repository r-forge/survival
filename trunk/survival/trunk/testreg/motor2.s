#
# Reproduce example 1 in the SAS lifereg documentation
#

# this fit doesn't give the same log-lik that they claim
fit1 <- survreg(Surv(time, status) ~ I(1000/(273.2+temp)), motor,
		subset=(temp>150), dist='lognormal')
summary(fit1)

# This one, with the loglik on the transformed scale (the inappropriate
#   scale, Ripley & Venables would argue) does agree.
# All coefs are of course identical.
fit2 <- survreg(Surv(log(time), status) ~ I(1000/(273.2+temp)), motor,
		subset=(temp>150), dist='gaussian')


# Give the quantile estimates

pp1 <- predict(fit1, newdata=list(temp=c(130,150)), p=c(.1, .5, .9),
		      type='quantile', se=T)
pp2 <- predict(fit1, newdata=list(temp=c(130,150)), p=c(.1, .5, .9),
		      type='uquantile', se=T)
pp1

temp130 <- matrix(0, nrow=3, ncol=6)
temp130[,1] <- pp1$fit[1,]
temp130[,2] <- pp1$se.fit[1,]
temp130[,3] <- pp2$fit[1,]
temp130[,4] <- pp2$se.fit[1,]
temp130[,5] <- exp(pp2$fit[1,] - 1.64*pp2$se.fit[1,])
temp130[,6] <- exp(pp2$fit[1,] + 1.64*pp2$se.fit[1,])
dimnames(temp130) <- list(c("p=.1", "p=.2", "p=.3"),
     c("Time", "se(time)", "log(time)", "se[log(time)]", 
       "lower 90", "upper 90"))
print(temp130)
