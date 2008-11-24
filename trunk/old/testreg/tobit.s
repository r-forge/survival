#
# Another example from the SAS manual, claiming to be Tobit regression
#
#  Data is from J. Tobin, Estimation of relationships for limited dependent
#   variables, Econometrica, v26, 24-36, 1958.
#
#  Durable Goods Purchase, Age in Years, Liquidity Ratio * 1000
#
#  The data set is copied from the SAS manual.  The answers below are not
#    exactly the ones printed in the SAS manual, but they are the ones that
#    you get using SAS and this data.  (Typo in the data set?)

tobin <- read.table('data.tobin', col.names=c('durable', 'age', 'quant'))

tfit <- survreg(Surv(durable, durable>0, type='left') ~age + quant,
		data=tobin, dist='gaussian')

#
# Repeat as interval censored data
#
temp <- Surv(ifelse(tobin$durable==0, NA, tobin$durable),
	     tobin$durable, type='interval2')
tfit2 <- survreg(temp ~ age + quant, data=tobin, dist='gaussian')

# Make uncensored interval data that mimics the above
temp1 <- ifelse(tobin$durable==0, -300, tobin$durable)
tfit3 <- survreg(Surv(temp1, durable, type='interval2') ~age + quant,
		 dist='gaussian', data=tobin)


tgauss <- survreg.distributions[['gaussian']]
tgauss$name <- 'gauss2'
tfit4 <-survreg(Surv(durable, durable>0, type='left') ~age + quant,
		data=tobin, dist=tgauss) 

# test equality
# tfit3 may differ somewhat, due to a different starting estimate
all.equal(tfit$coef, tfit2$coef)
all.equal(tfit$coef, tfit4$coef)
all.equal(tfit$coef, tfit3$coef)  

