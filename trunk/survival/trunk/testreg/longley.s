#
# Do a ridge regression/ ridge trace on the longley data set
#
longley <- read.table("data.longley", 
		      col.names=c("employed", "deflator", "gnp", "unemployed",
			          "military", "pop", "year"))
lfit0 <- lm(employed ~ deflator + gnp + unemployed + military + pop + year,
	    data=longley)

lfit1 <- survreg(Surv(employed) ~ deflator + gnp + unemployed + military + 
		      pop + year, dist='gaussian', data=longley,
		      control=survreg.control(toler.chol=1e-13))

lfit2 <-  survreg(Surv(employed) ~ ridge(deflator, gnp, unemployed, military,
					 pop, year, theta=1e-6),
		    dist='gaussian', data=longley,
		    control=survreg.control(toler.chol=1e-13))

all.equal(lfit0$coef, lfit1$coef)
print(lfit2)


