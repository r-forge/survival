#
# The Stanford data from 1980 is used in Escobar and Meeker
#	t5 = T5 mismatch score
#  Their case numbers correspond to a data set sorted by age
#
stanford2 <- read.table('data.stanford', 
			col.names=c('id', 'time', 'status', 'age', 't5'))
 
stanford2$t5 <- ifelse(stanford2$t5 <0, NA, stanford2$t5)
stanford2 <- stanford2[order(stanford2$age, stanford2$time),]
stanford2$time <- ifelse(stanford2$time==0, .5, stanford2$time)

cage <- stanford2$age - mean(stanford2$age)
fit1 <- survreg(Surv(time, status) ~ cage + cage^2, stanford2,
		dist='lognormal')
fit1
ldcase <- resid(fit1, type='ldcase')
ldresp <- resid(fit1, type='ldresp')
print(ldresp)
# The ldcase and ldresp should be compared to table 1 in Escobar and 
#  Meeker, Biometrics 1992, p519; the colum they label as (1/2) A_{ii}

plot1 <- function() {
    # make their figure 1, 2, and 6
    plot(stanford2$age, stanford2$time, log='y', xlab="Age", ylab="Days",
	 ylim=c(.01, 10^6))
    temp <- predict(fit1, type='response', se.fit=T) 
    matlines(stanford2$age, cbind(temp$fit, temp$fit-1.96*temp$se.fit,
				            temp$fit+1.96*temp$se.fit),
	     lty=c(1,2,2))
    # these are the wrong CI lines, he plotted std dev, I plotted std err
    # here are the right ones
    #  Using uncentered age gives different coefs, but makes prediction over an
    #    extended range somewhat simpler 
    refit <- survreg(Surv(time,status)~ age + age^2, stanford2,
		     dist='lognormal')
    plot(stanford2$age, stanford2$time, log='y', xlab="Age", ylab="Days",
	 ylim=c(.01, 10^6), xlim=c(0,75))
    temp2 <- predict(refit, list(age=1:75), type='quantile', p=c(.05, .5, .95))
    matlines(1:75, temp2, lty=c(1,2,2), col=2)

    tsplot(ldcase, xlab="Case Number", ylab="(1/2) A")
    title (main="Case weight pertubations")
    tsplot(ldresp, xlab="Case Number", ylab="(1/2) A")
    title(main="Response pertubations")
    }

postscript('meekerplot.ps')
plot1()
dev.off()
