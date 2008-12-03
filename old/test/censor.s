# Fit the Fleming data with censorReg

temp1 <- ifelse(is.na(flem2$ltime), 2,
		ifelse(is.na(flem2$rtime),0, 1+ 2*(flem2$ltime < flem2$rtime)))
temp2 <- ifelse(temp1==2, flem2$rtime, flem2$ltime)

# Should be the same as fitfw2
xfit <- censorReg(censor(temp2, flem2$rtime, temp1, type='interval') ~
		  age + ecog.ps, ovarian, dist='weibull')
