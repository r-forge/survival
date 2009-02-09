options(na.action=na.exclude, contrasts=c(contr.treatment, contr.poly))  #preserve length of missings
library(survival)

capacitor <- read.table('data.capacitor', row.names=1,
			col.names=c('', 'days', 'event', 'voltage'))

fitig <- survreg(Surv(days, event)~voltage, 
	dist = "gaussian", data = capacitor)
summary(fitig)

fitix <- survreg(Surv(days, event)~voltage, 
	dist = "extreme", data = capacitor)
summary(fitix)

fitil <- survreg(Surv(days, event)~voltage, 
	dist = "logistic", data = capacitor)
summary(fitil)
