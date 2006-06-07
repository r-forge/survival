#
# An example to break the maximizer, suggested by Bill Meeker
#

hard <- data.frame(time=c(5,25,88,120, 245), status=c(1,1,1,1,0),
		   wt=c(1,1,1,1,500000))

hfit <- survreg(Surv(time, status) ~1, data=hard, weight=wt,
	control=survreg.control(iter=50))

hard2 <- data.frame(time=c(5,25,88,120, 245), status=c(1,1,1,1,0),
		   wt=c(1,1,1,1,50000))
hfit2 <- survreg(Surv(time, status) ~1, data=hard2, weight=wt,
	control=survreg.control(iter=50))
