#
# Do a contour plot of the donnell data
#
npt <- 20
beta0  <- seq(.4, 3.4, length=npt)
logsig <- seq(-1.4, 0.41, length=npt)
donlog <- matrix(0,npt, npt)

for (i in 1:npt) {
    for (j in 1:npt) {
	fit <- survreg(Surv(time1, time2, status, type='interval') ~1,
			donnell, init=c(beta0[i],logsig[j]),
		        maxiter=0)
	donlog[i,j] <- fit$log[1]
	}
    }

clev <- -c(51, 51.5, 52:60, 65, 75, 85, 100, 150)
#clev <-  seq(-51, -50, length=10)

contour(beta0, logsig, pmax(donlog, -200), levels=clev, xlab="Intercept",
	ylab="Log(sigma)")
points(2.39, log(.7885), pch=1, col=2)
title("Donnell data")
#
# Compute the path of the iteration
#
niter <- 14
donpath <- matrix(0,niter+1,2)
for (i in 0:niter){
    fit <- survreg(Surv(time1, time2, status, type='interval') ~1,
		    donnell, maxiter=i)
    donpath[i+1,] <- c(fit$coef, log(fit$scale))
    }
points(donpath[,1], donpath[,2])
lines(donpath[,1], donpath[,2], col=4)
