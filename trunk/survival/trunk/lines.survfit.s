# SCCS $Id: lines.survfit.s,v 4.5 1994-01-06 12:21:39 therneau Exp $
lines.survfit <- function(x, type='s', mark=3, col=1, lty=1, lwd=1,
		       mark.time =T, ...) {
    if (inherits(x, 'survexp')) {
	if (missing(type)) type <- 'l'
	if (!is.numeric(mark.time)) mark.time <- F
	}
    if (is.numeric(mark.time)) mark.time <- sort(unique(mark.time[mark.time>0]))

    if (is.null(x$strata)) {
	if (is.matrix(x$surv)) ncurv <- ncol(x$surv)
	else ncurve <- 1
	nstrat <- 1
	strata <- length(x$surv)/nstrat
	}
    else {
	strata <- x$strata
	nstrat <- length(strata)
	ncurve <- nstrat * length(x$surv)/ sum(strata)
	}

    mark <- rep(mark, length=ncurve)
    col  <- rep(col , length=ncurve)
    lty  <- rep(lty , length=ncurve)
    lwd  <- rep(lwd , length=ncurve)
    time <- rep(x$time, length=length(x$surv))
    j <- 1
    for (i in 1:ncurve) {
	n <- strata[1+(j-1)%%nstrat]
	who <- seq(from=j, length=n)
	xx <- c(0, time[who])
	yy <- c(1, x$surv[who])
	lines(xx, yy, type=type, col=col[i], lty=lty[i], lwd=lwd[i], ...)

	if (is.numeric(mark.time)) {
	    nn <- length(xx)
	    indx <- mark.time
	    for (k in seq(along=mark.time))
		indx[k] <- sum(mark.time[k] > xx)
	    points(mark.time[indx<nn], yy[indx[indx<nn]],
		   pch=mark[i],col=col[i], ...)
	    }
	else if (mark.time==T) {
	    deaths <- c(-1, surv$n.event[who])
	    if ( any(deaths==0))
		points(xx[deaths==0], yy[deaths==0],
			      pch=mark[i],col=col[i], ...)
	    }
	}
    }
