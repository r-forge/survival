# SCCS $Id: lines.survfit.s,v 4.8 1994-12-14 14:51:24 therneau Exp $
lines.survfit <- function(x, type='s', mark=3, col=1, lty=1, lwd=1,
		       mark.time =T, xscale=1, yscale=1,  ...) {
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
	n <- strata[1+(i-1)%%nstrat]
	who <- seq(from=j, length=n)
	j <-  j+n
	xx <- c(0, time[who])/xscale
	yy <- c(1, x$surv[who])*yscale
	nn <- length(xx)
	#
	# This 'whom' addition is to replace verbose horizonal sequences
	#  like (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
	#  with (1, .2), (3, .1) if type='s' and (1, .2), (2.9, .2), (3, .1)
	#  otherwise.  They are slow, and can smear the looks of a line type
	#
	whom <- c(match(unique(yy[-nn]), yy), nn)
	if (type!='s') whom <- sort(unique(c(whom, whom[-1]-1)))
	lines(xx[whom], yy[whom], type=type, col=col[i], lty=lty[i], lwd=lwd[i], ...)

	if (is.numeric(mark.time)) {
	    indx <- mark.time
	    for (k in seq(along=mark.time))
		indx[k] <- sum(mark.time[k] > xx)
	    points(mark.time[indx<nn], yy[indx[indx<nn]],
		   pch=mark[i],col=col[i], ...)
	    }
	else if (mark.time==T) {
	    deaths <- c(-1, x$n.event[who])
	    if ( any(deaths==0))
		points(xx[deaths==0], yy[deaths==0],
			      pch=mark[i],col=col[i], ...)
	    }
	}
    invisible()
    }
