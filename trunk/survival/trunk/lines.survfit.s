# SCCS $Id: lines.survfit.s,v 4.10 1997-04-23 16:29:03 therneau Exp $
lines.survfit <- function(x, type='s', mark=3, col=1, lty=1, lwd=1,
		       mark.time =T, xscale=1, yscale=1, 
		        ptype=c('surv', 'event', 'cumhaz'),  ...) {

    ptype <- match.arg(ptype)
    if (inherits(x, 'survexp')) {
	if (missing(type)) type <- 'l'
	if (!is.numeric(mark.time)) mark.time <- F
	}
   
    if (is.numeric(mark.time)) mark.time<- sort(unique(mark.time[mark.time>0]))

    if (is.matrix(x$surv)) {
	ncol.per.strat <- ncol(x$surv)
	ncurve <- ncol(x$surv)
	coffset <- nrow(x$surv)*(1:ncurve -1)     #within matrix offset
	strata <- nrow(x$surv)                    #length of a curve
        }
    else {
	ncol.per.strat <- 1
	ncurve <- 1
	coffset <- 0
	strata <- length(x$time)
        }

    if (is.null(x$strata)) {
	nstrat <- 1
	soffset <- 0
	}
    else {
	strata <- x$strata			#actual length of curves
	nstrat <- length(strata)
	ncurve <- ncurve * nstrat
	soffset<- ncol.per.strat * c(0, cumsum(strata))
	}

    mark <- rep(mark, length=ncurve)
    col  <- rep(col , length=ncurve)
    lty  <- rep(lty , length=ncurve)
    lwd  <- rep(lwd , length=ncurve)
    time <- rep(x$time, ncol.per.strat)

    for (i in 1:nstrat) {
      for (j in 1:ncol.per.strat) {
	who <- seq(soffset[i]+ coffset[j]+1, length=strata[i])  
	xx <- c(0, time[who])/xscale
	nn <- length(xx)

	if (ptype=='event')   yy <- c(0, 1-x$surv[who]) * yscale
	else if(ptype=='surv') yy <- c(1, x$surv[who])*yscale
	else                   yy <- c(0, -log(x$surv[who])) *yscale     
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
      }
    invisible()
    }
