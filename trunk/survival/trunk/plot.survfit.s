#SCCS $Date: 1997-08-04 09:42:10 $ $Id: plot.survfit.s,v 4.13 1997-08-04 09:42:10 atkinson Exp $
plot.survfit<- function(surv, conf.int,  mark.time=T,
		 mark=3,col=1,lty=1, lwd=1, cex=1,log.=F, yscale=1,
		 xscale=1,  ptype=c('surv', 'event', 'cumhaz'),
		 xlab="", ylab="", xaxs='i', ...) {

    if (!inherits(surv, 'survfit'))
	  stop("First arg must be the result of survfit")

    ptype <- match.arg(ptype)
    if (missing(conf.int)) {
	if (is.null(surv$strata) && !is.matrix(surv$surv)) conf.int <-T
	else conf.int <- F
	}

    if (log. & ptype!='surv') 
	    stop("The log option is valid only for ptype='surv'")
    if (ptype=='event') {
	firsty <- 0
	ssurv <- 1-surv$surv
	ymax <- 1
	if (!is.null(surv$upper)) {
	    surv$upper <- 1-surv$upper
	    surv$lower <- 1-surv$lower
	    }
        }
    else if (ptype=='cumhaz') {
	firsty <- 0
	ssurv <- -log(surv$surv)
	ymax <- max(ssurv[!is.inf(ssurv)])
	if (!is.null(surv$upper)) {
	    surv$upper <- -log(surv$upper)
	    surv$lower <- -log(surv$lower)
	    if (conf.int) ymax <- max(surv$lower[!is.inf(surv$lower)], ymax,
				      na.rm=T)
	    }
        }
    else {
	firsty <- 1
	ssurv <- surv$surv
	ymax <- 1
        }
    stime <- surv$time / xscale

    if (is.null(surv$strata)) {
	nstrat <- 1
	stemp <- rep(1, length(surv$time))
	}
    else {
	nstrat <- length(surv$strata)
	stemp <- rep(1:nstrat,surv$strata)
	}
    if (is.null(surv$n.event)) mark.time <- F   #expected survival curve

    # set default values for missing parameters
    if (is.matrix(ssurv)) ncurve <- nstrat * ncol(ssurv)
    else                  ncurve <- nstrat
    mark <- rep(mark, length=ncurve)
    col  <- rep(col, length=ncurve)
    lty  <- rep(lty, length=ncurve)
    lwd  <- rep(lwd, length=ncurve)

    if (is.numeric(mark.time)) mark.time <- sort(mark.time[mark.time>0])
    if (missing(xaxs)) temp <- 1.04*max(stime)
    else               temp <- max(stime)
    #
    # for log plots we have to be tricky about the y axis scaling
    #
    if   (log.) {
	    ymin <- min(.1,ssurv[!is.na(ssurv) &ssurv>0])
	    ssurv[!is.na(ssurv) &ssurv==0] <- ymin
	    plot(c(0, temp),
	       yscale*c(.99, ymin),
	       type ='n', log='y', xlab=xlab, ylab=ylab, xaxs=xaxs,...)
	    }
     else
	 plot(c(0, temp), yscale*c(0,ymax),
	      type='n', xlab=xlab, ylab=ylab, xaxs=xaxs, ...)

    if (yscale !=1) par(usr=par("usr")/ c(1,1,yscale, yscale))
    #
    # put up the curves one by one
    #   survfit has already put them into the "right" order
    i _ 0
    xend _ NULL
    yend _ NULL

    #
    # The 'whom' addition is to replace verbose horizonal sequences
    #  like (1, .2), (1.4, .2), (1.8, .2), (2.3, .2), (2.9, .2), (3, .1)
    #  with (1, .2), (3, .1) -- remember that type='s'.  They are slow,
    #  and can smear the looks of a line type
    #
    for (j in unique(stemp)) {
	who _ (stemp==j)
	xx _ c(0,stime[who])
	deaths <- c(-1, surv$n.event[who])
	if (is.matrix(ssurv)) {
	    for (k in 1:ncol(ssurv)) {
		i _ i+1
		yy _ c(firsty, ssurv[who,k])
		nn <- length(xx)
		whom <- c(match(unique(yy[-nn]), yy), nn)
		lines(xx[whom], yy[whom], lty=lty[i], col=col[i], lwd=lwd[i], type='s')

		if (is.numeric(mark.time)) {
		    indx <- mark.time
		    for (k in seq(along=mark.time))
			indx[k] <- sum(mark.time[k] > xx)
		    points(mark.time[indx<nn], yy[indx[indx<nn]],
			   pch=mark[i],col=col[i],cex=cex)
		    }
		else if (mark.time==T && any(deaths==0))
		    points(xx[deaths==0], yy[deaths==0],
			   pch=mark[i],col=col[i],cex=cex)
		xend _ c(xend,max(xx))
		yend _ c(yend,min(yy))

		if (conf.int==T && !is.null(surv$upper)) {
		    if (ncurve==1) lty[i] <- lty[i] +1
		    yy _ c(firsty, surv$upper[who,k])
		    lines(xx,yy, lty=lty[i], col=col[i], lwd=lwd[i], type='s')
		    yy _ c(firsty, surv$lower[who,k])
		    lines(xx,yy, lty=lty[i], col=col[i], lwd=lwd[i], type='s')
		    }
		}
	    }

	else {
	    i <- i+1
	    yy _ c(firsty, ssurv[who])
	    nn <- length(xx)
	    whom <- c(match(unique(yy[-nn]), yy), nn)
	    lines(xx[whom], yy[whom], lty=lty[i], col=col[i], lwd=lwd[i], type='s')

	    if (is.numeric(mark.time)) {
		indx <- mark.time
		for (k in seq(along=mark.time))
		    indx[k] <- sum(mark.time[k] > xx)
		points(mark.time[indx<nn], yy[indx[indx<nn]],
		       pch=mark[i],col=col[i],cex=cex)
		}
	    else if (mark.time==T && any(deaths==0))
		points(xx[deaths==0], yy[deaths==0],
		       pch=mark[i],col=col[i],cex=cex)
	    xend _ c(xend,max(xx))
	    yend _ c(yend,min(yy))

	    if (conf.int==T && !is.null(surv$upper)) {
		if (ncurve==1) lty[i] <- lty[i] +1
		yy _ c(firsty, surv$upper[who])
		lines(xx,yy, lty=lty[i], col=col[i], lwd=lwd[i], type='s')
		yy _ c(firsty, surv$lower[who])
		lines(xx,yy, lty=lty[i], col=col[i], lwd=lwd[i], type='s')
		}
	    }
	}
    invisible(list(x=xend, y=yend))
}

