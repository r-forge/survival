#SCCS $Date: 1992-03-04 16:48:11 $ $Id: plot.survfit.s,v 4.1 1992-03-04 16:48:11 therneau Exp $
plot.surv.fit<- function(surv, conf.int,  mark.time=T,
		 mark=3,col=1,lty=1,na.strata=T,cex=1,log=F,yscale=1,
		 xlab="", ylab="", xaxs='i', ...) {

    if (!inherits(surv, 'surv.fit'))
	  stop("First arg must be the result of surv.fit")

    if (missing(conf.int)) {
	if (is.null(surv$strata)) conf.int <-T
	else conf.int <- F
	}

    if (is.null(surv$strata)) {
	nstrat <- 1
	stemp <- rep(1, length(surv$time))
	}
    else {
	nstrat <- length(surv$strata)
	stemp <- rep(1:nstrat,surv$strata)
	}

    # set default values for missing parameters
    mark <- rep(mark, length=nstrat)
    col  <- rep(col, length=nstrat)
    lty  <- rep(lty, length=nstrat)

    if (is.numeric(mark.time)) mark.time <- sort(mark.time[mark.time>0])
    #
    # for log plots we have to be tricky about the y axis scaling
    #
    if   (log) {
	    plot(c(0, max(surv$time)),
	       yscale*c(.99,min(.1,surv$surv[surv$surv>0],na.rm=T)),
	       type ='n', log='y', xlab=xlab, ylab=ylab, xaxs=xaxs,...)
	    }
     else
	 plot(c(0,max(surv$time)), yscale*c(0,1),
	      type='n', xlab=xlab, ylab=ylab, xaxs=xaxs, ...)

    if (yscale !=1) par(usr=par("usr")/ c(1,1,yscale, yscale))
    #
    # put up the curves one by one
    #   surv.fit has already put them into the "right" order
    i _ 0
    xend _ NULL
    yend _ NULL

    for (j in unique(stemp)) {
	i _ i+1
	who _ (stemp==j)
	# next line identifies all of the 'step downs' or 'last point'
	drops _ (surv$n.event>0 | surv$time==max(surv$time[who]))
	xx _ c(0,surv$time[who & drops])
	yy _ c(1,surv$surv[who &drops])
	lines(stepfun(xx,yy), lty=lty[i], col=col[i])

	if (is.numeric(mark.time)) {
	    nn <- length(xx)
	    indx <- mark.time
	    for (k in seq(along=mark.time))
		indx[k] <- sum(mark.time[k] > xx)
	    points(mark.time[indx<nn], yy[indx[indx<nn]],
		   pch=mark[i],col=col[i],cex=cex)
	    }
	else if (mark.time==T & any(surv$n.event[who]==0))
	    points(surv$time[who & surv$n.event==0],
		   surv$surv[who & surv$n.event==0],
		   pch=mark[i],col=col[i],cex=cex)
	xend _ c(xend,max(xx))
	yend _ c(yend,min(yy))

	if (conf.int==T && !is.null(surv$upper)) {
	    if (nstrat==1) lty[i] <- lty[i] +1
	    yy _ c(1,surv$upper[who &drops])
	    lines(stepfun(xx,yy), lty=lty[i], col=col[i])
	    yy _ c(1,surv$lower[who &drops])
	    lines(stepfun(xx,yy), lty=lty[i], col=col[i])
	    }

	}
    invisible(list(x=xend, y=yend))
}

