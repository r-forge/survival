#SCCS $Date: 1992-09-20 23:26:34 $ $Id: summary.survfit.s,v 1.5 1992-09-20 23:26:34 therneau Exp $
summary.survfit <- function(fit, times, censored=F, scale=1, ...) {
    if (!inherits(fit, 'survfit'))
	    stop("Invalid data")

    n <- length(fit$surv)
    stime <- fit$time/scale
    if (is.null(fit$strata)) {
	stemp <- rep(1,n)
	nstrat <- 1
	}
    else {
	nstrat <- length(fit$strata)
	stemp <- rep(1:nstrat,fit$strata)
	}

    surv <- as.matrix(fit$surv)
    if (is.null(fit$std.err)) std.err <- NULL
    else                      std.err <- fit$std.err * surv

    if (!is.null(fit$lower)) {
	lower <- as.matrix(fit$lower)
	upper <- as.matrix(fit$upper)
	}

    if (missing(times)) {
	if (censored) {
	    times <- stime
	    n.risk<- fit$n.risk
	    n.event <- fit$n.event
	    }
	else {
	    who    <- (fit$n.event >0)
	    times  <-  stime[who]
	    n.risk <-  fit$n.risk[who]
	    n.event <- fit$n.event[who]
	    stemp <- stemp[who]
	    surv <- surv[who,,drop=F]
	    if (!is.null(std.err)) std.err <- std.err[who,,drop=F]
	    if (!is.null(fit$lower)) {
		lower <- lower[who,,drop=F]
		upper <- upper[who,,drop=F]
		}
	    }
	}

    else {  #this case is much harder
	if (any(times<0)) stop("Invalid time point requested")
        if(max(fit$time) < min(times))
            stop("Requested times are all beyond the end of the survival curve")
	if (length(times) >1 )
	    if (any(diff(times)<0)) stop("Times must be in increasing order")

	temp <- .C("survindex2", as.integer(n),
				  as.double(stime),
				  as.integer(stemp),
				  as.integer(length(times)),
				  as.double(times),
				  as.integer(nstrat),
				  indx = integer(nstrat*length(times)),
				  indx2= integer(nstrat*length(times)) )
	keep <- temp$indx >=0
	indx <- temp$indx[keep]
	ones <- (temp$indx2==1)[keep]
	ties <- (temp$indx2==2)[keep]  #data set time === requested time

	times <- rep(times, nstrat)[keep]
	n.risk <- fit$n.risk[indx+1 - (ties+ones)]
	surv   <- surv[indx,,drop=F];   surv[ones,] <- 1
	if (!is.null(std.err)) {
	    std.err<- std.err[indx,,drop=F]
	    std.err[ones,] <-0
	    }
	fit$n.event[stime>max(times)] <- 0
	n.event <- (cumsum(c(0,fit$n.event)))[ifelse(ones, indx, indx+1)]
	n.event<-  diff(c(0, n.event))

	if (!is.null(fit$lower)) {
	    lower <- lower[indx,,drop=F];  lower[ones,] <- 1;
	    upper <- upper[indx,,drop=F];  upper[ones,] <- 1;
	    }

	stemp <- stemp[indx]
	}

    ncurve <- ncol(surv)
    temp <- list(surv=surv, time=times, n.risk=n.risk, n.event=n.event)

    if (ncurve==1) {
	temp$surv <- drop(temp$surv)
	if (!is.null(std.err)) temp$std.err <- drop(std.err)
	if (!is.null(fit$lower)) {
	    temp$lower <- drop(lower)
	    temp$upper <- drop(upper)
	    }
	}
    else {
	if (!is.null(std.err)) temp$std.err <- std.err
	if (!is.null(fit$lower)) {
	    temp$lower <- lower
	    temp$upper <- upper
	    }
	}
    if (!is.null(fit$strata))
	temp$strata <- factor(stemp,
	    labels = names(fit$strata)[sort(unique(stemp))])
    temp$call <- fit$call
    if (!is.null(fit$na.action)) temp$na.action <- fit$na.action
    class(temp) <- 'summary.survfit'
    temp
    }
