#SCCS $Date: 1992-03-30 02:46:47 $ $Id: print.survfit.s,v 4.2 1992-03-30 02:46:47 therneau Exp $
print.surv.fit <- function(fit.list, times, censored=F,
		       print.it=T,  digits=3, ... ) {
    fit <- fit.list
    if (!inherits(fit.list, 'surv.fit'))
	    stop("Invalid data")

    if (!is.null(cl<- fit$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
	}

    savedig <- options(digits=digits)
    on.exit(options(savedig))

    n <- length(fit$surv)
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
	    times <- fit$time
	    n.risk<- fit$n.risk
	    n.event <- fit$n.event
	    }
	else {
	    who    <- (fit$n.event >0)
	    times  <-  fit$time[who]
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
	if (length(times) >1 )
	    if (any(diff(times)<0)) stop("Times must be in increasing order")

	temp <- .C("survindex", as.integer(n),
				  as.double(fit$time),
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
	fit$n.event[fit$time>max(times)] <- 0
	n.event <- (cumsum(c(0,fit$n.event)))[ifelse(ones, indx, indx+1)]
	n.event<-  diff(c(0, n.event))

	if (!is.null(fit$lower)) {
	    lower <- lower[indx,,drop=F];  lower[ones,] <- 1;
	    upper <- upper[indx,,drop=F];  upper[ones,] <- 1;
	    }

	stemp <- stemp[indx]
	}

    #
    # Now for the printout logic, which is based on width of paper worries
    #  If "print.it=F" return everything (easy)
    #  else if ncurve==1, print all vars
    #          ncurve> 1, skip std error and CI limits
    #
    ncurve <- ncol(surv)
    if (print.it==F) {
	temp <- list(time=times, n.risk=n.risk, n.event=n.event)
	if (ncurve==1) {
	    temp$surv <- drop(temp$surv)
	    if (!is.null(std.err)) temp$std.err <- drop(std.err)
	    if (!is.null(fit$lower)) {
		temp$lower <- drop(lower)
		temp$upper <- drop(upper)
		}
	    }
	else {
	    temp$surv <- temp$surv
	    if (!is.null(std.err)) temp$std.err <- std.err
	    if (!is.null(fit$lower)) {
		temp$lower <- lower
		temp$upper <- upper
		}
	    }
	if (!is.null(fit$strata))
	    temp$strata <- factor(stemp, labels=names(fit$strata))
	return(invisible(temp))
	}

    mat <- cbind(times, n.risk, n.event, surv)
    cnames <- c("time", "n.risk", "n.event")
    if (ncurve==1) {
	cnames <- c(cnames, "survival")
	if (!is.null(std.err)) {
	    if (is.null(fit$lower)) {
		mat <- cbind(mat, std.err)
		cnames <- c(cnames, "std.err")
		}
	    else {
		mat <- cbind(mat, std.err, lower, upper)
		cnames <- c(cnames, 'std.err',
			  paste("lower ", fit$conf.int*100, "% CI", sep=''),
			  paste("upper ", fit$conf.int*100, "% CI", sep=''))
		}
	    }
	}
    else cnames <- c(cnames, paste("survival", seq(ncurve), sep=''))

    dimnames(mat) <- list(NULL, cnames)
    if (nstrat==1) {
	prmatrix(mat, rowlab=rep("", nrow(mat)))
	}
    else  { #print it out one strata at a time
	for (i in unique(stemp)) {
	    who <- (stemp==i)
	    cat("               ", names(fit$strata)[i], "\n")
	    if (sum(who) ==1) print(mat[who,])
	    else    prmatrix(mat[who,], rowlab=rep("", sum(who)))
	    cat("\n")
	    }
	dimnames(mat) <- list(names(fit$strata)[stemp], cnames)
	}
    invisible(mat)
    }
