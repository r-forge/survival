#SCCS $Id: print.survfit.s,v 4.12 1999-02-06 23:27:10 therneau Exp $
print.survfit <- function(x, scale=1, 
			  digits = max(options()$digits - 4, 3), ...) {

    if (!is.null(cl<- x$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
	}
    omit <- x$na.action
    if (length(omit)) cat("  ", naprint(omit), "\n")

    savedig <- options(digits=digits)
    on.exit(options(savedig))
    pfun <- function(stime, surv, n.risk, n.event, lower, upper) {
	#compute the mean, median, se(mean), and ci(median)
	minmin <- function(y, xx) {
	     if (any(!is.na(y) & y==.5)) {
	       if (any(!is.na(y) & y <.5))
		 .5*( min(xx[!is.na(y) & y==.5]) + min(xx[!is.na(y) & y<.5]))
	       else
		 .5*( min(xx[!is.na(y) & y==.5]) + max(xx[!is.na(y) & y==.5]))
	       }
	     else  min(xx[!is.na(y) & y<=.5])
	     }
	n <- length(stime)
	hh <- c(n.event[-n]/(n.risk[-n]*(n.risk[-n]-n.event[-n])), 0)
	nused <- n.risk[1]
	ndead<- sum(n.event)
	dif.time <- c(diff(c(0, stime)), 0)
	if (is.matrix(surv)) {
	    n <- nrow(surv)
	    mean <- dif.time * rbind(1, surv)
	    if (n==1) temp <- mean[2,,drop=F]
	    else temp <- (apply(mean[(n+1):2,,drop=F], 2, cumsum))[n:1,,drop=F]
	    varmean <- c(hh %*% temp^2)
	    med <- apply(surv, 2, minmin, stime)
	    if (!is.null(upper)) {
		upper <- apply(upper, 2, minmin, stime)
		lower <- apply(lower, 2, minmin, stime)
		cbind(nused, ndead, apply(mean, 2, sum),
			  sqrt(varmean), med, lower, upper)
		}
	    else {
		cbind(nused, ndead, apply(mean, 2, sum),
			   sqrt(varmean), med)
		}
	    }
	else {
	    mean <- dif.time*c(1, surv)
	    varmean <- sum(rev(cumsum(rev(mean))^2)[-1] * hh)
	    med <- minmin(surv, stime)
	    if (!is.null(upper)) {
		upper <- minmin(upper, stime)
		lower <- minmin(lower, stime)
		c(nused, ndead, sum(mean), sqrt(varmean), med, lower, upper)
		}
	    else {
		c(nused, ndead, sum(mean), sqrt(varmean), med)
		}
	    }
	}

    stime <- x$time/scale
    surv <- x$surv
    plab <- c("n", "events", "mean", "se(mean)", "median")
    if (!is.null(x$conf.int))
	plab2<- paste(x$conf.int, c("LCL", "UCL"), sep='')

    #Four cases: strata Y/N  by  ncol(surv)>1 Y/N
    #  Repeat the code, with minor variations, for each one
    if (is.null(x$strata)) {
	x1 <- pfun(stime, surv, x$n.risk, x$n.event, x$lower, x$upper)
	if (is.matrix(x1)) {
	    if (is.null(x$lower)) dimnames(x1) <- list(NULL, plab)
	    else                    dimnames(x1) <- list(NULL, c(plab, plab2))
	    }
	else {
	    if (is.null(x$lower)) names(x1) <- plab
	    else                    names(x1) <- c(plab, plab2)
	    }
	print(x1)
	}
    else {   #strata case
	nstrat <- length(x$strata)
	stemp <- rep(1:nstrat,x$strata)
	x1 <- NULL
	for (i in unique(stemp)) {
	    who <- (stemp==i)
	    if (is.matrix(surv)) {
		temp <- pfun(stime[who], surv[who,,drop=F],
			  x$n.risk[who], x$n.event[who],
			  x$lower[who,,drop=F], x$upper[who,,drop=F])
		x1 <- rbind(x1, temp)
		}
	    else  {
		temp <- pfun(stime[who], surv[who], x$n.risk[who],
			  x$n.event[who], x$lower[who], x$upper[who])
		x1 <- rbind(x1, temp)
		}
	    }
	temp <- names(x$strata)
	if (nrow(x1) > length(temp)) {
	    nrep <- nrow(x1)/length(temp)
	    temp <- rep(temp, rep(nrep, length(temp)))
	    }
	if (is.null(x$lower)) dimnames(x1) <- list(temp, plab)
	else                    dimnames(x1) <- list(temp, c(plab, plab2))
	print(x1)
	}
    invisible(x)
    }
