#SCCS $Date: 1993-04-12 11:06:51 $ $Id: print.survfit.s,v 4.7 1993-04-12 11:06:51 therneau Exp $
print.survfit <- function(fit, scale=1, digits=3, ...) {

    if (!is.null(cl<- fit$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
	}
    omit <- fit$na.action
    if (length(omit)) cat("  ", naprint(omit), "\n")

    savedig <- options(digits=digits)
    on.exit(options(savedig))
    pfun <- function(stime, surv, n.risk, n.event, lower, upper) {
	#compute the mean, median, se(mean), and ci(median)
	minmin <- function(y, x) {
	     if (any(!is.na(y) & y==.5)) {
	       if (any(!is.na(y) & y <.5))
		 .5*( min(x[!is.na(y) & y==.5]) + min(x[!is.na(y) & y<.5]))
	       else
		 .5*( min(x[!is.na(y) & y==.5]) + max(x[!is.na(y) & y==.5]))
	       }
	     else  min(x[!is.na(y) & y<=.5])
	     }
	n <- length(stime)
	hh <- c(n.event[-n]/(n.risk[-n]*(n.risk[-n]-n.event[-n])), 0)
	nused <- n.risk[1]
	ndead<- sum(n.event)
	dif.time <- c(diff(c(0, stime)), 0)
	if (is.matrix(surv)) {
	    n <- nrow(surv)
	    mean <- dif.time * rbind(1, surv)
	    temp <- (apply(mean[(n+1):2,,drop=F], 2, cumsum))[n:1,,drop=F]
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

    stime <- fit$time/scale
    surv <- fit$surv
    plab <- c("n", "events", "mean", "se(mean)", "median")
    if (!is.null(fit$conf.int))
	plab2<- paste(fit$conf.int, c("CI", "CI"), sep='')

    #Four cases: strata Y/N  by  ncol(surv)>1 Y/N
    #  Repeat the code, with minor variations, for each one
    if (is.null(fit$strata)) {
	x <- pfun(stime, surv, fit$n.risk, fit$n.event, fit$lower, fit$upper)
	if (is.matrix(x)) {
	    if (is.null(fit$lower)) dimnames(x) <- list(NULL, plab)
	    else                    dimnames(x) <- list(NULL, c(plab, plab2))
	    }
	else {
	    if (is.null(fit$lower)) names(x) <- plab
	    else                    names(x) <- c(plab, plab2)
	    }
	print(x)
	}
    else {   #strata case
	nstrat <- length(fit$strata)
	stemp <- rep(1:nstrat,fit$strata)
	x <- NULL
	for (i in unique(stemp)) {
	    who <- (stemp==i)
	    if (is.matrix(surv)) {
		temp <- pfun(stime[who], surv[who,,drop=F],
			  fit$n.risk[who], fit$n.event[who],
			  fit$lower[who,,drop=F], fit$upper[who,,drop=F])
		x <- rbind(x, temp)
		}
	    else  {
		temp <- pfun(stime[who], surv[who], fit$n.risk[who],
			  fit$n.event[who], fit$lower[who], fit$upper[who])
		x <- rbind(x, temp)
		}
	    }
	temp <- names(fit$strata)
	if (nrow(x) > length(temp)) {
	    nrep <- nrow(x)/length(temp)
	    temp <- rep(temp, rep(nrep, length(temp)))
	    }
	if (is.null(fit$lower)) dimnames(x) <- list(temp, plab)
	else                    dimnames(x) <- list(temp, c(plab, plab2))
	print(x)
	}
    invisible(fit)
    }
