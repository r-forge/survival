#SCCS $Id: print.survexp.s,v 4.6 1994-01-04 14:56:59 therneau Exp $
print.survexp <- function(fit, scale=1, digits=3, ...) {
    if (!inherits(fit, 'survexp'))
	    stop("Invalid data")
    savedig <- options(digits=digits)
    on.exit(options(savedig))

    if (!is.null(cl<- fit$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
	}

    if (!is.null(fit$summ)) cat(fit$summ)
    omit <- fit$na.action
    if (length(omit))
	cat(naprint(omit), "\n")
    else cat("\n")

    mat <- cbind(fit$time, fit$n.risk, fit$surv)
    cnames <- c("time", "n.risk")
    if (is.matrix(fit$surv)) {
	ncurve <- ncol(fit$surv)
	cnames <- c(cnames, paste("survival", seq(ncurve), sep=''))
	}
    else {
	ncurve <- 1
	cnames <- c(cnames, "survival")
	}

    dimnames(mat) <- list(NULL, cnames)
    if (is.null(fit$strata)) {
	prmatrix(mat, rowlab=rep("", nrow(mat)))
	}
    else  { #print it out one strata at a time
	nstrat <- length(fit$strata)
	stemp <- rep(1:nstrat,fit$strata)
	strata <- factor(stemp,
	    labels = names(fit$strata)[sort(unique(stemp))])
	for (i in levels(strata)) {
	    who <- (strata==i)
	    cat("               ", i, "\n")
	    if (sum(who) ==1) print(mat[who,])
	    else    prmatrix(mat[who,], rowlab=rep("", sum(who)))
	    cat("\n")
	    }
	}
    invisible(fit)
    }
