#SCCS $Id: print.summary.survfit.s,v 4.1 1992-09-20 23:26:55 therneau Exp $
print.summary.survfit <- function(fit, digits=3, ...) {
    savedig <- options(digits=digits)
    on.exit(options(savedig))

    if (!is.null(cl<- fit$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
	}

    omit <- fit$na.action
    if (length(omit))
	cat(naprint(omit), "\n")

    mat <- cbind(fit$time, fit$n.risk, fit$n.event, fit$surv)
    cnames <- c("time", "n.risk", "n.event")
    if (is.matrix(fit$surv)) ncurve <- ncol(fit$surv)
    else ncurve <- 1

    if (ncurve==1) {                 #only 1 curve
	cnames <- c(cnames, "survival")
	if (!is.null(fit$std.err)) {
	    if (is.null(fit$lower)) {
		mat <- cbind(mat, fit$std.err)
		cnames <- c(cnames, "std.err")
		}
	    else {
		mat <- cbind(mat, fit$std.err, fit$lower, fit$upper)
		cnames <- c(cnames, 'std.err',
			  paste("lower ", fit$conf.int*100, "% CI", sep=''),
			  paste("upper ", fit$conf.int*100, "% CI", sep=''))
		}
	    }
	}
    else cnames <- c(cnames, paste("survival", seq(ncurve), sep=''))

    dimnames(mat) <- list(NULL, cnames)
    if (is.null(fit$strata)) {
	prmatrix(mat, rowlab=rep("", nrow(mat)))
	}
    else  { #print it out one strata at a time
	for (i in levels(fit$strata)) {
	    who <- (fit$strata==i)
	    cat("               ", i, "\n")
	    if (sum(who) ==1) print(mat[who,])
	    else    prmatrix(mat[who,], rowlab=rep("", sum(who)))
	    cat("\n")
	    }
	}
    invisible(fit)
    }
