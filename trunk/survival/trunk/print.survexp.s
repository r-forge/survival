#SCCS $Id: print.survexp.s,v 4.7 1994-01-06 11:02:46 therneau Exp $
print.survexp <- function(fit, scale=1, digits=3, naprint=F, ...) {
    if (!inherits(fit, 'survexp'))
	    stop("Invalid data")
    savedig <- options(digits=digits)
    on.exit(options(savedig))

    if (!is.null(cl<- fit$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}

    if (!is.null(fit$summ)) cat(fit$summ)
    omit <- fit$na.action
    if (length(omit))
	cat(naprint(omit), "\n")
    else cat("\n")

    if (is.null(fit$strata))  { #print it as a matrix
	mat <- cbind(fit$time/scale, fit$n.risk, fit$surv)
	if (!naprint) {
	    miss <- (is.na(mat)) %*% rep(1,ncol(mat))
	    mat <- mat[miss<(ncol(mat)-2),,drop=F]
	    }
	prmatrix(mat, rowlab=rep("", nrow(mat)),
		   collab=c("Time", "n.risk", dimnames(fit$surv)[[2]]))
	}
    else  { #print it out one strata at a time, since n's differ
	nstrat <- length(fit$strata)
	levs <- names(fit$strata)
	if (nrow(fit$surv)==1) {
	    mat <- cbind(c(fit$n.risk), c(fit$surv))
	    dimnames(mat) <- list(levs, c("n.risk", "survival"))
	    cat(" Survival at time", fit$time, "\n")
	    prmatrix(mat)
	    }
	else {
	    for (i in 1:nstrat) {
		cat("       ", levs[i], "\n")
		mat <- cbind(fit$time/scale, fit$n.risk[,i], fit$surv[,i])
		if (!naprint) mat <- mat[!is.na(mat[,3]),,drop=F]
		prmatrix(mat, rowlab=rep("",nrow(mat)),
				collab=c("Time", "n.risk", "survival"))
		cat("\n")
		}
	    }
	}
    invisible(fit)
    }
