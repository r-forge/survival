#SCCS $Date: 1992-03-04 16:48:22 $ $Id: print.survexp.s,v 4.1 1992-03-04 16:48:22 therneau Exp $
print.surv.exp <- function(object, ...) {
    if (!is.null(cl<- object$call)) {
	cat("Call:  ")
	dput(cl)
	cat("\n")
	}

    omit <- attr(object$n, 'omit')
    if (length(omit))
	cat("  n=", object$n, " (", length(omit), " deleted due to missing)\n",
				sep="")
    else cat("  n=", object$n, "\n")

    temp <- as.matrix(object$surv)
    if (ncol(temp) ==1)
	dimnames(temp) <- list(object$time, "Expected")
    else
	dimnames(temp) <- list(object$time,
			       paste("Expected", seq(ncol(temp)), sep=''))
    print.matrix(temp)
    invisible(temp)
    }
