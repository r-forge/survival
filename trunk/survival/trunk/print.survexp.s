#SCCS $Date: 1992-04-15 10:21:35 $ $Id: print.survexp.s,v 4.4 1992-04-15 10:21:35 therneau Exp $
print.survexp <- function(object, ...) {
    if (!is.null(cl<- object$call)) {
	cat("Call:  ")
	dput(cl)
	cat("\n")
	}

    omit <- object$na.action
    if (length(omit))
	cat("  n=", object$n, " (", naprint(omit), ")\n", sep="")
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
