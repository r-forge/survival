#SCCS $Date: 1992-03-30 02:45:56 $ $Id: print.survexp.s,v 4.2 1992-03-30 02:45:56 therneau Exp $
print.surv.exp <- function(object, ...) {
    if (!is.null(cl<- object$call)) {
	cat("Call:  ")
	dput(cl)
	cat("\n")
	}

    omit <- object$na.action
    if (length(omit))
	cat("  n=", object$n, " (", printna(omit), ")\n", sep="")
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
