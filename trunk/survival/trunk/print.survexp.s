#SCCS $Date: 1992-04-14 18:07:32 $ $Id: print.survexp.s,v 4.3 1992-04-14 18:07:32 grill Exp $
print.survexp <- function(object, ...) {
    if (!is.null(cl<- object$call)) {
	cat("Call:  ")
	dput(cl)
	cat("\n")
	}

    omit <- object$na.action
    if (length(omit))
	cat("  n=", object$n, " (", printnaomit), ")\n", sep="")
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
