#SCCS $Date: 1997-12-29 12:31:23 $ $Id: print.coxph.null.s,v 4.7 1997-12-29 12:31:23 boos Exp $
print.coxph.null <-
 function(x, digits=max(options()$digits - 4, 3), ...)
    {
    if (!is.null(cl<- x$call)) {
	cat("Call:  ")
	dput(cl)
	cat("\n")
	}

    cat("Null model\n  log likelihood=", format(x$loglik), "\n")
    omit <- x$na.action
    if (length(omit))
	cat("  n=", x$n, " (", naprint(omit), ")\n",
				sep="")
    else cat("  n=", x$n, "\n")
    }
