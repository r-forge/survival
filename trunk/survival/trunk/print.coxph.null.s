#SCCS $Date: 1993-05-31 00:15:24 $ $Id: print.coxph.null.s,v 4.4 1993-05-31 00:15:24 therneau Exp $
print.coxph.null <-
 function(cox, digits=3, ...)
    {
    if (!is.null(cl<- cox$call)) {
	cat("Call:  ")
	dput(cl)
	cat("\n")
	}

    cat("Null model\n  log likelihood=", format(cox$loglik), "\n")
    omit <- cox$na.action
    if (length(omit))
	cat("  n=", cox$n, " (", naprint(omit), ")\n",
				sep="")
    else cat("  n=", cox$n, "\n")
    }
