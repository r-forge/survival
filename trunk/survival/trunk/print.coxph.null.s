#SCCS $Date: 1996-09-27 11:02:48 $ $Id: print.coxph.null.s,v 4.6 1996-09-27 11:02:48 boos Exp $
print.coxph.null <-
 function(cox, digits=max(options()$digits - 4, 3), ...)
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
