#SCCS $Date: 1992-03-30 02:43:03 $ $Id: print.coxph.null.s,v 4.2 1992-03-30 02:43:03 therneau Exp $
print.coxreg.null <-
 function(cox, digits=3, ...)
    {
    if (!is.null(cl<- cox$call)) {
	cat("Call:  ")
	dput(cl)
	cat("\n")
	}

    cat("Null model\n  log likelihood=", format(cox$loglik),
		  "  (Breslow approx)\n")
    omit <- cox$na.action
    if (length(omit))
	cat("  n=", cox$n, " (", naprint(omit), ")\n",
				sep="")
    else cat("  n=", cox$n, "\n")
    }
