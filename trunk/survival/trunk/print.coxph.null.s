#SCCS $Date: 1992-04-14 18:07:21 $ $Id: print.coxph.null.s,v 4.3 1992-04-14 18:07:21 grill Exp $
print.coxph.null <-
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
