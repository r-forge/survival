#SCCS $Date: 1992-03-04 16:48:16 $ $Id: print.coxph.null.s,v 4.1 1992-03-04 16:48:16 therneau Exp $
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
    omit <- attr(cox$n, 'omit')
    if (length(omit))
	cat("  n=", cox$n, " (", length(omit), " deleted due to missing)\n",
				sep="")
    else cat("  n=", cox$n, "\n")
    }
