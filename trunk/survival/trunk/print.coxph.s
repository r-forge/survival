# SCCS $Id: print.coxph.s,v 4.2 1992-08-05 13:45:12 therneau Exp $
print.coxph <-
 function(cox, digits=3, ...)
    {
    if (!is.null(cl<- cox$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
    if (!is.null(cox$fail)) {
	cat(" Coxrph failed.", cox$fail, "\n")
	return()
	}
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    coef <- cox$coef
    se <- sqrt(diag(cox$var))
    if(is.null(coef) | is.null(se))
        stop("Input is not valid")
    tmp <- cbind(coef, exp(coef), se, coef/se, 1 - pchisq((coef/
	se)^2, 1))
    dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
	"se(coef)", "z", "p"))
    prmatrix(tmp)

    logtest <- -2 * (cox$loglik[1] - cox$loglik[2])
    df <- length(coef)
    cat("\n")
    cat("Likelihood ratio test=", format(round(logtest, 2)), "  on ",
	df, " df,", " p=", format(1 - pchisq(logtest, df)),  sep="")
    omit <- cox$na.action
    if (length(omit))
	cat("  n=", cox$n, " (", naprint(omit), ")\n", sep="")
    else cat("  n=", cox$n, "\n")
    invisible()
    }
