# SCCS $Id: print.coxph.s,v 4.6 1995-12-22 17:05:18 therneau Exp $
print.coxph <-
 function(cox, digits=.Options$digits -4, ...)
    {
    if (!is.null(cl<- cox$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
    if (!is.null(cox$fail)) {
	cat(" Coxph failed.", cox$fail, "\n")
	return()
	}
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    coef <- cox$coef
    se <- sqrt(diag(cox$var))
    if(is.null(coef) | is.null(se))
        stop("Input is not valid")

    if (is.null(cox$naive.var)) {
	tmp <- cbind(coef, exp(coef), se, coef/se,
	       signif(1 - pchisq((coef/ se)^2, 1), digits -1))
	dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
	    "se(coef)", "z", "p"))
	}
    else {
	nse <- sqrt(diag(cox$naive.var))
	tmp <- cbind(coef, exp(coef), nse, se, coef/se,
	       signif(1 - pchisq((coef/se)^2, 1), digits -1))
	dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
	    "se(coef)", "robust se", "z", "p"))
	}
    cat("\n")
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
    if (length(cox$icc))
	cat("   number of clusters=", cox$icc[1],
	    "    ICC=", format(cox$icc[2:3]), "\n")
    invisible()
    }
