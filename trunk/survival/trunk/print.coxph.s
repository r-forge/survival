# SCCS $Id: print.coxph.s,v 4.4 1993-07-04 16:08:54 therneau Exp $
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

    if (is.null(cox$robust.var)) {
	tmp <- cbind(coef, exp(coef), se, coef/se,
	       signif(1 - pchisq((coef/ se)^2, 1), digits -1))
	dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
	    "se(coef)", "z", "p"))
	}
    else {
	rse <- sqrt(diag(cox$robust.var))
	tmp <- cbind(coef, exp(coef), se, rse, coef/rse,
	       signif(1 - pchisq((coef/rse)^2, 1), digits -1))
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
    invisible()
    }
