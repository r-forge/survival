#SCCS $Date: 1992-03-04 16:48:24 $ $Id: print.survreg.s,v 4.1 1992-03-04 16:48:24 therneau Exp $
print.surv.reg <-
 function(cox, digits=3, ...)
    {
    if (!is.null(cl<- cox$call)) {
	cat("Call:  ")
	dput(cl)
	cat("\n")
	}
    if (!is.null(cox$fail)) {
	cat(" Coxreg failed.", cox$fail, "\n")
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
    omit <- attr(cox$n, 'omit')
    cat("\n")
    cat("Likelihood ratio test=", format(round(logtest, 2)), "  on ",
	df, " df,", " p=", format(1 - pchisq(logtest, df)),  sep="")
    if (length(omit))
	cat("\nn=", cox$n, " (", length(omit), " deleted due to missing)\n",
				sep="")
    else cat("    n=", cox$n, "\n")
    invisible()
    }
