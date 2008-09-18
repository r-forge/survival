# $Id$
print.summary.coxph <-
 function(x, coef = TRUE, conf.int=T,
			digits = max(options()$digits - 4, 3)) {

    if (!is.null(cl<- x$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
    if (!is.null(x$fail)) {
	cat(" Coxreg failed.", x$fail, "\n")
	return()
	}
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    omit <- x$na.action
    if (length(omit))
	cat("  n=", x$n, " (", naprint(omit), ")\n", sep="")
    else cat("  n=", x$n, "\n")
    if (length(x$icc))
	cat("  robust variance based on", x$icc[1],
	    "groups, intra-class correlation =", format(x$icc[2:3]), "\n")
    if (is.null(x$coef)) {   # Null model
	cat ("   Null model\n")
	return()
	}

    if(coef) {
        cat("\n")
        prmatrix(x$coefficient)
        }
    if(conf.int) {
        cat("\n")
        prmatrix(x$cimat)
        }
 
    logtest <- -2 * (x$loglik[1] - x$loglik[2])
    sctest <- x$score
    df <- sum(!is.na(x$coefficient[,1]))

    cat("\n")
    cat("Rsquare=", format(round(1-exp(-logtest/x$n),3)),
	"  (max possible=", format(round(1-exp(2*x$loglik[1]/x$n),3)),
	")\n" )
    cat("Likelihood ratio test= ", format(round(logtest, 2)), "  on ",
	df, " df,", "   p=", format(1 - pchisq(logtest, df)),
	"\n", sep = "")
    cat("Wald test            = ", format(round(x$wald.test, 2)), "  on ",
	df, " df,", "   p=", format(1 - pchisq(x$wald.test, df)),
	"\n", sep = "")
    cat("Score (logrank) test = ", format(round(sctest, 2)), "  on ", df,
        " df,", "   p=", format(1 - pchisq(sctest, df)), sep ="") 
    if (is.null(x$rscore)) cat("\n\n")
    else cat(",   Robust = ", format(round(x$rscore, 2)), 
	   "  p=", format(1 - pchisq(x$rscore, df)), "\n\n", sep="")   

    if (!is.null(x$naive.var))
	cat("  (Note: the likelihood ratio and score tests",
	  "assume independence of\n     observations within a cluster,",
	    "the Wald and robust score tests do not).\n")
    invisible()
    }
