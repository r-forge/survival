#SCCS $Date: 1992-06-24 09:35:11 $ $Id: summary.coxph.s,v 4.1 1992-06-24 09:35:11 therneau Exp $
summary.coxph <-
 function(cox, table = T, coef = T, conf.int = 0.95, scale = 1, digits=3)
    {
    if (!is.null(cl<- cox$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
    if (!is.null(cox$fail)) {
	cat(" Coxreg failed.", cox$fail, "\n")
	return()
	}
    omit <- cox$na.action
    if (length(omit))
	cat("  n=", cox$n, " (", naprint(omit), ")\n", sep="")
    else cat("  n=", cox$n, "\n")

    if (is.null(cox$coef)) {   # Null model
	cat ("   Null model\n")
	return()
	}

    savedig <- options(digits = digits)
    on.exit(options(savedig))
    beta <- cox$coef
    se <- sqrt(diag(cox$var))
    if(is.null(beta) | is.null(se))
        stop("Input is not valid")
    if(coef) {
        tmp <- cbind(beta, exp(beta), se, beta/se, 1 - pchisq((beta/
            se)^2, 1))
        dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)", 
            "se(coef)", "z", "p"))
        cat("\n")
        prmatrix(tmp)
        }
    if(conf.int) {
        z <- qnorm((1 + conf.int)/2, 0, 1)
        beta <- beta * scale
        se <- se * scale
        tmp <- cbind(exp(beta), exp(.Uminus(beta)), exp(beta - z * se),
            exp(beta + z * se))
        dimnames(tmp) <- list(names(beta), c("exp(coef)", "exp(-coef)",
            paste("lower .", round(100 * conf.int, 2), sep = ""),
            paste("upper .", round(100 * conf.int, 2), sep = "")))
        cat("\n")
        prmatrix(tmp)
        }
    logtest <- -2 * (cox$loglik[1] - cox$loglik[2])
    sctest <- cox$score
    df <- length(beta)
    cat("\n")
    cat("Rsquare=", format(round(1-exp(-logtest/cox$n),3)),
	"  (max possible=", format(round(1-exp(2*cox$loglik[1]/cox$n),3)),
	")\n" )
    cat("Likelihood ratio test= ", format(round(logtest, 2)), "  on ",
	df, " df,", "   p=", format(1 - pchisq(logtest, df)),
	"\n", sep = "")
    cat("Efficient score test = ", format(round(sctest, 2)), "  on ", df,
        " df,", "   p=", format(1 - pchisq(sctest, df)), "\n\n", sep = 
        "")
    invisible()
    }
