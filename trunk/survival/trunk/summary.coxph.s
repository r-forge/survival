#SCCS $Date: 1995-12-22 17:05:19 $ $Id: summary.coxph.s,v 4.4 1995-12-22 17:05:19 therneau Exp $
summary.coxph <-
 function(cox, table = T, coef = T, conf.int = 0.95, scale = 1,
			digits = max(options()$digits-4,3) )
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
    savedig <- options(digits = digits)
    on.exit(options(savedig))

    omit <- cox$na.action
    if (length(omit))
	cat("  n=", cox$n, " (", naprint(omit), ")\n", sep="")
    else cat("  n=", cox$n, "\n")
    if (length(cox$icc))
	cat("  robust variance based on", cox$icc[1],
	    "groups, intra-class correlation =", format(cox$icc[2:3]), "\n")
    if (is.null(cox$coef)) {   # Null model
	cat ("   Null model\n")
	return()
	}

    beta <- cox$coef
    if(is.null(beta) | is.null(cox$var))
        stop("Input is not valid")

    if (is.null(cox$naive.var)) {
	se <- sqrt(diag(cox$var))
	wald.test <-  sum(beta * solve(cox$var, beta))
	}
    else {
	nse <- sqrt(diag(cox$naive.var))        #naive se
	se <- sqrt(diag(cox$var))
	wald.test <-  sum(beta * solve(cox$var, beta))
	}
    if(coef) {
	if (is.null(cox$naive.var)) {
	    tmp <- cbind(beta, exp(beta), se, beta/se,
		   signif(1 - pchisq((beta/ se)^2, 1), digits -1))
	    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
		"se(coef)", "z", "p"))
	    }
	else {
	    tmp <- cbind(beta, exp(beta), nse, se, beta/se,
		   signif(1 - pchisq((beta/ se)^2, 1), digits -1))
	    dimnames(tmp) <- list(names(beta), c("coef", "exp(coef)",
		"se(coef)", "robust se", "z", "p"))
	    }
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
    cat("Wald test            = ", format(round(wald.test, 2)), "  on ",
	df, " df,", "   p=", format(1 - pchisq(wald.test, df)),
	"\n", sep = "")
    cat("Efficient score test = ", format(round(sctest, 2)), "  on ", df,
        " df,", "   p=", format(1 - pchisq(sctest, df)), "\n\n", sep = 
        "")
    invisible()
    }
