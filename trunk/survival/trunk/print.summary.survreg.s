# SCCS $Id: print.summary.survreg.s,v 4.13 1999-01-06 07:21:40 therneau Exp $
print.summary.survreg <- function(x, digits = max(options()$digits - 4, 3), 
				  quote = T, prefix = "") {
    correl <- x$correl
    n <- x$n

    if(is.null(digits))
        digits <- options()$digits
    cat("\nCall:\n")
    dput(x$call)

    print(x$table, digits = digits)
    if (nrow(x$var)==length(x$coefficients)) 
	    cat("\nScale fixed at",format(x$scale, digits=digits),"\n") 
    else if (length(x$scale)==1) 
	    cat ("\nScale=", format(x$scale, digits=digits), "\n")
    else {
	cat("\nScale:\n")
	print(x$scale, digits=digits, ...)
	}

    cat("\n", x$parms, "\n", sep='')
    df  <- sum(x$df) - x$idf   # The sum is for penalized models
    cat("Loglik(model)=", format(round(x$loglik[2],1)),
	"  Loglik(intercept only)=", format(round(x$loglik[1],1)))
    if (df > 0)
	    cat("\n\tChisq=", format(round(x$chi,2)), "on", round(df,1),
		"degrees of freedom, p=", 
		format(signif(1-pchisq(x$chi, df),2)), "\n")
    else cat("\n")
    cat("Number of Newton-Raphson Iterations:", format(trunc(x$iter)),
        "\n")
    omit <- x$na.action
    if (length(omit))
	cat("n=", x$n, " (", naprint(omit), ")\n", sep="")
    else cat("n=", x$n, "\n")

    if(!is.null(correl)) {
        p <- dim(correl)[2]
        if(p > 1) {
            cat("\nCorrelation of Coefficients:\n")
            ll <- lower.tri(correl)
            correl[ll] <- format(round(correl[ll], digits=digits))
            correl[!ll] <- ""
            print(correl[-1,  - p, drop = F], quote = F)
            }
        }
    cat("\n")
    invisible(NULL)
    }
