# SCCS $Id: print.summary.survreg.s,v 4.12 1998-11-30 08:29:14 therneau Exp $
print.summary.survreg <- function(x, digits = max(options()$digits - 4, 3), 
				  quote = T, prefix = "") {
    correl <- x$correl
    n <- x$n

    if(is.null(digits))
        digits <- options()$digits
    else options(digits = digits)
    cat("\nCall:\n")
    dput(x$call)

    print(x$table, digits = digits)
    if (nrow(x$var)==length(x$coefficients)) 
	    cat("\nScale fixed at",format(x$scale),"\n") 
    else if (length(x$scale)==1) cat ("\nScale=", format(x$scale), "\n")
    else {
	cat("\nScale:\n")
	print(x$scale, ...)
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
            correl[ll] <- format(round(correl[ll], digits))
            correl[!ll] <- ""
            print(correl[-1,  - p, drop = F], quote = F, digits = 
                digits)
            }
        }
    cat("\n")
    invisible(NULL)
    }
