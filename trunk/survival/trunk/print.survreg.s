#SCCS %#% 7/13/92
print.survreg <- function(x, ...)
{
    if(!is.null(cl <- x$call)) {
        cat("Call:\n")
        dput(cl)
        }
    if (!is.null(x$fail)) {
	cat(" Survreg failed.", x$fail, "\n")
	return(invisible(x))
	}
    coef <- x$coef
    if(any(nas <- is.na(coef))) {
        if(is.null(names(coef))) names(coef) <- paste("b", 1:length(
                coef), sep = "")        #               coef <- coef[!nas]
        cat("\nCoefficients: (", sum(nas), 
            " not defined because of singularities)\n", sep = "")
        }
    else cat("\nCoefficients:\n")
    print(coef, ...)
    rank <- x$rank
    if(is.null(rank))
        rank <- sum(!nas)
    nobs <- length(x$residuals)
    rdf <- x$df.resid
    if(is.null(rdf))
        rdf <- nobs - rank
    omit <- x$na.action
    if (length(omit))
	cat("  n=", nobs, " (", naprint(omit), ")\n", sep="")
    sd <- survreg.distributions[[x$family[1]]]
    cat("\n", sd$print(x$parms, x$fixed), "\n", sep='')
    cat("Degrees of Freedom:", nobs, "Total;", rdf, "Residual\n")
    cat("Residual Deviance:", format(x$deviance), "\n")
    invisible(x)
    }
