# SCCS $Id: print.summary.survreg.s,v 4.11 1998-10-27 17:46:59 therneau Exp $
print.summary.survreg <- function(x, digits = max(options()$digits - 4, 3), quote = T, prefix = "")
{
    nas <- x$nas
    coef <- x$coef
    correl <- x$correl
    if(any(nas)) {
        nc <- length(nas)
        cnames <- names(nas)
        coef1 <- array(NA, c(nc, 3), list(cnames, dimnames(coef)[[2]]))
            
        coef1[!nas,  ] <- coef
        coef <- coef1
        if(!is.null(correl)) {
            correl1 <- matrix(NA, nc, nc, dimnames = list(cnames,
                cnames))
            correl1[!nas, !nas] <- correl
            correl <- correl1
            }
        }
    if(is.null(digits))
        digits <- options()$digits
    else options(digits = digits)
    cat("\nCall:\n")
    dput(x$call)
    dresid <- x$deviance.resid
    n <- length(dresid)
    rdf <- x$df[2]
    if(rdf > 5) {
        cat("Deviance Residuals:\n")
        rq <- quantile(as.vector(dresid))
        names(rq) <- c("Min", "1Q", "Median", "3Q", "Max")
        print(rq, digits = digits)
        }
    else if(rdf > 0) {
        cat("Deviance Residuals:\n")
        print(dresid, digits = digits)
        }
    if(any(nas))
        cat("\nCoefficients: (", sum(nas), 
            " not defined because of singularities)\n", sep = "")
        
    else cat("\nCoefficients:\n")
    print(coef, digits = digits)
    omit <- x$na.action
    if (length(omit))
	cat("  n=", n, " (", naprint(omit), ")\n", sep="")

    cat("\n", x$parms, "\n", sep='')
    int <- attr(x$terms, "intercept")
    if(is.null(int))
        int <- 1
    temp <- format(round(c(x$null.deviance, x$deviance), digits))
    cat("\n    Null Deviance:", temp[1], "on",
		     n - int, "degrees of freedom\n")
    cat("Residual Deviance:", temp[2], "on",
	   round(rdf, digits), "degrees of freedom  (LL=",
		format(x$loglik), ")\n")
    cat("Number of Newton-Raphson Iterations:", format(trunc(x$iter)),
        "\n")
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
