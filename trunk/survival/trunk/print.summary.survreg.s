# SCCS @(#)print.summary.survreg.s	4.2 7/13/92
print.summary.survreg <- function(x, digits = 3, quote = T, prefix = "")
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

    cat(paste("\n(Scale Parameter for", names(x$scale),
	"family taken to be", format(round(x$scale, digits)),
        ")\n"))
    int <- attr(x$terms, "intercept")
    if(is.null(int))
        int <- 1
    cat("\n    Null Deviance:", format(round(x$null.deviance, digits)),
        "on", n - int, "degrees of freedom\n")
    cat("Residual Deviance:", format(round(x$deviance, digits)), "on",
        round(rdf, digits), "degrees of freedom\n")
    cat("Number of Newton-Rhaphson Iterations:", format(trunc(x$iter)),
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
    invisible(NULL)
    }
