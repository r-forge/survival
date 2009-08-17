# $Id: print.coxme.s,v 1.6 2003/08/09 21:46:42 therneau Exp $
print.coxme <- function(x, digits=options()$digits, ...) {
    cat("Cox mixed-effects model fit by maximum likelihood\n")
    if (!is.null(x$call$data)) 
        cat("  Data:", deparse(x$call$data))
    if(!is.null(x$call$subset)) {
        cat(";  Subset:", deparse(x$call$subset), "\n")
	}
    else cat("\n")

    beta <- x$coefficients
    nvar <- length(beta$fixed)
    nfrail<- nrow(x$var) - nvar

    omit <- x$na.action
    if(length(omit))
        cat("  n=", x$n, " (", naprint(omit), ")\n", sep = "")
    else cat("  n=", x$n, "\n")
    temp <- matrix(x$loglik, nrow=1)
    cat("  Iterations=", x$iter, "\n")
    dimnames(temp) <- list("Log-likelihood", 
                           c("NULL", "Integrated", "Penalized"))
    print(temp)
    chi <- 2*diff(x$loglik[c(1,3)]) 
    cat("\n  Penalized loglik: chisq=", format(round(chi,2)), 
        "on", format(round(x$df[2],2)), "degrees of freedom, p=",
        format(signif(1- pchisq(chi,x$df[2]),2)),"\n")
    chi <- 2*diff(x$loglik[1:2]) 
    cat(" Integrated loglik: chisq=", format(round(chi,2)), 
        "on", format(round(x$df[1],2)), "degrees of freedom, p=",
        format(signif(1- pchisq(chi,x$df[1]),2)),"\n\n")

    cat ("Model: ", deparse(x$call$formula), "\n")
    if (nvar > 0)  { # Not a ~1 model
        cat("Fixed effects \n")
        coef <- beta$fixed
        se <- sqrt(diag(x$var)[nfrail+1:nvar])
        tmp <- cbind(coef, exp(coef), se, round(coef/se,2),
               signif(1 - pchisq((coef/ se)^2, 1), 2))
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
            "se(coef)", "z", "p"))
        print(tmp)
        }

    cat("\nRandom effects\n")
    random <- beta$random
    rlength <- unlist(lapply(random, length))
    rrow <- ceiling(rlength+1)/2  #Assume 1 element=1 row, 3= 2 rows, 5=3,...
    temp <- matrix("", nrow=sum(rrow), ncol=5)
    temp[cumsum(rrow),1] <- names(beta$random)
    nocor <- function(x) x[seq(1, length=length(x), by=2)]
    temp[,2] <- unlist(lapply(random, function(x) nocor(names(x))))
    temp[,3] <- format(unlist(lapply(random, nocor)))
    temp[,4] <- format(sqrt(unlist(lapply(random, nocor))))
    dimnames(temp) <- list(rep("", sum(rrow)), 
                           c("Group", "Variable", "Variance", "Std Dev", "Corr"))
    print(temp, quote=F)
    invisible(x)
    }
