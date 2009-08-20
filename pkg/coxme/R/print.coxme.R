# $Id: print.coxme.s,v 1.6 2003/08/09 21:46:42 therneau Exp $
print.coxme <- function(x, rcoef=FALSE, digits=options()$digits, ...) {
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
        coef <- beta$fixed
        se <- sqrt(diag(x$var)[nfrail+1:nvar])
        tmp <- cbind(coef, exp(coef), se, round(coef/se,2),
               signif(1 - pchisq((coef/ se)^2, 1), 2))
        dimnames(tmp) <- list(names(coef), c("coef", "exp(coef)",
            "se(coef)", "z", "p"))
        }
    if (rcoef) { # print the random coefs
        coef <- unlist(x$frail)
        se <- sqrt(diag(x$var)[1:nfrail])
        rtmp <- cbind(coef, exp(coef), se, round(coef/se,2),
               signif(1 - pchisq((coef/ se)^2, 1), 2))
        dimnames(rtmp) <- list(names(coef), c("coef", "exp(coef)",
            "se(coef)", "z", "p"))
        }
        
    if (nvar > 0 & rcoef) {
        cat("Fixed and penalized coefficients \n") 
        print(rbind(tmp, rtmp))
        }
    else if (nvar>0) {
        cat("Fixed coefficients\n")
        print(tmp)
        }
    else if (rcoef) {
        cat("Penalized coefficients\n")
        print(rtmp)
        }

    cat("\nRandom effects\n")

    random <- x$coefficients$random
    gname <- names(random)
    nrow <-  sapply(random, 
                    function(x) if (is.matrix(x)) nrow(x) else length(x))
    maxcol <-max(sapply(random,
                        function(x) if (is.matrix(x)) 1+ncol(x) else 2))
    temp1 <- matrix(NA, nrow=sum(nrow), ncol=maxcol)
    indx <- 0
    for (i in  random) {
        if (is.matrix(i)) {
            k <- ncol(i)
            temp1[1:k + indx, 1] <- sqrt(diag(i))  #std
            temp1[1:k + indx, 2] <- diag(i)        #variancw
            for (j in 1:(k-1)) temp1[j+indx, 1+(j+1):k ] <- i[j, (j+1):k]  #corr
            }
        else {
            k <- length(i)
            temp1[1:k + indx,1] <- sqrt(i)
            temp1[1:k + indx,2] <- i
            }
        indx <- indx + k
        }
        
    indx <- cumsum(c(1, nrow))   # starting row of each effect
    temp3 <- rep("", nrow(temp1))
    temp3[indx[-length(indx)]] <- names(random)
    xname <- unlist(lapply(random, 
                  function(x) if (is.matrix(x)) dimnames(x)[[1]] else names(x)))
    temp <- cbind(temp3, xname, format(temp1))
    if (maxcol == 2)
        temp4 <- c("Group", "Variable", "Std Dev", "Variance")
    else 
        temp4 <- c("Group","Variable", "Std Dev", "Variance", "Corr", 
                   rep("", maxcol-2))
    dimnames(temp) <- list(rep("", nrow(temp)), temp4)
    print(temp, quote=F)
    invisible(x)
    }
