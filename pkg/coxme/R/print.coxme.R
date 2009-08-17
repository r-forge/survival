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
    if (rcoef) { #random coefs
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
    # for each term: grouping name and nrow/ncol for each group
    tfun <- function(x) {
        gname <- as.list(names(x))
        nrow <- ncol <- xname <- vector('list', length(x))
        for (i in 1:length(x)) {
            if (is.list(x[[i]])) {
                temp <- tfun(x[[i]])
                gname[[i]] <- temp$gname
                nrow[[i]] <- temp$nrow
                ncol[[i]] <- temp$ncol
                xname[[i]]<- temp$xname
                }
            else {
                if (is.matrix(x[[i]])) {
                    nrow[[i]] <- nrow(x[[i]])
                    ncol[[i]] <- ncol(x[[i]])
                    xname[[i]] <- dimnames(x[[i]])[[1]]
                    }
                else {
                    nrow[[i]] <- length(x[[i]])
                    ncol[[i]] <- 1
                    xname <- names(x[[i]])
                    }
                }
            }
        list(gname=unlist(gname), nrow=unlist(nrow), ncol=unlist(ncol),
             xname=unlist(xname))
        }

    rtype <- tfun(beta$random)
    maxcol <- max(rtype$ncol)
    temp1 <- matrix(NA, nrow=sum(rtype$nrow), ncol=1+ maxcol)
    indx <- 0
    for (i in  beta$random) {
        if (!is.list(i)) i <- list(i)
        for (j in i) {
            if (is.matrix(j)) {
                k <- ncol(j)
                temp1[1:k + indx, 1] <- sqrt(diag(j))
                for (i in 1:k) temp1[i+indx, i:k +1] <- j[i, i:k]
                }
            else {
                k <- length(j)
                temp1[1:k + indx,1] <- sqrt(j)
                temp1[1:k + indx,2] <- j
                }
            }
        }

    indx <- cumsum(c(1, rtype$nrow))
    temp3 <- rep("", nrow(temp1))
    temp3[indx[-length(indx)]] <- rtype$gname
    temp <- cbind(temp3, rtype$xname, format(temp1))
    if (maxcol == 1)
        temp4 <- c("Group", "Variable", "Std Dev", "Variance")
    else 
        temp4 <- c("Group","Variable", "Std Dev", "Var/Cov", rep("", maxcol-2))
    dimnames(temp) <- list(rep("", nrow(temp)), temp4)
    print(temp, quote=F)
    invisible(x)
    }
