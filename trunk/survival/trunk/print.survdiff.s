#SCCS $Date: 1992-03-04 16:48:21 $ $Id: print.survdiff.s,v 4.1 1992-03-04 16:48:21 therneau Exp $
print.surv.diff <- function(diff.list, digits=4, ...) {

    fit <- diff.list
    saveopt <-options(digits=digits)
    on.exit(options(saveopt))

    if (!inherits(diff.list, 'surv.diff'))
	stop("Object is not the result of surv.diff")
    if (!is.null(cl<- fit$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}

    temp <- cbind(fit$n, fit$obs, fit$exp, ((fit$obs-fit$exp)^2)/ fit$exp)
    dimnames(temp) <- list(names(fit$n), c("N", "Observed", "Expected",
					   "(O-E)^2/E"))
    print(temp)
    df <- length(fit$n) -1
    cat("\n Chisq=", format(round(fit$chisq,1)),
	     " on", df, "degrees of freedom, p=",
	     format(signif(1-pchisq(fit$chisq, df),digits)), "\n")
    invisible()
    }
