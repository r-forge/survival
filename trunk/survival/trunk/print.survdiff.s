#SCCS $Date: 1992-03-08 20:13:48 $ $Id: print.survdiff.s,v 4.2 1992-03-08 20:13:48 therneau Exp $
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

    omit <- attr(fit$n, 'omit')
    if (length(omit)) cat("n=", sum(fit$n), ", ", length(omit),
					  "deleted due to missing.\n\n")

    temp <- cbind(fit$n, fit$obs, fit$exp, ((fit$obs-fit$exp)^2)/ fit$exp)
    if (length(fit$n)==1)  {
	df <- fit$n   #one sample test
	dimnames(temp) <- list("One sample test",
				c("N", "Observed", "Expected",  "(O-E)^2/E"))
	}
    else {
	df <- length(fit$n) -1
	dimnames(temp) <- list(names(fit$n), c("N", "Observed", "Expected",
					   "(O-E)^2/E"))
       }
    print(temp)

    cat("\n Chisq=", format(round(fit$chisq,1)),
	     " on", df, "degrees of freedom, p=",
	     format(signif(1-pchisq(fit$chisq, df),digits)), "\n")
    invisible(fit)
    }
