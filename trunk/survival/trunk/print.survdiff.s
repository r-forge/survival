#SCCS $Date: 1996-01-07 01:35:51 $ $Id: print.survdiff.s,v 4.9 1996-01-07 01:35:51 therneau Exp $
print.survdiff <- function(fit, digits=max(options()$digits-4,3), ...) {

    saveopt <-options(digits=digits)
    on.exit(options(saveopt))

    if (!inherits(fit, 'survdiff'))
	stop("Object is not the result of survdiff")
    if (!is.null(cl<- fit$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}

    omit <- fit$na.action
    if (length(omit)) cat("n=", sum(fit$n), ", ", naprint(omit),
					  ".\n\n", sep='')

    if (length(fit$n)==1)  {
	z <- sign(fit$exp - fit$obs) * sqrt(fit$chisq)
	temp <- c(fit$obs, fit$exp, z, signif(1-pchisq(fit$chisq, 1),digits))
	names(temp) <- c("Observed", "Expected", "Z", "p")
	print(temp)
	}
    else {
	if (is.matrix(fit$obs)){
	    otmp <- apply(fit$obs,1,sum)
	    etmp <- apply(fit$exp,1,sum)
	    }
	else {
	    otmp <- fit$obs
	    etmp <- fit$exp
	    }
	df <- (sum(1*(etmp>0))) -1
	temp <- cbind(fit$n, otmp, etmp, ((otmp-etmp)^2)/ etmp,
					 ((otmp-etmp)^2)/ diag(fit$var))
	dimnames(temp) <- list(names(fit$n), c("N", "Observed", "Expected",
				  "(O-E)^2/E", "(O-E)^2/V"))
	print(temp)
	cat("\n Chisq=", format(round(fit$chisq,1)),
		 " on", df, "degrees of freedom, p=",
		 format(signif(1-pchisq(fit$chisq, df),digits)), "\n")
       }
    invisible(fit)
    }
