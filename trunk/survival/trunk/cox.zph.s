# SCCS $Id: cox.zph.s,v 1.8 1993-01-13 01:01:00 therneau Exp $
#  Test proportional hazards
#
cox.zph <- function(fit, transform='rank', global=F) {
    if (!inherits(fit, 'coxph')) stop ("Argument must be the result of coxph")
    if (inherits(fit, 'coxph.null'))
	stop("The are no score residuals for a Null model")

    sresid <- as.matrix(resid(fit, 'scaledsch'))
    varnames <- names(fit$coef)

    times <- as.numeric(dimnames(sresid)[[1]])
    if (is.character(transform)) {
	times <- switch(transform,
			       'identity'= times,
			       'rank'    = rank(times),
			       stop("Unrecognized transform"))
	}
    else times <- transform(times)

    ndead<- sum(fit$y[,ncol(fit$y)])
    test <- apply((times/sum(times))*sresid, 2, sum)

    if (global) {
	z <- (test%*% fit$var) %*%test * ndead^2 / (ndead-1)
	z <- c(chisq=z, p=pchisq(z, ncol(sresid)))
	z
	}
    else {
	corel <- cor(sresid, times)
	z <- test^2 * diag(fit$var)* ndead^2/(ndead-1)
	Z.ph <- cbind(corel, z, pchisq(z,1))
	dimnames(Z.ph) <- list(varnames, c("rho", "chisq", "p"))
	Z.ph
	}
    }
