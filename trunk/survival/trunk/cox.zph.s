# SCCS $Id: cox.zph.s,v 1.6 1992-05-19 20:53:25 therneau Exp $
#  Do the Z:PH test on a Cox model fit
#
cox.zph <- function(fit, ranks=T, global=T) {
    if (!inherits(fit, 'coxph')) stop ("Argument must be the result of coxph")
    if (inherits(fit, 'coxph.null'))
	stop("The are no score residuals for a Null model")

    sresid <- as.matrix(resid(fit, 'scho'))
    varnames <- names(fit$coef)

    if (global && ncol(sresid)>1) {
	sresid <- cbind(sresid, sresid %*% fit$coef)
	varnames <- c(varnames, "GLOBAL")
	}
    times <- as.numeric(dimnames(sresid)[[1]])
    if (ranks) times <- rank(times)
    corel <- cor(times, sresid)
    n <- length(times)
    Z.ph <- .5*log((1+corel)/(1-corel))*sqrt(n-3)
    Z.ph <- cbind(corel, Z.ph, 2*pnorm(-abs(Z.ph)))
    dimnames(Z.ph) <- list(varnames, c("rho", "Z:ph", "p"))
    Z.ph
    }
