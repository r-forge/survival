# SCCS $Id: cox.zph.s,v 1.5 1992-04-30 16:37:46 therneau Exp $
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
    Z.ph <- rbind(corel, Z.ph, 2*pnorm(-abs(Z.ph)))
    dimnames(Z.ph) <- list(c("rho", "Z:ph", "p"), varnames)
    Z.ph
    }
