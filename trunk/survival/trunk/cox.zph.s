# SCCS $Id: cox.zph.s,v 1.3 1992-04-14 18:33:55 splus Exp $
#  Do the Z:PH test on a Cox model fit
#
cox.zph <- function(fit, ranks=T, global=T) {
    if (!inherits(fit, 'coxph')) stop ("Argument must be the result of coxph")
    if (inherits(fit, 'coxph.null'))
	stop("The are no score residuals for a Null model")

    sresid <- resid(fit, 'scho')
    names <- dimnames(sresid)[[2]]
    if (global) {
	sresid <- cbind(sresid, sresid %*% fit$coef)
	names <- c(names, "GLOBAL")
	}
    times <- as.numeric(dimnames(sresid)[[1]])
    if (ranks) times <- rank(times)
    corel <- cor(times, sresid)
    n <- length(times)
    Z.ph <- .5*log((1+corel)/(1-corel))*sqrt(n-3)
    Z.ph <- rbind(corel, Z.ph, 2*pnorm(-abs(Z.ph)))
    dimnames(Z.ph) <- list(c("rho", "Z:ph", "p"), names)
    Z.ph
    }
