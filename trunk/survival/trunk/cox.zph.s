# SCCS $Id: cox.zph.s,v 1.1 1992-03-24 09:27:06 therneau Exp $
#  Do the Z:PH test on a Cox model fit
#
coxreg.zph <- function(fit, ranks=T, global=T) {
    if (!inherits(fit, 'coxreg')) stop ("Argument must be the result of coxreg")
    if (inherits(fit, 'coxreg.null'))
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
