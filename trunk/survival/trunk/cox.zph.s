# SCCS $Id: cox.zph.s,v 1.9 1993-03-04 11:57:38 therneau Exp $
#  Test proportional hazards
#
cox.zph <- function(fit, transform='km', global=T) {
    if (is.character(transform)) tname<- transform
    else        tname <- deparse(substitute(transform))
    if (!inherits(fit, 'coxph')) stop ("Argument must be the result of coxph")
    if (inherits(fit, 'coxph.null'))
	stop("The are no score residuals for a Null model")

    sresid <- resid(fit, 'schoenfeld')
    varnames <- names(fit$coef)
    nvar <- length(varnames)
    ndead<- length(sresid)/nvar
    if (nvar==1) times <- as.numeric(names(sresid))
    else         times <- as.numeric(dimnames(sresid)[[1]])

    if (is.character(transform)) {
	times <- switch(transform,
			       'identity'= times,
			       'rank'    = rank(times),
			       'km' = {
				    temp <- survfit.km(factor(rep(1,nrow(fit$y))),
							fit$y, se.fit=F)
				    1-(c(1,temp$surv))[match(times,temp$time)]
				    },
			       stop("Unrecognized transform"))
	}
    else times <- transform(times)
    xx <- times - mean(times)

    r2 <- sresid %*% fit$var * ndead
    test <- xx %*% r2        # time weighted col sums
    corel <- c(cor(xx, r2))
    z <- c(test^2 /(diag(fit$var)*ndead* sum(xx^2)))
    Z.ph <- cbind(corel, z, 1- pchisq(z,1))

    if (global && nvar>1) {
	test <- c(xx %*% sresid)
	z    <- c(test %*% fit$var %*% test) * ndead / sum(xx^2)
	Z.ph <- rbind(Z.ph, c(NA, z, 1-pchisq(z, ncol(sresid))))
	dimnames(Z.ph) <- list(c(varnames, "GLOBAL"), c("rho", "chisq", "p"))
	}
    else dimnames(Z.ph) <- list(varnames, c("rho", "chisq", "p"))

    dimnames(r2) <- dimnames(sresid)
    temp <-list(table=Z.ph, x=times, y=r2, var=fit$var, transform=tname)
    class(temp) <- "cox.zph"
    temp
    }
