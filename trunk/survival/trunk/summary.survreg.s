# SCCS $Id: summary.survreg.s,v 4.9 1992-12-30 14:21:45 therneau Exp $
summary.survreg<- function(object, correlation = T)
{
    if (!is.null(object$fail)) {
	warning(" Survreg failed.", x$fail, "   No summary provided\n")
	return(invisible(object))
	}
    wt <- object$weights
    fparms <- object$fixed
    coef <- c(object$coef, object$parms[!fparms])
    resid <- object$residuals
    dresid <- object$dresiduals
    n <- length(resid)
    p <- sum(!is.na(coef))
    if(!p) {
        warning("This model has zero rank --- no summary is provided")
        return(object)
        }
    nsingular <- length(coef) - p
    rdf <- object$df.resid
    if(is.null(rdf))
        rdf <- n - p
    R <- object$R   #check for rank deficiencies
    if(p < max(dim(R)))
        R <- R[1:p,     #coded by pivoting
        1:p]
    if(!is.null(wt)) {
        wt <- wt^0.5
        resid <- resid * wt
        excl <- wt == 0
        if(any(excl)) {
            warning(paste(sum(excl), 
                "rows with zero weights not counted"))
            resid <- resid[!excl]
            if(is.null(object$df.residual))
                rdf <- rdf - sum(excl)
            }
        }
    famname <- object$family["name"]
    if(is.null(famname))
	famname <- "gaussian"
    nas <- is.na(coef)
    cnames <- names(coef[!nas])
    coef <- matrix(rep(coef[!nas], 4), ncol = 4)
    dimnames(coef) <- list(cnames, c("Value", "Std. Error", "z value", "p"))
    stds <- sqrt(diag(object$var[!nas,!nas,drop=F]))
    coef[, 2] <- stds
    coef[, 3] <- coef[, 1]/stds
    coef[, 4] <- 2*pnorm(-abs(coef[,3]))
    if(correlation && sum(!nas)>1 ) {
	correl <- diag(1/stds) %*% object$var[!nas, !nas] %*% diag(1/stds)
        dimnames(correl) <- list(cnames, cnames)
        }
    else correl <- NULL
    ocall <- object$call
    if(!is.null(form <- object$formula)) {
        if(is.null(ocall$formula))
	    ocall <- match.call(get("survreg"), ocall)
        ocall$formula <- form
        }
    sd <- survreg.distributions[[famname]]
    pprint<- paste(sd$name, 'distribution:', sd$print(object$parms, fparms))
    structure(list(call = ocall, terms = object$terms, coefficients = coef,
	scale= scale, df = c(p, rdf), deviance.resid = dresid,
	var=object$var, correlation = correl, deviance = deviance(object),
	null.deviance = object$null.deviance, iter = object$iter,
	nas = nas, parms=pprint, loglik=object$loglik[2]),
	class = "summary.survreg")
    }
