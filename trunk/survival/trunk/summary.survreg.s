# SCCS $Id: summary.survreg.s,v 4.11 1998-11-30 08:29:17 therneau Exp $
summary.survreg<- function(object, correlation = T)
{
    if (!is.null(object$fail)) {
	warning(" Survreg failed.", x$fail, "   No summary provided\n")
	return(invisible(object))
	}
    wt <- object$weights
    
    nvar0 <- length(object$coef)
    nvar <- nrow(object$var)
    if (nvar > nvar0) {
	coef <- c(object$coef, log(object$scale))
	if ( (nvar-nvar0)==1) cname <- c(names(object$coef), "Log(scale)")
	else cname <- c(names(object$coef), names(object$scale))
	}
    else {
	coef <- object$coef
	cname <- names(object$coef)
	}

    n <- length(object$linear.predictors)
    p <- sum(!is.na(coef))
    if(!p) {
        warning("This model has zero rank --- no summary is provided")
        return(invisible(object))
        }

    nsingular <- nvar - p
    table <- matrix(rep(coef, 4), ncol = 4)
    dimnames(table) <- list(cname, c("Value", "Std. Error", "z", "p"))
    stds <- sqrt(diag(object$var))
    table[, 2] <- stds
    table[, 3] <- table[, 1]/stds
    table[, 4] <- 2*pnorm(-abs(table[,3]))
    if(correlation) {
	nas <- is.na(coef)
	stds <- stds[!nas]
	correl <- diag(1/stds) %*% object$var[!nas, !nas] %*% diag(1/stds)
        dimnames(correl) <- list(cname, cname)
        }
    else correl <- NULL

    dist <- object$dist
    if (is.character(dist)) sd <- survreg.distributions[[dist]]
    else sd <- dist

    if (length(object$parms)) 
	    pprint<- paste(sd$name, 'distribution: parmameters=', object$parms)
    else    pprint<- paste(sd$name, 'distribution')

    x <- object[c('call', 'df', 'loglik', 'iter', 'na.action', 'idf', 'scale',
		  'coefficients', 'var')]
    x$table <- table
    x$correlation  <- correl
    x$parms <- pprint
    x$n <- n
    x$chi <- 2*diff(object$loglik)

    oldClass(x) <- 'summary.survreg'
    x
    }
