# 
#  SCCS $Id: survreg.fit.s,v 5.4 1998-11-30 08:34:08 therneau Exp $
#
survreg.fit<- function(x, y, weights, offset, init, controlvals, dist, 
		       scale=0, nstrat=1, strata) {

    iter.max <- controlvals$iter.max
    eps <- controlvals$rel.tol
    toler.chol <- controlvals$toler.chol
    debug <- controlvals$debug

    if (!is.matrix(x)) stop("Invalid X matrix ")
    n <- nrow(x)
    nvar <- ncol(x)
    ny <- ncol(y)
    if (is.null(offset)) offset <- rep(0,n)
    if (missing(weights)|| is.null(weights)) weights<- rep(1.0,n)
    else if (any(weights<=0)) stop("Invalid weights, must be >0")

    if (scale <0) stop("Invalid scale")
    if (scale >0 && nstrat >1) 
	    stop("Cannot have both a fixed scale and strata")
    if (nstrat>1 && (missing(strata) || length(strata)!= n))
	    stop("Invalid strata variable")
    if (nstrat==1) strata <- rep(1,n)
    if (scale >0) nstrat2 <- 0; else nstrat2 <- nstrat

    if (is.character(dist)) {
	sd <- survreg.distributions[[dist]]
	if (is.null(sd)) stop ("Unrecognized distribution")
	}
    else sd <- dist
    dnum <- match(sd$name, c("Extreme value", "Logistic", "Gaussian"))
    if (is.na(dnum)) {
	# Not one of the "built-in distributions
	dnum <- 4
	if (!is.function(sd$density)) 
		stop("Missing density function in user defined distribution")
	stop ('function not yet finished')
	}
    #
    # Fit the model with just a mean and scale
    #    assume initial values don't apply here
    # Unless, of course, someone is fitting a mean only model!
    #
    meanonly <- (nvar==1 && all(x==1))
    if (!meanonly) {
	yy <- ifelse(y[,ny]!=3, y[,1], (y[,1]+y[,2])/2 )
	coef <- sd$init(yy, weights)
	if (scale >0) coef[2] <- scale
	variance <- log(coef[2])/2   # init returns \sigma^2, I need log(sigma)
	coef <- c(coef[1], rep(variance, nstrat))
	# get a better initial value for the mean using the "glim" trick
	deriv <- .C("survreg3",
		    as.integer(n),
		    as.double(y),
		    as.integer(ny),
		    as.double(yy),
		    nstrat = as.integer(nstrat2),
		    strata = as.integer(strata),
		    vars= as.double(coef[-1]),
		    deriv = matrix(double(n * 3),nrow=n),
		    as.integer(3),
		    as.integer(dnum))$deriv
	coef[1] <- coef[1] - sum(weights*deriv[,2])/sum(weights*deriv[,3])

	# Now the fit proper (intercept only)
	nvar2 <- 1 +nstrat2
	fit0 <- .C("survreg2",
		       iter = as.integer(iter.max),
		       as.integer(n),
		       as.integer(1),
		       as.double(y),
		       as.integer(ny),
		       rep(1.0, n),
		       as.double(weights),
		       as.double(offset),
		       coef= as.double(coef),
		       as.integer(nstrat2),
		       as.integer(strata),
		       u = double(3*(nvar2) + nvar2^2),
		       var = matrix(0.0, nvar2, nvar2),
		       loglik=double(1),
		       flag=integer(1),
		       as.double(eps),
		       as.double(toler.chol), 
		       as.integer(dnum),
		       debug = as.integer(floor(debug/2)))
	}

    #
    # Fit the model with all covariates
    #
    nvar2 <- nvar + nstrat2
    if (is.numeric(init)) {
	if (length(init) != nvar2) stop("Wrong length for initial parameters")
	}
    else  {
	# Do the 'glim' method of finding an initial value of coef
	if (meanonly) {
	    yy <- ifelse(y[,ny]!=3, y[,1], (y[,1]+y[,2])/2 )
	    coef <- sd$init(yy, weights)
	    if (scale >0) coef[2] <- scale
	    vars  <- rep(log(coef[2])/2, nstrat)  
	    }
	else vars <- fit0$coef[-1]
	eta <- yy - offset     #what would be true for a 'perfect' model

	deriv <- .C("survreg3",
		       as.integer(n),
		       as.double(y),
		       as.integer(ny),
		       as.double(yy),
		       nstrat = as.integer(nstrat2),
		       strata = as.integer(strata),
		       vars= as.double(vars),
		       deriv = matrix(double(n * 3),nrow=n),
		       as.integer(3),
		       as.integer(dnum))$deriv

	wt <-  -1*deriv[,3]*weights
	coef <- solve(t(x)%*% (wt*x), c((wt*eta + weights*deriv[,2])%*% x))
	init <- c(coef, vars)
	}

    # Now for the fit in earnest
    fit <- .C("survreg2",
		   iter = as.integer(iter.max),
		   as.integer(n),
		   as.integer(nvar),
		   as.double(y),
		   as.integer(ny),
		   as.double(x),
	           as.double(weights),
		   as.double(offset),
		   coef= as.double(init),
	           as.integer(nstrat2),
	           as.integer(strata),
		   u = double(3*(nvar2) + nvar2^2),
		   var = matrix(0.0, nvar2, nvar2),
		   loglik=double(1),
		   flag=integer(1),
		   as.double(eps),
	           as.double(toler.chol), 
		   as.integer(dnum),
	           debug = as.integer(debug))

    if (debug>0) browser()
    if (iter.max >1 && fit$flag <nvar) {
	if (controlvals$failure==1)
	       warning("Ran out of iterations and did not converge")
	else if (controlvals$failure==2)
	       return("Ran out of iterations and did not converge")
	}

    cname <- dimnames(x)[[2]]
    if (is.null(cname)) cname <- paste("x", 1:ncol(x))
    if (scale==0) cname <- c(cname, rep("Log(scale)", nstrat))
    dimnames(fit$var) <- list(cname, cname)
    if (scale>0) fit$coef <- fit$coef[1:nvar2]
    names(fit$coef) <- cname

    if (meanonly) {
	coef0 <- fit$coef
	loglik <- rep(fit$loglik,2)
	}
    else {
	coef0 <- fit0$coef
	names(coef0) <- c("Intercept", rep("Log(scale)", nstrat))
	loglik <- c(fit0$loglik, fit$loglik)
	}
    temp <- list(coefficients   = fit$coef,
		 icoef  = coef0, 
		 var    = fit$var,
		 loglik = loglik, 
		 iter   = fit$iter,
		 linear.predictors = c(x %*% fit$coef[1:nvar]),	
		 df     = length(fit$coef)
		 )
    if (debug>0) {
	temp$u <- fit$u[1:nvar2]
	JJ     <- matrix(fit$u[-seq(1, 3*nvar2)], nvar2, nvar2)
	temp$JJ <- JJ
	temp$var2 <- fit$var %*% JJ %*% fit$var
	}

    temp
    }
