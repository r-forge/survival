#SCCS $Id: survreg.fit.s,v 4.3 1992-07-10 09:18:14 therneau Exp $
#
# This handles the one parameter distributions-- extreme, logistic, and
#       gaussian
survreg.fit<- function(x, y, offset, init, iter.max,
			eps, fun, scale, singular.ok=F)
    {
    if (!is.matrix(x)) stop("Invalid formula for fitting function")
    n <- nrow(x)
    nvar <- ncol(x)
    ny <- ncol(y)
    if (is.null(offset)) offset <- rep(0,n)
    if (missing(scale) || is.null(scale)) {
	scale <- 0   #needs to be estimated
	nvar2 <- nvar +1
	}
    else {
	nvar2 <- nvar
	if (scale <=0)  stop("Invalid value for the scale parameter")
	}

    if (!missing(init) &&  !is.null(init)) {
	if (length(init) != nvar2) stop("Wrong length for inital values")
	if (any(is.na(init))) stop("Init contains missing values")
	if (scale==0) {
	    if (init[nvar2] <=0)
		stop("Invalid initial value for scale")
	    else init[nvar2] <- log(init[nvar2])
	    }
	}
    mysig <- log(var(y[,1]))

    if (scale==0 && any(y[,ny]==3)) ncol.deriv <- 6
    else                            ncol.deriv <- 3
    #
    # Fit the null model --- only an intercept and sigma
    #
    nfit <- .C(fun, iter= as.integer(iter.max),
		   as.integer(n),
		   as.integer(1),
		   y,
		   as.integer(ny),
		   rep(1,n),
		   as.double(offset),
		   coef= as.double(c(0,mysig)),
		   as.double(scale),
		   double(6),
		   imat= matrix(0, 2,2),
		   loglik=double(2),
		   flag=integer(1),
		   as.double(eps),
		   deriv = double(n * ncol.deriv),
		   as.integer(1))

    if (missing(init) || is.null(init)) {
	init <- c(rep(0,nvar), nfit$coef[2])[1:nvar2]
	doinit <- 1
	}
    else doinit <- 0

    fit <- .C(fun, iter=as.integer(iter.max),
		   as.integer(n),
		   as.integer(nvar),
		   y,
		   as.integer(ny),
		   x,
		   as.double(offset),
		   coef= as.double(init),
		   as.double(scale),
		   double(3*(nvar2) + 2*n),
		   imat= matrix(0, nvar2, nvar2),
		   loglik=double(2),
		   flag=integer(1),
		   as.double(eps),
		   deriv = matrix(double(ncol.deriv*n),ncol=ncol.deriv),
		   as.integer(doinit))

    temp <- dimnames(x)[[2]]
    if (is.null(temp)) temp <- paste("x", 1:ncol(x))
    if (scale==0)  temp <- c(temp, "log(scale)")
    dimnames(fit$imat) <- list(temp, temp)
    names(fit$coef) <- temp

    c(fit[c("iter", "coef", "imat", "loglik", "flag", "deriv")],
		list(null=nfit$loglik))
    }
