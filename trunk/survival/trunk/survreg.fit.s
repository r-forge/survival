#SCCS $Date: 1992-06-26 23:27:46 $ $Id: survreg.fit.s,v 4.1 1992-06-26 23:27:46 therneau Exp $
#
# This handles the one parameter distributions-- extreme, logistic, and
#       gaussian
survreg.fit<- function(x, y, offset, init, iter.max,
			eps, fun, scale, singular.ok=F)
    {
    debug<-1
    if (!is.matrix(x)) stop("Invalid formula for fitting function")
    n <- nrow(x)
    nvar <- ncol(x)
    ny <- ncol(y)
    if (is.null(offset)) offset <- rep(0,n)
    if (is.null(scale)) {
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
	if (scale >0 && init[nvar2] <=0)
		stop("Invalid initial value for sigma")
	}

    #
    # Fit the null model --- only an intercept and sigma
    #
    nfit <- CC(fun, iter= as.integer(iter.max),
		   as.integer(n),
		   as.integer(1),
		   y,
		   as.integer(ny),
		   rep(1,n),
		   as.double(offset),
		   coef= as.double(c(0,1)),
		   as.double(scale),
		   double(6 + 2*n),
		   imat= matrix(0, 2,2),
		   loglik=double(2),
		   flag=integer(1),
		   as.double(eps),
		   deriv = matrix(double(3*n),ncol=3),
		   as.integer(debug),
		   as.integer(1))
    # If necessary, get starting values for the regression
    #   It seems to be more stable if I leave sigma out of the iteration
    #
    if (missing(init) || is.null(init) ) {
	if (scale==0) sigma <- nfit$coef[2]
	else          sigma <- scale
	fit <- .C(fun, iter=as.integer(1),
		       as.integer(n),
		       as.integer(nvar),
		       y,
		       as.integer(ny),
		       x,
		       as.double(offset),
		       coef= as.double(rep(0,nvar)),
		       as.double(sigma),
		       double(3*(nvar) + 2*n),
		       imat= matrix(0, nvar, nvar),
		       loglik=double(2),
		       flag=integer(1),
		       as.double(eps),
		       deriv = matrix(double(3*n),ncol=3),
		       as.integer(debug),
		       as.integer(0))
#       xfit <- .C(fun, iter=as.integer(0),
#                      as.integer(n),
#                      as.integer(nvar),
#                      y,
#                      as.integer(ny),
#                      x,
#                      as.double(offset),
#                      coef= as.double(rep(0,nvar)),
#                      as.double(1),
#                      double(3*(nvar) + 2*n),
#                      imat= matrix(0, nvar, nvar),
#                      loglik=double(2),
#                      flag=integer(1),
#                      as.double(eps),
#                      deriv = matrix(double(3*n),ncol=3),
#                      as.integer(debug))
#       xlm <- lm.wfit(x, xfit$deriv[,2]/xfit$deriv[,3], -xfit$deriv[,3], "qr")
#       yfit <- .C(fun, iter=as.integer(0),
#                      as.integer(n),
#                      as.integer(nvar),
#                      y,
#                      as.integer(ny),
#                      x,
#                      as.double(offset),
#                      coef= as.double(rep(0,nvar)),
#                      as.double(sigma),
#                      double(3*(nvar) + 2*n),
#                      imat= matrix(0, nvar, nvar),
#                      loglik=double(2),
#                      flag=integer(1),
#                      as.double(eps),
#                      deriv = matrix(double(3*n),ncol=3),
#                      as.integer(debug))
#       ylm <- lm.wfit(x, yfit$deriv[,2]/yfit$deriv[,3], -yfit$deriv[,3], "qr")
#
#       browser()
	if (scale==0) {
	    #Estimate sigma from the residuals  - I'm a little fuzzy on this
	    init <- c(fit$coef, sigma)
	    }
	else init <- fit$coef
	}

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
		   deriv = matrix(double(3*n),ncol=3),
		   as.integer(debug),
		   as.integer(1))
browser()
    c(fit[c("iter", "coef", "imat", "loglik", "flag", "deriv")],
		list(null=nfit$loglik))
    }
