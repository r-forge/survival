#SCCS $Id: survreg.fit.s,v 4.8 1992-12-29 15:14:57 therneau Exp $
#
# This handles the one parameter distributions-- extreme, logistic,
#       gaussian, and cauchy.
# The parent routine survreg can't allow Cauchy.  It gives negative weights
#   when we try to fit it into the IRLS model, as Cauchy is not log-covex.
#
survreg.fit<- function(x, y, offset, init, controlvals, dist, fixed,
				nullmodel=T) {

    iter.max <- controlvals$maxiter
    eps <- controlvals$rel.tol

    if (!is.matrix(x)) stop("Invalid X matrix ")
    n <- nrow(x)
    nvar <- ncol(x)
    ny <- ncol(y)
    if (is.null(offset)) offset <- rep(0,n)

    sd <- survreg.distributions[[dist]]
    if (is.null(sd)) stop ("Unrecognized distribution")
    dnum <- match(dist, c("extreme", "logistic", "gaussian", "cauchy"))
    if (is.na(dnum)) stop ("Unknown distribution given")

    ncol.deriv <- if (any(y[,ny]==3)) 6  else 3
    nvar2 <- nvar + length(sd$parms) - length(fixed)
    yy <- ifelse(y[,ny]!=3, y[,1], (y[,1]+y[,2])/2 )

    if (is.numeric(init)) {
	if (length(init) != nvar2) stop("Wrong length for initial parameters")
	eta <- x %*% init[1:nvar]
	tfix <- sd$init(yy - eta, fixed, init)
	}
    else {
	if (is.null(eta <- init$eta))  eta <- mean(yy)
	else if (length(eta) != n) stop ("Wrong length for init$eta")

	# Do the 'glim' method of finding an initial value of coef,
	tfix <- sd$init(yy - (eta+offset), init, fixed)
	deriv <- .C("survreg_g",
		       as.integer(n),
		       y,
		       as.integer(ny),
		       as.double(eta + offset),
		       coef= as.double(c(0,tfix[,1])),
		       deriv = matrix(double(n * 3),nrow=n),
		       as.integer(3),
		       as.integer(dnum))$deriv
	wt <-  -1*deriv[,3]
	coef <- solve(t(x)%*% (wt*x), c((wt*eta + deriv[,2])%*% x))
	eta <- x %*% coef  + offset
	tfix <- sd$init(yy-eta, fixed, init)
	init <- c(coef, tfix[tfix[,2]==0,1])
	}

    fit <- .C("survreg",
		   iter = as.integer(iter.max),
		   as.integer(n),
		   as.integer(nvar),
		   y,
		   as.integer(ny),
		   x,
		   as.double(offset),
		   coef= as.double(init),
		   as.integer(nrow(tfix)),
		   tfix,
		   double(3*(nvar2) + 2*n),
		   imat= matrix(0, nvar2, nvar2),
		   loglik=double(2),
		   flag=integer(1),
		   as.double(eps),
		   deriv = matrix(double(ncol.deriv*n),nrow=n),
		   as.integer(dnum))

    if (iter.max >1 && fit$flag== -1) {
	if (controlvals$failure==1)
	       warning("Ran out of iterations and did not converge")
	else if (controlvals$failure==2)
	       stop("Ran out of iterations and did not converge")
	}

    temp <- dimnames(x)[[2]]
    if (is.null(temp)) temp <- paste("x", 1:ncol(x))
    temp <- c(temp, (dimnames(tfix)[[1]])[tfix[,2]==0])
    dimnames(fit$imat) <- list(temp, temp)
    names(fit$coef) <- temp
    parms <- tfix[,1]
    parms[tfix[,2]==0] <- fit$coef[-(1:nvar)]

    if (!nullmodel)
	c(fit[c("iter", "coef", "imat", "loglik", "flag", "deriv")],
		list(parms=parms, fixed=tfix[,2]==1))
    else {
	init <- c(mean(x%*% fit$coef[1:nvar]), fit$coef[-(1:nvar)])
	temp <- cbind(parms, 1)     # "nail down" extras
	nfit <- .C("survreg",
		   iter = as.integer(iter.max),
		   as.integer(n),
		   nvar= as.integer(1),
		   y,
		   as.integer(ny),
		   x= rep(1,n),
		   as.double(offset),
		   coef= as.double(init),
		   as.integer(nrow(tfix)),
		   temp,
		   double(3*(nvar) + 2*n),
		   imat= matrix(0, nvar, nvar),
		   loglik=double(2),
		   flag=integer(1),
		   as.double(eps),
		   deriv = matrix(double(ncol.deriv*n),nrow=n),
		   as.integer(dnum))

	c(fit[c("iter", "coef", "imat", "loglik", "flag", "deriv")],
	       list(ndev=nfit$loglik, ncoef=nfit$coef, parms=parms,
		    fixed=tfix[,2]==1))
	}
    }
