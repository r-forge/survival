#SCCS 8/6/92 @(#)agexact.fit.s	4.7
agexact.fit <- function(x, y, strata, offset, iter.max,
			eps, weights, init, method, rownames)
    {
    if (!is.matrix(x)) stop("Invalid formula for cox fitting function")
    if (!is.null(weights) && any(weights!=1))
	  stop("Case weights are not supported for the exact method")
    n <- nrow(x)
    nvar <- ncol(x)
    if (ncol(y)==3) {
	start <- y[,1]
	stopp <- y[,2]
	event <- y[,3]
	}
    else {
	start <- rep(0,n)
	stopp <- y[,1]
	event <- y[,2]
	}

    # Sort the data (or rather, get a list of sorted indices)
    if (is.null(strata)) {
	sorted <- order(stopp, -event)
	newstrat <- as.integer(rep(0,n))
	}
    else {
	sorted <- order(strata, stopp, -event)
	strata <- (as.numeric(strata))[sorted]
	newstrat <- as.integer(c(1*(diff(strata)!=0), 1))
	}
    if (is.null(offset)) offset <- rep(0,n)

    sstart <- as.double(start[sorted])
    sstop <- as.double(stopp[sorted])
    sstat <- as.integer(event[sorted])

    if (is.null(nvar)) {
	# A special case: Null model.  Not worth coding up
	stop("Cannot handle a null model + exact calculation (yet)")
	}

    if (!is.null(init)) {
	if (length(init) != nvar) stop("Wrong length for inital values")
	}
    else init <- rep(0,nvar)

    agfit <- .C("agexact", iter= as.integer(iter.max),
		   as.integer(n),
		   as.integer(nvar), sstart, sstop,
		   sstat,
		   x= x[sorted,],
		   as.double(offset[sorted] - mean(offset)),
		   newstrat,
		   means = double(nvar),
		   coef= as.double(init),
		   u = double(nvar),
		   imat= double(nvar*nvar), loglik=double(2),
		   flag=integer(1),
		   double(2*nvar*nvar +nvar*4 + n),
		   integer(2*n),
		   as.double(eps),
		   sctest=double(1) )

    var <- matrix(agfit$imat,nvar,nvar)
    coef <- agfit$coef
    if (agfit$flag < nvar) which.sing <- diag(var)==0
    else which.sing <- rep(F,nvar)

    infs <- abs(agfit$u %*% var)
    if (iter.max >1) {
	if (agfit$flag == 1000)
	       warning("Ran out of iterations and did not converge")
	else if (any((infs > eps) & (infs > sqrt(eps)*abs(coef))))
	    warning(paste("Loglik converged before variable ",
		      paste((1:nvar)[(infs>eps)]),
		      ", beta may be infinite. "))
	}

    names(coef) <- dimnames(x)[[2]]
    lp  <- x %*% coef + offset - sum(coef *agfit$means)
    score <- as.double(exp(lp[sorted]))
    agres <- .C("agmart",
		   as.integer(n),
		   as.integer(0),
		   sstart, sstop,
		   sstat,
		   score,
		   rep(1,n),
		   newstrat,
		   resid=double(n))
    resid _ double(n)
    resid[sorted] <- agres$resid
    names(resid) <- rownames
    coef[which.sing] <- NA

    list(coefficients  = coef,
		var    = var,
		loglik = agfit$loglik,
		score  = agfit$sctest,
		iter   = agfit$iter,
		linear.predictors = lp,
		residuals = resid,
		means = agfit$means,
		method= 'coxph')
    }
