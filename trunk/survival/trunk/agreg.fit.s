#SCCS $Id: agreg.fit.s,v 4.18 1998-12-22 13:17:53 therneau Exp $
agreg.fit <- function(x, y, strata, offset, init, iter.max,
			eps, toler.chol, weights, method, rownames)
    {
    n <- nrow(y)
    nvar <- ncol(x)
    start <- y[,1]
    stopp <- y[,2]
    event <- y[,3]

    # Sort the data (or rather, get a list of sorted indices)
    if (length(strata)==0) {
	sorted <- order(stopp, -event)
	newstrat <- as.integer(rep(0,n))
	}
    else {
	sorted <- order(strata, stopp, -event)
	strata <- (as.numeric(strata))[sorted]
	newstrat <- as.integer(c(1*(diff(strata)!=0), 1))
	}
    if (missing(offset) || is.null(offset)) offset <- rep(0.0, n)
    if (missing(weights)|| is.null(weights))weights<- rep(1.0, n)
    else {
	if (any(weights<=0)) stop("Invalid weights, must be >0")
	weights <- weights[sorted]
	}
    sstart <- as.double(start[sorted])
    sstop <- as.double(stopp[sorted])
    sstat <- as.integer(event[sorted])

    if (is.null(nvar)) {
	# A special case: Null model.  Just return obvious stuff
	score <- as.double(exp(offset[sorted]))
	agfit <- .C("agfit_null",
		       as.integer(n),
		       as.integer(method=='efron'),
		       sstart, sstop,
		       sstat,
		       as.double(offset[sorted]),
		       as.double(weights),
		       newstrat,
		       loglik=double(1))

	agres <- .C("agmart",
		       as.integer(n),
		       as.integer(method=='efron'),
		       sstart, sstop,
		       sstat,
		       score,
		       as.double(weights),
		       newstrat,
		       resid=double(n))

	resid _ double(n)
	resid[sorted] <- agres$resid
	names(resid) <- rownames

	list(loglik=agfit$loglik,
	     linear.predictors = offset,
	     residuals = resid,
	     method= c("coxph.null", 'coxph') )
	}

    else {
	if (!is.null(init)) {
	    if (length(init) != nvar) stop("Wrong length for inital values")
	    }
	else init <- rep(0,nvar)

	agfit <- .C("agfit2", iter= as.integer(iter.max),
		       as.integer(n),
		       as.integer(nvar), sstart, sstop,
		       sstat,
		       x= x[sorted,],
		       as.double(offset[sorted] - mean(offset)),
		       as.double(weights),
		       newstrat,
		       means = double(nvar),
		       coef= as.double(init),
		       u = double(nvar),
		       imat= double(nvar*nvar), loglik=double(2),
		       flag=integer(1),
		       double(2*nvar*nvar +nvar*3 + n),
		       integer(n),
		       as.double(eps),
		       as.double(toler.chol),
		       sctest=as.double(method=='efron') )

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
			  paste((1:nvar)[(infs>eps)],collapse=","),
			  "; beta may be infinite. "))
	    }

	names(coef) <- dimnames(x)[[2]]
	lp  <- x %*% coef + offset - sum(coef *agfit$means)
	score <- as.double(exp(lp[sorted]))
	agres <- .C("agmart",
		       as.integer(n),
		       as.integer(method=='efron'),
		       sstart, sstop,
		       sstat,
		       score,
		       as.double(weights),
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
		    linear.predictors = as.vector(lp),
		    residuals = resid,
		    means = agfit$means,
		    method= 'coxph')
	}
    }
