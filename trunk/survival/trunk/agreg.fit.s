#SCCS $Date: 1992-08-06 16:49:25 $ $Id: agreg.fit.s,v 4.5 1992-08-06 16:49:25 therneau Exp $
agreg.fit <- function(x, y, strata, offset, init, iter.max,
			eps, method, rownames)
    {
    n <- nrow(y)
    nvar <- ncol(x)
    start <- y[,1]
    stopp <- y[,2]
    event <- y[,3]

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
	# A special case: Null model.  Just return obvious stuff
	score <- as.double(exp(offset[sorted]))
	agfit <- .C("agfit_null",
		       as.integer(n),
		       sstart, sstop,
		       sstat,
		       offset[sorted],
		       newstrat,
		       loglik=double(1),
		       hazard=double(n), cumhaz=double(n))

	aghaz <- .C("aghaz",
		       as.integer(n),
		       sstart, sstop,
		       sstat,
		       score,
		       newstrat,
		       hazard=double(n), cumhaz=double(n))

	resid _ double(n)
	resid[sorted] <- sstat - score*aghaz$cumhaz
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

	agfit <- .C("agfit", iter= as.integer(iter.max),
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
		       double(2*nvar*nvar +nvar*3 + n),
		       as.double(eps),
		       sctest=as.double(method=='efron') )

	infs <- abs(agfit$u %*% matrix(agfit$imat,nvar))
	if (iter.max >1) {
	    if (agfit$flag == 1000)
		   warning("Ran out of iterations and did not converge")
	    else if (any((infs > eps) & (infs > eps*abs(agfit$coef))))
		warning(paste("Loglik converged before variable ",
			  (1:nvar)[(infs>eps)], ", beta may be infinite. ",
			   collapse=''))
	    }
	if (agfit$flag < 0)
	      return(paste("X matrix deemed to be singular; variable",
			-agfit$flag, 'at iteration', agfit$iter))

	names(agfit$coef) <- dimnames(x)[[2]]
	lp  <- x %*% agfit$coef + offset - sum(agfit$coef *agfit$means)
	score <- as.double(exp(lp[sorted]))
	aghaz <- .C("aghaz",
		       as.integer(n),
		       sstart, sstop,
		       sstat,
		       score,
		       newstrat,
		       hazard=double(n), cumhaz=double(n))

	resid _ double(n)
	resid[sorted] <- sstat - score*aghaz$cumhaz
	names(resid) <- rownames

	list(coefficients  = agfit$coef,
		    var    = matrix(agfit$imat, ncol=nvar),
		    loglik = agfit$loglik,
		    score  = agfit$sctest,
		    iter   = agfit$iter,
		    linear.predictors = as.vector(lp),
		    residuals = resid,
		    means = agfit$means,
		    method= 'coxph')
	}
    }
