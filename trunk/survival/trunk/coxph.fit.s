#SCCS $Date: 1993-06-17 12:26:16 $ $Id: coxph.fit.s,v 4.12 1993-06-17 12:26:16 therneau Exp $
coxph.fit <- function(x, y, strata, offset, init, iter.max,
			eps, weights, method, rownames)
    {
    n <-  nrow(y)
    nvar <- ncol(x)
    time <- y[,1]
    status <- y[,2]

    # Sort the data (or rather, get a list of sorted indices)
    if (is.null(strata)) {
	sorted <- order(time)
	newstrat <- as.integer(rep(0,n))
	}
    else {
	sorted <- order(strata, time)
	strata <- (as.numeric(strata))[sorted]
	newstrat <- as.integer(c(1*(diff(strata)!=0), 1))
	}
    if (missing(offset) || is.null(offset)) offset <- rep(0,n)
    if (missing(weights)|| is.null(weights))weights<- rep(1,n)
    else {
	if (any(weights<=0)) stop("Invalid weights, must be >0")
	weights <- weights[sorted]
	}
    stime <- as.double(time[sorted])
    sstat <- as.integer(status[sorted])

    if (is.null(nvar)) {
	# A special case: Null model.
	#  (This is why I need the rownames arg- can't use x' names)
	score <- exp(offset[sorted])
	coxfit <- .C("coxfit_null", as.integer(n),
				    as.integer(method=='efron'),
				    stime,
				    sstat,
				    exp(offset[sorted]),
				    weights,
				    newstrat,
				    loglik=double(1),
				    resid = double(n) )
	resid <- double(n)
	resid[sorted] <- coxfit$resid
	names(resid) <- rownames

	list( loglik = coxfit$loglik,
	      linear.predictors = offset,
	      residuals = resid,
	      method= c('coxph.null', 'coxph') )
	}

    else {
	if (!missing(init) && !is.null(init)) {
	    if (length(init) != nvar) stop("Wrong length for inital values")
	    }
	else init <- rep(0,nvar)

	coxfit <- .C("coxfit2", iter=as.integer(iter.max),
		       as.integer(n),
		       as.integer(nvar), stime,
		       sstat,
		       x= x[sorted,] ,
		       as.double(offset[sorted] - mean(offset)),
		       weights,
		       newstrat,
		       means= double(nvar),
		       coef= as.double(init),
		       u = double(nvar),
		       imat= double(nvar*nvar), loglik=double(2),
		       flag=integer(1),
		       double(2*n + 2*nvar*nvar + 3*nvar),
		       as.double(eps),
		       sctest=as.double(method=="efron") )

	var <- matrix(coxfit$imat,nvar,nvar)
	coef <- coxfit$coef
	if (coxfit$flag < nvar) which.sing <- diag(var)==0
	else which.sing <- rep(F,nvar)

	infs <- abs(coxfit$u %*% var)
	if (iter.max >1) {
	    if (coxfit$flag == 1000)
		   warning("Ran out of iterations and did not converge")
	    else if (any((infs > eps) & (infs > sqrt(eps)*abs(coef))))
		warning(paste("Loglik converged before variable ",
			  paste((1:nvar)[(infs>eps)]),
			  ", beta may be infinite. "))
	    }

	names(coef) <- dimnames(x)[[2]]
	lp <- c(x %*% coef) + offset - sum(coef*coxfit$means)
	score <- exp(lp[sorted])
	coxres <- .C("coxmart", as.integer(n),
				as.integer(method=='efron'),
				stime,
				sstat,
				newstrat,
				score,
				weights,
				resid=double(n))
	resid <- double(n)
	resid[sorted] <- coxres$resid
	names(resid) <- rownames
	coef[which.sing] <- NA

	list(coefficients  = coef,
		    var    = var,
		    loglik = coxfit$loglik,
		    score  = coxfit$sctest,
		    iter   = coxfit$iter,
		    linear.predictors = as.vector(lp),
		    residuals = resid,
		    means = coxfit$means,
		    method='coxph')
	}
    }
