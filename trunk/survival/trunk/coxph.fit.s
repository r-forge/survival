#SCCS $Date: 1992-03-30 09:44:23 $ $Id: coxph.fit.s,v 4.2 1992-03-30 09:44:23 therneau Exp $
coxreg.fit <- function(x, y, strata, offset, init, iter.max,
			eps, inf.ratio, method, rownames)
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
    if (is.null(offset)) offset <- rep(0,n)

    stime <- as.double(time[sorted])
    sstat <- as.integer(status[sorted])

    if (is.null(nvar)) {
	# A special case: Null model.
	#  (This is why I need the rownames arg- can't use x' names)
	score <- exp(offset[sorted])
	coxfit <- .C("coxfit_null", as.integer(n),
				    stime,
				    sstat,
				    offset[sorted],
				    newstrat,
				    loglik=double(1),
				    cumhaz = double(n) )
	resid <- double(n)
	resid[sorted] <- sstat - score*coxfit$cumhaz
	names(resid) <- rownames

	list( loglik = coxfit$loglik,
	      linear.predictors = offset,
	      residuals = resid,
	      method= c('coxreg.null', 'coxreg') )
	}

    else {
	if (!is.null(init)) {
	    if (length(init) != nvar) stop("Wrong length for inital values")
	    }
	else init <- rep(0,nvar)

	coxfit <- .C("coxfit", iter=as.integer(iter.max),
		       as.integer(n),
		       as.integer(nvar), stime,
		       sstat,
		       x= x[sorted,] ,
		       as.double(offset[sorted] - mean(offset)),
		       newstrat,
		       means= double(nvar),
		       coef= as.double(init),
		       double(nvar),
		       rep(log(inf.ratio), nvar),
		       imat= double(nvar*nvar), loglik=double(2),
		       flag=integer(1),
		       mark = integer(n),
		       double(2*nvar*nvar + 3*nvar),
		       as.double(eps),
		       sctest=as.double(method=="efron") )


	if (coxfit$flag == 1000 && iter.max>1)
	    warning("Ran out of iterations and did not converge")
	if (coxfit$flag < 0)
	      return(paste("X matrix deemed to be singular; variable",-coxfit$flag))
	if (coxfit$flag>0 & coxfit$flag<=nvar){
	      return(paste("Variable ", coxfit$flag,"is becoming infinite:",
			    coxfit$coef[coxfit$flag]))
	      }

	names(coxfit$coef) <- dimnames(x)[[2]]
	lp <- c(x %*% coxfit$coef) + offset - sum(coxfit$coef*coxfit$means)
	score <- exp(lp[sorted])
	coxhaz <- .C("coxhaz", as.integer(n), score, coxfit$mark, newstrat,
		      hazard=double(n), cumhaz=double(n))
	resid <- double(n)
	resid[sorted] <- sstat - score*coxhaz$cumhaz
	names(resid) <- rownames

	list(coefficients  = coxfit$coef,
		    var    = matrix(coxfit$imat, ncol=nvar),
		    loglik = coxfit$loglik,
		    score  = coxfit$sctest,
		    iter   = coxfit$iter,
		    linear.predictors = as.vector(lp),
		    residuals = resid,
		    means = coxfit$means,
		    method='coxreg')
	}
    }
