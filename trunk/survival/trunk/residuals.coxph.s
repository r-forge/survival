# SCCS $Id: residuals.coxph.s,v 4.27 1997-05-08 09:03:50 therneau Exp $
residuals.coxph <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld",
			  "dfbeta", "dfbetas", "scaledsch"),
	    collapse=F, weighted=F)
    {
    type <- match.arg(type)
    otype <- type
    if (type=='dfbeta' || type=='dfbetas') type <- 'score'
    if (type=='scaledsch') type<-'schoenfeld'
    n <- length(object$residuals)
    rr <- object$residual
    y <- object$y
    x <- object$x
    vv <- object$naive.var
    if (is.null(vv)) vv <- object$var
    weights <- object$weights
    strat <- object$strata
    method <- object$method
    if (method=='exact' && (type=='score' || type=='schoenfeld'))
	stop(paste(type, 'residuals are not available for the exact method'))

    if (type == 'martingale') rr <- object$residual
    else {
	# I need Y, and perhaps the X matrix (and strata)
	Terms <- object$terms
	if (!inherits(Terms, 'terms'))
		stop("invalid terms component of object")
	strats <- attr(Terms, "specials")$strata
	if (is.null(y)  ||  (is.null(x) && type!= 'deviance')) {
	    temp <- coxph.getdata(object, y=T, x=T, strata=T)
	    y <- temp$y
	    x <- temp$x
	    if (length(strats)) strat <- temp$strata
	    }

	ny <- ncol(y)
	status <- y[,ny,drop=T]

	if (type != 'deviance') {
	    nstrat <- as.numeric(strat)
	    nvar <- ncol(x)
	    if (is.null(strat)) {
		ord <- order(y[,ny-1], -status)
		newstrat <- rep(0,n)
		}
	    else {
		ord <- order(nstrat, y[,ny-1], -status)
		newstrat <- c(diff(as.numeric(nstrat[ord]))!=0 ,1)
		}
	    newstrat[n] <- 1

	    # sort the data
	    x <- x[ord,]
	    y <- y[ord,]
	    score <- exp(object$linear.predictor)[ord]
	    if (is.null(weights)) {weights <- rep(1,n); weighted <- F}
	    else                  weights <- weights[ord]
	    }
	}

    #
    # Now I have gotton the data that I need-- do the work
    #
    if (type=='schoenfeld') {
	if (ny==2)  y <- cbind(-1,y)
	temp <- .C("coxscho", n=as.integer(n),
			    as.integer(nvar),
			    as.double(y),
			    resid= x,
			    score * weights,
			    as.integer(newstrat),
			    as.integer(method=='efron'),
			    double(3*nvar))

	deaths <- y[,3]==1

	if (nvar==1) rr <- temp$resid[deaths]
	else rr <- matrix(temp$resid[deaths,], ncol=nvar) #pick rows, and kill attr
	if (length(strats)) attr(rr, "strata")  <- table((strat[ord])[deaths])
	time <- c(y[deaths,2])  # 'c' kills all of the attributes
	if (is.matrix(rr)) dimnames(rr)<- list(time, names(object$coef))
	else               names(rr) <- time

	if (otype=='scaledsch') {
	    ndead <- sum(deaths)
	    coef <- ifelse(is.na(object$coef), 0, object$coef)
	    if (nvar==1) rr <- rr*vv *ndead + coef
	    else         rr <- rr %*%vv * ndead +
						outer(rep(1,nrow(rr)),coef)
	    }
	return(rr)
	}

    if (type=='score') {
	if (ny==2) {
	    resid <- .C("coxscore", as.integer(n),
				as.integer(nvar),
				as.double(y),
				x=x,
				as.integer(newstrat),
				score,
				weights,
				as.integer(method=='efron'),
				resid= double(n*nvar),
				double(2*nvar))$resid
	    }
	else {
	    resid<- .C("agscore",
				as.integer(n),
				as.integer(nvar),
				as.double(y),
				x,
				as.integer(newstrat),
				score,
				weights,
				as.integer(method=='efron'),
				resid=double(n*nvar),
				double(nvar*6))$resid
	    }
	if (nvar >1) {
	    rr <- matrix(0, n, nvar)
	    rr[ord,] <- matrix(resid, ncol=nvar)
	    dimnames(rr) <- list(names(object$resid), names(object$coef))
	    }
	else rr[ord] <- resid

	if      (otype=='dfbeta') {
	    if (is.matrix(rr)) rr <- rr %*% vv
	    else               rr <- rr * vv
	    }
	else if (otype=='dfbetas') {
	    if (is.matrix(rr))  rr <- (rr %*% vv) %*% diag(sqrt(1/diag(vv)))
	    else                rr <- rr * sqrt(vv)
	    }
	}

    #
    # Multiply up by case weights, if requested
    #
    if (!is.null(weights) & weighted) {
	weights[ord] <- weights
	rr <- rr * weights
	}

    #Expand out the missing values in the result
    if (!is.null(object$na.action)) {
	rr <- naresid(object$na.action, rr)
	if (is.matrix(rr)) n <- nrow(rr)
	else               n <- length(rr)
	if (type=='deviance') status <- naresid(object$na.action, status)
	}

    # Collapse if desired
    if (!missing(collapse)) {
	if (length(collapse) !=n) stop("Wrong length for 'collapse'")
	rr <- rowsum(rr, collapse)
	status <- rowsum(status, collapse)
	}

    # Deviance residuals are computed after collapsing occurs
    if (type=='deviance')
	sign(rr) *sqrt(-2* (rr+
			      ifelse(status==0, 0, status*log(status-rr))))
    else rr
    }
