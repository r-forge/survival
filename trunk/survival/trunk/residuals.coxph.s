#SCCS $Id: residuals.coxph.s,v 4.13 1993-01-19 14:23:55 therneau Exp $
residuals.coxph <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld",
			  "dbeta", "dfbetas", "scaledsch"),
	    collapse=F)
    {
    type <- match.arg(type)
    otype <- type
    if (type=='dbeta' || type=='dfbetas') type <- 'score'
    if (type=='scaledsch') type<-'schoenfeld'
    n <- length(object$residuals)
    rr <- object$residual
    y <- object$y
    x <- object$x
    strat <- object$strata
    method <- object$method
    if (method=='exact' && (type=='score' || type=='schoenfeld'))
	stop(paste(type, 'residuals are not available for the exact method'))

    if (type != 'martingale') {
	# I need Y, and perhaps the X matrix, score, and strata
	Terms <- object$terms
	if (!inherits(Terms, 'terms'))
		stop("invalid terms component of object")
	strats <- attr(Terms, "specials")$strata
	if (is.null(y)  ||  (is.null(x) && type!= 'deviance')) {
	    m <- object$model
	    if (is.null(m)) m <- model.frame(object)

	    if (is.null(x) && type!= 'deviance') {
		if (length(strats)) {
		    temp <- untangle.specials(Terms, 'strata', 1)
		    x <- model.matrix(Terms[-temp$terms], m)[,-1,drop=F]
		    strat <- as.numeric(strata(m[temp$vars]))
		    }
		else x <- model.matrix(Terms, m)[,-1,drop=F]   #remove column of 1's though
		}
	    if (is.null(y)) y <- model.extract(m, 'response')
	    }

	ny <- ncol(y)
	status <- y[,ny,drop=T]

	if (type != 'deviance') {
	    nvar <- ncol(x)
	    if (is.null(strat)) {
		ord <- order(y[,ny-1], -status)
		newstrat <- rep(0,n)
		}
	    else {
		ord <- order(strat, y[,ny-1], -status)
		newstrat <- c(diff(as.numeric(strat[ord]))!=0 ,1)
		}
	    newstrat[n] <- 1

	    # sort the data
	    x <- x[ord,]
	    y <- y[ord,]
	    score <- exp(object$linear.predictor)[ord]
	    }
	}

    #
    # Now I have gotton the data that I need-- do the work
    #
    if (type=='schoenfeld') {
	if (ny==2)  y <- cbind(0,y)
	temp <- .C("coxscho", n=as.integer(n),
			    as.integer(nvar),
			    as.double(y),
			    resid= x,
			    score,
			    as.integer(newstrat),
			    as.integer(method=='efron'),
			    double(2*nvar))

	deaths <- y[,3]==1
	rr <- matrix(temp$resid[deaths,], ncol=nvar) #pick rows, and kill attr
	if (length(strats)) attr(rr, "strata")  <- strat[deaths]
	time <- c(y[deaths,2])  # 'c' kills all of the attributes
	if (is.matrix(rr)) dimnames(rr)<- list(time, names(object$coef))
	else               names(rr) <- time

	if (otype=='scaledsch') {
	    if (nvar==1) rr <- rr*object$var + object$coef
	    else         rr <- rr %*%object$var + outer(rep(1,nrow(rr)),
							     object$coef)
	    }
	return(rr)
	}

    if (type=='score') {
	if (ny==2) {
	    temp <- .C("coxres12", as.integer(n),
				as.integer(nvar),
				as.double(y),
				x=x,
				as.integer(newstrat),
				score,
				as.integer(method=='efron'),
				double(2*nvar))
	    resid <- (rep(rr[ord],nvar) * x) - temp$x
	    }
	else {
	    temp <- .C("agres12", as.integer(n),
				as.integer(nvar),
				as.double(y),
				x,
				as.integer(newstrat),
				score,
				as.integer(method=='efron'),
				resid=double(n*nvar),
				double(nvar*2))
	    resid <- temp$resid
	    }
	if (nvar >1) {
	    rr <- matrix(0, n, nvar)
	    rr[ord,] <- matrix(resid, ncol=nvar)
	    dimnames(rr) <- list(names(object$resid), names(object$coef))
	    }
	else rr[ord] <- resid
	}

    #Expand out the missing values in the result
    if (!is.null(object$na.action)) {
	rr <- naresid(object$na.action, rr)
	if (is.matrix(rr)) n <- nrow(rr)
	else               n <- length(rr)
	}

    # Collapse if desired
    if (!missing(collapse)) {
	if (length(collapse) !=n) stop("Wrong length for 'collapse'")
	if (!is.matrix(rr))  rr <- tapply(rr, list(collapse), "sum")
	else  rr<- tapply(rr, list(collapse[row(rr)], col(rr)), sum)
	}

    # Deviance residuals are computed after collapsing occurs
    if (type=='deviance')
	rr <- sign(rr) *sqrt(-2* (rr+
			      ifelse(status==0, 0, status*log(status-rr))))

    if      (otype=='dbeta') rr %*% object$var
    else if (otype=='dfbetas') (rr %*% object$var) %*% diag(sqrt(1/diag(object$var)))
    else  rr
    }
