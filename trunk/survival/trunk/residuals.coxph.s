#SCCS $Id: residuals.coxph.s,v 4.5 1992-07-14 00:01:24 therneau Exp $
residuals.coxph <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld"),
	    collapse)
    {
    type <- match.arg(type)
    n <- length(object$residuals)
    rr <- object$residual
    y <- object$y
    x <- object$x
    strat <- object$strata

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

	    if (ny ==3) subs <- paste("agres", 1:2, sep='')
	    else        subs <- paste("coxres",1:2, sep='')
	    }
	}

    #
    # Now I have gotton the data that I need-- do the work
    #
    if (type=='schoenfeld') {
	temp <- .C(subs[2], n=as.integer(n),
			    as.integer(nvar),
			    indx = as.double(y),
			    x,
			    as.integer(newstrat),
			    score,
			    resid=double(n*nvar),
			    double(n*nvar))
	if (ny==2) indx <- temp$indx[1:(temp$n)] #unique death times
	else       indx <- temp$indx[(n+1):(n+temp$n)]  #col 2 of y
	if (length(strats)) {
	    rr <- matrix(temp$resid, ncol=nvar)[1:(temp$n),,drop=F]
	    }
	else  {
	    #put the resids in time order, rather than time within strata
	    ord2  <- order(y[indx, ny-1])
	    indx  <- indx[ord2]
	    rr   <- matrix(temp$resid, ncol=nvar)[ord2,,drop=F]
	    attr(rr, "strata")  <- strats[ord][indx]
	    }
	time <- c(y[indx, ny-1])  # 'c' kills all of the attributes
	dimnames(rr)<- list(time, names(object$coef))
	return(drop(rr))
	}

    if (type=='score') {
	# I need the hazard
	if (ny==2) { #I can get the hazard from the residuals
	    cumhaz <- (status[ord] - rr[ord])/ score
	    temp <- c(1, newstrat[-n])  #marks first obs of new strata
	    hazard <- ifelse(temp==1, cumhaz, diff(c(0, cumhaz)))
	    }
	else {  # I need to call a routine
	    aghaz <- .C("aghaz",
			   as.integer(n),
			   as.double(y[,1]),
			   as.double(y[,2]),
			   as.integer(y[,3]),
			   score,
			   as.integer(newstrat),
			   hazard=double(n), cumhaz=double(n))
	    hazard <- aghaz$hazard
	    cumhaz <- aghaz$cumhaz
	    }
	temp <- .C(subs[1], as.integer(n),
			    as.integer(nvar),
			    as.double(y),
			    x,
			    as.integer(newstrat),
			    score,
			    hazard,
			    cumhaz,
			    resid=double(n*nvar),
			    wmean= double(n*nvar))
	if (nvar >1) {
	    rr <- matrix(0, n, nvar)
	    rr[ord,] <- matrix(temp$resid, ncol=nvar)
	    dimnames(rr) <- list(names(object$resid), names(object$coef))
	    }
	else rr[ord] <- temp$resid
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
	if (ncol(rr)==1)  rr <- tapply(rr, list(collapse), "sum")
	else  rr<- tapply(rr, list(collapse[row(rr)], col(rr)), sum)
	}

    # Deviance residuals are computed after collapsing occurs
    if (type=='deviance')
	rr <- sign(rr) *sqrt(-2* (rr+
			      ifelse(status==0, 0, status*log(status-rr))))

    rr
    }
