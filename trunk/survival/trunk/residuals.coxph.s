#SCCS $Date: 1992-03-04 16:48:26 $ $Id: residuals.coxph.s,v 4.1 1992-03-04 16:48:26 therneau Exp $
residuals.coxreg <-
  function(object, type=c("martingale", "deviance", "score", "schoenfeld"),
	    miss.expand=T, collapse)
    {
    type <- match.arg(type)
    if (type=='martingale' || type=='deviance') return(NextMethod())

    n <- length(object$residuals)
    rr <- object$residual

    #Get y
    y <- object$y
    ny <- ncol(y)
    status <- y[,ny,drop=T]

    # I need the X matrix, score, and strata
    Terms <- object$terms
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of object")
    x <- object$x
    strata <- attr(Terms, "specials")$strata
    if (length(strata)) {
	m <- model.frame(object)
	if (is.null(x)) x <- model.matrix(Terms[-strata], m)
	if (length(strata)>1) stop("Only one strata() expression allowed")
	strata <- m[[(as.character(Terms))[strata]]]
	ord <- order(strata, y[,ny-1], -status)
	newstrat <- c(diff(as.numeric(strata[ord]))!=0 ,1)
	}
    else {
	if (is.null(x)) x<- model.matrix(Terms, model.frame(object))
	ord <- order(y[,ny-1], -status)
	newstrat <- rep(0,n)
	}
    newstrat[n] <- 1
    nvar <- ncol(x)

    # sort the data
    x <- x[ord,]
    y <- y[ord,]
    score <- exp(object$linear.predictor)[ord]

    if (ny ==3) subs <- paste("agres", 1:2, sep='')
    else        subs <- paste("coxres",1:2, sep='')

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
	if (length(strata)) {
	    rr <- matrix(temp$resid, ncol=nvar)[1:(temp$n),,drop=F]
	    }
	else  {
	    #put the resids in time order, rather than time within strata
	    ord2  <- order(y[indx, ny-1])
	    indx  <- indx[ord2]
	    rr   <- matrix(temp$resid, ncol=nvar)[ord2,,drop=F]
	    attr(rr, "strata")  <- strata[ord][indx]
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
    if (miss.expand && !is.null(omit <- attr(object$n, 'omit'))) {
	rr <- na.expand(rr, omit)
	n  <- n + length(omit)
	}

    # Collapse if desired
    if (!missing(collapse)) {
	if (length(collapse) !=n) stop("Wrong length for 'collapse'")
	if (ncol(rr)==1)  rr <- tapply(rr, list(collapse), "sum")
	else  rr<- tapply(rr, list(collapse[row(rr)], col(rr)), sum)
	}

    rr
    }
