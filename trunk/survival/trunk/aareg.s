# SCCS $Id: aareg.s,v 1.1 2001-04-17 08:23:40 therneau Exp $
# Aalen's additive regression model
#
aareg <- function(formula, data=sys.parent(), ..., qrtol=1e-7, nmin) {
    call <- match.call()

    # Use a Null coxph call, followed by coxph.detail, to add up
    #   all of the things we really need
    tfit <- coxph(formula, data, iter=0, method='breslow', x=T, ...)
    temp <- attr(tfit$terms, 'specials')
    if (!is.null(temp) && !is.null(temp$strata))
	    stop("Model cannot contain a strata() term")
    if (!is.null(temp) && !is.null(temp$cluster))
	    stop("Model cannot contain a cluster() term")
    
    dt   <- coxph.detail(tfit)
    nvar  <- ncol(tfit$x)	 # number of covariates
    if (!is.null(dt$strata))
	    stop("Model cannot contain a strata() or cluster() term")

    if (missing(nmin)) nmin <- 3*nvar
    n     <- length(dt$time[dt$nrisk> nmin]) # number of unique death times
    if (n<=1) stop("The threshold 'nmin' is too high, no model can be fit")

    # If (start, stop] data was used, time/status will be cols 2 and 3 of y,
    #   otherwise cols 1 and 2.  In both cases, these are the last 2 cols
    times <- tfit$y[,ncol(tfit$y)-1]
    status<- tfit$y[,ncol(tfit$y)]

    # This matches the death times in the data set (unsorted) to the
    #  sorted list of unique death times.  "0" = not a death 
    index <- match(times, dt$time, nomatch=0) * status

    # Increment is the step in Aalen's plots
    # Twt is the set of weights that he uses in his test

    if (ncol(tfit$x)==1)  { # special case of only 1 covariate
	xx <- c(tfit$x)
	xx <- tapply(xx, index, mean) 
	if (any(index==0)) xx <- xx[-1]
	xx <- xx - dt$means
	keep <- (dt$imat > 0 & dt$nrisk > nmin)
	twt  <- dt$nevent[keep]/dt$imat[keep]  # V^{-1}
	increment <- xx[keep] *twt * dt$nevent[keep]/ dt$nrisk[keep]
	n2 <- sum(keep)
	b0 <- dt$nevent[1:n2]/dt$nrisk[1:n2] - dt$means[1:n2] * increment
	}
    else {
	increment <- matrix(0,n, nvar)
	twt <- increment
	n2 <- 0   #actual number of elements used
	for (i in 1:n) {
	    # Get the x-vectors of the deaths, and subtract the mean at 
	    #    this time
	    #  If there are multiple deaths, add up the (x - mean) vectors
	    xx <- tfit$x[index==i,]
	    if (is.matrix(xx)){
		xx <- t(xx) - as.vector(dt$means[i,])
		xx <- xx %*% rep(1, ncol(xx))
		}
	    else xx <- xx - dt$means[i,]

	    #solve, and check for singularity
	    qri <- qr(dt$imat[,,i], tol=qrtol)
	    if (qri$rank < nvar) break
	    increment[i,] <- qr.coef(qri, xx)/dt$nrisk[i]	
	    twt[i,] <- 1/diag(qr.coef(qri, diag(nvar)))
	    n2 <- i    
	    }
	if (n2 < n) increment <- increment[1:n2,]
	b0 <- dt$nevent[1:n2]/dt$nrisk[1:n2] - 
		 apply(dt$means[1:n2,]*increment, 1, sum)
	}

    increment <- cbind(b0,increment)
    dimnames(increment) <- list(dt$time[1:n2], 
				c("Intercept", names(tfit$coef)))

    ans <- list(ntime=c(n2,length(dt$time)), times=dt$time[1:n2], 
		  nrisk=dt$nrisk[1:n2], 
		nevent=dt$nevent[1:n2], increment=increment, 
		tweight = cbind(1,twt), call=call) 
    oldClass(ans) <- 'aareg'
    ans
    }

"[.aareg" <- function(x, ..., drop=F) {
    if (!inherits(x, 'aareg')) stop ("Must be an aareg object")

    # There is a bug in Splus6, the "=F" on drop is ignored as a
    #   default.  Add temporary hack
    drop <- F

    i <- ..1
    if (is.matrix(x$increment)) {
	x$increment <- x$increment[,i, drop=drop]
	x$tweight   <- x$tweight[,i,drop=drop]
	}
    else stop("Subsripting impossible, increment component not a matrix")

    x
    }

