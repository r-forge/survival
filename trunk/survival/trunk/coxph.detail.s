#SCCS  $Id: coxph.detail.s,v 4.5 1993-06-17 16:03:04 therneau Exp $
coxph.detail <-  function(object) {
    method <- object$method
    if (method!='breslow' && method!='efron')
	stop(paste("Detailed output is not available for the", method,
			"method"))
    n <- length(object$residuals)
    rr <- object$residual
    y <- object$y
    x <- object$x
    weights <- object$weights
    strat <- object$strata
    Terms <- object$terms
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of object")
    strats <- attr(Terms, "specials")$strata

    if (is.null(y)  ||  is.null(x)) {
	m <- object$model
	if (is.null(m)) m <- model.frame(object)

	if (is.null(x) ) {
	    if (length(strats)) {
		temp <- untangle.specials(Terms, 'strata', 1)
		x <- model.matrix(Terms[-temp$terms], m)[,-1,drop=F]
		strat <- as.numeric(strata(m[temp$vars]))
		}
	    else x <- model.matrix(Terms, m)[,-1,drop=F]   #remove column of 1's though
	    }
	if (is.null(y)) y <- model.extract(m, 'response')
	weights <- model.extract(m, 'weights')
	}

    nvar <- ncol(x)
    if (ncol(y)==2) y <- cbind(0,y)
    if (is.null(strat)) {
	ord <- order(y[,2], -y[,3])
	newstrat <- rep(0,n)
	}
    else {
	ord <- order(strat, y[,2], -y[,3])
	newstrat <- c(diff(as.numeric(strat[ord]))!=0 ,1)
	}
    newstrat[n] <- 1

    # sort the data
    x <- x[ord,]
    y <- y[ord,]
    storage.mode(y) <- 'double'
    score <- exp(object$linear.predictor)[ord]
    if (is.null(weights)) weights <- rep(1,n)
    else                  weights <- weights[ord]

    ndeath <- sum(y[,3])
    ff <- .C("coxdetail", as.integer(n),
			  as.integer(nvar),
			  ndeath= as.integer(ndeath),
			  y = y,
			  as.double(x),
			  as.integer(newstrat),
			  index =as.double(score),
			  as.double(weights),
			  means= c(method=='efron', double(ndeath*nvar)),
			  u = double(ndeath*nvar),
			  i = double(ndeath*nvar*nvar),
			  double(nvar*(3 + 2*nvar)) )
    keep <- 1:ff$ndeath
    vname<- dimnames(x)[[2]]
    time <- y[ff$index[keep],2]
    names(time) <- NULL
    means<- (matrix(ff$means,ndeath, nvar))[keep,]
    score<-  matrix(ff$u, ndeath, nvar)[keep,]
    var <- array(ff$i, c(nvar, nvar, ndeath))[,,keep]
    if (nvar>1) {
	dimnames(means) <- list(time, vname)
	dimnames(score) <- list(time, vname)
	dimnames(var) <- list(vname, vname, time)
	}
    else {
	names(means) <- time
	names(score) <- time
	names(var) <- time
	}

    dimnames(ff$y) <- NULL
    temp <- list(time = time, means=means, nevent=ff$y[keep,1],
	 nrisk = ff$y[keep,2], hazard= ff$y[keep,3], score= score,  imat=var,
	 y=y, x=x)
    if (length(strats)) temp$strata <- strat[keep]
    if (!all(weights==1)) temp$weights <- weights
    temp
    }
