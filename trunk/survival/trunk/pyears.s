#SCCS  $Id: pyears.s,v 4.1 1993-12-02 21:37:07 therneau Exp $
pyears <- function(formula=formula(data), data=sys.parent(),
	weights, subset, na.action,
	ratetable=survexp.uswhite, scale=365.25,
	model=F, x=F, y=F) {

    call <- match.call()
    m <- match.call(expand=F)
    m$ratetable <- m$model <- m$x <- m$y <- m$scale<- NULL

    Terms <- if(missing(data)) terms(formula, 'ratetable')
	     else              terms(formula, 'ratetable',data=data)
    if (any(attr(Terms, 'order') >1))
	    stop("Pyears cannot have interaction terms")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    Y <- model.extract(m, 'response')
    if (is.null(Y)) stop ("Follow-up time must appear in the formula")
    if (!is.Surv(Y)){
	if (any(Y <0)) stop ("Negative follow up time")
	Y <- as.matrix(Y)
	if (ncol(Y) >2) stop("Y has too many columns")
	if (ncol(Y)==2 && any(Y[,2] <= Y[,1]))
	    stop("stop time must be > start time")
	}
    n <- nrow(Y)

    weights <- model.extract(m, 'weights')

    rate <- attr(Terms, "specials")$ratetable
    if (length(rate) >1 )
	stop ("Can have only 1 ratetable() call in a formula")
    else if (length(rate)==1) {
	ovars <- (dimnames(attr(Terms, 'factors'))[[1]])[-c(1, rate)]
	temp <- match.ratetable(m[,rate], ratetable)
	R <- temp$R
	if (!is.null(temp$call)) {  #need to dop some dimensions from ratetable
	    ratetable <- eval(parse(text=temp$call))
	    }
	}
    else {
	ovars <- (dimnames(attr(Terms, 'factors'))[[1]])[-1]
	}

    # Now process the other (non-ratetable) variables
    if (length(ovars)==0)  {
	# no categories!
	X <- rep(1,n)
	ofac <- odim <- odims <- ocut <- 1
	}
    else {
	odim <- length(ovars)
	ocut <- NULL
	odims <- ofac <- double(odim)
	X <- matrix(0, n, odim)
	outdname <- vector("list", odim)
	for (i in 1:odim) {
	    temp <- m[[ovars[i]]]
	    ctemp <- class(temp)
	    if (!is.null(ctemp) && ctemp=='tcut') {
		X[,i] <- temp
		temp2 <- attr(temp, 'cutpoints')
		odims[i] <- length(temp2) -1
		ocut <- c(ocut, temp2)
		ofac[i] <- 0
		outdname[[i]] <- attr(temp, 'labels')
		}
	    else {
		temp2 <- factor(temp)
		X[,i] <- temp2
		temp3 <- levels(temp2)
		odims[i] <- length(temp3)
		ofac[i] <- 1
		outdname[[i]] <- temp3
		}
	    }
	}

    # Now do the computations
    ocut <-c(ocut,0)   #just in case it were of length 0
    osize <- prod(odims)
    if (length(rate)) {  #include expected
	temp <- attributes(ratetable)
	temp <- .C("pyears1",
			as.integer(n),
			as.integer(ncol(Y)),
			as.integer(is.Surv(Y)),
			as.double(Y),
			as.integer(length(temp$dim)),
			as.integer(temp$factors),
			as.integer(temp$dim),
			as.double(unlist(temp$cutpoints)),
			ratetable,
			as.double(R),
			as.integer(odim),
			as.integer(ofac),
			as.integer(odims),
			as.double(ocut),
			X,
			pyears=double(osize),
			pn    =double(osize),
			pcount=double(if(is.Surv(Y)) osize else 1),
			pexpect=double(osize),
			offtable=double(1))[16:20]
	}
    else {
	temp <- .C('pyears2',
			as.integer(n),
			as.integer(ncol(Y)),
			as.integer(is.Surv(Y)),
			as.double(Y),
			as.integer(odim),
			as.integer(ofac),
			as.integer(odims),
			as.double(ocut),
			X,
			pyears=double(osize),
			pn    =double(osize),
			pcount=double(if(is.Surv(Y)) osize else 1),
			offtable=double(1)) [10:13]
	}

    out <- list(call = call,
		pyears= array(temp$pyears/scale, dim=odims, dimnames=outdname),
		n     = array(temp$pn,     dim=odims, dimnames=outdname),
		offtable = temp$offtable/scale)
    if (length(rate)) {
	out$expected <- array(temp$pexpect, dim=odims, dimnames=outdname)
	out$unitcheck <- (attr(ratetable, "Rfun"))(R)
	}
    if (is.Surv(Y))
	out$event    <- array(temp$pcount, dim=odims, dimnames=outdname)
    na.action <- attr(m, "na.action")
    if (length(na.action))  out$na.action <- na.action
    if (model) out$model <- m
    else {
	if (x) out$x <- cbind(X, R)
	if (y) out$y <- Y
	}
    out
    }
