residuals.surv.reg <-
  function(object, type=c("martingale", "deviance"), collapse )
    {
    type <- match.arg(type)
    n <- length(object$residuals)
    rr <- object$residual

    if (type != 'martingale') { # I need the Y matrix
	y <- object$y
	if (is.null(y)){
	    m <- model.frame(object)
	    y <- model.extract(m, 'response')
	    }
	ny <- ncol(y)
	status <- y[,ny,drop=T]
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
	rr <- tapply(rr, list(collapse), "sum")
	if (type=='deviance' )
	    status <- tapply(status, list(collapse), 'sum')
	}

    # Deviance residuals are computed after collapsing occurs
    if (type=='deviance')
	rr <- sign(rr) *sqrt(-2* (rr+
			      ifelse(status==0, 0, status*log(status-rr))))

    rr
    }
