#SCCS $Date: 1992-03-04 16:48:27 $ $Id: residuals.survreg.s,v 4.1 1992-03-04 16:48:27 therneau Exp $
residuals.surv.reg <-
  function(object, type=c("martingale", "deviance"),
	    miss.expand=T, collapse)
    {
    type <- match.arg(type)
    n <- length(object$residuals)
    rr <- object$residual

    if (type != 'martingale') { # I need the Y matrix
	y <- object$y
	ny <- ncol(y)
	status <- y[,ny,drop=T]
	}

    #Expand out the missing values in the result
    if (miss.expand && !is.null(omit <- attr(object$n, 'omit'))) {
	rr <- na.expand(rr, omit)
	n  <- n + length(omit)
	if (type=='deviance') status <- na.expand(status, omit)
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
