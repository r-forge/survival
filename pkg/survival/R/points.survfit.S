# $Id$
points.survfit <- function(x, xscale, xmax, fun, ...) {
    ssurv <- x$surv
    stime <- x$time


    if (is.null(x$strata)) {
	nstrat <- 1
	stemp <- rep(1, length(x$time))
	}
    else {
	nstrat <- length(x$strata)
	stemp <- rep(1:nstrat, x$strata)
	}
    if (!missing(xmax) && any(x$time>xmax)) {
	# prune back the survival curves
	# I need to replace x's over the limit with xmax, and y's over the
	#  limit with either the prior y value or firsty
	keepx <- keepy <- NULL  # lines to keep
	yzero <- NULL           # if all points on a curve are < xmax
	tempn <- table(stemp)
	offset <- cumsum(c(0, tempn))
	for (i in 1:nstrat) {
	    ttime <-stime[stemp==i]
	    if (all(ttime <= xmax)) {
		keepx <- c(keepx, 1:tempn[i] + offset[i])
		keepy <- c(keepy, 1:tempn[i] + offset[i])
		}
	    else {
		bad <- min((1:tempn[i])[ttime>xmax])
		if (bad==1)  {
		    keepy <- c(keepy, 1+offset[i])
		    yzero <- c(yzero, 1+offset[i])
		    }
		else  keepy<- c(keepy, c(1:(bad-1), bad-1) + offset[i])
		keepx <- c(keepx, (1:bad)+offset[i])
		stime[bad+offset[i]] <- xmax
		x$n.event[bad+offset[i]] <- 1   #don't plot a tick mark
		}
	    }

	# ok, now actually prune it
	stime <- stime[keepx]
	stemp <- stemp[keepx]
	x$n.event <- x$n.event[keepx]
	if (is.matrix(ssurv)) {
	    if (length(yzero)) ssurv[yzero,] <- firsty
	    ssurv <- ssurv[keepy,,drop=FALSE]
	    }
	else {
	    if (length(yzero)) ssurv[yzero] <- firsty
	    ssurv <- ssurv[keepy]
	    }
	}

    stime <- stime/xscale



    if (!missing(fun)) {
        if (is.character(fun)) {
            tfun <- switch(fun,
                            'log' = function(x) x,
                            'event'=function(x) 1-x,
                            'cumhaz'=function(x) -log(x),
                            'cloglog'=function(x) log(-log(x)),
                            'pct' = function(x) x*100,
                            'logpct'= function(x) 100*x,
                            stop("Unrecognized function argument")
                            )
            }
        else if (is.function(fun)) tfun <- fun
        else stop("Invalid 'fun' argument")

        ssurv <- tfun(ssurv)
        }


    if (!is.matrix(ssurv))
	    points(stime, ssurv, ...)
    else
	    matpoints(stime, ssurv, ...)
    }
