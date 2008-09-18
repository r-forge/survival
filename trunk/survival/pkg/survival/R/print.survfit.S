# $Id$
print.survfit <- function(x, scale=1, 
			  digits = max(options()$digits - 4, 3), ...) {

    if (!is.null(cl<- x$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
    }	
    omit <- x$na.action
    if (length(omit)) cat("  ", naprint(omit), "\n")

    savedig <- options(digits=digits)
    on.exit(options(savedig))

    print(survmean(x, scale=scale))
    x
    }


#
# The function that does all of the actual work --create a matrix
#
survmean <- function(x, scale=1) {

    # The starting point for the integration of the AUC
    if (!is.null(x$start.time)) start.time <- x$start.time
    else                        start.time <- min(0, x$time)

    #
    # The function below is called once for each line of output,
    #  i.e., once per curve.  It creates the line of output
    #
    pfun <- function(nused, time, surv, n.risk, n.event, lower, upper, 
		      start.time) {
	# compute the mean of the curve, with "start.time" as 0
	#   start by drawing rectangles under the curve
        # Lining up the terms for "varmean" is tricky -- the easiest 
        #   check is to look at the homework solution on page 195-196
        #   of Miller, Survival Analysis, Wiley, 1981.
	n <- length(time)
	delta <- diff(c(start.time, time))     #width of rectangles
	rectangles <- delta * c(1, surv[-n])   #area of rectangles
	hh <- ifelse((n.risk-n.event)==0, 0, 
		       n.event /(n.risk *(n.risk -n.event)))
	varmean <- sum( cumsum(rev(rectangles[-1]))^2 * rev(hh)[-1])
	mean <- sum(rectangles) + start.time

        #compute the median  and ci(median)
	minmin <- function(y, xx) {
	     if (any(!is.na(y) & y==.5)) {	
		if (any(!is.na(y) & y <.5))
 		  .5*(min(xx[!is.na(y) & y==.5]) + min(xx[!is.na(y) & y<.5]))
		else
		  .5*(min(xx[!is.na(y) & y==.5]) + max(xx[!is.na(y) & y==.5]))
	        }
	     else   min(xx[!is.na(y) & y<=.5])
	     }

	med <- minmin(surv, time)
	if (!is.null(upper)) {
	    upper <- minmin(upper, time)
	    lower <- minmin(lower, time)
	    c(nused, max(n.risk), n.risk[1], 
              sum(n.event), sum(mean), sqrt(varmean), med, lower, upper)
	    }
	else
		c(nused, max(n.risk), n.risk[1], sum(n.event), 
                  sum(mean), sqrt(varmean), med)
	}

    stime <- x$time/scale
    surv <- x$surv
    plab <- c("n.obs", "n.max", "n.first", "events", 
              "mean", "se(mean)", "median")  #col labels
    ncols <- 7     #number of columns in the output
    if (!is.null(x$conf.int)) {
	plab <- c(plab, paste(x$conf.int, c("LCL", "UCL"), sep=''))
	ncols <- 9
	}

    #Four cases: strata Y/N  by  ncol(surv)>1 Y/N
    #  Repeat the code, with minor variations, for each one
    if (is.null(x$strata)) {
	if (is.matrix(surv)) {
	    out <- matrix(0, ncol(surv), ncols)
	    for (i in 1:ncol(surv)) {
		if (is.null(x$conf.int))
		     out[i,] <- pfun(x$n, stime, surv[,i], x$n.risk, x$n.event,
				      NULL, NULL, start.time)
		else out[i,] <- pfun(x$n, stime, surv[,i], x$n.risk, x$n.event,
				    x$lower[,i], x$upper[,i], start.time)
		}
	    dimnames(out) <- list(dimnames(surv)[[2]], plab)
	    }
	else {
	    out <- pfun(x$n, stime, surv, x$n.risk, x$n.event, x$lower, 
			x$upper, start.time)
	    names(out) <- plab
 	    }
        }
    else {   #strata case
	nstrat <- length(x$strata)
	stemp <- rep(1:nstrat,x$strata)  # the index vector for strata1, 2, etc

	if (is.matrix(surv)) {
	    ns <- ncol(surv)
	    out <- matrix(0, nstrat*ns, ncols)
            if (is.null(dimnames(surv)[[2]]))
                dimnames(out) <- list(rep(names(x$strata), rep(ns,nstrat)), 
                                      plab)
            else {
                cname <- outer(dimnames(surv)[[2]], names(x$strata), paste,
                               sep=", ")
                dimnames(out) <- list(c(cname), plab)
                }
	    k <- 0
	    for (i in 1:nstrat) {
		who <- (stemp==i)
		for (j in 1:ns) {
		    k <- k+1
		    if (is.null(x$lower))
		         out[k,] <- pfun(x$n[i], stime[who], surv[who,j],
					 x$n.risk[who], x$n.event[who],
					 NULL, NULL, start.time)
		    else out[k,] <- pfun(x$n[i], stime[who], surv[who,j],
					 x$n.risk[who], x$n.event[who],
					 x$lower[who,j], x$upper[who,j], 
					 start.time)
		    }
		}
	    }
	else {
	    out <- matrix(0, nstrat, ncols)
	    dimnames(out) <- list(names(x$strata), plab)
	    for (i in 1:nstrat) {
		who <- (stemp==i)
		if (is.null(x$lower))
		     out[i,] <- pfun(x$n[i], stime[who], surv[who], 
				     x$n.risk[who], x$n.event[who], 
				     NULL, NULL, start.time)
		else out[i,] <- pfun(x$n[i], stime[who], surv[who], 
				     x$n.risk[who], x$n.event[who], 
				     x$lower[who], x$upper[who], start.time)
		}
	    }
	}
    out
    }


