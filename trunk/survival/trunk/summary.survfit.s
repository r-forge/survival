#SCCS $Id: summary.survfit.s,v 5.9 2001-12-31 14:57:52 therneau Exp $
#
# Version with no C code, using approx() to do the subscript
#  calculations

summary.survfit <- function(object, times, censored=F, 
			    scale=1, extend=F, ...) {
    fit <- object
    if (!inherits(fit, 'survfit'))
	    stop("Invalid data")

    table <- survmean(fit, scale=scale)  #for inclusion in the output list

    # The fit$surv object is sometimes a vector and sometimes a matrix.
    #  Make a copy of it that is always a matrix, to simplify the number of
    #  cases for our subscripting work below.  At the end of the routine
    #  we'll turn it back into a vector if needed.  Similar treatment is
    #  given to the standard error and confidence limits.
    surv <- as.matrix(fit$surv)
    if (is.null(fit$strata)) {
	nstrat <- 1
	stemp <- rep(1, length(surv))
        strata.names <- ""
	}
    else   {
	nstrat <- length(fit$strata)
	stemp <- rep(1:nstrat, fit$strata)
        strata.names <- names(fit$strata)
	}

    if (is.null(fit$std.err)) std.err <- NULL
    else 		      std.err <- fit$std.err * surv

    if (!is.null(fit$lower)) {
	lower <- as.matrix(fit$lower)
	upper <- as.matrix(fit$upper)
        }

    if (missing(times)) {
	# just pick off the appropriate rows of the output
	if (censored) {
	    # No subscripting at all
	    times <- fit$time
	    n.risk<- fit$n.risk
	    n.event <- fit$n.event
	    n.enter <- fit$n.enter
	    n.censor  <- fit$n.censor
	    strata <- factor(stemp, labels=strata.names)
	    }
	else {
	    # select off rows with at least one event
	    who <- (fit$n.event > 0)
	    times <- fit$time[who]
	    n.risk <- fit$n.risk[who]
	    n.event <- fit$n.event[who]
	    n.enter <- fit$n.enter[who]
	    n.censor <- fit$n.censor[who]
	    surv <- surv[who,,drop=F]
	    if (!is.null(std.err)) std.err <- std.err[who,,drop=F]
	    if (!is.null(fit$lower)) {
		lower <- lower[who,,drop=F]
		upper <- upper[who,,drop=F]
	        }
	    strata <- (factor(stemp, labels=strata.names))[who, drop=T]
	    }
        }

    else {  
	#this case is harder, since it involves "in between" points
	times <- sort(times)   #just in case the user didn't

	# The one line function below might be opaque (even to me) --
	# For n.event, we want to know the number since the last chosen
	#  printout time point.  Start with the curve of cumulative
	#  events at c(0, stime) (the input time points), which is
	#  the cumsum below; pluck off the values corresponding to our
	#  time points, the [x] below; then get the difference since the
	#  last chosen time point (or from 0, for the first chosen point).
	cfun <- function(x, data) diff(c(0, cumsum(c(0,data))[x]))

	# Now to work
	# The basic idea is to process the curves one at a time,
	#   adding the results for that curve onto a list, so the
	#   number of events will be n.enter[[1]], n.enter[[2]], etc.
	# For the survival, stderr, and confidence limits it suffices
	#   to create a single list 'indx1' containing a subscripting vector
	indx1 <- n.risk <- n.event <- newtimes <- vector('list', nstrat)
	if (!is.null(fit$n.enter))  n.enter <- vector('list', nstrat)
	if (!is.null(fit$n.censor)) n.censor<- vector('list', nstrat)
	n <- length(stemp)
	for (i in 1:nstrat) {
	    who <- (1:n)[stemp==i]
	    stime <- fit$time[who]

	    # First, toss any printing times that are outside our range
	    if (is.null(fit$start.time)) mintime <- min(stime, 0)
	    else                         mintime <- fit$start.time
	    ptimes <- times[times >= mintime]

	    if (!extend) {
		maxtime <- max(stime)
		ptimes <- ptimes[ptimes <= maxtime]
		}

	    newtimes[[i]] <- ptimes

	    # If we tack a 0 onto the front of the vector of survival
	    #  times, then indx1 is the subscript for that vector
	    #  corresponding to the list of "ptimes".  If the input
	    #  data had stime=c(10,20) and ptimes was c(5,10,15,20),
	    #  the result would be 1,2,2,3.
	    # For n.risk we want a slightly different index: 2,2,2,3.
	    #  The number at risk at time 15 = number at risk at time 10,
	    #  but the number at risk before time 15 is not 0 (like the
	    #  survival), it is the number at risk at time 10 - #entered at 10
	    #
	    temp1 <- approx(c(mintime, stime), 0:length(stime), xout=ptimes,
			    method='constant', f=0, rule=2)$y
	    indx1[[i]] <- ifelse(temp1==0, 1, 1+ who[pmax(1,temp1)])
	    n.event[[i]] <- cfun(temp1+1, fit$n.event[who])

	    if (!is.null(fit$n.censor))
		    n.censor[[i]] <- cfun(temp1+1, fit$n.censor[who])
	    if (is.null(fit$n.enter)) temp2 <- 0
	    else  {
		temp2 <- fit$n.enter[who[1]]
		n.enter[[i]] <- cfun(temp1+1, fit$n.enter[who])
		}

	    n.risk[[i]] <- ifelse(temp1==0, fit$n.risk[who[1]] - temp2,
				            fit$n.risk[who[pmax(1,temp1)]])
            # Why not just "who[temp1]" instead of who[pmax(1,temp1)] in the
            #  line just above?  When temp1 has zeros, the first expression
            #  gives a vector that is shorter than temp1, and the ifelse
            #  doesn't work right due to mismatched lengths.  
	    }

	# Now create the output list
	times  <- unlist(newtimes)
	n.risk <-  unlist(n.risk)
	n.event <- unlist(n.event)
	if (!is.null(fit$n.enter))  n.enter <- unlist(n.enter)
	if (!is.null(fit$n.censor)) n.censor<- unlist(n.censor)
	indx1 <- unlist(indx1)
	surv <- (rbind(1.,surv))[indx1,,drop=F]
	if (!is.null(std.err)) std.err <- rbind(0.,std.err)[indx1,,drop=F]
	if (!is.null(fit$lower)) {
	    lower <- rbind(1.,lower)[indx1,,drop=F]
	    upper <- rbind(1.,upper)[indx1,,drop=F]
	    }
	if (!is.null(fit$strata)) {
	    scount <- unlist(lapply(newtimes, length))
	    strata <- factor(rep(1:nstrat, scount), labels=names(fit$strata))
	    }
	}

    #
    # Final part of the routine: paste the material together into
    #  the correct output structure
    #
    if (fit$type == 'right' || inherits(fit, 'survfit.cox')) 
	    temp <- list(surv=surv, time=times, n.risk=n.risk, n.event=n.event,
			 conf.int=fit$conf.int, type=fit$type, table=table)
    
    if (fit$type == 'counting')
	    temp <- list(surv=surv, time=times, n.risk=n.risk, n.event=n.event,
			 conf.int=fit$conf.int, n.enter=n.enter,
			 n.censor=n.censor, type=fit$type, 
			 table=table)

    if (!is.null(fit$start.time)) temp$start.time <- fit$start.time

    if (ncol(surv)==1) {
	# Make surve & etc vectors again
	temp$surv <- drop(temp$surv)
	if (!is.null(std.err)) temp$std.err <- drop(std.err)
	if (!is.null(fit$lower)) {
	    temp$lower <- drop(lower)
	    temp$upper <- drop(upper)
	    }
	}
    else {
	if (!is.null(std.err)) temp$std.err <- std.err
	if (!is.null(fit$lower)) {
	    temp$lower <- lower
	    temp$upper <- upper
	    }
	}

    if (!is.null(fit$strata)) {
	temp$strata <- strata
	}
    temp$call <- fit$call
    if (!is.null(fit$na.action)) temp$na.action <- fit$na.action
  
    oldClass(temp) <- 'summary.survfit'
    temp
    }









