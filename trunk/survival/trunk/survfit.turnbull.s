# SCCS $Id: survfit.turnbull.s,v 1.1 2001-12-31 09:35:47 therneau Exp $
# Compute the K-M for left/right/interval censored data via Turnbull's
#      slow EM calculation
#
survfit.turnbull <- function(x, y, casewt=rep(1,n),
                       type=c('kaplan-meier', 'fleming-harrington', 'fh2'),
                       error=c('greenwood', "tsiatis"), se.fit=T,
                       conf.int= .95,
                       conf.type=c('log',  'log-log',  'plain', 'none'),
                       conf.lower=c('usual', 'peto', 'modified'),
		       start.time) {
			     
    type <- match.arg(type)
    method <- match(type, c("kaplan-meier", "fleming-harrington", "fh2"))

    error <- match.arg(error)
    error.int <- match(error, c("greenwood", "tsiatis"))
    conf.type <- match.arg(conf.type)
    conf.lower<- match.arg(conf.lower)

    if (!is.Surv(y)) stop("y must be a Surv object")
    if (!is.factor(x)) stop("x must be a factor")
    xlev <- levels(x)   # Will supply names for the curves
    x <- as.numeric(x)  # keep only the levels

    if (!missing(start.time)) { 
        # The use has requested that survival be "survival given that they
        #  made it to start.time".  We do this by just tossing those who
        #  are known to end before start.time.  Now if one of the times were
        #  interval censored (15,42) and start.time were 20, perhaps it should
        #  be modified too, but we don't.  I really don't know what the 
        #  correct action would be, actually.
        #
	ny <- ncol(y)      
	n.all <- c(table(x))   # remember the original data size
	# remove any obs whose end time is <= start.time
	keep <- (y[,ny-1] >= start.time)
	if (all(keep==F))
		stop(paste("start.time =", start.time,
			   "is greater than all time points."))
	x <- x[keep]
	y <- y[keep,,drop=F]  #make sure y remains a matrix
	casewt <- casewt[keep]
        }
    n.used <- as.vector(table(x))    # This is for the printout
    nstrat <- length(n.used)

    # Make sure that the time variable is not "counting" type,
    #  and convert "left" to "interval" style.
    stype <- attr(y, 'type')
    if (stype=='counting') 
	    stop("survfit.turnbull not appropriate for counting process data")
    if (stype=='interval') status <- y[,3]
    if (stype=='left') status <- ifelse(y[,2]==0,2,1)
    if (stype=='right')status <- y[,2]


    # the code below actually does the estimate, one curve at a time
    doit <- function(y,status, wt, ...) {
	n <- length(status)
	# Find all of the jump points for the KM in the data set, which are
	#  the exact times, plus any right-followed-by-left pairs
	# Interval censoreds count as two obs in this calculation, (right,left)
        # The variables time2, stat2, n2 are never needed after this point
	if (any(status==3)) { #interval censored
	    stat2 <- c(ifelse(status==3, 1,status), rep(2, sum(status==3)))
	    time2 <- c(y[,1], y[status==3,2])
	    }
	else {
	    stat2 <- status
	    time2 <- y[,1]
	    }
	ord <- order(time2, stat2)
	time2 <- time2[ord]
	stat2 <- stat2[ord]
	n2 <- length(time2)
	pairs <- (stat2[-n2]==0 & stat2[-1]==2)
	jtimes <- c(time2[stat2==1], .5*(time2[-n2] + time2[-1])[pairs])

	#
	# If any of the left censored times are < min(jtime), then treat
	#  them as though they were exact (for now).  The formal MLE
        #  algebra puts all their mass at an arbitray point between the
        #  smallest of such times and -infinity.  
	#
	mintime <- min(jtimes)
	who <- (status==2 & y[,1] < mintime)
	if (any(who)) {
	    status[who] <- 1
	    jtimes <- c(y[who,1], jtimes)
	    }

	# The initial "starter" KM is evenly spaced on these times
	jtimes <- sort(unique(jtimes))
	njump <- length(jtimes)
	tfit  <- list(time=jtimes, surv= 1- (1:njump)/njump)

	# The KM is computed on a fake data set with njump points
	#  standing in for the left and interval censored observations
        # So tempy contains the exact and right censored y data, followed
        #  by the fakes
	nreal <- sum(status<2)
	tempx <- factor(rep(1, njump + nreal))  #dummy x var for survfit.km
	tempy <- Surv(c(y[status<2, 1], jtimes),
		      c(status[status<2], rep(1, njump)))

	# wtmat marks, for each left/interval obs, which jump points are in it
        # A column is a "fake" time point, a row is an observation
        # For a left censored obs, we assume that the true event time is 
        #   <= the time recorded, and for an interval one that (a, b] contains
        #   the true event time.  This is motivated by data that would come
        #   from repeated visits, and agrees with Turnbull's paper.
	temp <- matrix(jtimes, nrow=sum(status>1), ncol=njump, byrow=T)
	indx <- (1:n)[status>1] #these are the interval and lc rows of the data
        temp1 <- (temp <= y[indx,1])   # logical matrix for the left censored
        temp2 <- ( temp >  y[indx,1] & temp <= y[indx,2]) # same for interval
        temp3 <- rep(status[indx]==2, njump)
        wtmat <- matrix(as.numeric((temp3&temp1) | (!temp3 & temp2)),
                        ncol=njump)

	lwt <- wt[indx]  # the input vector of case weights, for these
	eps <- 1
	old <- tfit$surv

	while (eps > .0001) {
	    # partition each left/interval person out over the jumps
	    jumps <- diff(c(1, tfit$surv[match(jtimes, tfit$time)])) #KM jumps
	    wt2 <- wtmat %*% diag(-jumps)
	    wt2 <- (lwt/(apply(wt2,1,sum))) * wt2 
	    wt2 <- apply(wt2, 2, sum)
	    tfit <- survfit.km(tempx, tempy, casewt=c(wt[status<2], wt2), ...)
	    if (F) {
                # these lines are in for debugging: change the above to 
                #  " if (T)" to turn on the printing
		cat("survival=",
		    format(round(tfit$surv[tfit$n.event>0],3)),  "\n")
		cat(" weights=", format(round(wt2,3)), "\n")
		}
            stemp <- tfit$surv[match(jtimes, tfit$time)] 
	    eps <- mean(abs(old-stemp))
	    old <- stemp
	    }	
	#
	# Now, fix up the "cheating" I did for any left censoreds which were
	#  less than the smallest jump time
	who <- (tfit$time < mintime)
	if (any(who)) {
	    indx <- min((1:length(tfit$time))[!who])  #first "real" point
	    tfit$surv[who] <- tfit$surv[indx]
	    tfit$n.event[who] <- 0
	    if (!is.null(tfit$std.err)) {
		tfit$std.err[who] <- tfit$std.err[indx]
		tfit$lower[who]   <- tfit$lower[indx]
		tfit$upper[who]   <- tfit$upper[indx]
		}
	    #remove any duplicate times
	    }
	tfit
	}		
    #
    # Now to work, one curve at a time
    #
    time   <- vector('list', nstrat)
    n.risk <- vector('list', nstrat)
    surv   <- vector('list', nstrat)
    n.cens <- vector('list', nstrat)
    n.event<- vector('list', nstrat)
    strata <- integer(nstrat)

    uniquex <- sort(unique(x))
    for (i in 1:nstrat) {
	who <- (x== uniquex[i])
	tfit <- doit(y[who,,drop=F], status[who], casewt[who])
	time[[i]]   <- tfit$time
	n.risk[[i]] <- tfit$n.risk
	surv[[i]]   <- tfit$surv
	n.cens[[i]] <- tfit$n.cens
	n.event[[i]]<- tfit$n.event
	if (i==1) {
	    if (!is.null(tfit$std.err)) {
		std.err <- vector('list', nstrat)
		conf.lower <- vector('list', nstrat)
		conf.upper <- vector('list', nstrat)
		se.fit <- T
		}
	    else se <- F
	    }
	if (se.fit) {
	    std.err[[i]]    <- tfit$std.err
	    conf.lower[[i]] <- tfit$lower
	    conf.upper[[i]] <- tfit$upper
	    }
	}

    temp <- list(n=n.used,
		 time = unlist(time),
		 n.risk = unlist(n.risk),
		 n.event= unlist(n.event),
		 n.censor = unlist(n.cens),
		 surv = unlist(surv),
		 type='right')
    
    if (nstrat >1) {
	names(strata) <- xlev[sort(unique(x))]
	temp$strata <- strata
	}

    if (se.fit) {
	temp$std.err <- unlist(std.err)
	temp$lower <- unlist(conf.lower)
	temp$upper <- unlist(conf.upper)
	temp$conf.type <- tfit$conf.type
	temp$conf.int  <- tfit$conf.int
	}
    temp
    }
