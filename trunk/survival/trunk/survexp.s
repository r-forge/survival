#SCCS  $Id: survexp.s,v 4.14 1994-01-04 14:56:58 therneau Exp $
survexp <- function(formula=formula(data), data=sys.parent(),
	weights, subset, na.action,
	times,  cohort=T,  conditional=T,
	ratetable=survexp.us, scale=365.25,
	model=F, x=F, y=F) {

    call <- match.call()
    m <- match.call(expand=F)
    m$ratetable <- m$model <- m$x <- m$y <- m$scale<- m$cohort <- NULL
    m$times <- m$conditional <- NULL

    Terms <- if(missing(data)) terms(formula, 'ratetable')
	     else              terms(formula, 'ratetable',data=data)
    if (any(attr(Terms, 'order') >1))
	    stop("Survexp cannot have interaction terms")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    n <- nrow(m)

    Y <- model.extract(m, 'response')
    no.Y <- is.null(Y)
    if (!no.Y) {
	if (is.matrix(Y) && ncol(Y)>1) {
	    if (cohort || ncol(Y)!=2)
		stop("Illegal response value, for cohort survival")
	    if (Y[,1] >= Y[,2]) stop("Stop time must be > start time")
	    if (Y[,1] <0) stop("Start time must be >=0")
	    }
	else if (any(Y <0)) stop ("Negative follow up time")
	}
    else {
	conditional <- F
	if (missing(times))
	   stop("There is no times argument, and no follow-up times are given in the formula")
	Y <- rep(max(times), n)
	}
    if (!missing(times)) newtime <- sort(unique(c(Y, times)))
    else                 newtime <- sort(unique(Y))
    weights <- model.extract(m, 'weights')

    rate <- attr(Terms, "specials")$ratetable
    if (length(rate)==0)
	stop("Must have a ratetable() call in the formula")
    if (length(rate) >1 )
	stop ("Can have only 1 ratetable() call in a formula")

    ovars <- (dimnames(attr(Terms, 'factors'))[[1]])[-c(1, rate)]
    rtemp <- match.ratetable(m[,rate], ratetable)
    R <- rtemp$R
    if (!is.null(rtemp$call)) {  #need to dop some dimensions from ratetable
	ratetable <- eval(parse(text=rtemp$call))
	}

    if (cohort) {
	# Now process the other (non-ratetable) variables
	if (length(ovars)==0)  X <- rep(1,n)  #no categories
	else {
	    odim <- length(ovars)
	    for (i in 1:odim) {
		temp <- m[[ovars[i]]]
		ctemp <- class(temp)
		if (!is.null(ctemp) && ctemp=='tcut')
		    stop("Can't use tcut variables in expected survival")
		}
	    X <- strata(m[ovars])
	    }
	temp <- survexp.fit(cbind(as.numeric(X),R), Y, newtime,
			       conditional, ratetable)

	if (missing(times)) {
	    n.risk <- temp$n
	    if (is.matrix(temp$surv)) surv <- apply(temp$surv, 2,cumprod)
	    else                      surv <- cumprod(temp$surv)
	    }
	else {
	    keep <- match(times, newtime)
	    if (is.matrix(temp$surv)) {
		surv <- apply(temp$surv, 2, cumprod)
		surv <- surv[keep,]
		n.risk <- temp$n[keep,]
		}
	    else {
		surv <- cumprod(temp$surv)
		surv <- surv[keep]
		n.risk <- temp$n[keep]
		}
	    newtime <- times
	    }
	newtime <- newtime/scale
	if (length(ovars)) {    #matrix output
	    if (no.Y) { #use matrix formulation
		out <- list(call=call, surv=surv, n.risk=c(n.risk[,1]),
			    time=newtime)
		}
	    else {
		#I can't use the "matrix" form of surv, due to different n's
		out <- list(call=call, surv=c(surv), n.risk=c(n.risk),
				time = rep(newtime, ncol(surv)))
		tstrat <- rep(nrow(surv), ncol(surv))
		names(tstrat) <- levels(X)
		out$strata <- tstrat
		}
	    }
	else  out <- list(call=call, surv=c(surv), n.risk=c(n.risk),
			   time=newtime)

	na.action <- attr(m, "na.action")
	if (length(na.action))  out$na.action <- na.action
	if (model) out$model <- m
	else {
	    if (x) out$x <- structure(cbind(X, R),
		dimnames=list(row.names(m), c("group", dimid)))
	    if (y) out$y <- Y
	    }
	out$summ <- rtemp$summ
	class(out) <- c("survexp", "survfit")
	out
	}

    else { #individual survival
	temp <- survexp.fit(cbind(1:n,R), Y, max(Y), conditional, ratetable)
	xx <- temp$surv
	names(xx) <- row.names(m)
	na.action <- attr(m, "na.action")
	if (length(na.action)) naresid(na.action, xx)
	else xx
	}
    }
