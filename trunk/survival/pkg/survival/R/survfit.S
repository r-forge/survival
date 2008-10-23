# $Id$
survfit <- function (formula, data, weights, subset, na.action, ...) {
    Call <- match.call()
    # Real tricky -- find out if the first arg is "Surv(...)" without
    #  evaluating it. 
    # At one time, we allowed this to act like a formula, by adding the ~1
    #    on for the user.  Very non-standard, and a bad idea.
    # We removed all refernce to this useage from the documentation several
    #    years ago.  The line below is a temporary, "nice" failure message
    #    for this no-longer supported use.  In a couple more years, axe this
    #    section.
    #
    if ((mode(Call[[2]]) == 'call' &&  Call[[2]][[1]] == as.name('Surv'))
		|| inherits(formula, 'Surv'))  {
	stop(paste("Survfit requires a formula or a coxph fit as the",
		   "first argument, a Surv object is not supported"))
	}

    # Okay, now to the real work.
    # if the first object is not a formula, go to the next step
    if (!inherits(formula, 'formula')) {
	temp <- UseMethod("survfit")
	temp$call <- Call
	}
    else {
	# Ok, I have a formula
	# grab the data and process it
	temp <- call('model.frame', formula=formula)
	for (i in c('data', 'na.action', 'weights', 'subset')) 
		if (!is.null(Call[[i]])) temp[[i]] <- Call[[i]]
	if (is.R()) m <- eval(temp, parent.frame())
	    else        m <- eval(temp, sys.parent())

	Terms <- terms(formula, 'strata')
	ord <- attr(Terms, 'order')
	if (length(ord) & any(ord !=1))
	    stop("Interaction terms are not valid for this function")

	n <- nrow(m)
	Y <- model.extract(m, 'response')
	if (!is.Surv(Y)) stop("Response must be a survival object")

	casewt <- model.extract(m, "weights")
	# The second line below works around a bug in Splus 3.0.1, which later
	#    went away, i.e., casewt is returned as an unevaluated arg.
	if (is.null(casewt)) casewt <- rep(1,n)
	else if (mode(casewt)=='argument') casewt <- eval(casewt[[1]])

	if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

	ll <- attr(Terms, 'term.labels')
	if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
	else X <- strata(m[ll])

	if (!is.Surv(Y)) stop("y must be a Surv object")
	if (attr(Y, 'type') != 'right' && attr(Y, 'type') != 'counting')
		temp <- survfit.turnbull(X, Y, casewt, ...)
	else    temp <- survfit.km(X, Y, casewt, ...)
	if (!is.null(attr(m, 'na.action'))) 
		temp$na.action <- attr(m, 'na.action')

	temp$call <- Call
	if (is.R()) class(temp) <- 'survfit'
	else        oldClass(temp) <- 'survfit'
	}

    temp
    }

# The subscript function is bundled in here, although used most
#  often in plotting

"[.survfit" <- function(x, ..., drop=FALSE) {
    if (missing(..1)) i<- NULL  else i <- ..1
    if (missing(..2)) j<- NULL  else j <- ..2
    if (is.null(x$strata)) {
	if (is.matrix(x$surv)) {
	    x$surv <- x$surv[,i,drop=drop]
	    if (!is.null(x$std.err)) x$std.err <- x$std.err[,i,drop=drop]
	    if (!is.null(x$upper)) x$upper <- x$upper[,i,drop=drop]
	    if (!is.null(x$lower)) x$lower <- x$lower[,i,drop=drop]
	    }
	else warning("Survfit object has only a single survival curve")
	}
    else {
	if (is.null(i)) keep <- seq(along.with=x$time)
	else {
            if (is.character(i)) strat <- rep(names(x$strata), x$strata)
	    else                 strat <- rep(1:length(x$strata), x$strata)
            indx <- match(strat, i, nomatch=0) # match each row to "i"
            # the "order" part of the next line is needed if someone
            #  writes surv[2:1] to get the curves re-ordered
            # It assumes that the sort is order preserving (lines in the same
            #  strata stay in the same relative order to each other).
            keep <- (which(indx>0))[order(indx[indx>0])]
	    if (length(i) <=1) x$strata <- NULL
	    else               x$strata  <- x$strata[i]

	    x$n       <- x$n[i]
	    if (!is.null(x$n.all)) x$n.all <- x$n.all[i]
	    x$time    <- x$time[keep]
	    x$n.risk  <- x$n.risk[keep]
	    x$n.event <- x$n.event[keep]
	    x$n.censor<- x$n.censor[keep]
	    if (!is.null(x$n.enter)) x$n.enter <- x$n.enter[keep]
	    }
	if (is.matrix(x$surv)) {
	    if (is.null(j)) {
		x$surv <- x$surv[keep,,drop=drop]
		if (!is.null(x$std.err)) 
			x$std.err <- x$std.err[keep,,drop=drop]
		if (!is.null(x$upper)) x$upper <-x$upper[keep,,drop=drop]
		if (!is.null(x$lower)) x$lower <-x$lower[keep,,drop=drop]
		}
	    else {
		x$surv <- x$surv[keep,j]
		if (!is.null(x$std.err)) x$std.err <- x$std.err[keep,j]
		if (!is.null(x$upper)) x$upper <- x$upper[keep,j]
		if (!is.null(x$lower)) x$lower <- x$lower[keep,j]
		}
	    }
	else {
	    x$surv <- x$surv[keep]
	    if (!is.null(x$std.err)) x$std.err <- x$std.err[keep]
	    if (!is.null(x$upper)) x$upper <- x$upper[keep]
	    if (!is.null(x$lower)) x$lower <- x$lower[keep]
	    }
	}
    x
    }


