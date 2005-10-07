#SCCS $Id: survfit.s,v 4.21 2005-10-07 22:48:20 lunde Exp $
survfit <- function (formula, data, weights, subset, na.action, ...) {
    call <- match.call()
    # Real tricky -- find out if the first arg is "Surv(...)" without
    #  evaluating it.  If this is so, or it is a survival object, turn it
    #  into a formula
    # (This allows people to leave off the "~1" from a formula.  This
    # is an option I allowed early on, now wish I never had done it,
    # and have removed all references to this option from all documentation.
    # I a couple of more years I'll axe it.)
    if ((mode(call[[2]]) == 'call' &&  call[[2]][[1]] == as.name('Surv'))
		|| inherits(formula, 'Surv'))  {
	# The dummy function stops an annoying warning message "Looking for
	#  'formula' of mode function, ignored one of mode ..."
	xx <- function(x) formula(x)
	formula <- xx(paste(deparse(call[[2]]), 1, sep="~"))
	}

    # if the first object is a Cox model, call survfit.coxph
    if (!inherits(formula, 'formula')) temp <- UseMethod("survfit")
    else {
	# Ok, I have a formula
        # grab the data and process it
	m <- match.call(expand=F)
	m$... <- NULL

	Terms <- terms(formula, 'strata')
	ord <- attr(Terms, 'order')
	if (length(ord) & any(ord !=1))
	    stop("Interaction terms are not valid for this function")
	m$formula <- Terms
	m[[1]] <- as.name("model.frame")
	m <- eval(m, sys.parent())

	n <- nrow(m)
	Y <- model.extract(m, response)
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
	oldClass(temp) <- "survfit"
	if (!is.null(attr(m, 'na.action'))) 
		temp$na.action <- attr(m, 'na.action')
	}
    temp$call <- call
    # added to change INFs to NAs for C program - cmb 8/25/2000
    if (any(is.inf(temp$std.err))) temp$std.err[is.inf(temp$std.err)] <- NA
    temp
    }

# The subscript function is bundled in here, although used most
#  often in plotting

"[.survfit" <- function(fit, ..., drop=F) {
    if (missing(..1)) i<- NULL  else i <- ..1
    if (missing(..2)) j<- NULL  else j <- ..2
    if (is.null(fit$strata)) {
	if (is.matrix(fit$surv)) {
	    fit$surv <- fit$surv[,i,drop=drop]
	    if (!is.null(fit$std.err)) fit$std.err <- fit$std.err[,i,drop=drop]
	    if (!is.null(fit$upper)) fit$upper <- fit$upper[,i,drop=drop]
	    if (!is.null(fit$lower)) fit$lower <- fit$lower[,i,drop=drop]
	    }
	else warning("Survfit object has only a single survival curve")
	}
    else {
	if (is.null(i)) keep <- seq(along=fit$time)
	else {
            if (is.character(i)) strat <- rep(names(fit$strata), fit$strata)
	    else                 strat <- rep(1:length(fit$strata), fit$strata)
	    keep <- seq(along=strat)[match(strat, i, nomatch=0)>0]
	    if (length(i) <=1) fit$strata <- NULL
	    else               fit$strata  <- fit$strata[i]

	    fit$n       <- fit$n[i]
	    if (!is.null(fit$n.all)) fit$n.all <- fit$n.all[i]
	    fit$time    <- fit$time[keep]
	    fit$n.risk  <- fit$n.risk[keep]
	    fit$n.event <- fit$n.event[keep]
	    fit$n.censor<- fit$n.censor[keep]
	    if (!is.null(fit$n.enter)) fit$n.enter <- fit$n.enter[keep]
	    }
	if (is.matrix(fit$surv)) {
	    if (is.null(j)) {
		fit$surv <- fit$surv[keep,,drop=drop]
		if (!is.null(fit$std.err)) 
			fit$std.err <- fit$std.err[keep,,drop=drop]
		if (!is.null(fit$upper)) fit$upper <-fit$upper[keep,,drop=drop]
		if (!is.null(fit$lower)) fit$lower <-fit$lower[keep,,drop=drop]
		}
	    else {
		fit$surv <- fit$surv[keep,j]
		if (!is.null(fit$std.err)) fit$std.err <- fit$std.err[keep,j]
		if (!is.null(fit$upper)) fit$upper <- fit$upper[keep,j]
		if (!is.null(fit$lower)) fit$lower <- fit$lower[keep,j]
		}
	    }
	else {
	    fit$surv <- fit$surv[keep]
	    if (!is.null(fit$std.err)) fit$std.err <- fit$std.err[keep]
	    if (!is.null(fit$upper)) fit$upper <- fit$upper[keep]
	    if (!is.null(fit$lower)) fit$lower <- fit$lower[keep]
	    }
	}
    fit
    }


