#SCCS $Date: 2000-03-02 21:48:10 $ $Id: survfit.s,v 4.15 2000-03-02 21:48:10 boos Exp $
survfit <- function(formula, data, ...) {
    call <- match.call()  #make the "call" returned object correct
    # Real tricky -- for backwards compatability I want to allow the
    #   expression "survfit(Surv(time, status), data=...", that is, no ~1
    #   in the formula.  In order to make the data= part of this work, I have
    #   to detect this without evaluating the first arg, i.e., I can't 
    #   ask for it's class, even indirectly through UseMethod(). 
    # In order to evaluate class(some typed expression) S has to evaluate the
    #   expression.  This is ok for a formula, because one can go from
    #   character string to formula without reference to the data set itself.
    #   Not so if the string is "Surv(x,y)".
    # This trick is not perfect.  If a data set "mine" contains a variable
    #   "x" which is already a Surv object, there is no way to find that out
    #    without fetching "x" so survfit(x, data=mine) will fail at the
    #    "inherits" line below when it can't find 'x'.
    # If the first arg is "Surv", then force a call to survfit.Surv,
    #   otherwise use the class of the arg as is usual with methods.

    ## other code that used to be here is now in:
    ##      survfit.formula.s
    if ((mode(call[[2]]) == 'call' && call[[2]][[1]] == as.name('Surv'))
	    || inherits(formula, 'Surv')) {
	# The dummy function stops an annoying warning message "Looking for
	#  'formula' of mode function, ignored one of mode ..."
	xx <- function(x) formula(x)
	formula <- xx(paste(deparse(call[[2]]), 1, sep="~"))
        }

    UseMethod("survfit", formula, data=data, ..., call=call)
    }

# The subscript function is bundled in here, although used most
#  often in plotting

"[.survfit" <- function(fit, ..., drop=F) {
    if (missing(..1)) i<- NULL else i <- ..1
    if (missing(..2)) j<- NULL else j <- ..2
    if (is.null(fit$strata)) {
	if (is.matrix(fit$surv)) {
	    fit$surv <- fit$surv[,i,drop=drop]
	    if (!is.null(fit$std.err)) 
		    fit$std.err <- fit$std.err[,i,drop=drop]
	    if (!is.null(fit$upper)) fit$upper <- fit$upper[,i,drop=drop]
	    if (!is.null(fit$lower)) fit$lower <- fit$lower[,i,drop=drop]
	    } 
	else warning("Survfit object has only a single survival curve")
        } 
    else {
	if (is.null(i)) keep <- seq(along=fit$time)
	else {
	    if (is.character(i)) 
		    strat <- rep(names(fit$strata), fit$ntimes.strata)
	    else
		    strat <- rep(1:length(fit$strata), fit$ntimes.strata)
	    keep <- seq(along=strat)[match(strat, i, nomatch=0)>0]
	    if (length(i) <=1) fit$strata <- NULL
	    else 	       fit$strata <- fit$strata[i] 
	    fit$strata.all <- fit$strata.all[i]
	    fit$ntimes.strata <- fit$ntimes.strata[i]
	    fit$time    <- fit$time[keep]
	    fit$n.risk  <- fit$n.risk[keep]
	    fit$n.event <- fit$n.event[keep]
	    }
	if (is.matrix(fit$surv)) {
	    if (is.null(j)) {
		fit$surv <- fit$surv[keep,,drop=drop]
		if (!is.null(fit$std.err)) 
			fit$std.err <- fit$std.err[keep,,drop=drop]
		if (!is.null(fit$upper))
			fit$upper <- fit$upper[keep,,drop=drop]
		if (!is.null(fit$lower)) 
			fit$lower <- fit$lower[keep,,drop=drop]
	        }
	    else {
		fit$surv <- fit$surv[keep,j]
		if (!is.null(fit$std.err)) fit$std.err <- fit$std.err[keep,j]
		if (!is.null(fit$upper)) fit$upper <- fit$upper[keep,j]
		if (!is.null(fit$lower)) fit$lower <- fit$lower[keep,j]
	        }
	    }
	else { # added information for start and stop time data
	    fit$surv <- fit$surv[keep]
	    if (!is.null(fit$enter)) fit$enter <- fit$enter[keep]
	    if (!is.null(fit$exit.censored)) 
		    fit$exit.censored <- fit$exit.censored[keep]
	    if (!is.null(fit$std.err)) fit$std.err <- fit$std.err[keep]
	    if (!is.null(fit$upper)) fit$upper <- fit$upper[keep]
	    if (!is.null(fit$lower)) fit$lower <- fit$lower[keep]
	    }
        }
    fit
    }




