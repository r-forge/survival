# SCCS $Id: coxph.getdata.s,v 4.1 1994-04-08 15:23:38 therneau Exp $
#
# Reconstruct the Cox model data.  This is done in so many routines
#  that I extracted it out.
#
coxph.getdata <- function(fit, y=T, x=T, strata=T, offset=F) {
    ty <- fit$y
    tx <- fit$x
    strat <- fit$strata
    Terms <- fit$terms
    if (!inherits(Terms, 'terms'))
	    stop("invalid terms component of fit")
    strats <- attr(Terms, "specials")$strata

    if ( (y && is.null(ty)) || (x && is.null(tx)) ||
	     (strata && length(strats) && is.null(strat)) ||
		(offset && is.null(offset))) {
	# get the model frame
	m <- fit$model
	if (is.null(m)) m <- model.frame(fit)

	# Pull things out
	if (y && is.null(ty)) ty <- model.extract(m, 'response')

	if (offset) toff <- model.extract(m, 'response')

	# strata was saved in the fit if and only if x was
	if (x && is.null(tx)) {
	    if (length(strats)) {
		temp <- untangle.specials(Terms, 'strata', 1)
		tx <- model.matrix(Terms[-temp$terms], m)[,-1,drop=F]
		if (strata)
		    strat <- (get("strata",mode="function"))(m[temp$vars], shortlabel=T)
		}
	    else tx <- model.matrix(Terms, m)[,-1,drop=F]   #remove column of 1's though
	    }
	}

    temp <- NULL
    if (y) temp <- c(temp, "y=ty")
    if (x) temp <- c(temp, "x=tx")
    if (strata && !is.null(strat)) temp <- c(temp, "strata=strat")
    if (offset && !is.null(toff))  temp <- c(temp, "offset=toff")

    eval(parse(text=paste("list(", paste(temp, collapse=','), ")")))
    }
