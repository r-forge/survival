#SCCS $Date: 1992-08-06 16:32:16 $ $Id: survfit.s,v 4.4 1992-08-06 16:32:16 therneau Exp $
survfit <- function (formula, data, weights, subset, na.action, ...) {
    call <- match.call()
    # Real tricky -- find out if the first arg is "Surv(...)" without
    #  evaluating it.  If this is so, or it is a survival object, turn it
    #  into a formula
    if ((mode(call[[2]]) == 'call' &&  call[[2]][[1]] == as.name('Surv'))
		|| inherits(formula, 'Surv'))  {
	# The dummy function stops an annoying warning message "Looking for
	#  'formula' of mode function, ignored one of mode ..."
	xx <- function(x) formula(x)
	formula <- xx(paste(deparse(call[[2]]), 1, sep="~"))
	}

    if (!inherits(formula, 'formula')) temp <- UseMethod("survfit")
    else {
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
	casewt <- model.extract(m, "weights")
	# The second line below works around a bug in Splus 3.0.1, which later
	#    went away, i.e., casewt is returned as an unevaluated arg.
	if (is.null(casewt)) casewt <- rep(1,n)
	else if (mode(casewt)=='argument') casewt <- eval(casewt[[1]])

	ll <- attr(Terms, 'term.labels')
	if (length(ll) == 0) X <- factor(rep(1,n))
	else {
	    temp <-  rep(1, length(ll))
	    strat <- untangle.specials(Terms, 'strata',1)$terms
	    if (length(strat)) temp[strat] <- 0
	    lname <- ifelse(temp==1, paste(ll,'=', sep=''), "")
	    temp <- paste("'",lname, "', ","as.character(m[['", ll, "']])", sep='')
	    temp <- paste(temp, collapse=", ', ',")
	    temp <- paste("paste(", temp, ",sep='')")
	    X <- factor(eval(parse(text=temp)))
	    }

	temp <- survfit.km(X, Y, casewt, ...)
	attr(temp, "class") <- "survfit"
	if (!is.null(attr(m, 'na.action'))) temp$na.action <- attr(m, 'na.action')
	}
    temp$call <- call
    temp
    }
