# SCCS $Id: predict.survreg.s,v 4.8 1998-11-30 08:30:55 therneau Exp $
# Type "p" removed till I understand it better
predict.survreg <-
    function(object, newdata, type=c("lp", 'response', 'terms', 'quantile'),
				se.fit=F,  terms=labels.lm(object),
	                        p=c(.1, .9))
    {
#
# What do I need to do predictions ?
#   
#  linear predictor: exists
#           +se    : X matrix
#          newdata : new X matrix
#
#  response -- same as lp, +transform, from distribution
#  
#  p --  density function from distribution
#          scale(s) -- if multiple I need the strata
#          +se : disallowed
#	   newdata: new X
#
    type <-match.arg(type)
    if (type=='quantile' && se.fit) 
	    stop("Standard errors not available for probability predictions")
    n <- length(object$linear.predictors)
    Terms <- object$terms
    if(!inherits(Terms, "terms"))
	    stop("invalid terms component of  object")

    strata <- attr(Terms, 'specials')$strata
    Terms <- delete.response(Terms)
    coef <- object$coefficients
    intercept <- attr(Terms, "intercept")
    

    if (missing(newdata) && (type=='terms' || se.fit)) need.x <- T
    else  need.x <- F

    if (length(strata) && type=='quantile') {
	if (is.null(object$m)) m <- model.frame(object)
	else m <- object$m
	temp <- untangle.specials(Terms, 'strata', 1)
	dropx <- temp$terms
	if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
	else strata.keep <- strata(m[,temp$vars], shortlabel=T)
	strata <- as.numeric(strata.keep)
	nstrata <- max(strata)
	    
	if (missing(newdata) && need.x){
	    x <- object$x
	    if (is.null(x)) {
		if (is.null(object$m)) 
			x <- model.matrix(Terms[-dropx], model.frame(object))
		else    x <- model.matrix(Terms[-dropx], object$m)
		}
	    }

	else if (!missing(newdata)) {
	    if (length(temp$vars)==1) newstrat <- newdata[[temp$vars]]
	    else newstrat <- strata(newdata[,temp$vars], shortlabel=T)
	    strata <- match(newstrat, strata.keep)
	    x <- model.matrix(Terms[-dropx], newdata)
	    offset <- model.extract(newdata, 'offset')
	    }
	}

    else {  # per subject strata not needed
	strata <- rep(1,n)
	if (missing(newdata) && need.x) {
	    x <- object$x
	    if (is.null(x)) {
		if (is.null(object$model)) 
			x <- model.matrix(Terms, model.frame(object))
		else    x <- model.matrix(Terms, object$model)
		}
	    }

	else if (!missing(newdata)) {
	    x <- model.matrix(Terms, newdata)
	    offset <- 0
	    strata <- rep(1, nrow(x))
	    }
	}
    scale <- object$scale[strata]
    #center x if terms are to be computed
    if(type=='p' || (type == "terms" && intercept)) 
	    x <- sweep(x, 2, object$means)

    #
    # Grab the distribution
    #
    if (is.character(object$dist)) dd <- survreg.distributions[[object$dist]]
    else dd <- object$dist
    if (is.null(dd$itrans)) {
	itrans <- dtrans <-function(x)x
	}
    else {
	itrans <- dd$itrans
	dtrans <- dd$dtrans
	}
    if (!is.null(dd$dist)) dd <- survreg.distributions[[dd$dist]]
    nvar <- length(object$coef)
    vv <- object$var[1:nvar, 1:nvar]

    #
    # Now, lay out the code one case at a time.
    #  There is some repetition this way, but otherwise the code just gets
    #    too complicated.
    #
    if (type=='lp' || type=='response') {
	if (missing(newdata)) {
 	    pred <- object$linear.predictors
#	    names(pred) <- names(object$residuals)
	    }
	else  pred <- x %*% coef  + offset
	if (se.fit) se <- sqrt(diag(x %*% vv %*% t(x)))

	if (type=='response') {
	    pred <- itrans(pred)
	    if (se.fit) se <- se/ sqrt(dtrans(pred))
	    }
	}
    else if (type=='quantile') {
	if (missing(newdata)) pred <- object$linear.predictors
	else  pred <- x %*% coef 
	# "pred" is the mean of the distribution,
	#   now add quantiles and then invert
	qq <- dd$quantile(p, dd$parm)
	if (length(qq)==1) pred <- pred + qq*scale
	else pred <- c(pred) + outer(scale, qq)
	pred <- itrans(pred)
	}

    else {  #terms
	asgn <- attr(x, 'assign')
	attr(x, 'constant') <- object$means
	terms <- match.arg(terms, labels.lm(object))
	asgn <- asgn[terms]

	if (se.fit) {
	    temp <- Build.terms(x, coef, vv, asgn, F)
	    pred <- temp$fit
	    se   <- temp$se.fit
	    }
	else pred<- Build.terms(x, coef, NULL, asgn, F)
	const <- attr(pred, 'constant')
	}

    #Expand out the missing values in the result
    # But only if operating on the original dataset
    if (missing(newdata) && !is.null(object$na.action)) {
	pred <- naresid(object$na.action, pred)
	if(se.fit) se <- naresid(object$na.action, se)

	if (type=='terms') attr(pred, 'constant') <- const
	}
    if (se.fit) list(fit=pred, se.fit=se)
    else pred
    }
