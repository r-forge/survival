#SCCS $Date: 1992-03-04 16:48:14 $ $Id: predict.survreg.s,v 4.1 1992-03-04 16:48:14 therneau Exp $
#What do I need to do predictions --
#
#linear predictor:  exists
#        +se     :  X matrix
#        +newdata:  means of old X matrix, new X matrix, new offset
#
#risk -- same as lp
#
#expected --    cumulative hazard for subject= baseline haz + time + risk
#        +se :  sqrt(expected)
#      +new  :  baseline hazard function, new time, new x, means of old X,
#                        new offset, new strata
#
#terms -- : X matrix and the means
#    +se  : I matrix
#   +new  : new X matrix and the old means
predict.surv.reg <-
function(object, newdata, type=c("lp", "risk", "expected", "terms"),
		se.fit=F,
		terms=labels(object), miss.expand=T, collapse, ...)

    {
    type <- match.arg(type)
    n <- object$n
    omit <- attr(n, 'omit')


    # Toss out some simple cases
    if (!se.fit && missing(newdata)) {
	if (type=='lp') return(object$linear.predictors)
	if (type=='risk') return(exp(object$linear.predictors))
	}

    if (type=='expected') {
	if (missing(newdata)) {
	    y <- object$y
	    if (is.null(y)) {
		y <- eval(attr(Terms, "variables")[attr(Terms, "response")])
		if (!is.null(omit)) y <- y[omit,]
		}
	    pred <- y[,ncol(y)] - object$residual
	    se   <- sqrt(pred)
	    }
	else {
	    pred <- coxreg(object, newdata, se.fit=F, at.failtime=T)
	    se <- sqrt(pred)
	    }
	}

    else {
	# Ok, its going to be harder
	Terms <- object$terms
	if (!inherits(Terms, 'terms'))
	    stop("Invalid terms component for object")
	offset <- attr(Terms, 'offset')

	if (!missing(newdata)) {
	    m <- model.newframe(object, newdata)
	    x <- model.matrix(Terms, m)
	    if (!is.null(offset))  riskwt <- m[[riskwt]]
	    else riskwt <- 0
	    }
	else {
	    x <- object$x
	    if (is.null(x)) {
		x <- model.matrix(Terms, model.frame(object))
		}
	    }
	x <- sweep(x, 2, object$means)

	if (type == 'terms') {
	    attr(x, "constant") <- rep(0, ncol(x))
	    asgn <- object$assign
	    terms <- match.arg(terms, labels(object))
	    asgn <- asgn[terms]
	    if (se.fit) {
		temp <- Build.terms(x, object$coef, object$var, asgn, F)
		pred <- temp$fit
		se   <- temp$se.fit
		}
	    else pred<- Build.terms(x, object$coef, NULL, asgn, F)
	    }

	else {
	    if (length(offset)) offset <- m[[offset]]
	    else                offset <- 0
	    pred <- x %*% object$coef
	    se   <- diag(x %*% object$var %*% t(x))
	    if (type=='risk') {
		pred <- exp(pred)
		se  <- se * pred
		}
	    }
	}

    if (se.fit) se <- drop(se)
    pred <- drop(pred)
    #Expand out the missing values in the result
    if (miss.expand && !is.null(omit <- attr(object$n, 'omit'))) {
	pred <- na.expand(pred, omit)
	if(se.fit) se <- na.expand(se, omit)
	n  <- n + length(omit)
	}

    # Collapse over subjects, if requested
    if (!missing(collapse)) {
	if (length(collapse) != n) stop("Collapse vector is the wrong length")
	if (!is.matrix(pred)) pred <- tapply(pred, list(collapse), "sum")
	else
	    pred<- tapply(pred, list(collapse[row(pred)], col(pred)), sum)
	if (se.fit){
	    if (!is.matrix(pred)) se <- sqrt(tapply(se^2, list(collapse), "sum"))
	    else
		se<- sqrt(tapply(se^2, list(collapse[row(pred)], col(pred)), sum))
	    }
	}

    if (se.fit) list(fit=pred, se.fit=se)
    else pred
    }
