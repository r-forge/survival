#SCCS $Date: 1992-07-14 00:01:23 $ $Id: predict.coxph.s,v 4.5 1992-07-14 00:01:23 therneau Exp $
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
#    +se  :  ""  + I matrix
#   +new  : new X matrix and the old means + I matrix
predict.coxph <-
function(object, newdata, type=c("lp", "risk", "expected", "terms"),
		se.fit=F,
		terms=labels.lm(object), collapse, safe=F, ...)

    {
    type <-match.arg(type)
    n <- object$n
    Terms <- object$terms
    strata <- attr(Terms, 'specials')$strata
    if (length(strata)) {
	   temp <- untangle.specials(Terms, 'strata', 1)
	   Terms2 <- Terms[-temp$terms]
	   }
    else  Terms2 <- Terms
    offset <- attr(Terms, "offset")
    resp <- attr(Terms, "variables")[attr(Terms, "response")]

    if (missing(newdata)) {
	if (type=='terms' || (se.fit && (type=='lp' || type=='risk'))) {
	    x <- object$x
	    if (is.null(x)) {
		x <- model.matrix(Terms2, model.frame(object))[,-1,drop=F]
		}
	    x <- sweep(x, 2, object$means)
	    }
	else if (type=='expected') {
	    y <- object$y
	    if (missing(y)) {
		m <- model.frame(object)
		y <- model.extract(m, 'response')
		}
	    }
	}
    else {
	if (type=='expected')
	     m <- model.newframe(Terms, newdata, response=T)
	else m <- model.newframe(Terms2, newdata)

	x <- model.matrix(Terms2, m)[,-1,drop=F]
	x <- sweep(x, 2, object$means)
	if (length(offset)) {
	    if (type=='expected') offset <- as.numeric(m[[offset]])
	    else {
		offset <- attr(Terms2, 'offset')
		offset <- as.numeric(m[[offset]])
		}
	    }
	else offset <- 0
	}

    #
    # Now, lay out the code one case at a time.
    #  There is some repetition this way, but otherwise the code just gets
    #    too complicated.
    if (type=='lp' || type=='risk') {
	if (missing(newdata)) {
	    pred <- object$linear.predictors
	    names(pred) <- names(object$residuals)
	    }
	else                  pred <- x %*% object$coef  + offset
	if (se.fit) se <- sqrt(diag(x %*% object$var %*% t(x)))

	if (type=='risk') {
	    pred <- exp(pred)
	    if (se.fit) se <- se * sqrt(pred)
	    }
	}

    else if (type=='expected') {
	if (missing(newdata)) pred <- y[,ncol(y)] - object$residual
	else  stop("Method not yet finished")
	se   <- sqrt(pred)
	}

    else {  #terms
	attr(x, "constant") <- rep(0, ncol(x))
	asgn <- object$assign
	terms <- match.arg(Terms2, labels.lm(object))
	asgn <- asgn[terms]
	if (se.fit) {
	    temp <- Build.terms(x, object$coef, object$var, asgn, F)
	    pred <- temp$fit
	    se   <- temp$se.fit
	    }
	else pred<- Build.terms(x, object$coef, NULL, asgn, F)
	}

    if (se.fit) se <- drop(se)
    pred <- drop(pred)
    #Expand out the missing values in the result
    if (!is.null(object$na.action)) {
	pred <- naresid(object$na.action, pred)
	if (is.matrix(pred)) n <- nrow(pred)
	else               n <- length(pred)
	if(se.fit) se <- naresid(object$na.action, se)
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
