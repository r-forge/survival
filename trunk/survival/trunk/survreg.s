#SCCS $Id: survreg.s,v 4.9 1992-08-06 16:33:01 therneau Exp $
survreg <- function(formula=formula(data), data=sys.parent(),
	subset, na.action,
	link=c('log', 'identity'),
	dist=c("extreme", "logistic", "gaussian", "exponential"),
	scale,
	eps=.0001, init, iter.max=10,
	model=F, x=F, y=F, ...) {

    call <- match.call()
    m <- match.call(expand=F)
    m$dist <- m$link <- m$model <- m$x <- m$y <- m$... <-  NULL
    m$eps <- m$init <- m$iter.max <- NULL
    m$scale <- NULL

    Terms <- if(missing(data)) terms(formula, 'strata')
	     else              terms(formula, 'strata',data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    dist <- match.arg(dist)
    link <- match.arg(link)
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    X <- model.matrix(Terms, m)
    n <- nrow(X)
    nvar <- ncol(X)
    offset<- attr(Terms, "offset")
    if (!is.null(offset)) offset <- as.numeric(m[[offset]])
    else                  offset <- rep(0, n)

    type <- attr(Y, "type")
    linkfun <- glm.links["link", link][[1]]
    if (type== 'counting') stop ("Invalid survival type")
    else if (type=='interval') {
	if (any(Y[,3]==3))
	     Y <- cbind(linkfun(Y[,1:2]), Y[,3])
	else Y <- cbind(linkfun(Y[,1]), Y[,3])
	}
    else if (type=='left')
	     Y <- cbind(linkfun(Y[,1]), 2-Y[,2])
    else     Y <- cbind(linkfun(Y[,1]), Y[,2])

    if (missing(init)) init <- NULL
    if (missing(scale)) scale <- NULL

    if( dist=="extreme")
	sfit <- survreg.fit(X, Y, offset, init=init, iter.max=iter.max,
			    eps = eps, fun='exvalue', scale=scale)
    else if (dist=="logistic")
	sfit <- survreg.fit(X, Y, offset, init=init, iter.max=iter.max,
			    eps = eps, fun='logisticfit', scale=scale)
    else if (dist=="exponential") {
	scale <- 1
	sfit <- survreg.fit(X, Y, offset, init=init, iter.max=iter.max,
			    eps = eps, fun='exvalue', scale=1)
	}
    else stop(paste ("Unknown distribution", dist))

    if (is.character(sfit))  fit <- list(fail=sfit)  #error message
    else {
	# There may be more clever ways to do this, but ....
	#  In order to make it look like IRLS, and get all the components
	#  that I need for glm inheritance, do one step of weighted least
	#  squares.
	eta <- c(X %*% sfit$coef[1:nvar]) + offset
	wt<- -sfit$deriv[,3]
	fit <- lm.wfit(X, eta + sfit$deriv[,2]/wt, wt, "qr", ...)

	if (link=='log') fit$fitted.values <- exp(fit$fitted.values)
	fit$family <- c(name=dist, link=link, "")
	fit$linear.predictors <- eta
	fit$deviance <- sfit$loglik[1]
	fit$null.deviance <- sfit$null[1]
	fit$iter <- sfit$iter

	# If singular.ok=T, there may be NA coefs.  The var matrix should
	#   be an inversion of the "non NA" portion.
	if (is.null(scale)) {
	    scale <- exp(as.vector(sfit$coef[nvar+1]))
	    good <- c(!is.na(fit$coef),T)
	    var <- matrix(0, nvar+1, nvar+1)
	    }
	else {
	    good <- !is.na(fit$coef)
	    var <- matrix(0, nvar, nvar)
	    }
	var[good,good] <- solve(sfit$imat[good,good])
	fit$scale<- scale
	fit$var <- var

	fit$flag <- sfit$flag
	fit$dresiduals <- sign(fit$residuals)*sqrt(sfit$deriv[,1])
	}

    na.action <- attr(m, "na.action")
    if (length(na.action)) fit$na.action <- na.action
    attr(fit, "class") <-  c("survreg", "glm", "lm")
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    if (model) fit$model <- m
    if (x)     fit$x <- X
    if (y)     fit$y <- Y
    fit
    }
