#SCCS $Date: 1992-06-26 23:27:25 $ $Id: survreg.s,v 4.5 1992-06-26 23:27:25 therneau Exp $
survreg <- function(formula=formula(data), data=sys.parent(),
	subset, na.action,
	link=c('log', 'identity'), dist=c("extreme", "logistic", "exp"),
	scale,
	eps=.0001, init, iter.max=10,
	model=F, x=F, y=F, ...) {

    call <- match.call()
    m <- match.call(expand=F)
    m$dist <- m$link <- m$model <- m$x <- m$y <- m$... <-  NULL
    m$eps <- m$inf.ratio <- m$init <- m$iter.max <- m$n.table <- NULL

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
	if (any(Y[,3]==4))
	     Y <- cbind(linkfun(Y[,1:2]), Y[,3])
	else Y <- cbind(linkfun(Y[,1]), Y[,3])
	}
    else     Y <- cbind(linkfun(Y[,1]), Y[,2])

    if (missing(init)) init <- NULL
    if (missing(scale)) scale <- NULL

    if( dist=="extreme")
	sfit <- survreg.fit(X, Y, offset, init=init, iter.max=iter.max,
			    eps = eps, fun='exvalue', scale=scale, ...)
    else if (dist=="logistic")
	sfit <- survreg.fit(X, Y, offset, init=init, iter.max=iter.max,
			    eps = eps, fun='logisticfit', scale=scale, ...)
    else if (dist=="exp")
	sfit <- survreg.fit(X, Y, offset, init=init, iter.max=iter.max,
			    eps = eps, fun='exvalue', scale=1, ...)

    else stop(paste ("Unknown method", dist))

    if (is.character(sfit))  fit <- list(fail=sfit)  #error message
    else {
	# There are certainly more clever ways to do this, but ....
	#  In order to make it look like IRLS, and get all the components
	#  that I need for glm inheritance, do one step of weighted least
	#  squares.
	if (missing(scale)) scale <- sfit$coef[nvar+1]
	eta <- c(X %*% sfit$coef[1:nvar]) + offset
	wt<- -sfit$deriv[,3]/scale^2
	fit <- lm.wfit(X, eta + sfit$deriv[,2]/(scale*wt), wt, "qr", ...)

	na.action <- attr(m, "na.action")
	if (length(na.action)) fit$na.action <- na.action
	if (link=='log') fit$fitted.values <- exp(fit$fitted.values)
	fit$family <- c(name=dist, link=link, "")
	fit$linear.predictors <- eta
	fit$deviance <- -2*sfit$loglik[1]
	fit$null.deviance <- -2*sfit$null[1]
	fit$iter <- sfit$iter
	fit$var  <- solve(sfit$imat)
	fit$flag <- sfit$flag
	fit$dresiduals <- sign(fit$residuals)*sqrt(-2*sfit$deriv[,1])
	}

    attr(fit, "class") <-  c("survreg", "glm", "lm")
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    if (model) fit$model <- m
    if (x)     fit$x <- X
    if (y)     fit$y <- Y
    fit
    }
