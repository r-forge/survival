#SCCS $Id: survreg.s,v 4.14 1993-03-29 14:27:25 therneau Exp $
survreg <- function(formula=formula(data), data=sys.parent(),
	subset, na.action,
	link='log',
	dist=c("extreme", "logistic", "gaussian", "exponential",
	       "rayleigh"),
	init=NULL,  fixed=list(), control,
	model=F, x=F, y=T, ...) {

    call <- match.call()
    m <- match.call(expand=F)
    m$dist <- m$link <- m$model <- m$x <- m$y <- m$... <-  NULL
    m$start <- m$fixed <- m$control <- NULL
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    Terms <- attr(m, 'terms')

    dist <- match.arg(dist)
    lnames <- dimnames(glm.links)[[2]]
    link <- pmatch(link, lnames, 0)
    if (link==0) stop("Invalid link function")
    else link <- lnames[link]

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

    controlvals <- survreg.control()
    if (!missing(control)) 
	controlvals[names(control)] <- control

    if( dist=="exponential") {
	fixed$scale <- 1
	dist <- 'extreme'
	}
    else if (dist=="rayleigh") {
	fixed$scale <- .5
	dist <- 'extreme'
	}

    sd <- survreg.distributions[[dist]]
    if (length(fixed)>0) {
	ifix <- match(names(fixed), names(sd$parms), nomatch=0)
	if (any(ifix==0))
	    stop (paste("Parameter(s)", paste(names(fixed)[ifix==0]),
			"in the fixed list not valid for this dist"))
	}
    if (is.list(init) && length(init)>0) {
	ifix <- match(names(init), c('eta',names(sd$parms)), nomatch=0)
	if (any(ifix==0))
	    stop (paste("Parameter(s)", paste(names(init)[ifix==0]),
			"in the init list not valid for this dist"))
	}


    sfit <- survreg.fit(X, Y, offset, init=init, controlvals=controlvals,
			dist= dist, fixed=fixed)
    if (is.character(sfit))  fit <- list(fail=sfit)  #error message
    else {
	# There may be more clever ways to do this, but ....
	#  In order to make it look like IRLS, and get all the components
	#  that I need for glm inheritance, do one step of weighted least
	#  squares.
	eta <- c(X %*% sfit$coef[1:nvar]) + offset
	wt<- -sfit$deriv[,3]
	fit <- lm.wfit(X, eta + sfit$deriv[,2]/wt, wt, "qr", ...)

	ifun <- glm.links['inverse',link][[1]]
	fit$fitted.values <- ifun(fit$fitted.values)
	fit$family <- c(name=dist, link=link, "")
	fit$linear.predictors <- eta
	fit$iter <- sfit$iter
	fit$parms <- sfit$parms
	fit$df.residual <- fit$df.residual - sum(!sfit$fixed)

	# If singular.ok=T, there may be NA coefs.  The var matrix should
	#   be an inversion of the "non NA" portion.
	var <- 0*sfit$imat
	good <- c(!is.na(fit$coef), rep(T, ncol(var)-nvar))
	var[good,good] <- solve(qr(sfit$imat[good,good], tol=1e-12))
	fit$var <- var
	fit$fixed <- sfit$fixed
	dev <- sd$deviance(Y, fit$parms, sfit$deriv[,1])
	fit$dresiduals <- sign(fit$residuals)*sqrt(dev)
	fit$deviance <- sum(dev)
	fit$null.deviance <- fit$deviance +2*(sfit$loglik[2]- sfit$ndev[2])
	fit$loglik <- c(sfit$ndev[2], sfit$loglik[2])
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
