#SCCS $Date: 1992-04-13 23:05:51 $ $Id: survreg.s,v 4.3 1992-04-13 23:05:51 therneau Exp $
surv.reg <- function(formula=formula(data), data=sys.parent(),
	subset, na.action,
	eps=.0001, init, iter.max=10,
	method= c("weibull"),
	model=F, x=F, y=F, ...) {

    method <- match.arg(method)
    call <- match.call()
    m <- match.call(expand=F)
    m$method <- m$model <- m$x <- m$y <- m$... <-  NULL
    m$eps <- m$inf.ratio <- m$init <- m$iter.max <- m$n.table <- NULL

    Terms <- if(missing(data)) terms(formula, 'strata')
	     else              terms(formula, 'strata',data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    if (method== 'model.frame') return(m)

    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    offset<- attr(Terms, "offset")
    if (!is.null(offset)) offset <- as.numeric(m[[offset]])
    X <- model.matrix(Terms, m)

    type <- attr(Y, "type")
    if (type!='right' && type!='counting') stop ("Invalid survival type")

    if( method=="weibull") fitter <- weibull.fit
    else stop(paste ("Unknown method", method))

    casewt <- rep(1,ncol(X))
    if (missing(init)) init <- NULL
    fit <- fitter(Y, X, casewt, offset, iter.max=iter.max,
			eps=eps, init=init, ...)

    if (is.character(fit))  fit <- list(fail=fit)
    else {
	fit$n <- nrow(Y)
	na.action <- attr(m, "na.action")
	if (length(na.action)) fit$na.action <- na.action
	}

    attr(fit, "class") <-  c(method, "surv.reg", "lm")
    fit$terms <- Terms
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    fit$assign <- attr(X, 'assign')
    if (model) fit$model <- m
    if (x)     fit$x <- X
    if (y)     fit$y <- Y
    fit
    }
