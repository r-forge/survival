#SCCS 3/4/92 @(#)coxreg.s	4.1
coxreg <- function(formula=formula(data), data=sys.parent(),
	weights,  subset, na.action,
	eps=.0001, inf.ratio=200, init, iter.max=10,
	method= c("cox", "cox.efron", "cox.exact", "model.frame"),
	model=F, x=F, y=T) {

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
    casewt<- model.extract(m, "weights")
    offset<- attr(Terms, "offset")
    if (!is.null(offset)) offset <- as.numeric(m[[offset]])

    attr(Terms,"intercept")<- 0  #Cox model has no explicit intercept term
    strata <- attr(Terms, "specials")$strata
    if (length(strata)>1) stop("Only one strata() expression allowed")
    else if (length(strata)==1) {
	X <- model.matrix(Terms[-strata], m)
	strata <- m[[(as.character(Terms))[strata]]]
	}
    else X <- model.matrix(Terms, m)

    type <- attr(Y, "type")
    if( method=="cox" || method =="cox.efron") {
	if (type== 'right')  fitter <- get("coxreg.fit")
	else if (type=='counting') fitter <- get("agreg.fit")
	else stop(paste("Cox model doesn't support \"", type,
			  "\" survival data", sep=''))
	}
    else if (method=='cox.exact') fitter <- get("agexact.fit")
    else stop(paste ("Unknown method", method))

    if (missing(init)) init <- NULL
    fit <- fitter(Y, X, strata, casewt, offset, iter.max=iter.max,
			eps=eps, inf.ratio=inf.ratio, init=init,
			method=method, row.names(m) )

    if (is.character(fit)) {
	fit <- list(fail=fit)
	attr(fit, 'class') <- 'coxreg'
	}
    else {
	fit$n <- nrow(Y)
	omit <- attr(m, "omit")
	if (length(omit)) attr(fit$n, "omit") <- omit
	attr(fit, "class") <-  c(fit$method, "surv.reg")
	fit$method <- NULL
	fit$terms <- Terms
	fit$assign <- attr(X, 'assign')
	if (model) fit$model <- m
	if (x)  {
	    fit$x <- X
	    if (length(strata)) fit$strata <- strata
	    }
	if (y)     fit$y <- Y
	}
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    fit
    }
