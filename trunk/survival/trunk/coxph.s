#SCCS 3/4/92 @(#)coxph.s	4.1
coxph <- function(formula=formula(data), data=sys.parent(),
	subset, na.action,
	eps=.0001, inf.ratio=200, init, iter.max=10,
	method= c("breslow", "efron", "exact"),
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

    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
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
    if( method=="breslow" || method =="efron") {
	if (type== 'right')  fitter <- get("coxph.fit")
	else if (type=='counting') fitter <- get("agreg.fit")
	else stop(paste("Cox model doesn't support \"", type,
			  "\" survival data", sep=''))
	}
    else if (method=='exact') fitter <- get("agexact.fit")
    else stop(paste ("Unknown method", method))

    if (missing(init)) init <- NULL
    fit <- fitter(X, Y, strata, offset, iter.max=iter.max,
			eps=eps, inf.ratio=inf.ratio, init=init,
			method=method, row.names(m) )

    if (is.character(fit)) {
	fit <- list(fail=fit)
	attr(fit, 'class') <- 'coxph'
	}
    else {
	fit$n <- nrow(Y)
	na.action <- attr(m, "na.action")
	if (length(na.action)) fit$na.action <- na.action
	attr(fit, "class") <-  c(fit$method, "survreg")
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
