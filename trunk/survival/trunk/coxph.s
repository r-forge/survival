#SCCS  $Id: coxph.s,v 4.10 1993-06-17 12:26:18 therneau Exp $
coxph <- function(formula=formula(data), data=sys.parent(),
	weights, subset, na.action,
	eps=.0001, init, iter.max=10,
	method= c("efron", "breslow", "exact"),
	singular.ok =T,
	model=F, x=F, y=T) {

    method <- match.arg(method)
    call <- match.call()
    m <- match.call(expand=F)
    m$method <- m$model <- m$x <- m$y <- m$... <-  NULL
    m$eps <- m$init <- m$iter.max <- m$n.table <- NULL

    Terms <- if(missing(data)) terms(formula, 'strata')
	     else              terms(formula, 'strata',data=data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")
    weights <- model.extract(m, 'weights')
    offset<- attr(Terms, "offset")
    tt <- length(offset)
    offset <- if(tt == 0)
		    rep(0, nrow(Y))
	      else if(tt == 1)
		      m[[offset]]
	      else {
		    ff <- frame[[offset[1]]]
		    for(i in 2:tt)
			    ff <- ff + frame[[offset[i]]]
		    ff
		    }

    attr(Terms,"intercept")<- 1  #Cox model always has \Lambda_0
    strats <- attr(Terms, "specials")$strata
    if (length(strats)) {
	temp <- untangle.specials(Terms, 'strata', 1)
	X <- model.matrix(Terms[-temp$terms], m)[,-1,drop=F]
	if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
	else strata.keep <- strata(m[temp$vars], shortlabel=T)
	strats <- as.numeric(strata.keep)
	}
    else X <- model.matrix(Terms, m)[,-1,drop=F]   #remove column of 1's though

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
    fit <- fitter(X, Y, strats, offset, init=init, iter.max=iter.max,
			eps=eps, weights=weights,
			method=method, row.names(m))

    if (is.character(fit)) {
	fit <- list(fail=fit)
	attr(fit, 'class') <- 'coxph'
	}
    else {
	if (any(is.na(fit$coef))) {
	   vars <- (1:length(fit$coef))[is.na(fit$coef)]
	   msg <-paste("X matrix deemed to be singular; variable",
			   paste(vars, collapse=" "))
	   if (singular.ok) warning(msg)
	   else             stop(msg)
	   }
	fit$n <- nrow(Y)
	na.action <- attr(m, "na.action")
	if (length(na.action)) fit$na.action <- na.action
	attr(fit, "class") <-  fit$method
	fit$method <- NULL
	fit$terms <- Terms
	fit$assign <- attr(X, 'assign')
	if (model) fit$model <- m
	else {
	    if (x)  {
		fit$x <- X
		if (length(strats)) fit$strata <- strata.keep
		if (length(weights))fit$weights<- weights
		}
	    if (y)     fit$y <- Y
	    if (!is.null(weights) && any(weights!=1)) fit$weights <- weights
	    }
	}
    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    fit$method <- method
    fit
    }
