#SCCS  $Id: coxph.s,v 5.7 1998-11-02 19:52:56 therneau Exp $
# Version with general penalized likelihoods
setOldClass(c('coxph.penal', 'coxph'))

coxph <- function(formula=formula(data), data=sys.parent(),
	weights, subset, na.action,
	eps=.0001, init, iter.max=10, toler.chol=1e-9,
	method= c("efron", "breslow", "exact"),
	singular.ok =T, robust=F, outer.max=10,
	model=F, x=F, y=T) {

    method <- match.arg(method)
    call <- match.call()
    m <- match.call(expand=F)
    temp <- c("", "formula", "data", "weights", "subset", "na.action")
    m <- m[ match(temp, names(m), nomatch=0)]
    special <- c("strata", "cluster")
    Terms <- if(missing(data)) terms(formula, special)
	     else              terms(formula, special, data=data)
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
		    ff <- m[[offset[1]]]
		    for(i in 2:tt)
			    ff <- ff + m[[offset[i]]]
		    ff
		    }

    attr(Terms,"intercept")<- 1  #Cox model always has \Lambda_0
    strats <- attr(Terms, "specials")$strata
    cluster<- attr(Terms, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
	if (missing(robust)) robust <- T
	tempc <- untangle.specials(Terms, 'cluster', 1:10)
	ord <- attr(Terms, 'order')[tempc$terms]
	if (any(ord>1)) stop ("Cluster can not be used in an interaction")
	cluster <- strata(m[,tempc$vars], shortlabel=T)  #allow multiples
	dropx <- tempc$terms
	}
    if (length(strats)) {
	temp <- untangle.specials(Terms, 'strata', 1)
	dropx <- c(dropx, temp$terms)
	if (length(temp$vars)==1) strata.keep <- m[[temp$vars]]
	else strata.keep <- strata(m[,temp$vars], shortlabel=T)
	strats <- as.numeric(strata.keep)
	}

    if (length(dropx)) X <- model.matrix(Terms[-dropx], m)[,-1,drop=F]
    else               X <- model.matrix(Terms, m)[,-1,drop=F]
	
    type <- attr(Y, "type")
    if (type!='right' && type!='counting')
	stop(paste("Cox model doesn't support \"", type,
			  "\" survival data", sep=''))
    if (missing(init)) init <- NULL

    # Check for penalized terms
    pterms <- sapply(m, inherits, 'coxph.penalty')
    if (any(pterms)) {
	pattr <- lapply(m[pterms], attributes)
	# 
	# the 'order' attribute has the same components as 'term.labels'
	#   pterms always has 1 more (response), sometimes 2 (offset)
	# drop the extra parts from pterms
	temp <- c(attr(Terms, 'response'), attr(Terms, 'offset'))
	if (length(dropx)) temp <- c(temp, dropx+1)
	pterms <- pterms[-temp]
	temp <- match((names(pterms))[pterms], attr(Terms, 'term.labels'))
	ord <- attr(Terms, 'order')[temp]
	if (any(ord>1)) stop ('Penalty terms cannot be in an interaction')
	pcols <- (attr(X, 'assign')[-1])[pterms]  
  
	if (missing(eps)) eps <- 1e-7
	if (missing(iter.max)) iter.max <- 20  #penalized are hard sometimes
        fit <- coxpenal.fit(X, Y, strats, offset, init=init,
				iter.max=iter.max, outer.max=outer.max, 
			        eps=eps, toler.chol=toler.chol,
				weights=weights, method=method,
				row.names(m), pcols, pattr)
	}
    else {
	if( method=="breslow" || method =="efron") {
	    if (type== 'right')  fitter <- get("coxph.fit")
	    else                 fitter <- get("agreg.fit")
	    }
	else if (method=='exact') fitter <- get("agexact.fit")
	else stop(paste ("Unknown method", method))

	fit <- fitter(X, Y, strats, offset, init=init, iter.max=iter.max,
			    eps=eps, toler.chol=toler.chol, weights=weights,
			    method=method, row.names(m))
	}

    if (is.character(fit)) {
	fit <- list(fail=fit)
	oldClass(fit) <- 'coxph'
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
	oldClass(fit) <-  fit$method[1]
	fit$terms <- Terms
	fit$assign <- attr(X, 'assign')
	if (robust) {
	    fit$naive.var <- fit$var
	    fit$method    <- method
	    # a little sneaky here: by calling resid before adding the
	    #   na.action method, I avoid having missings re-inserted
	    # I also make sure that it doesn't have to reconstruct X and Y
	    fit2 <- c(fit, list(x=X, y=Y, weights=weights))
	    if (length(strats)) fit2$strata <- strata.keep
	    if (length(cluster)) {
		temp <- residuals.coxph(fit2, type='dfbeta', collapse=cluster,
					  weighted=T)
		# get score for null model
		if (is.null(init))
			fit2$linear.predictors <- 0*fit$linear.predictors
		else fit2$linear.predictors <- c(X %*% init)
		temp0 <- residuals.coxph(fit2, type='score', collapse=cluster,
					 weighted=T)
		}
	    else {
		temp <- residuals.coxph(fit2, type='dfbeta', weighted=T)
		fit2$linear.predictors <- 0*fit$linear.predictors
		temp0 <- residuals.coxph(fit2, type='score', weighted=T)
	        }
	    fit$var <- t(temp) %*% temp
	    u <- apply(as.matrix(temp0), 2, sum)
	    fit$rscore <- coxph.wtest(t(temp0)%*%temp0, u, toler.chol)
	    }
	#Wald test
	if (length(fit$coef) && is.null(fit$wald.test)) {  
	    #not for intercept only models, or if test is already done
	    nabeta <- !is.na(fit$coef)
	    if (is.null(init)) temp <- fit$coef[nabeta]
	    else temp <- (fit$coef - init)[nabeta]
	    fit$wald.test <-  coxph.wtest(fit$var[nabeta,nabeta], temp,
					  toler.chol)$test
	    }
	na.action <- attr(m, "na.action")
	if (length(na.action)) fit$na.action <- na.action
	if (model) fit$model <- m
	else {
	    if (x)  {
		fit$x <- X
		if (length(strats)) fit$strata <- strata.keep
		}
	    if (y)     fit$y <- Y
	    }
	}
    if (!is.null(weights) && any(weights!=1)) fit$weights <- weights

    fit$formula <- as.vector(attr(Terms, "formula"))
    fit$call <- call
    fit$method <- method
    fit
    }
