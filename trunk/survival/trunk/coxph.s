# SCCS @(#)coxph.s	5.2 09/01/98
setOldClass(c("coxph.null", "coxph"))
coxph <- function(formula=formula(data), data=sys.parent(),
	weights, subset, na.action,
	eps=1e-07, init, iter.max=10, toler.chol=eps/1000, 
	method= c("efron", "breslow", "exact"),
	singular.ok =T, robust=F,
	model=F, x=F, y=T) {

    method <- match.arg(method)
    call <- match.call()
    m <- match.call(expand=F)
    m$method <- m$model <- m$x <- m$y <- m$... <-  NULL
    m$eps <- m$init <- m$iter.max <- m$robust <- m$singular.ok <- NULL

    Terms <- if(missing(data)) terms(formula, c('strata', 'cluster'))
	     else              terms(formula, c('strata', 'cluster'),data=data)
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
    if( method=="breslow" || method =="efron") {
	if (type== 'right')  fitter <- get("coxph.fit")
	else if (type=='counting') fitter <- get("agreg.fit")
	}
    else if (method=='exact') fitter <- get("agexact.fit")
    else stop(paste ("Unknown method", method))

    if (missing(init)) init <- NULL
    fit <- fitter(X, Y, strats, offset, init=init, iter.max=iter.max,
			eps=eps, toler.chol=toler.chol, weights=weights,
			method=method, row.names(m))

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
	if (robust & length(fit$coef)) {
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
	    #fit$rscore <- c(u %*% solve(t(temp0)%*%temp0, u))
	    tsvd <- svd(temp0)
	    u <- c(u %*% tsvd$v)
	    fit$rscore <- sum((u/tsvd$d)^2)
	    }

	#Wald test
	if (length(fit$coef)) {  #not for intercept only models
	    nabeta <- !is.na(fit$coef)
	    if (is.null(init)) temp <- fit$coef[nabeta]
	    else               temp <- (fit$coef - init)[nabeta]
	    #
	    # The solve function has proven to be too inaccurate
	    #  some data sets not singular in coxfit2 it considers singular
	    #fit$wald.test <-  sum(temp * solve(fit$var[nabeta,nabeta], temp))
	    tsvd <- svd(fit$var[nabeta, nabeta])
	    temp <- c(temp %*% tsvd$u)
	    fit$wald.test <- sum(temp^2/tsvd$d)
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
