#SCCS  $Id: coxph.s,v 4.20 1996-08-30 13:38:37 therneau Exp $
coxph <- function(formula=formula(data), data=sys.parent(),
	weights, subset, na.action,
	eps=.0001, init, iter.max=10,
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
	attr(fit, "class") <-  fit$method
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
		indx <- match(cluster, unique(cluster))
		k    <- as.vector(table(indx))
		if (any(k>1)) {
		    #compute the ICC for the m-residuals
		    N    <- sum(k * (k-1))
		    mu   <- sum(fit$resid * (k-1)[indx])/N   #grand mean
		    mu2  <- tapply(fit$resid, indx, mean)    # group means
		    sig  <- tapply((fit$resid - mu)^2, indx, sum)  #SS
		    icc1 <- sum( (k*(mu2-mu))^2 - sig) / sum((k-1)*sig)
		    #rank residuals
		    rr <- rank(fit$resid)
		    mu   <- sum(rr * (k-1)[indx])/N   #grand mean
		    mu2  <- tapply(rr, indx, mean)    # group means
		    sig  <- tapply((rr - mu)^2, indx, sum)  #SS
		    icc2 <- sum( (k*(mu2-mu))^2 - sig) / sum((k-1)*sig)

		    fit$icc <- c(length(k), icc1, icc2)
		    names(fit$icc) <- c("nclust", "icc(resid)",
						    "icc(rank(resid))")
		    }
		}
	    else temp <- residuals.coxph(fit2, type='dfbeta', weighted=T)
	    fit$var <- t(temp) %*% temp
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
