#SCCS $Id: survfit.coxph.s,v 4.10 1993-03-19 15:24:55 therneau Exp $
survfit.coxph <-
  function(object, newdata, se.fit=T, conf.int=.95, individual=F,
	    type=c('tsiatis', 'kaplan-meier'),
	    conf.type=c('log', 'log-log', 'plain', 'none'))
		 {

    call <- match.call()
    Terms <- terms(object)
    strat <- attr(Terms, "specials")$strata
    resp <-  attr(Terms, "variables")[attr(Terms, "response")]
    n <- object$n
    nvar <- length(object$coef)
    score <- exp(object$linear.predictor)
    method <- match.arg(type)
    coxmethod <- object$method
    if (!se.fit) conf.type <- 'none'
    else conf.type <- match.arg(conf.type)

    x <- object$x
    y <- object$y
    if (is.null(x) && (length(strat) || se.fit)) {  # I need both X and Y
	stratum <- object$strata
	m <- model.frame(object)
	if (is.null(x)) {   #Both strata and X will be null, or neither
	    if (length(strat)) {
		temp <- untangle.specials(Terms, 'strata', 1)
		x <- model.matrix(Terms[-temp$terms], m)[,-1,drop=F]
		stratum <- strata(m[temp$vars])
		}
	    else {
		x <- model.matrix(Terms, m)[,-1,drop=F]
		stratum <- rep(1,n)
		}
	    }
	if (is.null(y)) y <- model.extract(m, 'response')
	}
    else {
	y <- object$y
	if (is.null(y)) {
	    m <- model.frame(object)
	    y <- model.extract(m, 'response')
	    }
	stratum <- rep(1,n)
	}

    ny <- ncol(y)
    if (nrow(y) != n) stop ("Mismatched lengths: logic error")
    type <- attr(y, 'type')
    if (type=='counting') {
	ord <- order(stratum, y[,2], -y[,3])
	if (method=='kaplan-meier')
	      stop ("KM method not valid for counting type data")
	}
    else if (type=='right') {
	ord <- order(stratum, y[,1], -y[,2])
	y <- cbind(-1, y)
	}
    else stop("Cannot handle \"", type, "\" type survival data")

    if (length(strat)) {
	newstrat <- (as.numeric(stratum))[ord]
	newstrat <- as.integer(c(1*(diff(newstrat)!=0), 1))
	}
    else newstrat <- as.integer(c(rep(0,n-1),1))

    if (individual && !missing(newdata)) stype <- 1
    else {
	stype <- 2
	if (length(strat)) Terms <- Terms[-temp$terms]  #don't need it
	}
    offset2 <- mean(object$linear.predictors)
    if (!missing(newdata)) {
	m2 <- model.newframe(Terms, newdata, response=(stype==1))
	if (!inherits(m2, 'data.frame'))  {
	    x2 <- as.matrix(m2)
	    if (ncol(x2) != nvar) stop ("Wrong # of variables in new data")
	    n2 <- nrow(x2)
	    if (stype==1) stop("Program error #3")
	    }

	else  {
	    x2 <- model.matrix(Terms, m2)[,-1,drop=F]
	    n2 <- nrow(x2)
	    offset2 <- model.extract(m2, 'offset')
	    if (is.null(offset2)) offset2 <- 0
	    if (stype==1) {
		#
		# The case of an agreg, with a multiple line newdata
		#
		if (length(strat)) {
		    strata2 <- factor(x2[,strat], levels=levels(stratum))
		    x2 <- x2[, -strat, drop=F]
		    }
		else strata2 <- rep(1, nrow(x2))
		y2 <- model.extract(m2, 'response')
		if (attr(y2,'type') != type)
		    stop("Survival type of newdata does not match the fitted model")
		if (nrow(y2) != n2) stop("Wrong # of rows for Y")
		}
	    }
	}
    else x2 <- matrix(object$means, nrow=1)
    n2 <- nrow(x2)
    newrisk <- exp(c(x2 %*% object$coef) + offset2 - sum(object$coef*object$means))

    dimnames(y) <- NULL   #I only use part of Y, so names become invalid
    if (stype==1) {
	surv <- .C("agsurv1", as.integer(n),
			     as.integer(nvar),
			     y[ord,],
			     score,
			     strata=newstrat,
			     surv=double(n*n2),
			     varh=double(n*n2),
			     nsurv=as.integer(2+ 1*coxmethod=='efron'),
			     x[ord,],
			     double(3*nvar),
			     object$var,
			     y = double(3*n*n2),
			     as.integer(n2),
			     y2,
			     x2,
			     newrisk,
			     as.integer(strata2) )
	ntime <- 1:surv$nsurv
	temp <- (matrix(surv$y, ncol=3))[ntime,]
	temp <- list(time = temp[,1],
		     n.risk= temp[,2],
		     n.event=temp[,3],
		     surv = surv$surv[ntime])
	if (se.fit) temp$std.err <- sqrt(surv$varh[ntime])
	}
    else {
	temp <- ifelse(method=='kaplan-meier', 1,
					2+as.integer(coxmethod=='efron'))
	surv <- .C('agsurv2', as.integer(n),
			      as.integer(nvar* se.fit),
			      y = y[ord,],
			      score[ord],
			      strata = newstrat,
			      surv = double(n*n2),
			      varhaz = double(n*n2),
			      x[ord,],
			      object$var,
			      nsurv = as.integer(temp),
			      double(3*nvar),
			      as.integer(n2),
			      x2,
			      newrisk)
	nsurv <- surv$nsurv
	ntime <- 1:nsurv
	if (n2>1) {
	    tsurv <- matrix(surv$surv[1:(nsurv*n2)], ncol=n2)
	    tvar  <- matrix(surv$varhaz[1:(nsurv*n2)], ncol=n2)
	    dimnames(tsurv) <- list(NULL, dimnames(x2)[[1]])
	    }
	else {
	    tsurv <- surv$surv[ntime]
	    tvar  <- surv$varhaz[ntime]
	    }
	if (surv$strata[1] <=1)
	    temp _ list(time=surv$y[ntime,1],
		     n.risk=surv$y[ntime,2],
		     n.event=surv$y[ntime,3],
		     surv=tsurv )
	else {
	    temp <- surv$strata[1:(1+surv$strata[1])]
	    tstrat <- diff(c(0, temp[-1])) #n in each strata
	    names(tstrat) <- levels(stratum)
	    temp _ list(time=surv$y[ntime,1],
		     n.risk=surv$y[ntime,2],
		     n.event=surv$y[ntime,3],
		     surv=tsurv,
		     strata= tstrat)
	    }
	if (se.fit) temp$std.err <- sqrt(tvar)
	}

    zval _ qnorm(1- (1-conf.int)/2, 0,1)
    if (conf.type=='plain') {
	temp1 <- temp$surv + zval* temp$std * temp$surv
	temp2 <- temp$surv - zval* temp$std * temp$surv
	temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
			conf.type='plain', conf.int=conf.int))
	}
    if (conf.type=='log') {
	xx <- ifelse(temp$surv==0,1,temp$surv)  #avoid some "log(0)" messages
	temp1 <- ifelse(temp$surv==0, 0*temp$std, exp(log(xx) + zval* temp$std))
	temp2 <- ifelse(temp$surv==0, 0*temp$std, exp(log(xx) - zval* temp$std))
	temp <- c(temp, list(upper=pmin(temp1,1), lower=temp2,
			conf.type='log', conf.int=conf.int))
	}
    if (conf.type=='log-log') {
	who <- (temp$surv==0 | temp$surv==1) #special cases
	xx <- ifelse(who, .1,temp$surv)  #avoid some "log(0)" messages
	temp1 <- exp(-exp(log(-log(xx)) + zval*temp$std/log(xx)))
	temp1 <- ifelse(who, temp$surv + 0*temp$std, temp1)
	temp2 <- exp(-exp(log(-log(xx)) - zval*temp$std/log(xx)))
	temp2 <- ifelse(who, temp$surv + 0*temp$std, temp2)
	temp <- c(temp, list(upper=temp1, lower=temp2,
			conf.type='log-log', conf.int=conf.int))
	}

    temp$call <- call
    attr(temp, 'class') <- c("survfit.cox", "survfit")
    temp
    }
