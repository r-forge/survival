#SCCS  $Id: survfit.coxph.null.s,v 4.5 1992-07-14 00:02:18 therneau Exp $ % G%
survfit.coxph.null <-
  function(object, newdata, se.fit=T, conf.int=.95, individual=F,
	    type=c('tsiatis', 'kaplan-meier'),
	    conf.type=c('log', 'log-log', 'plain', 'none'), ...) {
    # May have strata and/or offset terms, linear predictor = offset
    #  newdata doesn't make any sense
    #  This is survfit.coxph with lots of lines removed

    call <- match.call()
    Terms <- terms(object)
    strat <- attr(Terms, "specials")$strata
    n <- object$n
    score <- exp(object$linear.predictor)
    method <- match.arg(type)
    if (!se.fit) conf.type <- 'none'
    else conf.type <- match.arg(conf.type)

    y <- object$y
    strata <- object$strata
    if (is.null(y) || (length(strat) && is.null(strata))) {
	# I need the model frame
	m <- model.frame(object)
	if (is.null(strata)) {
	    temp <- untangle.specials(Terms, 'strata', 1)
	    strata <- strata(m[temp$vars])
	    }
	if (is.null(y)) y <- model.extract(m, 'response')
	}
    if (is.null(strata)) strata <- rep(1,n)
    ny <- ncol(y)
    if (nrow(y) != n) stop ("Mismatched lengths: logic error")

    type <- attr(y, 'type')
    if (type=='counting') {
	ord <- order(strata, y[,2], -y[,3])
	if (method=='kaplan-meier')
	      stop ("KM method not valid for counting type data")
	}
    else if (type=='right') {
	ord <- order(strata, y[,1], -y[,2])
	y <- cbind(0, y)
	}
    else stop("Cannot handle \"", type, "\" type survival data")

    if (length(strat)) {
	newstrat <- (as.numeric(strata))[ord]
	newstrat <- as.integer(c(1*(diff(newstrat)!=0), 1))
	}
    else newstrat <- as.integer(c(rep(0,n-1),1))

    if ( !missing(newdata))
	stop("A newdata argument does not make sense for a null model")

    dimnames(y) <- NULL   #I only use part of Y, so names become invalid
    surv <- .C('agsurv2', as.integer(n),
			  as.integer(0),
			  y = y[ord,],
			  score[ord],
			  strata = newstrat,
			  surv = double(n),
			  varhaz = double(n),
			  double(1),
			  double(0),
			  nsurv = as.integer(method=='kaplan-meier'),
			  double(2),
			  as.integer(1),
			  double(1),
			  newrisk= as.double(1))
    nsurv <- surv$nsurv
    ntime <- 1:nsurv
    tsurv <- surv$surv[ntime]
    tvar  <- surv$varhaz[ntime]
    if (surv$strata[1] <=1)
	temp _ list(time=surv$y[ntime,1],
		 n.risk=surv$y[ntime,2],
		 n.event=surv$y[ntime,3],
		 surv=tsurv )
    else {
	temp <- surv$strata[1:(1+surv$strata[1])]
	tstrat <- diff(c(0, temp[-1])) #n in each strata
	names(tstrat) <- levels(strata)
	temp _ list(time=surv$y[ntime,1],
		 n.risk=surv$y[ntime,2],
		 n.event=surv$y[ntime,3],
		 surv=tsurv,
		 strata= tstrat)
	}
    if (se.fit) temp$std.err <- sqrt(tvar)

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
