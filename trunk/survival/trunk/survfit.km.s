#SCCS $Id: survfit.km.s,v 4.15 2000-03-03 10:22:19 boos Exp $
survfit.km <- function(x, y, casewt=rep(1,n),
	    type=c('kaplan-meier', 'fleming-harrington', 'fh2'),
	    error=c('greenwood', "tsiatis"), se.fit=T,
	    conf.int= .95,
	    conf.type=c('log',  'log-log',  'plain', 'none'),
	    conf.lower=c('usual', 'peto', 'modified'))
    {
    type <- match.arg(type)
    method <- match(type, c("kaplan-meier", "fleming-harrington", "fh2"))

    error <- match.arg(error)
    error.int <- match(error, c("greenwood", "tsiatis"))
    conf.type <- match.arg(conf.type)
    conf.lower<- match.arg(conf.lower)

    ny <- ncol(y)
    n <- nrow(y)

    if (!is.Surv(y)) stop("y must be a Surv object")
    if (!is.factor(x)) stop("x must be a factor")
    if (attr(y, 'type') != 'right') stop("Can only handle right censored data")

    sorted <- (1:n)[order(x, y[,ny-1])]
    y <- y[sorted,]
    newstrat <- as.numeric(x[sorted])
    newstrat <- as.integer(c(1*(diff(newstrat)!=0), 1))
    if (sum(newstrat) > n/2)
	stop("Number of strata > number of observations/2")
    if (method==3 && any(floor(casewt) != casewt))
	stop("The fh2 method is not valid for fractional case weights")

    storage.mode(y) <- "double"
    dimnames(y) <- NULL
    surv <- .C("survfit2", as.integer(n),
			  y = y,
			  as.integer(ny),
			  as.double(casewt[sorted]),
			  strata= as.integer(newstrat),
			  nstrat= as.integer(method),
			  as.integer(error.int),
			  mark=double(n),
			  surv=double(n),
			  varhaz=double(n),
			  risksum=double(n),
			  ntime = integer(1))
    ntime <- surv$ntime
    if (error.int==1) surv$varhaz[surv$surv==0] <- NA
    ntime <- 1:ntime
    if (surv$nstrat ==1)
	temp _ list(time=surv$y[ntime,1],
		 n.risk=surv$risksum[ntime],
		 n.event=surv$mark[ntime],
		 surv=surv$surv[ntime])
    else {
	temp <- surv$strata[1:surv$nstrat]
	tstrat <- diff(c(0, temp)) #n in each strata
	names(tstrat) <- levels(x)
	temp _ list(time=surv$y[ntime,1],
		 n.risk=surv$risksum[ntime],
		 n.event=surv$mark[ntime],
		 surv=surv$surv[ntime],
		 strata= tstrat)
	}

    if (se.fit) {
	std.err <- sqrt(surv$varhaz[ntime])
	temp$std.err <- std.err
	events <- temp$n.event >0
	n.lag <- rep(c(temp$n.risk[1], temp$n.risk[events]),
	              diff(c(ntime[1], ntime[events], 1+max(ntime))))
	std.low <- switch(conf.lower,
			'usual'   = std.err,
			'peto'    = sqrt((1-temp$surv)/ temp$n.risk),
			'modified'= std.err * sqrt(n.lag/temp$n.risk)
			)
	zval _ qnorm(1- (1-conf.int)/2, 0,1)
	if (conf.type=='plain') {
	    temp1 <- temp$surv + zval* std.err * temp$surv
	    temp2 <- temp$surv - zval* std.low * temp$surv
	    temp <- c(temp, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
			    conf.type='plain', conf.int=conf.int))
	    }
	if (conf.type=='log') {
	    xx <- ifelse(temp$surv==0,1,temp$surv)  #avoid some "log(0)" messages
	    temp1 <- ifelse(temp$surv==0, NA, exp(log(xx) + zval* std.err))
	    temp2 <- ifelse(temp$surv==0, NA, exp(log(xx) - zval* std.low))
	    temp <- c(temp, list(upper=pmin(temp1,1), lower=temp2,
			    conf.type='log', conf.int=conf.int))
	    }
	if (conf.type=='log-log') {
	    who <- (temp$surv==0 | temp$surv==1) #special cases
	    temp3 <- ifelse(temp$surv==0, NA, 1)
	    xx <- ifelse(who, .1,temp$surv)  #avoid some "log(0)" messages
	    temp1 <- exp(-exp(log(-log(xx)) + zval*std.err/log(xx)))
	    temp1 <- ifelse(who, temp3, temp1)
	    temp2 <- exp(-exp(log(-log(xx)) - zval*std.low/log(xx)))
	    temp2 <- ifelse(who, temp3, temp2)
	    temp <- c(temp, list(upper=temp1, lower=temp2,
			    conf.type='log-log', conf.int=conf.int))
	    }
	}
    temp
    }
