#SCCS $Id: survfit.km.s,v 4.2 1992-04-14 18:08:14 grill Exp $
survfit.km <- function(x, y, casewt=rep(1,n),
	    type=c('kaplan-meier', 'fleming-harrington'),
	    error=c('greenwood', "tsiatis"), se.fit=T,
	    conf.int= .95,
	    conf.type=c('log', 'log-log', 'plain', 'none'))
    {
    type <- match.arg(type)
    method <- match(type, c("kaplan-meier", "fleming-harrington"))

    error <- match.arg(error)
    error.int <- match(error, c("greenwood", "tsiatis"))
    conf.type <- match.arg(conf.type)

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

    storage.mode(y) <- "double"
    surv <- .C("survfit", as.integer(n),
			  y = y,
			  as.integer(ny),
			  as.double(casewt),
			  strata= as.integer(newstrat),
			  as.integer(method),
			  as.integer(error.int),
			  mark=integer(n),
			  surv=double(n),
			  varhaz=double(n),
			  risksum=double(n),
			  ntime = integer(1))
    ntime <- surv$ntime
    if (error.int==1) surv$varhaz[surv$surv==0] <- NA
    ntime <- 1:ntime
    if (surv$strata[1] ==1)
	temp _ list(time=surv$y[ntime,1],
		 n.risk=surv$risksum[ntime],
		 n.event=surv$mark[ntime],
		 surv=surv$surv[ntime])
    else {
	temp <- surv$strata[1:(1+surv$strata[1])]
	tstrat <- diff(c(0, temp[-1])) #n in each strata
	names(tstrat) <- levels(x)
	temp _ list(time=surv$y[ntime,1],
		 n.risk=surv$risksum[ntime],
		 n.event=surv$mark[ntime],
		 surv=surv$surv[ntime],
		 strata= tstrat)
	}

    if (se.fit) {
	temp$std.err <- sqrt(surv$varhaz[ntime])
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
	}
    temp
    }
