#
# SCCS $Id: survfit.ci.s,v 1.2 2003-05-08 09:20:31 therneau Exp $
#
# Compute the current incidence curve for a data set
#   A strata() statement identifies the outcomes
#   A cluster() statement identifies repeated subjects
#
setOldClass(c('survfit.ci', 'survfit'))

survfit.ci <- function(formula=formula(data), data=sys.parent(),
		       weights, subset, na.action,
                       absorb,
		       type=c("kaplan-meier", "fleming-harrington", "fh2"),
		       conf.int= .95, se.fit=T,
		       conf.type=c('log',  'log-log',  'plain', 'none'),
		       conf.lower=c('usual', 'peto', 'modified'),
		       model=F, x=F, y=F) {

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
    n <- nrow(Y)
    conf.type <- match.arg(conf.type)
    conf.lower<- match.arg(conf.lower)

    method <- match.arg(type)
    if (method != 'kaplan-meier')
        warning("type=kaplan-meier assumed for competing risk estimation")
    weights <- model.extract(m, 'weights')
    if (length(weights)==0) weights <- rep(1.0, n)

    strats <- attr(Terms, "specials")$strata
    cluster<- attr(Terms, "specials")$cluster
    dropx <- NULL
    if (length(cluster)) {
	tempc <- untangle.specials(Terms, 'cluster', 1:10)
	ord <- attr(Terms, 'order')[tempc$terms]
	if (any(ord>1)) stop ("Cluster can not be used in an interaction")
        if (length(tempc$vars) >1) stop("Only 1 cluster term allowed")
	id <- m[,tempc$vars]
        id <- match(id, unique(id))  # a label from 1 to n
	dropx <- tempc$terms
	}
    if (length(strats)) {
	temp <- untangle.specials(Terms, 'strata', 1)
	dropx <- c(dropx, temp$terms)
	if (length(temp$vars)==1) state <- m[[temp$vars]] #event type
        else stop("Only 1 strata term allowed")
        }
    else stop("A strata() term is required for competing risks estimates")

    if (length(dropx)) 
          ll <- attr(Terms, 'term.labels')[-dropx]
    else  ll <- attr(Terms, 'term.labels') 
    if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
    else X <- strata(m[ll])
    ncurve <- length(levels(X))
	
    type <- attr(Y, "type")
    if (type!='right' && type!='counting')
	stop(paste("Cumulative incidence function doesn't support \"", type,
			  "\" survival data", sep=''))

    # Now to work
    #       Surv(time, status) ~ sex + strata(event.type)
    # is the canonical call.  This will produce a matrix of survivals,
    # stratified by sex.  Optionally, there can be a cluster statement
    # to indicate a multiple transitions/subject data set.
    #
    #   If status==1, then a transition at time "time" is assumed (or at stop
    # time for (start, stop] data.  If it is 0, then the subject is
    # censored.  If censored, the value of the event type variable
    # is completely ignored.  However, don't use NA in that case, since
    # the na.action will remove that line of data!  Subjects are assumed
    # to start in a "null" state, which is not tabulated for survival.
    # To change this behavior, give subjects a transition at time 0.
    #   
    # If there is not a cluster term, then the result can be computed
    #  by post-processing the KM.  If there is a cluster term, then we
    #  must explicitly use the redistribut-to-the-right algorithm to 
    #  compute the estimate
    #
    time   <- Y[,ncol(Y) -1]
    status <- Y[,ncol(Y)]
    # make sure that the states for status==0 are not tabulated
    #  in our output
    if (is.factor(state)) {
        temp <- levels(state[status !=0, drop=T])
        state <- factor(ifelse(status==0, NA, state))
        levels(state) <- temp
        }
    else state <- factor(ifelse(status==0, NA, state))
    
    if (length(cluster)==0) {
        # this case implicitly assumes one transition/subject
        # use survfit.km to get the # at risk and etc
	kfit <- survfit.km(X, Y, se.fit=F)
        
        # partition out the event type
	if (is.null(kfit$strata)) {
	    temp <- table(time, state, status, exclude=NA)[,,2]
	    temp <- temp / ifelse(kfit$n.event==0, 1, kfit$n.event)  # percents
	    jumps <- -diff(c(1, kfit$surv))
	    newjump <- jumps * temp
	    newsurv <- apply(newjump, 2, function(x) 1 - cumsum(x))
	    }
	else {  # do it one curve at a time
	    nstate <- length(levels(state))
	    newsurv <- matrix(1.0, nrow=length(kfit$time), ncol=nstate)
	    newjump <- 0*newsurv
	    stemp <- names(kfit$strata)
	    istart <- cumsum(c(1,strata)) #starting points
	    for (i in 1:length(strata)) {
		# j will index the rows in the survival (kfit) object
		# who will index into the original data
		j <- seq(from=istart[i], length=kfit$strata[i])
		who <- (X==stemp[i] & status==2)
		if (any(who)) { # if any deaths
		    temp <- table(factor(time[who], levels=kfit$time[j]),
				  state[who], exclude=NA)
		    tevent <- kfit$n.event[j]
		    temp <- temp / ifelse(tevent, 1, tevent)  # percents
		    jumps <- -diff(c(1, kfit$surv[j]))
		    jumps <- jumps * temp
		    newsurv[j,] <- apply(jumps, 2, function(x) 1 - cumsum(x))
		    newjump[j,] <- jumps  # save for se computations
		    }
		}
	    }
	dimnames(newsurv) <- list(NULL, levels(state))
	kfit$surv <- newsurv
		
	# compute the Greenwood error
	if (se.fit) {
	    greenwood <- newjump /(kfit$n.risk*(1-newjump))
	    kfit$std.err <- sqrt(apply(greenwood, 2, cumsum))
	    }
        }
    
    else {
        # Have to do real work in this case
 	if (ncol(Y) !=2) {
	    stop("(start, stop] data for repeated events not yet allowed")
	    }
	time   <- Y[,1]
	status <- Y[,2]
		    
	# ensure consistency of the weights
        idlist <- sort(unique(id))  #should be 1,2,3,...
	firstone <- match(idlist, id) # first obs for each subject
	wt2 <- weights[firstone]
	if (any(wt2[id] != weights))
		stop("All case weights for a given subject must agree")

	if (ncurve >1) {
	    # Ensure that any given subject is in only one curve
	    xtemp <- X[firstone]
	    if (any(xtemp[as.numeric(X)] != X))
		stop("No subject can be in more than 1 curve")
	    }

	# Do some preprocessing, to make things nice integers (1,2, etc)
	state2 <- as.numeric(state)      # numeric states
        nstate <- length(levels(state))

        # If the absorb argument is present, it lists the states that are
        #  "absorbing" states.  If someone's last transition is at time
        #  t to a non-absorbing state, that is, to a state that allows
        #  further transitions out of it, then we have to censor them at
        #  time t+0.
        # If the argument is absent, consider all states prior to the last of
        #  each subject to be non-absorbing (a transition was from them).
        # To process the data, we need to figure out each subject's last
        #  observation. The status variable is set to 2 for
        #  "had an event, it's their last obs, new state not absorbing".
        temp <- order(id, time)
        id2 <- id[temp]                 #sorted by time within subject
        indx <- match(idlist, id2)      #first obs for each subject
        lastone <- c(indx[-1]-1, length(id2))  #last obs for each person
        lastone <- temp[lastone]   #put index back into data set order
        
        if (!missing(absorb)) {
            absorb <- match(absorb, levels(state))
            if (any(is.na(absorb)))
              stop("The absorb argument contains levels not found in the data")
            not.absorb <- c(0, (1:nstate)[-absorb])
            }
        else {
            # All states except the last are non-absorbing ones.
            not.absorb <- sort(unique(c(0, state[-lastone])))
            }

        who <-  (status==1) & (!is.na(match(state2, not.absorb)))
        status[lastone[who[lastone]]] <- 2

	# This function creates a single survival curve
	docurve <- function(time, status, state, weights, id, nstate) {
	    n <- length(time)
	    timelist <- sort(unique(time))   # unique event/censor times
	    ntime <- length(timelist)
	    surv   <- matrix(0., nrow=ntime, ncol=nstate)
	    varhaz <- surv
	    nrisk <- nevent <- double(ntime)

	    # create a unique weight for each subject
	    # cstate is the current state for each subject, initially 0 or
	    #  "no state at all".  The factor causes all survival matrices
	    #  to be the same size.
	    idlist <- unique(id)
	    wtindex <- match(id, idlist)
	    wt2 <- weights[match(idlist, id)]
	    cstate <- factor(rep(0, length(idlist)), levels=0:nstate)

            # For computing "number at risk", I will need the original
            #  weights (censorings not redistributed)
            wt.orig <- wt2

	    # Walk down the timelist, one at a time
	    for (i in 1:ntime) {
		xtime <- timelist[i]
                at.risk <- unique(wtindex[time>= xtime])
		nrisk[i] <- sum(wt.orig[at.risk]) #still at risk for transition
		event <- (time==xtime & status >0)
		nevent[i] <- sum(wt.orig[wtindex[event]])
            
		# reassign the state of those with an event and 
		#  recompute survival
		if (any(event)) {
		    cstate[wtindex[event]] <- state[event]
		    temp <- tapply(wt2, cstate, sum)/sum(wt2)
		    surv[i,] <- ifelse(is.na(temp), 0, temp)[-1]
		    varhaz[i,] <- (surv[i,]/ (1-surv[i,])) * sum(wt2^2)/
			(sum(wt2))^2
		    }
		else {
		    if (i>1) {
			surv[i,] <- surv[i-1,]
			varhaz[i,] <- varhaz[i-1,]
			}
		    }
		# Now, redistribute the weights if those who are censored
		redo <- (status!=1  & time==xtime)  #censored at this time
		if (any(redo)) {
		    who2 <- (1:n)[redo]  #list of those censored here
		    for (j in who2) {
			# list of others in the same state, at this time, with
			#   a positive weight
			k <- wtindex[j]
			child <- (cstate==cstate[k] & wt2 >0)
			nchild <- sum(child) -1  # since "child" counts me too
			# if no children, the weight is not redistributed
			if (nchild >0){
			    wt2[child] <- wt2[child] + wt2[k]/nchild
			    wt2[k] <- 0
			    }
			}
		    }
		}
	    list(time=timelist, surv= 1.0 - surv, n.risk=nrisk,
		 n.event=nevent, std.err= sqrt(varhaz))
	    }
	
        # Ok, now we're actually ready to do the work.  Copy the
	#  list/unlist trick of survfit.km, although the matrices have
	#  to be transposed in order to make it work.
	#
	if (ncurve==1) {
	    kfit <- docurve(time, status, state2, weights, id, nstate)
	    kfit$type <- 'right'
	    kfit$call <- 'call'
	    kfit$n <- n
	    }
	else {
	    timex  <- vector('list', ncurve)
	    n.risk <- vector('list', ncurve)
	    surv   <- vector('list', ncurve)
	    n.event<- vector('list', ncurve)
	    strata <- integer(ncurve)
	    if (se.fit) varhaz <- vector('list', ncurve)
	    uniquex <- levels(X)
	    for (i in 1:ncurve) {
		who <- (x== uniquex[i])
		temp <- docurve(time[who], status[who], state2[who],
				weights[who], id[who], nstate)
		timex[[i]]   <- temp$time
		n.risk[[i]]  <- temp$n.risk
		n.event[[i]] <- temp$n.event
		surv[[i]]    <- c(t(temp$surv))
		strata[i]    <- length(temp$time)
		if (se.fit) varhaz[[i]] <- c(t(temp$std.err))
		}
	    names(strata) <- uniquex
	    kfit <- list(n=table(X),
			 time= unlist(timex),
			 n.risk= unlist(n.risk), 
			 n.event=unlist(n.nevent),
			 surv = t(matrix(unlist(surv), nrow=nstate)),
			 type='right', call=call)
	    if (se.fit) kfit$std.err <- t(matrix(unlist(varhaz),nrow=nstate))
	    }
	dimnames(kfit$surv) <- list(NULL, levels(as.factor(state)))
        kfit$absorb <- levels(state)[-not.absorb]
        }	

    #	
    # Last bit: add in the confidence bands (also stolen from survfit.km)
    #
    if (se.fit) {
	std.err <- kfit$std.err
	#
	# n.lag = the # at risk the last time there was an event (or
	#   the first time of a strata)
	#
	events <- kfit$n.event >0
	if (ncurve==1) events[1] <- 1
	else           events[1 + cumsum(c(0, strata[-ncurve]))] <- 1
	zz <- 1:length(events)
	n.lag <- rep(kfit$n.risk[events], diff(c(zz[events], 1+max(zz))))
	std.low <- switch(conf.lower,
			  'usual' = std.err,
			  'peto' = sqrt((1-kfit$surv)/ kfit$n.risk),
			  'modified' = std.err * sqrt(n.lag/kfit$n.risk))
	zval <- qnorm(1- (1-conf.int)/2, 0,1)

	if (conf.type=='plain') {
	    temp1 <- kfit$surv + zval* std.err * kfit$surv
	    temp2 <- kfit$surv - zval* std.low * kfit$surv
	    kfit <- c(kfit, list(upper=pmin(temp1,1), lower=pmax(temp2,0),
				 conf.type='plain', conf.int=conf.int))
	    }

	if (conf.type=='log') {
	    #avoid some "log(0)" messages
	    xx <- ifelse(kfit$surv==0,1,kfit$surv)  

	    temp1 <- ifelse(kfit$surv==0, NA, exp(log(xx) + zval* std.err))
	    temp2 <- ifelse(kfit$surv==0, NA, exp(log(xx) - zval* std.low))
	    kfit <- c(kfit, list(upper=pmin(temp1,1), lower=temp2,
				 conf.type='log', conf.int=conf.int))
	    }

	if (conf.type=='log-log') {
	    who <- (kfit$surv==0 | kfit$surv==1) #special cases
	    temp3 <- ifelse(kfit$surv==0, NA, 1)
	    xx <- ifelse(who, .1,kfit$surv)  #avoid some "log(0)" messages
	    temp1 <- exp(-exp(log(-log(xx)) + zval*std.err/log(xx)))
	    temp1 <- ifelse(who, temp3, temp1)
	    temp2 <- exp(-exp(log(-log(xx)) - zval*std.low/log(xx)))
	    temp2 <- ifelse(who, temp3, temp2)
	    kfit <- c(kfit, list(upper=temp1, lower=temp2,
				 conf.type='log-log', conf.int=conf.int))
	    }
        }
    
    kfit$call <- call
    oldClass(kfit) <- 'survfit.ci'
    kfit
    }

