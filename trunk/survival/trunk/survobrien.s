#SCCS $Date: 1992-03-04 16:48:38 $ $Id: survobrien.s,v 4.1 1992-03-04 16:48:38 therneau Exp $
# SCCS  $Id: survobrien.s,v 4.1 1992-03-04 16:48:38 therneau Exp $
#
# The test for survival proposed by Peter O'Brien
#
surv.obrien <- function(formula) {
    m <- model.frame(formula, na.action= function(x) x )
    n <- nrow(m)


    num <- rep(0, nvar)
    denom <- num
    nline <- 0
    for (i in unique(time[status==1])) nline <- nline + sum(time >=i)
    start <- stop <- event <- double(nline)
    xx <- matrix(double(nline*nvar), ncol=nvar)
    ltime <- 0
    j<- 1

    for (i in unique(time[status==1])) {
	who <- (time >=i)
	nrisk <- sum(who)
	if (nrisk<2) break
	temp <- apply(x[who,], 2, rank)
	temp <- (2*temp -1)/ (2* nrisk)   #percentiles
	logit<- log(temp/(1-temp))           #logits
	deaths <- (status[who]==1 & time[who]==i)
     # The mean logit is deterministicly=0, hence the (n-1)/n adjustment
     #    (If there are ties it may not be zero).
	denom <- denom + sum(deaths)*apply(logit,2,var)*(nrisk-1)/nrisk
	if (sum(deaths) >1) num <- num + apply(logit[deaths,], 2, sum)
	else                num <- num + logit[deaths,]

	#Multivariate case
	k <- seq(from=j, length=nrisk)
	start[k] <- ltime
	stop[k] <-  i
	event[k] <- deaths
	xx[k,] <- logit
	j <- j + nrisk
	ltime <- i
	}
    dimnames(xx) <- list(NULL, dimnames(x)[[2]])
    temp <- agreg(start, stop, event, xx)
    temp <- list(num=num, denom=denom, agfit=temp)
    attr(temp, "class") <- 'surv.obrien'
    temp
    }
