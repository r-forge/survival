# SCCS  $Id: survobrien.s,v 4.2 1992-03-09 01:33:11 therneau Exp $
#
# The test for survival proposed by Peter O'Brien
#
surv.obrien <- function(formula, data= sys.parent()) {
    m <- model.frame(formula, data, na.action= function(x) x )
    n <- nrow(m)
    Terms <- attr(m, 'terms')

    y <- model.extract(m, 'response')
    if (!inherits(y, "Surv")) stop ("Response must be a survival object")
    if (attr(y, 'type') != 'right') stop("Can only handle right censored data")

    # Figure out which are the continuous predictor variables
    m.name <- names(m)
    temp <- match(attr(Terms, 'term.labels'), m.name)
    cont <- NULL
    for (i in temp) {if (!is.factor(m[[i]])) cont <- c(cont, i)}
    if (is.null(cont)) stop ("No continuous variables to modify")

    keepers <- rep(T, length(m))  #The ones to be kept "as is"
    keepers[cont] <- F
    keepers[attr(Terms, 'response')] <- F
    ord <- order(y[,1])
    x <- as.matrix(m[ord, cont])
    time <- y[ord,1]
    status <- y[ord,2]
    nvar <- length(cont)

    nline <- 0
    for (i in unique(time[status==1])) nline <- nline + sum(time >=i)
    start <- stop <- event <- double(nline)
    xx <- matrix(double(nline*nvar), ncol=nvar)
    ltime <- 0
    j<- 1
    keep.index <- NULL
    for (i in unique(time[status==1])) {
	who <- (time >=i)
	nrisk <- sum(who)
	if (nrisk<2) break
	temp <- apply(x[who,,drop=F], 2, rank)
	temp <- (2*temp -1)/ (2* nrisk)   #percentiles
	logit<- log(temp/(1-temp))           #logits
	deaths <- (status[who]==1 & time[who]==i)

	k <- seq(from=j, length=nrisk)
	start[k] <- ltime
	stop[k] <-  i
	event[k] <- deaths
	xx[k,] <- logit
	j <- j + nrisk
	ltime <- i
	keep.index <- c(keep.index, ord[who])
	}

    if (any(keepers)) {
	 temp <- list(m[keep.index, keepers], start, stop, event, xx)
	 names(temp) <- c(m.name[keepers], "start", "stop", "event",
				m.name[cont])
	 }
    else {
	temp <- list(start, stop, event, xx)
	names(temp) <- c(m.name[keepers], "start", "stop", "event",
				m.name[cont])
	}
    temp
    }
