# SCCS  $Id: summary.aareg.s,v 1.1 2001-04-17 08:23:45 therneau Exp $
summary.aareg <- function(x, maxtime, weight=c('var', 'nrisk')) {
    if (!inherits(x, 'aareg')) stop ("Must be an aareg object")

    weight <- match.arg(weight)
    if (!missing(maxtime)) ntime <- sum(x$time < maxtime)
    else 		   ntime <- x$ntime[1]

    # Compute the test statistic
    if (weight=='var') twt <- (as.matrix(x$tweight))[1:ntime,]
    else	       twt <- 1/x$nrisk[1:ntime]
    
    tx <- as.matrix(twt * x$increment[1:ntime,])
    test <- c(rep(1,ntime) %*% tx)  #sum of weighted increments
    var  <- t(tx) %*% tx	    #std Poisson, ind increments variance

    # Compute a fit (to get a "slope") 
    #  Since this is a single variable model, no intercept, I
    #  don't need to call lm.wfit!
    time <- x$time[1:ntime]
    nvar <- ncol(tx)

    ctx  <- apply(tx, 2, cumsum)

    if (is.matrix(twt) && ncol(twt) >1) tempwt <- apply(twt*time^2, 2, sum)
    else    			        tempwt <- sum(twt*time^2)
    if (ncol(ctx) >1)
	    slope<- apply(ctx* time, 2, sum)/ tempwt
    else    slope <- sum(ctx*time) / temp2
    se  <- sqrt(diag(var))
    mat <- cbind(slope, test, se, test/se, 2*pnorm(-abs(test/se)))
    dimnames(mat) <- list((dimnames(x$inc)[[2]]),
			  c("slope", "test", "se(test)", "z", "p"))
    chi <- test%*%solve(var,test) 
    temp <- list(table=mat, test=test, var=var, chisq=chi, 
		 ntime=c(ntime, x$ntime[2]))
    oldClass(temp) <- 'summary.aareg'
    temp
    }

print.summary.aareg <- function(x) {
    print(signif(x$table,3))
    chi <- x$chisq
    cat("\nChisq=", format(round(chi,2)), "on", ncol(x$increment), "df, p=",
	            signif(1- pchisq(chi, nrow(x$table)),2), "\n")
    invisible(x$table)
    }
