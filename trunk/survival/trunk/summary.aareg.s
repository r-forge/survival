# SCCS $Id: summary.aareg.s,v 1.2 2002-04-29 14:15:12 therneau Exp $
# The summary routine for aareg models.
# A lot of the work below relates to one particular issue: the coeffients
#  of an aareg model often get "wild" near the end (at the largest times).
#  So, a common case is to 
#       fit the model (very slow)
#       look at the printout -- Hmmm x1 is significant, x2 not, ...., why?
#       look at plot(fit) and
#           oh my gosh, I should have cut the time scale off at 520 days
#
# This routine allows one to do that.  If maxtime is given, the overall
#    test statistic is re-computed.  One consequence is that lots of the
#    intermediate material from the fit had to be included in the aareg
#    object.
# The "variance" based weighting for a test is not allowed, because it would
#    have meant an awful more stuff to pass, lots more work, for a test
#    that is rarely used.
#
summary.aareg <- function(x, maxtime, test=c('aalen', 'nrisk')) {
    if (!inherits(x, 'aareg')) stop ("Must be an aareg object")

    test <- match.arg(test)

    # Compute a "slope" for each line, using appropriate weighting
    #  Since this is a single variable model, no intercept, I
    #  don't need to call lm.wfit!
    if (!missing(maxtime)) ntime <- sum(x$time <= maxtime)
    else 		   ntime <- nrow(x$coefficient)
    times <- x$time[1:ntime]

    if (missing(test)) test <- x$test

    if (test=='aalen') twt <- (as.matrix(x$tweight))[1:ntime,]
    else	       twt <-  x$nrisk[1:ntime]
    tx <- as.matrix(twt * x$coefficient[1:ntime,])
    ctx  <- apply(tx, 2, cumsum)

    if (is.matrix(twt) && ncol(twt) >1) tempwt <- apply(twt*times^2, 2, sum)
    else    			        tempwt <- sum(twt*times^2)
    if (ncol(ctx) >1)
	    slope<- apply(ctx* times, 2, sum)/ tempwt
    else    slope <- sum(ctx*times) / temp2

    if (!missing(maxtime) || x$test != test) {
	# Compute the test statistic
	test.stat <- apply(tx, 2, sum)  #sum of tested coefficients
	test.var  <- t(tx) %*% tx	#std Poisson, ind coefficients variance
	
	if (!is.null(x$dfbeta)) {
	    dd <- dim(x$dfbeta)
	    indx <- match(unique(times), x$times)
	    influence <- matrix(0, dd[1], dd[2])
	    for (i in 1:length(indx)) {
		if (test=='aalen') 
		     influence <- influence + x$dfbeta[,,i] %*% 
			                            diag(twt[indx[i],])
		else influence <- influence + x$dfbeta[,,i]* x$nrisk[indx[i]]
		}
	    if (!is.null(x$cluster)) influence <- rowsum(influence, cluster)
	    test.var2 <- t(influence) %*% influence
	    }
	}
    else {  #use the value that was passed in
	test.stat <- x$test.statistic
	test.var  <- x$test.var
	test.var2 <- x$test.var2
	}

    # create the matrix for printing out
    #  The chisquare test does not include the intercept
    se1 <- sqrt(diag(test.var))
    if (is.null(test.var2)) {
	mat <- cbind(slope, test.stat, se1, test.stat/se1, 
		     2*pnorm(-abs(test.stat/se1)))
	dimnames(mat) <- list((dimnames(x$coefficient)[[2]]),
			      c("slope", "test", "se(test)", "z", "p"))
	chi <- test.stat[-1] %*% solve(test.var[-1,-1],test.stat[-1]) 
	}
    else {
	se2 <- sqrt(diag(test.var2))
	mat <- cbind(slope, test.stat, se1, se2, test.stat/se2, 
		                    2*pnorm(-abs(test.stat/se2)))
	dimnames(mat) <- list((dimnames(x$coefficient)[[2]]),
			 c("slope", "test", "se(test)", "robust se", "z", "p"))
	chi <- test.stat[-1] %*% solve(test.var2[-1,-1], test.stat[-1]) 
	}

    temp <- list(table=mat, test=test, test.statistic=test.stat, 
		 test.var=test.var, test.var2=test.var2, chisq=chi,
		 n = c(x$n[1], length(unique(times)), x$n[3]))

    
    oldClass(temp) <- 'summary.aareg'
    temp
    }

print.summary.aareg <- function(x) {
    print(signif(x$table,3))
    chi <- x$chisq
    df  <- length(x$test.statistic) -1
    cat("\nChisq=", format(round(chi,2)), " on ", df, " df, p=",
	            signif(1- pchisq(chi, df),2), 
	            "; test weights=", x$test, "\n", sep='')
    invisible(x$table)
    }
