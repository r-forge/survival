#SCCS $Id: print.summary.survfit.s,v 4.5 2000-03-02 20:06:49 boos Exp $
print.summary.survfit <- function(x, 
				  digits = max(options()$digits - 4, 3), ...) {
    savedig <- options(digits=digits)
    on.exit(options(savedig))

    if (!is.null(cl<- x$call)) {
	cat("Call: ")
	dput(cl)
	cat("\n")
	}

    omit <- x$na.action
    if (length(omit)) 
	    cat(naprint(omit), "\n")

    if (x$type == 'right') {
	mat <- cbind(x$time, x$n.risk, x$n.event, x$surv)
	if (x$cumm == 1) 
		cnames <- c("time", "n.risk", "Cumm. n.event")
	else
		cnames <- c("time", "n.risk", "n.event")
        }

    if (x$type == 'counting') {
	mat <- cbind(x$time, x$n.risk, x$n.event, x$n.entered,
		     x$n.exit.censored, x$surv)
	if (x$cumm == 1) { 
	    if (is.Surv(x)) {
		cnames <- c("time", "n.risk", "Cumm. n.event", 
			    "Cumm. n.entered", "Cumm. n.censored")
	        }	
	    else {
		cnames <- c("time", "n.risk", "Cumm. n.event") 
	        }
	    }	    
	else {
	    if (is.Surv(x)) {
		cnames <- c("time", "n.risk", "n.event", 
			    "n.entered", "n.censored")
	       }
	    else {
		cnames <- c("time", "n.risk", "n.event")
	        }
	    }
        }
    if (is.matrix(x$surv)) ncurve <- ncol(x$surv)
    else	           ncurve <- 1
    if (ncurve==1) {                 #only 1 curve
	cnames <- c(cnames, "survival")
	if (!is.null(x$std.err)) {
	    if (is.null(x$lower)) {
		mat <- cbind(mat, x$std.err)
		cnames <- c(cnames, "std.err")
	        }
	    else {
		mat <- cbind(mat, x$std.err, x$lower, x$upper)
		cnames <- c(cnames, 'std.err',
			    paste("lower ", x$conf.int*100, "% CI", sep=''),
			    paste("upper ", x$conf.int*100, "% CI", sep=''))
	        }	
	    }
        }
    else cnames <- c(cnames, paste("survival", seq(ncurve), sep=''))

    if (!is.null(x$new.start)) {
	mat.keep <- mat[,1] >= x$new.start
	mat <- mat[mat.keep,,drop=F]
	if (is.null(dim(mat)))
		stop(paste("No information available using new.start =", x$new.start, "."))
        }
    if (!is.matrix(mat)) mat <- matrix(mat, nrow=1)
    if (!is.null(mat)) {
	dimnames(mat) <- list(NULL, cnames)
	if (is.null(x$strata))
		prmatrix(mat, rowlab=rep("", nrow(mat)))
	else  { #print it out one strata at a time
	    if (!is.null(x$times.strata))
		    strata <- x$times.strata
	    else
		    strata <- x$strata
	   
	    if (!is.null(x$new.start))
		    strata <- strata[mat.keep]
	    for (i in levels(strata)) {
		who <- (strata==i)
		cat("               ", i, "\n")
		if (sum(who) ==1)
			print(mat[who,])
	        else
		    prmatrix(mat[who,], rowlab=rep("", sum(who)))

		cat("\n")
 	        }
	    }
        }
    else 
	stop("There are no events to print.  Please use the option ",
	    "censored=T with the summary function to see the censored ",
	    "observations.")
    invisible(x)
    }














