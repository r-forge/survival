# SCCS $Id: print.aareg.s,v 1.2 2002-04-29 14:15:11 therneau Exp $
print.aareg <- function(x, maxtime, test=c('aalen', 'nrisk')) {
    if (!inherits(x, 'aareg')) stop ("Must be an addreg object")
    if (!is.null(cl<- x$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}

    if (missing(test)) test <- x$test
    else test <- match.arg(test)

    if (missing(maxtime)) summ <- summary(x, test=test)
    else                  summ <- summary(x, maxtime=maxtime, test=test)

    omit <- x$na.action
    if (length(omit))
	cat("  n=", x$n[1], " (", naprint(omit), ")\n", sep="")
    else cat("  n=", x$n[1], "\n")
    cat("   ", summ$n[2], "out of", x$n[3], "unique event times used\n\n")
    print(signif(summ$table,3))
    chi <- summ$chisq
    df <- nrow(summ$table) -1
    cat("\nChisq=", format(round(chi,2)), " on ", df, " df, p=",
	            signif(1- pchisq(chi,df),2), 
	            "; test weights=", x$test, "\n", sep="")
    invisible(x)
    }

