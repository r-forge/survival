# SCCS $Id: print.aareg.s,v 1.1 2001-04-17 08:23:43 therneau Exp $
print.aareg <- function(x, maxtime, weight=c('var', 'nrisk')) {
    if (!inherits(x, 'aareg')) stop ("Must be an addreg object")
    if (!is.null(cl<- x$call)) {
	cat("Call:\n")
	dput(cl)
	cat("\n")
	}
    weight <- match.arg(weight)

    if (missing(maxtime)) summ <- summary(x, weight=weight)
    else                  summ <- summary(x, maxtime=maxtime, weight=weight)

    cat(x$ntime[1], "out of", x$ntime[2], "unique death times used\n\n")
    print(signif(summ$table,3))
    chi <- summ$chisq
    cat("\nChisq=", format(round(chi,2)), "on", nrow(summ$table), "df, p=",
	            signif(1- pchisq(chi, nrow(summ$table)),2), "\n")
    invisible(x)
    }

