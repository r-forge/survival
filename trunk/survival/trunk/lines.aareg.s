# SCCS $Id: lines.aareg.s,v 1.1 2002-01-29 12:47:04 therneau Exp $
lines.aareg <- function(x, se=F, maxtime, ...) {
    if (!inherits(x, 'aareg')) stop ("Must be an aareg object")

    if (missing(maxtime)) keep <- 1:length(x$time)
    else		  keep <- 1:sum(x$time <= maxtime)

    if (is.matrix(x$increment)) yylab <- dimnames(x$increment)[[2]]
    else                        yylab <- ""
    if (is.matrix(x$increment) && ncol(x$increment)>1) {
	cum <- apply(x$increment[keep,], 2,cumsum)
	cum <- rbind(0,cum)
	if (se) {
	    se.cum <- sqrt(apply(x$increment[keep,]^2, 2,cumsum)) 
	    se.cum <- rbind(0, se.cum)
	    }
	ncurve <- ncol(cum)
	}
    else {
	cum <- cumsum(c(0, x$increment[keep]))
	se.cum <- sqrt(cumsum(c(0, x$increment[keep]^2)))
	ncurve <- 1
	}
   
    xx <- c(0, x$time[keep])
    if (se) { 
	yy <- cbind(cum, cum + 1.96*se.cum,
			 cum - 1.96*se.cum)
	if (ncurve >1) {
	    for (i in 1:ncurve) {
		j <- c(i, i+ncurve, i+2*ncurve)
		matlines(xx, yy[,j], ..., type='s', col=1, lty=c(1,2,2))
		}
	    }
	else matlines(xx, yy, ..., type='s', col=1, lty=c(1,2,2)) 
	}
    else {
	matlines(xx, cum, ..., type='s')
	}
    }
