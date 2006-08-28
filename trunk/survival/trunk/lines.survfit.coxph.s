# $Id: lines.survfit.coxph.s,v 1.2 2006-08-28 14:01:23 m015733 Exp $
lines.survfit.coxph <- function(x, mark.time=FALSE, ...) {
    if (is.logical(mark.time) & mark.time)
	    stop("Invalid value for mark.time")
    invisible(NextMethod('lines', mark.time=mark.time))
    }
