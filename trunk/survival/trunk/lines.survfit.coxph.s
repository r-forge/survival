# $Id: lines.survfit.coxph.s,v 1.1 1997-04-23 16:28:16 therneau Exp $
lines.survfit.coxph <- function(x, mark.time=F, ...) {
    if (is.logical(mark.time) & mark.time)
	    stop("Invalid value for mark.time")
    invisible(NextMethod('lines', mark.time=mark.time))
    }
