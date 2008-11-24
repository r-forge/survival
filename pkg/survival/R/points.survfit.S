# $Id$
points.survfit <- function(x, ...) {
    if (!is.matrix(x$surv))
	    points(x$time, x$surv, ...)
    else
	    matpoints(x$time, x$surv, ...)
    }
