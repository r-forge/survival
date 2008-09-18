# $Date: 2006-08-28 14:18:14 $ $Id$
points.survfit <- function(object, ...) {
    if (!is.matrix(object$surv))
	    points(object$time, object$surv, ...)
    else
	    matpoints(object$time, object$surv, ...)
    }
