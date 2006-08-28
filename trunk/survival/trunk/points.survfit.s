# $Date: 2006-08-28 14:18:14 $ $Id: points.survfit.s,v 4.3 2006-08-28 14:18:14 m015733 Exp $
points.survfit <- function(object, ...) {
    if (!is.matrix(object$surv))
	    points(object$time, object$surv, ...)
    else
	    matpoints(object$time, object$surv, ...)
    }
