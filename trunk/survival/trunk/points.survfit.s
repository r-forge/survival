#SCCS $Date: 1992-04-14 18:07:12 $ $Id: points.survfit.s,v 4.2 1992-04-14 18:07:12 grill Exp $
points.survfit <- function(object, ...) {
    if (!is.matrix(object$surv))
	    points(object$time, object$surv, ...)
    else
	    matpoints(object$time, object$surv, ...)
    }
