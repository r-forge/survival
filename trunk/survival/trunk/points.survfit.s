#SCCS $Date: 1992-03-04 16:48:13 $ $Id: points.survfit.s,v 4.1 1992-03-04 16:48:13 therneau Exp $
points.surv.fit <- function(object, ...) {
    if (!is.matrix(object$surv))
	    points(object$time, object$surv, ...)
    else
	    matpoints(object$time, object$surv, ...)
    }
