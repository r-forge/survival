# SCCS $Id: lines.survfit.s,v 4.4 1992-04-14 18:06:43 grill Exp $
lines.survfit <- function(object, type='s', ...) {
    if (inherits(object, 'survexp') && missing(type)) type <- 'l'
    if (!is.matrix(object$surv))
	lines(object$time, object$surv, type=type, ...)
    else {
	for (i in 1:ncol(object$surv))
	    lines(object$time, object$surv[,i], type=type, ...)
	}
    }
