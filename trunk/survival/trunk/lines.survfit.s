# SCCS $Id: lines.survfit.s,v 4.3 1992-04-13 22:05:35 therneau Exp $
lines.surv.fit <- function(object, type='s', ...) {
    if (inherits(object, 'surv.exp') && missing(type)) type <- 'l'
    if (!is.matrix(object$surv))
	lines(object$time, object$surv, type=type, ...)
    else {
	for (i in 1:ncol(object$surv))
	    lines(object$time, object$surv[,i], type=type, ...)
	}
    }
