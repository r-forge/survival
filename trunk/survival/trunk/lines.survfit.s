# SCCS $Id: lines.survfit.s,v 4.2 1992-03-24 09:31:08 therneau Exp $
lines.surv.fit <- function(object, type='s', ...) {
    if (!is.matrix(object$surv))
	lines(object$time, object$surv, type=type, ...)
    else {
	for (i in 1:ncol(object$surv))
	    lines(object$time, object$surv[,i], type=type, ...)
	}
    }
