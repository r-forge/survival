#SCCS $Date: 1992-03-04 16:48:05 $ $Id: lines.survfit.s,v 4.1 1992-03-04 16:48:05 therneau Exp $
lines.surv.fit <- function(object, ...) {
    if (!is.matrix(object$surv))
	lines(stepfun(object$time, object$surv), ...)
    else {
	for (i in 1:ncol(object$surv))
	    lines(stepfun(object$time, object$surv[,i]), ...)
	}
    }
