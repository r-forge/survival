residuals.survreg <-
function(object, type = c("deviance", "pearson", "working"))
{
    type <- match.arg(type)
    rr <- switch(type,
	    working = object$residuals,
	    pearson = sqrt(object$weights) * object$residuals,
	    deviance = object$dresiduals)

    #Expand out the missing values in the result
    if (!is.null(object$na.action))
	 naresid(object$na.action, rr)
    else rr
    }
