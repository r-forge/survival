residuals.survreg <-
function(object, type = c("deviance", "pearson", "working", "matrix"))
{
    type <- match.arg(type)
    rr <- switch(type,
	    working = object$residuals,
	    pearson = sqrt(object$weights) * object$residuals,
	    deviance = object$dresiduals,
	    matrix= {
		eta <- object$linear.predictors
		n   <- length(eta)
		y   <- object$y
		dist<- match(object$dist, c("extreme", "logistic",
					    "gaussian", "cauchy"))
		temp <-.C("survreg_g", as.integer(n),
				as.double(y),
				as.integer(ncol(y)),
				eta,
				as.double(object$parms),
				deriv=matrix(double(n*6), ncol=6),
				as.integer(6),
				as.integer(dist))$deriv
		dimnames(temp) <- list(names(object$residuals),
				       c("loglik", "eta'", "eta''", "sig'",
					 "sig''", "eta'sig'"))
		temp
		}
	    )

    #Expand out the missing values in the result
    if (!is.null(object$na.action))
	 naresid(object$na.action, rr)
    else rr
    }
