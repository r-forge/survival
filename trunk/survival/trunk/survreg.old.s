#  SCCS $Id: survreg.old.s,v 1.2 1999-02-07 23:36:53 therneau Exp $
# Map the argument list of the old survreg to the new one
#
survreg.old <- function(formula, data=sys.parent(), ...,
        link='log',
        dist=c("extreme", "logistic", "gaussian", "exponential",
               "rayleigh"),
	fixed=list()) {
    
    if (link=='log') {
	if (dist=='extreme') dist <- 'weibull'
	else dist <- paste('log', dist, sep='')
	}
    if (is.null(fixed$scale)) scale <- 0
    else scale <- fixed$scale

    survreg(formula, data, ..., dist=dist, scale=scale)
    }
