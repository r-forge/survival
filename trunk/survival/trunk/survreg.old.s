# $Id: survreg.old.s,v 1.5 2006-08-28 18:23:13 m015733 Exp $
# Map the argument list of the old survreg to the new one
#
survreg.old <- function(formula, data=sys.parent(), ...,
        link='log',
        dist=c("extreme", "logistic", "gaussian", "exponential",
               "rayleigh", "weibull"),
	fixed=list()) {
    
    dist <- match.arg(dist)

    if ((dist!='weibull' && dist != 'rayleigh') && link=='log') {
	if (dist=='extreme') dist <- 'weibull'
	else dist <- paste('log', dist, sep='')
	}
    if (is.null(fixed$scale)) scale <- 0
    else scale <- fixed$scale

    survreg(formula, data, ..., dist=dist, scale=scale)
    }
