#SCCS $Id: is.ratetable.s,v 4.2 1995-01-31 15:00:01 therneau Exp $
is.ratetable <- function(x) {
    if (!inherits(x, 'ratetable')) return(F)
    att <- attributes(x)
    if (any(is.na(match(c("dim", "dimnames", "dimid", "factor", "cutpoints"),
			names(att))))) return(F)
    nd <- length(att$dim)
    if (length(x) != prod(att$dim)) return(F)
    if (!(is.list(att$dimnames) && is.list(att$cutpoints))) return(F)
    if (length(att$dimnames)!=nd || length(att$factor)!=nd ||
		     length(att$cutpoints)!=nd) return(F)
    fac <- as.numeric(att$factor)
    if (any(is.na(fac))) return(F)
    if (any(fac <0)) return(F)
    for (i in 1:nd) {
	n <- att$dim[i]
	if (length(att$dimnames[[i]]) !=n) return(F)
	if (fac[i]!=1 && length(att$cutpoints[[i]]) != n)  return(F)
	if (fac[i]!=1 && any(order(att$cutpoints[[i]])!= 1:n)) return(F)
	if (fac[i]==1 && !is.null(att$cutpoints[[i]]))  return(F)
	}
    T
    }
