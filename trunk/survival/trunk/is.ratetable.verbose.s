#
# SCCS $Id: is.ratetable.verbose.s,v 1.1 1995-10-27 12:17:47 therneau Exp $
#   A version of the function that tells you WHY it's not a ratetable
#

is.ratetable.verbose <- function(x) {
    if (!inherits(x, 'ratetable')) return("wrong class")
    att <- attributes(x)
    if (any(is.na(match(c("dim", "dimnames", "dimid", "factor", "cutpoints"),
			names(att))))) return('missing an attribute')
    nd <- length(att$dim)
    if (length(x) != prod(att$dim)) return('dims dont match length')
    if (!(is.list(att$dimnames) && is.list(att$cutpoints)))
	     return('dimnames or cutpoints not a list')
    if (length(att$dimnames)!=nd || length(att$factor)!=nd ||
		     length(att$cutpoints)!=nd) return('bad lengths')
    fac <- as.numeric(att$factor)
    if (any(is.na(fac))) return('a missing factor')
    if (any(fac <0)) return('factor <0')
    for (i in 1:nd) {
	n <- att$dim[i]
	if (length(att$dimnames[[i]]) !=n) return('dimnames wrong length')
	if (fac[i]!=1 && length(att$cutpoints[[i]])!=n) return('cutpnt missing')
	if (fac[i]!=1 && any(order(att$cutpoints[[i]])!= 1:n)) return('unsorted cutpoints')
	if (fac[i]==1 && !is.null(att$cutpoints[[i]]))  return('cutpnt should be null')
	if (fac[i]>1 && i<nd) return('only the last dim can be interpolated')
	}
    T
    }
