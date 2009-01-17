#
# $Id$
#
is.ratetable <- function(x, verbose=FALSE) {
    dlist <- c("dim", "dimnames", "dimid", "cutpoints")
    if (!verbose) {
	if (!inherits(x, 'ratetable')) return(FALSE)
	att <- attributes(x)
	if (any(is.na(match(dlist, names(att))))) return(FALSE)
	nd <- length(att$dim)
	if (length(x) != prod(att$dim)) return(FALSE)
	if (!(is.list(att$dimnames) && is.list(att$cutpoints)))
		 return(FALSE)
	if (length(att$dimnames)!=nd || length(att$factor)!=nd ||
			 length(att$cutpoints)!=nd) return(FALSE)
        # One of 'factor' (old style table) or 'type' (new style) should exist
        if (!is.null(att$factor)) {
            fac <- as.numeric(att$factor)
            if (any(is.na(fac))) return(FALSE)
            if (any(fac <0)) return(FALSE)
            }
        else if (!is.null(att$type)) {
            if (any(is.na(match(att$type, 0:3)))) return(FALSE)
            }
        else return(FALSE)

        if (length(att$dimid) != nd) return(FALSE)
	for (i in 1:nd) {
	    n <- att$dim[i]
	    if (length(att$dimnames[[i]]) !=n) return(FALSE)
	    if (fac[i]!=1 && length(att$cutpoints[[i]])!=n) return(FALSE)
	    if (fac[i]!=1 && any(order(att$cutpoints[[i]])!= 1:n)) return(FALSE)
	    if (fac[i]==1 && !is.null(att$cutpoints[[i]]))  return(FALSE)
	    if (fac[i]>1 && i<nd) return(FALSE)
	    }
	return(TRUE)
	}

    # verbose return messages, useful for debugging
    msg <- NULL
    if (!inherits(x, 'ratetable')) msg <- c(msg, "wrong class")
    att <- attributes(x)

    temp <- is.na(match(dlist, names(att)))
    if (any(temp)) 
        msg <- c(msg, paste("missing attribute:", dlist[temp]))

    # This next error is unlikely, since S itself squawks when you
    #   try to set a wrong dimension.  Ditto with dimnames issues.
    nd <- length(att$dim)
    if (length(x) != prod(att$dim)) 
        msg <- c(msg, 'length of the data does not match prod(dim)')

    if (!is.list(att$dimnames))
	     msg <- c(msg, 'dimnames is not a list')
    if (!is.list(att$cutpoints))
	     msg <- c(msg, 'cutpoints is not a list')

    if (length(att$dimnames)!=nd)
        msg <- c(msg, 'wrong length for dimnames')
    if (length(att$dimid)!=nd)
        msg <- c(msg, 'wrong length for dimid')
    if (length(att$factor)!=nd)
        msg <- c(msg, 'wrong length for factor')
    if (length(att$cutpoints)!=nd) 
        msg <- c(msg, 'wrong length for cutpoints')

    if (!is.na(att$factor)) {
        fac <- as.numeric(att$factor)
        if (any(is.na(fac))) msg <- c(msg, "illegal 'factor' level of NA")
        if (any(fac <0)) msg <- c(msg, 'factor <0')
        }
    else if (!is.na(att$type)) {
        if (any(is.na(match(att$type, 0:3))))
            msg <- c(msg, 'type attribute must be 0, 1, 2, or 3')
        }
    else msg <- c(msg, "missing the 'type' attribute")

    for (i in 1:nd) {
	n <- att$dim[i]
	if (length(att$dimnames[[i]]) !=n) 
		msg <- c(msg, paste('dimname', i, 'is the wrong length'))

	if (fac[i]!=1) {
            if (length(att$cutpoints[[i]]) != n)
                msg <- c(msg, paste('wrong length for cutpoints', i))
            else if (any(order(att$cutpoints[[i]])!= 1:n)) 
		msg <- c(msg, paste('unsorted cutpoints for dimension',i))
                }

	if (fac[i]==1 && !is.null(att$cutpoints[[i]]))  
		msg <- c(msg, paste('factor', i, 
                                    'is 0; cutpoint should be null'))
	if (fac[i]>1 && i<nd) 
		msg <- c(msg, 'only the last dimension can be interpolated')
	}
    if (length(msg)==0) TRUE
    else msg
    }
