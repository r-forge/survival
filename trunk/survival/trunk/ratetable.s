#SCCS $Id: ratetable.s,v 4.2 1993-12-15 13:31:22 therneau Exp $
#
# This is a 'specials' function for pyears
#   it is a stripped down version of as.matrix(data.frame(...))
# There is no function to create a ratetable.
# This function has a class, only so that data frame subscripting will work
#
ratetable <- function(...) {
    args <- list(...)
    nargs <- length(args)
    ll <- 1:nargs
    for (i in ll) ll[i] <- length(args[[i]])
    n <- max(ll)
    levlist <- vector("list", nargs)
    x <- matrix(0,n,nargs)
    dimnames(x) <- list(1:n, names(args))
    for (i in 1:nargs) {
	if (ll[i] ==n) {
	    if (!is.numeric(args[[i]])) args[[i]] <- factor(args[[i]])
	    if (is.factor(args[[i]])) {
		levlist[[i]] <- levels(args[[i]])
		x[,i] <- c(args[[i]])
		}
	    else x[,i] <- args[[i]]
	    }
	else if (ll[i] ==1) levlist[i] <- args[i]
	else stop("All arguments to ratetable() must be the same length")
	}
    attr(x, "constants") <- (ll==1) & (n>1)
    attr(x, "levlist")   <- levlist
    attr(x, "class")  <- "ratetable2"
    x
    }

# The function below should only be called internally, when missing
#   values cause model.frame to drop some rows
"[.ratetable2" <- function(x, rows, cols, drop=F) {
    if (!missing(cols)) {
       stop("This should never be called!")
       }
    aa <- attributes(x)
    attributes(x) <- aa[c("dim", "dimnames")]
    y <- x[rows,,drop=F]
    attr(y,'constants') <- aa$constants
    attr(y,'levlist')   <- aa$levlist
    class(y) <- aa$class
    y
    }

#
# Functions to manipulate rate tables
#

"[.ratetable" <- function(x, ..., drop=T) {
    aa <- attributes(x)
    attributes(x) <- aa[c("dim", "dimnames")]
    y <- NextMethod("[", drop=F)
    if (drop && any(dropped <- attr(y, 'dim')==1)){
	if (all(dropped)) as.numeric(y)   #single element
	else {
	    attributes(y) <- list( dim = dim(y)[!dropped],
				   dimnames = dimnames(y)[!dropped],
				   dimid = aa$dimid[!dropped],
				   factor = aa$factor[!dropped],
				   cutpoints =aa$cutpoints[!dropped],
				   class = aa$class)
	    y
	    }
	}
    else {
	attributes(y) <- c(attributes(y), aa[c('dimid','factor','cutpoints',
						'class')])
	y
	}
    }

Math.ratetable <- function(x, ...) {
    attr(x, 'dimid') <- attr(x, 'factor') <- attr(x,'cutpoints') <- NULL
    class(x) <- NULL
    NextMethod(.Generic)
    }

Ops.ratetable <- function(e1, e2) {
    #just treat it as an array
    if (nchar(.Method[1])) {
	attr(e1, 'dimid') <- attr(x, 'factor') <- attr(x,'cutpoints') <- NULL
	class(e1) <- NULL
	}
    if (nchar(.Method[2])) {
	attr(e2, 'dimid') <- attr(x, 'factor') <- attr(x,'cutpoints') <- NULL
	class(e2) <- NULL
	}
    NextMethod(.Generic)
    }

is.ratetable <- function(x) {
    if (!inherits(x, 'ratetable')) return(F)
    att <- attributes(x)
    if (any(is.na(match(c("dim", "dimnames", "dimid", "factor", "cutpoints"),
			names(att))))) return(F)
    nd <- length(att$dim)
    if (!(is.list(att$dimnames) && is.list(att$cutpoints))) return(F)
    if (length(att$dimnames)!=nd || length(att$factor)!=nd ||
		     length(att$cutpoints)!=nd) return(F)
    fac <- as.numeric(att$factor)
    if (any(is.na(fac))) return(F)
    if (any(fac <0)) return(F)
    for (i in 1:nd) {
	n <- att$dim[i]
	if (length(att$dimnames[[i]]) !=n) return(F)
	if (fac[i]!=1 && length(att$cutpoints[[i]])!=n) return(F)
	if (fac[i]==1 && !is.null(att$cutpoints[[i]]))  return(F)
	}
    T
    }
