# $Id$
#
# This source file has two distinct parts in it.  The first is the
#  ratetable(), which is used inside pyears and survexp only to allow
#  users to match the names of variables in their data set to the names
#  of the dimensions in a ratetable.  It returns a matrix with one
#  column for each argument; usually that argument will be a vector but
#  may also be a single constant.  The result has a class "ratetable2",
#  whose only purpose is to allow na.action functions to work properly.
#
# The second part of the file are the methods for actual rate tables, like
#  the table of US survival rates by age and sex (survexp.us).  Rate tables
#  have the "ratetable" class.  However, since each one is rather unique, 
#  there is no function to create a rate table.  Each consists of a multi-way
#  array of event rates along with a set of attributes. 
#
ratetable <- function(...) {
    args <- list(...)
    nargs <- length(args)
    ll <- sapply(args, length)
    n <- max(ll)  # We assume this is the dimension of the user's data frame
    levlist <- vector("list", nargs)
    x <- matrix(0,n,nargs)
    dimnames(x) <- list(1:n, names(args))
    for (i in 1:nargs) {
	if (ll[i] ==n) { # if the arg is a vector of the right length
	    # numerics or factors are allowed. The line below is
	    #  mostly for cut() and tcut(), which create class='category'
	    if (!is.numeric(args[[i]])) args[[i]] <- factor(args[[i]])
	    if (is.factor(args[[i]])) {
		levlist[[i]] <- levels(args[[i]])
		x[,i] <- as.numeric(args[[i]]) # the vector of levels
		}
	    else x[,i] <- args[[i]]
	    }
	else if (ll[i] ==1) {
	    # Treat it like a factor with one level.  The match.ratetable
	    #  function will sort out if it is actually numeric
	    levlist[i] <- args[i]
	    x[,i] <- 1 
	    }
	else stop(paste("Aguments do not all have the same length (arg ",
			i, ")", sep=''))
	}
    attr(x, "constants") <- (ll==1) & (n>1)
    attr(x, "levlist")   <- levlist
    if (is.R()) class(x) <- 'ratetable2'
    else        oldClass(x)  <- "ratetable2"
    x
    }

# The two functions below should only be called internally, when missing
#   values cause model.frame to drop some rows
is.na.ratetable2 <- function(x) {
    attributes(x) <- list(dim=dim(x))
    as.vector((1 * is.na(x)) %*% rep(1, ncol(x)) >0)
    }
"[.ratetable2" <- function(x, rows, cols, drop=FALSE) {
    if (!missing(cols)) {
       stop("This should never be called!")
       }
    aa <- attributes(x)
    attributes(x) <- aa[c("dim", "dimnames")]
    y <- x[rows,,drop=FALSE]
    attr(y,'constants') <- aa$constants
    attr(y,'levlist')   <- aa$levlist
    if (is.R()) class(y) <- 'ratetable2'
    else        oldClass(y) <- 'ratetable2'
    y
    }

#
# Functions to manipulate rate tables
#
"[.ratetable" <- function(x, ..., drop=TRUE) {
    aa <- attributes(x)
    attributes(x) <- aa[c("dim", "dimnames")]
    y <- NextMethod("[", drop=FALSE)
    newdim <- attr(y, 'dim')
    if (is.null(newdim)) stop("Invalid subscript")
    dropped <- (newdim==1)
    if (drop)  change <- (newdim!=aa$dim & !dropped)
    else       change <- (newdim!=aa$dim)

    if (any(change)) {  #dims that got smaller, but not dropped
	newcut <- aa$cutpoints
	for (i in (1:length(change))[change])
	    if (!is.null(newcut[[i]])) newcut[[i]] <-
		(newcut[[i]])[match(dimnames(y)[[i]], aa$dimnames[[i]])]
	aa$cutpoints <- newcut
	}
    if (drop && any(dropped)){
	if (all(dropped)) as.numeric(y)   #single element
	else {
	    #Note that we have to drop the summary function
	    attributes(y) <- list( dim = dim(y)[!dropped],
				   dimnames = dimnames(y)[!dropped],
				   dimid = aa$dimid[!dropped],
				   factor = aa$factor[!dropped],
				   cutpoints =aa$cutpoints[!dropped])
	    if (is.R()) class(y) <- 'ratetable'
	    else        oldClass(y) <- 'ratetable'
	    y
	    }
	}
    else {
	aa$dim <- aa$dimnames <- NULL
	attributes(y) <- c(attributes(y), aa)
	y
	}
    }

is.na.ratetable  <- function(x)
    structure(is.na(as.vector(x)), dim=dim(x), dimnames=dimnames(x))

Math.ratetable <- function(x, ...) {
    attributes(x) <- attributes(x)[c("dim", "dimnames")]
    NextMethod(.Generic)
    }

Ops.ratetable <- function(e1, e2) {
    #just treat it as an array
    if (nchar(.Method[1])) attributes(e1) <- attributes(e1)[c("dim","dimnames")]
    if (nchar(.Method[2])) attributes(e2) <- attributes(e2)[c("dim","dimnames")]
    NextMethod(.Generic)
    }

as.matrix.ratetable <- function(x, ...) {
    attributes(x) <- attributes(x)[c("dim", "dimnames")]
    x
    }
