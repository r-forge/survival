# $Date: 2006-08-28 14:12:46 $ $Id: naresid.omit.s,v 4.3 2006-08-28 14:12:46 m015733 Exp $
naresid.omit <- function(omit, x) {
    if (!length(omit) || !is.numeric(omit))
	stop("Invalid argument for 'omit'")

    if (is.matrix(x)) {
	n <- nrow(x)
	keep <- rep(NA,n+ length(omit))
	keep[-omit] <- 1:n
	x <- x[keep,,drop=FALSE]
	temp <- dimnames(x)[[1]]
	if (length(temp)) {
	    temp[omit] <- names(omit)
	    dimnames(x) <- list(temp, dimnames(x)[[2]])
	    }
	x
	}
    else {
	n <- length(x)
	keep <- rep(NA,n+ length(omit))
	keep[-omit] <- 1:n
	x <- x[keep]
	temp <- names(x)
	if (length(temp)) {
	    temp[omit] <- names(omit)
	    names(x) <- temp
	    }
	x
	}
    }
