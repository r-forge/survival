#SCCS $Date: 1992-03-30 02:37:45 $ $Id: na.omit.s,v 4.2 1992-03-30 02:37:45 therneau Exp $
na.omit <- function(frame)  {
    n <- length(frame)
    omit <- FALSE
    vars <- seq(length = n)
    for(j in vars) {
	x <- frame[[j]]
	if(!is.atomic(x)) next
    # variables are assumed to be either some sort of matrix, numeric or cat'y
	x <- is.na(x)
	d <- dim(x)
	if(is.null(d) || length(d) != 2)
		omit <- omit | x
	else {
	    for(ii in 1:d[2])
		    omit <- omit | x[, ii]
	    }
	}
    xx <- frame[!omit,  , drop = F]
    if (any(omit)) {
	temp <- seq(omit)[omit]
	names(temp) <- row.names(frame)[omit]
	attr(temp, 'class') <- 'omit'
	attr(xx, "na.action") <- temp
	}
    xx
    }

naprint.omit <- function(x)
    paste(length(x), "deleted due to missing")

# Put the missing values back into a vector.
#   And be careful about the labels too.
naresid.omit <- function(omit, x) {
    if (!length(omit) || !is.numeric(omit))
	stop("Invalid argument for 'omit'")

    if (is.matrix(x)) {
	n <- nrow(x)
	keep <- rep(NA,n+ length(omit))
	keep[-omit] <- 1:n
	x <- x[keep,,drop=F]
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
