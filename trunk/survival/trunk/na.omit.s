#SCCS $Date: 1992-03-04 16:48:10 $ $Id: na.omit.s,v 4.1 1992-03-04 16:48:10 therneau Exp $
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
	attr(xx, "omit") <- temp
	}
    xx
    }
