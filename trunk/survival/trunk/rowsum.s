#SCCS $Date: 1993-06-29 10:06:54 $ $Id: rowsum.s,v 4.1 1993-06-29 10:06:54 therneau Exp $
rowsum <- function(x, group) {
    if (!is.numeric(x)) stop("x must be numeric")
    if (is.matrix(x)) dd <- dim(x)
    else              dd <- c(length(x), 1)
    n <- dd[1]

    if (length(group) !=n)  stop("Incorrect length for 'group'")

    storage.mode(x) <- 'double'
    temp <- .C("rowsum", dd= as.integer(dd),
			 x = x,
			 as.double(group))
    new.n <- temp$dd[1]
    if (is.matrix(x)){
	new.x <- temp$x[1:new.n,]
	dimnames(new.x) <- list(unique(group), dimnames(x)[[2]])
	}
    else {
	new.x <- temp$x[1:new.n]
	names(new.x) <- unique(group)
	}
    new.x
    }
