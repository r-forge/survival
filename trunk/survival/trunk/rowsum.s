#SCCS $Date: 1998-07-22 08:31:34 $ $Id: rowsum.s,v 4.5 1998-07-22 08:31:34 therneau Exp $
rowsum <- function(x, group) {
    if (!is.numeric(x)) stop("x must be numeric")
    if (any(is.na(group)))  stop("Missing values for 'group'")

    if (is.matrix(x)) {
	if (length(group) != nrow(x)) stop("Incorrect length for 'group'")
	temp <- tapply(x, list(group[row(x)], col(x)), sum)
	dimnames(temp)<- list(dimnames(temp)[[1]], dimnames(x)[[2]])
	temp
	}
    else  {
        if (length(group) !=length(x))  stop("Incorrect length for 'group'")
	else tapply(x, group, sum)
	}
    }
