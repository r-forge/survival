# $Id: rowsum.s,v 5.2 2006-08-28 15:36:44 m015733 Exp $
rowsum <- function(x, group ) {
    if (!is.numeric(x)) stop("x must be numeric")
    if (any(is.na(group)))  stop("Missing values for 'group'")

    if (is.matrix(x)) {
	n <- nrow(x)
	if (length(group) != n) stop("Incorrect length for 'group'")
	tapply(x, list(rep(group, ncol(x)), col(x)), sum)
        }
    else {
	n <- length(x)
	if (length(group) !=n)  stop("Incorrect length for 'group'")
	tapply(x, group, sum)
        }
    }
