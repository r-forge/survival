# SCCS $Id: format.Surv.s,v 4.3 1995-05-30 16:47:04 therneau Exp $
#
# These two functions operate with the newer data.frame code, found on statlib
#   (Eventually part of S, I assume)
#
format.Surv <- function(x) format(as.character.Surv(x))

# The better definition for this is
#   "as.data.frame.Surv <- as.data.frame.model.matrix"
# but, since not everyone has installed the new data.frame software yet, this
# will fail when it can't find as.data.frame.model.matrix.
# So, for this release of survival, the function below is an exact copy of
# as.data.frame.model.matrix, as found on statlib 9/94.
#
#  Changed 5/30/95:  there is a bug in the code I copied: if there are no
#     row names then row.names == character(0), not NULL
as.data.frame.Surv <-
    function(x, row.names = NULL, optional = F)
    {
	d <- dim(x)
	nrows <- d[[1]]
	dn <- dimnames(x)
	row.names <- dn[[1]]
	value <- list(x)
#       if(!is.null(row.names)) {
	if(length(row.names)>0) {
		row.names <- as.character(row.names)
		if(length(row.names) != nrows)
			stop(paste("supplied", length(row.names), 
				"names for a data frame with", nrows, "rows"))
	}
	else if(optional)
		row.names <- character(nrows)
	else row.names <- as.character(seq(length = nrows))
	if(!optional)
		names(value) <- deparse(substitute(x))[[1]]
	attr(value, "row.names") <- row.names
	class(value) <- "data.frame"
	value
}
