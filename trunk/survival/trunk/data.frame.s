# SCCS $Id: data.frame.s,v 4.1 1993-03-14 19:33:06 therneau Exp $
data.frame <-
function(..., row.names = NULL, check.rows = F, check.names = T, na.strings = 
	"NA")
{
	if(check.rows && missing(row.names))
		data.row.names <- function(current, new, i)
		{
			new <- as.character(new)
			if(any(duplicated(new)))
				return(current)
			if(is.null(current))
				return(new)
			if(all(current == new) || all(current == ""))
				return(new)
			stop(paste(
				"mismatch of row names in elements of \"data.frame\", item",
				i))
		}
	else data.row.names <- function(current, new, i)
		{
			if(is.null(current) && !any(duplicated(new <- 
				as.character(new))))
				new
			else current
		}
	object <- as.expression(substitute(list(...)))[-1]
	x <- list(...)
	n <- length(x)
	value <- list()
	if(n < 1)
		return(structure(list(), class = "data.frame"))
	which <- 1
	vnames <- names(x)
	if(length(vnames) != n)
		vnames <- character(n)
	no.vn <- nchar(vnames) == 0
	vnames <- as.list(vnames)
	nrows <- numeric(n)
	for(i in 1:n) {
		xi <- x[[i]]
		if(inherits(xi, "factor"))
			classi <- "factor"
		else classi <- data.class(xi)[1]
		switch(classi,
			date = ,
			complex = ,
			numeric = {
				nrows[i] <- length(xi)
				value[[which]] <- xi
				if(no.vn[i])
				  vnames[[i]] <- deparse(object[[i]])
				if(missing(row.names) && length(ni <- names(xi)
				  ))
				  row.names <- data.row.names(row.names, ni, i)
				which <- which + 1
			}
			,
			factor = ,
			ordered = ,
			category = {
				if(is.null(class(xi)))
				  xi <- switch(classi,
				    factor = ,
				    category = factor(xi),
				    ordered = ordered(xi),
				    stop("invalid argument"))
				nrows[i] <- length(xi)
				value[[which]] <- xi
				if(no.vn[i])
				  vnames[[i]] <- deparse(object[[i]])
				if(missing(row.names) && length(ni <- names(xi)
				  ))
				  row.names <- data.row.names(row.names, ni, i)
				which <- which + 1
			}
			,
			AsIs = {
				di <- dim(xi)
				if(is.null(di))
				  nrows[i] <- length(xi)
				else nrows[i] <- di[1]
				value[[which]] <- xi
				if(no.vn[i])
				  vnames[[i]] <- deparse(object[[i]])
				which <- which + 1
			}
			,
			Surv = ,
			model.matrix = {
				nrows[i] <- dim(xi)[1]
				value[[which]] <- xi
				if(no.vn[i])
				  vnames[[i]] <- deparse(object[[i]])
				ni <- dimnames(xi)
				if(missing(row.names) && length(ni <- ni[[1]])
				  )
				  row.names <- data.row.names(row.names, ni, i)
				which <- which + 1
			}
			,
			matrix = {
				di <- dim(xi)
				nrows[i] <- di[1]
				ni <- dimnames(xi)
				labels <- ni[[2]]
				ni <- ni[[1]]
				if(length(ni))
				  if(missing(row.names))
				    row.names <- data.row.names(row.names, ni, 
				      i)
				ncol <- di[2]
				cols <- 1:ncol
				if(length(labels) != ncol) {
				  if(no.vn[i])
				    vnames[[i]] <- deparse(object[[i]])
				  labels <- (if(ncol > 1) paste(vnames[[i]], 
				      cols, sep = ".") else vnames[[i]])
				}
				else {
				  if(!no.vn[i])
				    labels <- paste(vnames[[i]], labels, sep = 
				      ".")
				}
				new <- which - 1 + cols
				which <- new[ncol] + 1
				vnames[[i]] <- labels
				length(value) <- new[ncol]
				for(i in cols)
				  value[[new[i]]] <- as.variable(xi[, i])
			}
			,
			design = ,
			data.frame = {
				ni <- attr(xi, "row.names")
				nrows[i] <- length(ni)
				if(missing(row.names) && length(ni) && !all(
				  as.character(1:length(ni)) == ni))
				  row.names <- data.row.names(row.names, ni, i)
				nnew <- length(xi)
				new <- seq(from = which, length = nnew)
				value[new] <- xi
				vnames[[i]] <- names(xi)
				which <- new[nnew] + 1
			}
			,
			list = {
				nnew <- length(xi)
				ni <- names(xi)
				if(any(nchar(ni)) == 0)
				  stop("lists for data frames must have names")
				call <- vector("call", nnew + 1)
				call[[1]] <- as.name("data.frame")
				for(i in 1:length(ni))
				  call[[i + 1]] <- as.name(ni[i])
				if(!check.names)
				  call$check.names <- F
				call$na.strings <- na.strings
				xi <- eval(call, xi)
				new <- seq(from = which, length = length(xi))
				value[new] <- xi
				if(no.vn[i])
				  vnames[[i]] <- names(xi)
				else vnames[[i]] <- paste(vnames[[i]], names(xi
				    ), sep = ".")
				which <- new[nnew] + 1
				ni <- attr(xi, "row.names")
				nrows[i] <- length(ni)
				if(missing(row.names) && length(ni) && !all(
				  as.character(1:length(ni)) == ni))
				  row.names <- data.row.names(row.names, ni, i)
			}
			,
			character = ,
			logical = {
				xi <- factor(xi, exclude = if(classi == 
				  "character") na.strings else NA)
				nrows[i] <- length(xi)
				value[[which]] <- xi
				if(no.vn[i])
				  vnames[[i]] <- deparse(object[[i]])
				which <- which + 1
			}
			,
			NULL = {
			}
			,
			stop(paste("class \"", classi, "\" not defined", sep = 
				"")))
	}
	nrows <- unique(nrows)
	if(length(nrows) > 1)
		stop(paste("arguments imply differing number of rows:", paste(
			nrows, collapse = ", ")))
	vnames <- unlist(vnames)
	noname <- nchar(vnames) == 0
	if(any(noname))
		vnames[noname] <- paste("Var", 1:length(vnames), sep = ".")[
			noname]
	if(check.names)
		vnames <- make.names(vnames)
	names(value) <- vnames
	n.row.names <- length(row.names)
	if(n.row.names == 0)
		row.names <- 1:nrows
	else if(n.row.names != nrows) {
		if(n.row.names != 1)
			stop("length of row.names should be 1 or the same as the number of rows "
				)
		if(is.character(row.names))
			row.names <- match(row.names, vnames, 0)
		if(row.names < 1 || row.names > length(vnames))
			stop("row.names should specify one of the variables")
		i <- row.names
		row.names <- value[[i]]
		value <- value[ - i]
	}
	row.names <- as.character(row.names)
	if(any(duplicated(row.names)))
		stop("row.names must be unique")
	attr(value, "row.names") <- row.names
	attr(value, "class") <- "data.frame"
	value
}
