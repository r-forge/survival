#SCCS $Date: 1992-04-14 18:06:06 $ $Id: print.data.frame.s,v 4.2 1992-04-14 18:06:06 grill Exp $
#
# In order to get objects with attributes to print correctly, I replace the
#   call to "as.matrix" with a copy of as.matrix.data.frame, one that knows
#   its output is character, and so calls the appropriate as.character routine
#
print.data.frame <-
function(x, ..., quote = F, right = T)
{
    if(!length(x))
        cat("NULL data frame with", length(row.names(x)), "rows\n")

    else if(!length(row.names(x))) {
        print.atomic(names(x), quote = F)
        cat("< 0 rows> \n")
        }
    else {
	DD <- dim(x)
	dn <- dimnames(x)
	collabs <- as.list(dn[[2]])
	class(x) <- NULL
	p <- DD[2]
	n <- DD[1]
	non.numeric <- non.atomic <- F
	for(j in 1:p) {
	    xj <- x[[j]]

	    # 4 Line below are the ones I added
	    if (inherits(xj, "data.frame"))
		xj <- x[[j]] <- as.matrix(x[[j]])
	    else if (!is.null(class(xj)))
		xj <- x[[j]] <- as.character(x[[j]])
	    if(length(dj <- dim(xj)) == 2 && dj[2] > 1) {
		dnj <- dimnames(xj)[[2]]
		collabs[[j]] <- paste(collabs[[j]], if(length(dnj)) dnj
		     else seq(1:dj[2]), sep = ".")
		}
	    if(length(levels(xj)) > 0 || !is.numeric(xj) )
		non.numeric <- TRUE
	    if(!is.atomic(xj))
		non.atomic <- TRUE
	    }
	if(non.atomic) {
	    for(j in 1:p) {
		xj <- x[[j]]
		if(is.recursive(xj)) {
		    }
		else x[[j]] <- as.list(as.vector(xj))
		}
	    }
	else if(non.numeric) {
	    for(j in 1:p) {
		xj <- x[[j]]
		if(length(ll <- levels(xj)))
		    x[[j]] <- ll[unclass(xj)]
		else x[[j]] <- format(xj)
		}
	    }
	x <- unlist(x, recursive = F)
	dim(x) <- c(n, length(x)/n)
	dimnames(x) <- list(dn[[1]], unlist(collabs))
	class(x) <- "matrix"
	print(x, quote=quote, right=right)
	}
    }
