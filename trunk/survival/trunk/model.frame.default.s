#SCCS $Date: 1992-06-10 20:47:36 $ $Id: model.frame.default.s,v 4.3 1992-06-10 20:47:36 therneau Exp $
# Only change -- look to options() for the default na.action
#
model.frame.default <- function(formula, data = NULL, na.action = na.fail, ...)
{
    if(missing(formula)) {
        if(!missing(data) && inherits(data, "data.frame") && length(
            attr(data, "terms")))
            return(data)
        formula <- as.formula(data)
        }
    else if(missing(data) && inherits(formula, "data.frame")) {
        if(length(attr(formula, "terms")))
            return(formula)
        data <- formula
        formula <- as.formula(data)
        }
    if(missing(na.action)) {
	if (!is.null(tj <- attr(data, "na.action")))  na.action <- tj
	else if (!is.null(tj <- options("na.action")[[1]])) na.action <- tj
	}

    if(!inherits(formula, "terms"))
        formula <- terms(formula, data = data)
    dots <- substitute(list(...))
# Work around a bug in Splus 3.0 -- doesn't check that the objects are
#   all the same length.  Later fixed in Bell S.  Someday replace all of
#   this by the original, commented line below.
#   .Internal(model.frame(formula, data, na.action, dots), "model_frame")
    temp <- .Internal(model.frame(formula, data, na.action, dots), "model_frame")
    if (is.matrix(temp[[1]])) n <- nrow(temp[[1]])
	else                  n <- length(temp[[1]])
    j <- 1
    for (i in temp[-1]) {
	j <- j+1
	mat <- is.matrix(i)
	if ((mat && n!=nrow(i)) ||  (!mat && n!=length(i)))
	    stop(paste("Length mismatch, term", j))
	}
    temp
    }
