#SCCS $Date: 1992-04-14 18:06:49 $ $Id: model.frame.default.s,v 4.2 1992-04-14 18:06:49 grill Exp $
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
    .Internal(model.frame(formula, data, na.action, dots), "model_frame")
    }
