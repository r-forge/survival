# SCCS $Id: model.frame.survreg.s,v 1.1 1998-11-25 21:08:04 therneau Exp $
model.frame.survreg <- function(object, ...) {
    Call <- object$call
    Call[[1]] <- as.name("model.frame")
    Call <- Call[match(c("", "formula", "data", "weights", "subset",
			   "na.action"), names(Call), 0)]
    eval(Call)
    }
