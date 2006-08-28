# $Id: model.frame.survreg.s,v 1.2 2006-08-28 14:09:21 m015733 Exp $
model.frame.survreg <- function(object, ...) {
    Call <- object$call
    Call[[1]] <- as.name("model.frame")
    Call <- Call[match(c("", "formula", "data", "weights", "subset",
			   "na.action"), names(Call), 0)]
    eval(Call)
    }
