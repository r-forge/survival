#  $Id: model.frame.coxph.s,v 4.5 2006-08-28 14:07:49 m015733 Exp $
model.frame.coxph <- function(object, ...) {
    Call <- object$call
    Call[[1]] <- as.name("model.frame")
    Call <- Call[match(c("", "formula", "data", "weights", "subset",
			   "na.action"), names(Call), 0)]
    eval(Call)
    }
