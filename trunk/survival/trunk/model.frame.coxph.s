#  SCCS  $Id: model.frame.coxph.s,v 4.4 1999-02-21 16:24:17 therneau Exp $
model.frame.coxph <- function(object, ...) {
    Call <- object$call
    Call[[1]] <- as.name("model.frame")
    Call <- Call[match(c("", "formula", "data", "weights", "subset",
			   "na.action"), names(Call), 0)]
    eval(Call)
    }
