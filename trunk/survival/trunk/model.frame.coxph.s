model.frame.coxreg <- function(object, ...) {
    Call <- object$call
    Call[[1]] <- as.name("model.frame")
    Call <- Call[match(c("", "formula", "data", "subset", "na.action"),
		     names(Call), 0)]
    eval(Call)
    }
