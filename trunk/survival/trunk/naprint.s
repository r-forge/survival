naprint <- function(x, ...)
    UseMethod("naprint")

naprint.default <- function(...)
    return("")
