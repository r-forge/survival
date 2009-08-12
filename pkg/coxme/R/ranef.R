ranef <- function(object, ...) {
    UseMethod("ranef")
    }
random.effects <- function(object,...) {
    UseMethood("ranef")
    }

fixef <- function(object, ...) {
    UseMethod("fixef")
    }
fixed.effects <- function(object, ...) {
    UseMethod("fixef")
    }

fixef.coxme <- function(object, ...)
    object$coefficients$fixed
