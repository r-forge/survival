# $Date: 2006-08-28 18:12:46 $ $Id$
survfit.formula <- function(formula, data, weights, subset, 
			    na.action, ...) {
    Call <- match.call()

    # create a copy of the call that has only the arguments we want,
    #  and use it to call model.frame()
    temp <- Call[c(1, match(c("formula", "data", "weights", "subset", 
                              "na.action"), names(Call), nomatch=0))]
    temp[[1]] <- as.name("model.frame")
    if (is.R()) m <- eval.parent(temp)
    else        m <- eval(temp, sys.parent())
    
    Terms <- terms(formula, 'strata')
    ord <- attr(Terms, 'order')
    if (length(ord) & any(ord !=1))
	    stop("Interaction terms are not valid for this function")

    n <- nrow(m)
    Y <- model.extract(m, response)
    if (!is.Surv(Y)) stop("Response must be a survival object")

    casewt <- model.extract(m, "weights")

    if (is.null(casewt)) casewt <- rep(1,n)

    if (!is.null(attr(Terms, 'offset'))) warning("Offset term ignored")

    ll <- attr(Terms, 'term.labels')
    if (length(ll) == 0) X <- factor(rep(1,n))  # ~1 on the right
    else X <- strata(m[ll])
    
    temp <- survfit.km(X, Y, casewt, ...)
    oldClass(temp) <- "survfit"
    if (!is.null(attr(m, 'na.action')))
	    temp$na.action <- attr(m, 'na.action')

    temp$call <- Call
    temp
    }









