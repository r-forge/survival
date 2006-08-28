# $Id: survfit.Surv.s,v 1.3 2006-08-28 18:08:31 m015733 Exp $
survfit.Surv <- function (formula, data, weights, subset,
			  na.action, call, ...) {
    #  turn the first arg into a real formula
    #avoid a waring about "looking for function formula, ignored..."
    tfun <- function(x)
            formula(paste(x, "~1"))
    formula <- tfun(deparse(substitute(formula)))

    # You might think that having modified the "Surv" object into a
    #  formula I could now call survfit.formula, but it isn't so.  One
    #  loses the odd matching of arguments necessary for a model.frame
    #  call, which depend on being only 1 call level below a UseMethod().  
    # And I can't call UseMethod() again, because it passes forward the
    #  unchanged argument list. So from here on down is a repeat 
    #  of survfit.formula, line by line.

    if (missing(call)) call <- match.call()

    m <- match.call(expand=TRUE)
    m <- m[match(c("", "data", "weights", "subset", "na.action"),
		 names(m), nomatch=0)]
    
    Terms <- terms(formula, 'strata')
    ord <- attr(Terms, 'order')
    if (length(ord) & any(ord !=1))
	    stop("Interaction terms are not valid for this function")
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())

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

    temp$call <- call
    temp
    }









