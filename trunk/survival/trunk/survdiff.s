#SCCS $Date: 1992-03-08 20:15:16 $ $Id: survdiff.s,v 4.2 1992-03-08 20:15:16 therneau Exp $
surv.diff <- function(formula, rho=0, riskwt,
		      weights, subset) {
    call <- match.call()
    m <- match.call(expand=F)
    m$... <- m$rho <- m$riskwt <- NULL

    if (!inherits(formula,"formula")) {
	# The dummy function stops an annoying warning message "Looking for
	#  'formula' of mode function, ignored one of mode ..."
	if (inherits(formula,"Surv")) {
	    xx <- function(x) formula(x)
	    m$formula <- xx(paste(deparse(substitute(formula)), 1, sep="~"))
	    }
	else stop("Invalid formula")
	}
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    Terms <- attr(m, 'terms')
    Y <- model.extract(m, response)
    ny <- ncol(Y)
    n <- nrow(Y)
    casewt <- model.extract(m, "weights")
    if (missing(riskwt)) riskwt <- rep(1,n)
    else if (length(riskwt) != n) stop("Wrong length for 'riskwt'")

    if (length(m) == 1) strats <- factor(rep(1,n))
    else if (length(m)==2) strats <- factor(m[[2]])
    else {
	i <- seq(2, length(m))
	temp <- paste("as.character(m[[" ,i,"]]), ", sep='', collapse="', ',")
	strats <- factor(eval(parse(text=paste("paste(", temp, "sep='')" ))))
	}
    strat2 <- as.numeric(strats)
    ngroup <- max(strat2)
    if (ngroup <2) stop ("There is only 1 group")

    if (ny!=2) stop("Surf.diff does not apply to interval time data (yet)")
    ord <- order(Y[,1])

    xx <- .C("survdiff", as.integer(n), as.integer(ngroup), as.double(rho),
		   as.double(Y[ord,1]),
		   as.integer(Y[ord,2]),
		   as.integer(strat2[ord]),
		   observed = double(ngroup),
		   expected = double(ngroup),
		   var.e    = double((ngroup-1) * (ngroup-1)),
		   double(ngroup), double(n))
    temp <- table(strats)
    temp2 <- (xx$observed - xx$expected) [-1]
    chi  <- sum(solve(matrix(xx$var.e, ncol=ngroup-1), temp2) * temp2)
    rval <-list(n= temp, obs = xx$observed, exp=xx$expected,
		    chisq= chi)
    attr(rval, "class") <- 'surv.diff'
    rval
    }
